"""Nodule test: do annotated lung nodules survive the bounded-error base layer?

For each annotated nodule (LUNA16 annotations on LIDC-IDRI scans):
  1. Build variants of the surrounding slices: original, random +/-1 control, and
     JPEG-LS near-lossless base layers at several bounds (same codec as splitdicom.py).
  2. Run the MONAI lung_nodule_ct_detection RetinaNet on one 128x128x64 patch
     (90 x 90 x 80 mm) centred on the nodule. Record the matched detection score.
  3. Measure the nodule with fixed-threshold volumetry (-400 HU) on the native grid, and with
     an adaptive threshold half-way between the nodule core and surrounding lung (needed for
     ground-glass nodules, which sit below -400 HU). The threshold is set on the original.

Extra columns in the annotation CSV are passed through to the output.

usage: nodule_test.py <annotations.csv> <dicom_zip_dir> <model.pt> <selected.txt> <out.jsonl>
"""
import sys, os, io, csv, json, hashlib, zipfile, collections, warnings, time
import numpy as np, pydicom, torch, imagecodecs as ic
from scipy import ndimage
warnings.filterwarnings("ignore")
from monai.apps.detection.networks.retinanet_network import RetinaNet, resnet_fpn_feature_extractor
from monai.apps.detection.networks.retinanet_detector import RetinaNetDetector
from monai.apps.detection.utils.anchor_utils import AnchorGeneratorWithAnchorShape
from monai.networks.nets import resnet

VARIANTS = ["orig", "rand1", "4", "8", "16"]
SP = np.array([0.703125, 0.703125, 1.25]); PATCH = (128, 128, 64)
THR = -400.0


def build_detector(weights):
    backbone = resnet.resnet50(spatial_dims=3, n_input_channels=1, conv1_t_stride=[2, 2, 1], conv1_t_size=[7, 7, 7])
    fe = resnet_fpn_feature_extractor(backbone, 3, False, [1, 2], None)
    net = RetinaNet(spatial_dims=3, num_classes=1, num_anchors=3, feature_extractor=fe, size_divisible=[16, 16, 8], use_list_output=False)
    net.load_state_dict(torch.load(weights, map_location="cpu", weights_only=True))
    ag = AnchorGeneratorWithAnchorShape(feature_map_scales=[1, 2, 4], base_anchor_shapes=[[6, 8, 4], [8, 6, 5], [10, 10, 6]])
    det = RetinaNetDetector(network=net, anchor_generator=ag, debug=False, spatial_dims=3, num_classes=1, size_divisible=[16, 16, 8])
    det.set_target_keys(box_key="box", label_key="label")
    det.set_box_selector_parameters(score_thresh=0.02, topk_candidates_per_level=1000, nms_thresh=0.22, detections_per_img=300)
    det.eval()
    return det


def load_series(path):
    z = zipfile.ZipFile(path); sl = []
    for n in z.namelist():
        if not n.lower().endswith(".dcm"): continue
        ds = pydicom.dcmread(io.BytesIO(z.read(n)))
        if ds.get("Modality") == "CT" and "ImagePositionPatient" in ds: sl.append(ds)
    iop = np.array(sl[0].ImageOrientationPatient, float)
    if not np.allclose(iop, [1, 0, 0, 0, 1, 0], atol=1e-3): raise ValueError("non-standard orientation")
    sl.sort(key=lambda s: float(s.ImagePositionPatient[2]))
    zs = np.array([float(s.ImagePositionPatient[2]) for s in sl]); dz = np.diff(zs)
    if len(dz) == 0 or np.abs(dz - np.median(dz)).max() > 0.05: raise ValueError("non-uniform slice spacing")
    return sl, zs, float(np.median(dz))


class Scan:
    def __init__(self, path, seed):
        self.sl, self.zs, self.dz = load_series(path)
        s = self.sl[0]
        self.ps_row, self.ps_col = [float(x) for x in s.PixelSpacing]
        self.x0, self.y0 = float(s.ImagePositionPatient[0]), float(s.ImagePositionPatient[1])
        self.slope, self.icpt = float(s.get("RescaleSlope", 1)), float(s.get("RescaleIntercept", 0))
        self.cache = {}; self.rng = np.random.default_rng(seed); self.bytes = collections.Counter()
        self.meta = {"slice_spacing": round(self.dz, 3), "pixel_spacing": round(self.ps_col, 4), "kernel": str(s.get("ConvolutionKernel", "")),
                     "manufacturer": str(s.get("Manufacturer", "")), "slope": self.slope, "n_slices": len(self.sl)}

    def slice_hu(self, kind, k):
        key = (kind, k)
        if key not in self.cache:
            a = self.sl[k].pixel_array.astype(np.int64)
            if kind == "orig": b = a
            elif kind == "rand1": b = a + np.random.default_rng(k).integers(-1, 2, a.shape)
            else:
                n = int(kind); off = int(a.min()) if a.min() < 0 else 0
                u = (a - off).astype(np.uint16); cs = ic.jpegls_encode(u, level=n)
                b = ic.jpegls_decode(cs).astype(np.int64) + off
                assert np.abs(b - a).max() <= n
                self.bytes[kind] += len(cs); self.bytes["raw" + kind] += a.size * 2
            self.cache[key] = (b * self.slope + self.icpt).astype(np.float32)
        return self.cache[key]

    def block(self, kind, k0, k1):
        return np.stack([self.slice_hu(kind, k) for k in range(k0, k1)])  # [k, row, col]

    def world_to_vox(self, x, y, z):
        return (x - self.x0) / self.ps_col, (y - self.y0) / self.ps_row, (z - self.zs[0]) / self.dz


def model_patch(scan, kind, c, r, k):
    """Resample a PATCH-sized block in RAS order at SP spacing, centred on voxel (c, r, k)."""
    half_z = PATCH[2] * SP[2] / 2 / scan.dz
    k0 = max(0, int(np.floor(k - half_z)) - 1); k1 = min(len(scan.sl), int(np.ceil(k + half_z)) + 2)
    blk = scan.block(kind, k0, k1)
    i = (np.arange(PATCH[0]) - (PATCH[0] - 1) / 2) * SP[0]; j = (np.arange(PATCH[1]) - (PATCH[1] - 1) / 2) * SP[1]
    m = (np.arange(PATCH[2]) - (PATCH[2] - 1) / 2) * SP[2]
    cc = c - i / scan.ps_col            # R axis: increasing index moves toward patient right = decreasing LPS x
    rr = r - j / scan.ps_row            # A axis: increasing index moves anterior = decreasing LPS y
    kk = k + m / scan.dz - k0           # S axis
    K, R, C = np.meshgrid(kk, rr, cc, indexing="ij")          # [z, y, x]
    out = ndimage.map_coordinates(blk, [K, R, C], order=1, mode="constant", cval=-1024.0)
    return np.ascontiguousarray(out.transpose(2, 1, 0))       # [x, y, z]


def detect(det, patch_hu):
    img = torch.from_numpy(np.clip((patch_hu + 1024.0) / 1324.0, 0, 1).astype(np.float32))[None]
    with torch.no_grad():
        o = det([img], use_inferer=False)[0]
    box = o["box"].numpy(); sc = o["label_scores"].numpy(); res = []
    for b, s in zip(box, sc):
        ctr = (b[:3] + b[3:]) / 2 - (np.array(PATCH) - 1) / 2
        size = (b[3:] - b[:3]) * SP
        res.append([float(ctr[0] * SP[0]), float(ctr[1] * SP[1]), float(ctr[2] * SP[2]), float(s), float((size[0] + size[1]) / 2)] + [float(v) for v in b])
    return res  # offsets in mm (R, A, S) from the nodule centre, score, in-plane box size mm, then the raw box (x1,y1,z1,x2,y2,z2 in patch voxels)


def volumetry(scan, kind, c, r, k, diam, ref_mask=None, thr=THR):
    rad = diam * 0.75 + 1.0
    kz = int(np.ceil(rad / scan.dz)) + 1; ry = int(np.ceil(rad / scan.ps_row)) + 1; rx = int(np.ceil(rad / scan.ps_col)) + 1
    kc, rc, cc = int(round(k)), int(round(r)), int(round(c))
    k0, k1 = max(0, kc - kz), min(len(scan.sl), kc + kz + 1)
    blk = scan.block(kind, k0, k1)[:, max(0, rc - ry):rc + ry + 1, max(0, cc - rx):cc + rx + 1]
    K, R, C = np.meshgrid(np.arange(k0, k1) - k, np.arange(max(0, rc - ry), max(0, rc - ry) + blk.shape[1]) - r,
                          np.arange(max(0, cc - rx), max(0, cc - rx) + blk.shape[2]) - c, indexing="ij")
    dist = np.sqrt((K * scan.dz) ** 2 + (R * scan.ps_row) ** 2 + (C * scan.ps_col) ** 2)
    cand = (blk > thr) & (dist <= rad)
    lab, n = ndimage.label(cand)
    vox = scan.dz * scan.ps_row * scan.ps_col
    if n == 0: return 0.0, None, None
    d = np.where(cand, dist, 1e9); seed = np.unravel_index(np.argmin(d), d.shape)
    mask = lab == lab[seed]
    mean_ref = float(blk[ref_mask].mean()) if ref_mask is not None and ref_mask.shape == blk.shape and ref_mask.any() else float(blk[mask].mean())
    return float(mask.sum() * vox), mask, mean_ref


def adaptive_threshold(scan, c, r, k, diam):
    """Half-way threshold between the nodule core and surrounding lung, measured on the original.

    Fixed -400 HU volumetry cannot see ground-glass nodules, which sit well below -400.
    Returns (threshold, core_mean, background_median) or None when there is no usable contrast.
    """
    rout = max(1.5 * diam, 8.0)
    kz = int(np.ceil(rout / scan.dz)) + 1; ry = int(np.ceil(rout / scan.ps_row)) + 1; rx = int(np.ceil(rout / scan.ps_col)) + 1
    kc, rc, cc = int(round(k)), int(round(r)), int(round(c))
    k0, k1 = max(0, kc - kz), min(len(scan.sl), kc + kz + 1); r0, c0 = max(0, rc - ry), max(0, cc - rx)
    blk = scan.block("orig", k0, k1)[:, r0:rc + ry + 1, c0:cc + rx + 1]
    K, R, C = np.meshgrid(np.arange(k0, k1) - k, np.arange(r0, r0 + blk.shape[1]) - r, np.arange(c0, c0 + blk.shape[2]) - c, indexing="ij")
    dist = np.sqrt((K * scan.dz) ** 2 + (R * scan.ps_row) ** 2 + (C * scan.ps_col) ** 2)
    core = blk[dist <= max(0.3 * diam, 1.0)]
    if core.size == 0: core = blk[dist <= dist.min() + 1e-6]
    shell = blk[(dist >= max(diam, 5.0)) & (dist <= rout) & (blk < -500)]
    if shell.size < 20: return None
    cm, bg = float(core.mean()), float(np.median(shell))
    if cm - bg < 60: return None
    return (cm + bg) / 2, cm, bg


def main():
    ann_csv, zdir, weights, sel, out = sys.argv[1:6]
    torch.set_num_threads(2)
    ann = collections.defaultdict(list)
    for r in csv.DictReader(open(ann_csv)):
        info = {a: b for a, b in r.items() if a not in ("seriesuid", "coordX", "coordY", "coordZ", "diameter_mm")}
        ann[r["seriesuid"]].append((float(r["coordX"]), float(r["coordY"]), float(r["coordZ"]), float(r["diameter_mm"]), info))
    # A run is identified by its code, weights, annotations, and settings. Resuming under a different
    # fingerprint is refused, so one output file never mixes configurations.
    h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
    run = {"code_sha256": h(os.path.abspath(__file__)), "weights_sha256": h(weights), "annotations_sha256": h(ann_csv),
           "variants": VARIANTS, "patch": list(PATCH), "spacing": [float(x) for x in SP], "threshold_hu": THR,
           "versions": {"torch": torch.__version__, "monai": __import__("monai").__version__, "numpy": np.__version__, "pydicom": pydicom.__version__}}
    fingerprint = hashlib.sha256(json.dumps(run, sort_keys=True).encode()).hexdigest()
    done = set()
    if os.path.exists(out) and os.path.getsize(out):
        recs = [json.loads(l) for l in open(out) if l.strip()]
        prior = [r["run"]["fingerprint"] for r in recs if "run" in r]
        if prior != [fingerprint]:
            sys.exit("error: output file was produced by a different run configuration (or has no run record); use a new output path")
        done = {r["scan"] for r in recs if "scan" in r}
    det = build_detector(weights); f = open(out, "a")
    if not done and not (os.path.exists(out) and os.path.getsize(out)):
        f.write(json.dumps({"run": dict(run, fingerprint=fingerprint)}) + "\n"); f.flush()
    for si, uid in enumerate([l.strip() for l in open(sel) if l.strip()]):
        if uid in done: continue
        p = os.path.join(zdir, uid + ".zip")
        if not os.path.exists(p):
            f.write(json.dumps({"scan": uid, "skipped": "scan file not found"}) + "\n"); f.flush(); continue
        t = time.time()
        try: scan = Scan(p, si)
        except Exception as e:
            f.write(json.dumps({"scan": uid, "skipped": str(e)}) + "\n"); f.flush(); continue
        nods = []
        for ni, (x, y, z, d, info) in enumerate(ann[uid]):
            c, r, k = scan.world_to_vox(x, y, z)
            others = [((x2 - x), (y2 - y), (z2 - z), d2) for j, (x2, y2, z2, d2, _) in enumerate(ann[uid]) if j != ni]
            rec = {"diameter_mm": d, "info": info, "others_lps_mm": others, "variants": {}}; ref = None; aref = None
            at = adaptive_threshold(scan, c, r, k, d)
            rec["adaptive"] = {"threshold": at[0], "core_mean": at[1], "background_median": at[2]} if at else None
            for kind in VARIANTS:
                dets = detect(det, model_patch(scan, kind, c, r, k))
                vol, mask, mean = volumetry(scan, kind, c, r, k, d, ref)
                if kind == "orig": ref = mask;
                rec["variants"][kind] = {"dets": [x_ for x_ in dets if x_[3] >= 0.05], "vol_mm3": vol, "mean_hu_in_orig_mask": mean}
                if at:
                    avol, amask, amean = volumetry(scan, kind, c, r, k, d, aref, thr=at[0])
                    if kind == "orig": aref = amask
                    rec["variants"][kind].update({"adaptive_vol_mm3": avol, "adaptive_mean_hu_in_orig_mask": amean})
            nods.append(rec)
        ratios = {k: round(scan.bytes["raw" + k] / scan.bytes[k], 2) for k in VARIANTS if scan.bytes[k]}
        f.write(json.dumps({"scan": uid, "meta": scan.meta, "base_ratio": ratios, "nodules": nods, "sec": round(time.time() - t, 1)}) + "\n"); f.flush()
        print(si, uid[-12:], len(nods), "nodules", round(time.time() - t), "s", flush=True)


if __name__ == "__main__":
    main()
