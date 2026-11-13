"""Nodule test at equal size: bounded-error coding vs lossy JPEG 2000 and lossy JPEG XL.

Same patch-centred detection and adaptive volumetry as nodule_test.py. Variants:
  orig   original
  q8     uniform quantizer (+/-8) + lossless JPEG XL
  j2k    lossy JPEG 2000 (9/7 wavelet), one quality setting per scan
  jxl    lossy JPEG XL (VarDCT), one quality setting per scan
For each scan the lossy settings are tuned on three slices through its first nodule so that their
bytes are at least the q8 bytes on those slices. Bytes on the tuning slices are recorded.

usage: nodule_test_lossy.py <annotations.csv> <dicom_zip_dir> <model.pt> <selected.txt> <out.jsonl>
"""
import sys, os, csv, json, hashlib, collections, time, warnings
import numpy as np, torch, imagecodecs as ic
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nodule_test as nt
from compare_lossy import tune

VARIANTS = ["orig", "q8", "j2k", "jxl"]
N = 8


class LossyScan(nt.Scan):
    def prepare(self, k_center):
        """Tune per-scan lossy settings on three slices around k_center."""
        ks = sorted({min(len(self.sl) - 1, max(0, int(round(k_center)) + d)) for d in (-3, 0, 3)})
        arrs = [self.sl[k].pixel_array.astype(np.int64) for k in ks]
        self.off = min(0, min(int(a.min()) for a in arrs))
        # bit depth must cover every slice that will be coded, so take it from the stored range
        s = self.sl[0]; self.bits = int(s.BitsStored) + (1 if self.off < 0 else 0)
        self.bits = min(16, max(self.bits, max(int((a - self.off).max()) for a in arrs).bit_length()))
        us = [(a - self.off).astype(np.uint16) for a in arrs]; q = 2 * N + 1
        target = sum(len(ic.jpegxl_encode(((u.astype(np.int64) + N) // q).astype(np.uint16), lossless=True, effort=7)) for u in us)
        self.enc = {"j2k": lambda u, p: ic.jpeg2k_encode(u, level=p, codecformat="J2K", bitspersample=self.bits, reversible=False),
                    "jxl": lambda u, p: ic.jpegxl_encode(u, distance=p, effort=7, bitspersample=self.bits, lossless=False)}
        self.dec = {"j2k": ic.jpeg2k_decode, "jxl": ic.jpegxl_decode}
        self.setting = {"j2k": tune(self.enc["j2k"], 20.0, 120.0, +1, us, target), "jxl": tune(self.enc["jxl"], 0.01, 25.0, -1, us, target)}
        self.tuning = {"slices": ks, "bits": self.bits, "q8_bytes": target, "settings": self.setting,
                       "bytes": {c: sum(len(self.enc[c](u, self.setting[c])) for u in us) for c in self.enc}}

    def slice_hu(self, kind, k):
        key = (kind, k)
        if key not in self.cache:
            a = self.sl[k].pixel_array.astype(np.int64)
            if kind == "orig": b = a
            else:
                lim = (1 << self.bits) - 1
                u = np.clip(a - self.off, 0, lim).astype(np.uint16)
                if kind == "q8":
                    q = 2 * N + 1; b = ((u.astype(np.int64) + N) // q) * q + self.off
                else:
                    b = np.asarray(self.dec[kind](self.enc[kind](u, self.setting[kind]))).astype(np.int64) + self.off
                if int(u.max()) == lim and int((a - self.off).max()) > lim: self.clipped = True
                err = np.abs(b - a); self.err[kind] = max(self.err[kind], int(err.max())); self.sq[kind] += float((err.astype(np.float64) ** 2).sum()); self.npx[kind] += a.size
            self.cache[key] = (b * self.slope + self.icpt).astype(np.float32)
        return self.cache[key]


def main():
    ann_csv, zdir, weights, sel, out = sys.argv[1:6]
    torch.set_num_threads(2)
    ann = collections.defaultdict(list)
    for r in csv.DictReader(open(ann_csv)):
        info = {a: b for a, b in r.items() if a not in ("seriesuid", "coordX", "coordY", "coordZ", "diameter_mm")}
        ann[r["seriesuid"]].append((float(r["coordX"]), float(r["coordY"]), float(r["coordZ"]), float(r["diameter_mm"]), info))
    h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
    run = {"code_sha256": [h(os.path.abspath(__file__)), h(nt.__file__)], "weights_sha256": h(weights), "annotations_sha256": h(ann_csv), "variants": VARIANTS, "bound": N,
           "patch": list(nt.PATCH), "versions": {"torch": torch.__version__, "imagecodecs": ic.__version__, "numpy": np.__version__}}
    fp = hashlib.sha256(json.dumps(run, sort_keys=True).encode()).hexdigest(); done = set()
    if os.path.exists(out) and os.path.getsize(out):
        recs = [json.loads(l) for l in open(out) if l.strip()]
        if [r["run"]["fingerprint"] for r in recs if "run" in r] != [fp]: sys.exit("error: output file belongs to a different run configuration")
        done = {r["scan"] for r in recs if "scan" in r}
    det = nt.build_detector(weights); f = open(out, "a")
    if not done and os.path.getsize(out) == 0: f.write(json.dumps({"run": dict(run, fingerprint=fp)}) + "\n"); f.flush()
    for si, uid in enumerate([l.strip() for l in open(sel) if l.strip()]):
        if uid in done or not ann[uid]: continue
        p = os.path.join(zdir, uid + ".zip")
        if not os.path.exists(p): f.write(json.dumps({"scan": uid, "skipped": "scan file not found"}) + "\n"); f.flush(); continue
        t = time.time()
        try: scan = LossyScan(p, si)
        except Exception as e:
            f.write(json.dumps({"scan": uid, "skipped": str(e)}) + "\n"); f.flush(); continue
        scan.err = collections.Counter(); scan.sq = collections.Counter(); scan.npx = collections.Counter(); scan.clipped = False
        scan.prepare(scan.world_to_vox(*ann[uid][0][:3])[2])
        nods = []
        for ni, (x, y, z, d, info) in enumerate(ann[uid]):
            c, r, k = scan.world_to_vox(x, y, z)
            others = [((x2 - x), (y2 - y), (z2 - z), d2) for j, (x2, y2, z2, d2, _) in enumerate(ann[uid]) if j != ni]
            rec = {"diameter_mm": d, "info": info, "others_lps_mm": others, "variants": {}}; ref = None; aref = None
            at = nt.adaptive_threshold(scan, c, r, k, d)
            rec["adaptive"] = {"threshold": at[0], "core_mean": at[1], "background_median": at[2]} if at else None
            for kind in VARIANTS:
                dets = nt.detect(det, nt.model_patch(scan, kind, c, r, k))
                vol, mask, mean = nt.volumetry(scan, kind, c, r, k, d, ref)
                if kind == "orig": ref = mask
                rec["variants"][kind] = {"dets": [x_ for x_ in dets if x_[3] >= 0.05], "vol_mm3": vol, "mean_hu_in_orig_mask": mean}
                if at:
                    avol, amask, amean = nt.volumetry(scan, kind, c, r, k, d, aref, thr=at[0])
                    if kind == "orig": aref = amask
                    rec["variants"][kind].update({"adaptive_vol_mm3": avol, "adaptive_mean_hu_in_orig_mask": amean})
            nods.append(rec)
        errs = {k: {"max_abs_err": scan.err[k], "rmse": round(float(np.sqrt(scan.sq[k] / scan.npx[k])), 3)} for k in VARIANTS[1:] if scan.npx[k]}
        f.write(json.dumps({"scan": uid, "meta": scan.meta, "tuning": scan.tuning, "pixel_error": errs, "clipped": scan.clipped, "nodules": nods, "sec": round(time.time() - t, 1)}) + "\n"); f.flush()
        print(si, uid[-12:], len(nods), "nodules", round(time.time() - t), "s", scan.tuning["bytes"], scan.tuning["q8_bytes"], errs, flush=True)


if __name__ == "__main__":
    main()
