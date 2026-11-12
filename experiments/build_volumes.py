"""Build the NIfTI volumes used in notes/03: originals, a random +/-1 control, and
bounded-error base layers (same JPEG-LS near-lossless coding as splitdicom.py).

usage: build_volumes.py <data_root> <out_dir>

Outputs <series>_<variant>.nii.gz for series chest2mm, wbtorso, petct375 and variants
orig, rand1, 2, 4, 8, 16. The whole-body series is cropped to slices 46..133 (88 torso
slices) because the full volume exceeds 8 GB of RAM in TotalSegmentator.
"""
import sys, os, collections, warnings, json
import numpy as np, pydicom, nibabel as nib, imagecodecs as ic
warnings.filterwarnings("ignore")
SER = {"chest2mm": "ohif/dcm/scoord3d-and-scoord/scoord-bounding-box", "wb5mm": "ohif/dcm/Juno", "petct375": "ohif/dcm/acrin"}
CROP = {"wb5mm": ("wbtorso", 46, 134)}       # series -> (output name, first slice, one past last slice)
NEARS = [2, 4, 8, 16]
RAND_SEED = 0                                 # one generator, consumed in SER order over full volumes


def main():
    root, out = sys.argv[1], sys.argv[2]; os.makedirs(out, exist_ok=True); info = {}
    rng = np.random.default_rng(RAND_SEED)
    for name, d in SER.items():
        groups = collections.defaultdict(list)
        for dp, _, fs in os.walk(os.path.join(root, d)):
            for f in fs:
                try: ds = pydicom.dcmread(os.path.join(dp, f))
                except Exception: continue
                if ds.get("Modality") == "CT" and "PixelData" in ds and "ImagePositionPatient" in ds: groups[ds.SeriesInstanceUID].append(ds)
        if not groups: sys.exit(f"error: no CT series found under {os.path.join(root, d)}")
        sl = max(groups.values(), key=len)
        iop = np.array(sl[0].ImageOrientationPatient, float); rc, cc = iop[:3], iop[3:]; nrm = np.cross(rc, cc)
        sl.sort(key=lambda s: float(np.dot(nrm, np.array(s.ImagePositionPatient, float))))
        ps = [float(x) for x in sl[0].PixelSpacing]; p0 = np.array(sl[0].ImagePositionPatient, float); p1 = np.array(sl[-1].ImagePositionPatient, float)
        aff = np.eye(4); aff[:3, 0] = rc * ps[1]; aff[:3, 1] = cc * ps[0]; aff[:3, 2] = (p1 - p0) / (len(sl) - 1); aff[:3, 3] = p0; aff[:2, :] *= -1
        vols = {k: [] for k in ["orig"] + NEARS}; raw = 0; size = collections.Counter()
        for s in sl:
            a = s.pixel_array; slope, icpt = float(s.get("RescaleSlope", 1)), float(s.get("RescaleIntercept", 0))
            off = int(a.min()) if a.min() < 0 else 0; u = (a.astype(np.int64) - off).astype(np.uint16)
            vols["orig"].append(np.round(a * slope + icpt).astype(np.int16)); raw += a.size * 2
            for k in NEARS:
                cs = ic.jpegls_encode(u, level=k); b = ic.jpegls_decode(cs).astype(np.int64) + off
                assert np.abs(b - a).max() <= k; size[k] += len(cs)
                vols[k].append(np.round(b * slope + icpt).astype(np.int16))
        full = {k: np.stack(v).transpose(2, 1, 0) for k, v in vols.items()}
        full["rand1"] = (full["orig"].astype(np.int32) + rng.integers(-1, 2, full["orig"].shape)).astype(np.int16)
        outname, z0, z1 = CROP.get(name, (name, 0, full["orig"].shape[2]))
        for k, vol in full.items():
            nib.save(nib.Nifti1Image(vol, aff).slicer[:, :, z0:z1], os.path.join(out, f"{outname}_{k}.nii.gz"))
        info[outname] = {"source_series": name, "slices": len(sl), "slices_used": [z0, z1], "spacing": [ps[1], ps[0], float(np.linalg.norm(aff[:3, 2]))],
                         "desc": str(sl[0].get("SeriesDescription")), "kernel": str(sl[0].get("ConvolutionKernel")), "slope": float(sl[0].get("RescaleSlope", 1)),
                         "base_codestream_ratio_full_series": {k: round(raw / size[k], 2) for k in NEARS}}
        print(outname, info[outname], flush=True)
    json.dump(info, open(os.path.join(out, "info.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
