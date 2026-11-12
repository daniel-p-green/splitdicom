"""How much smaller can the total get? Measures, on the CT sample series:

  1. Lossless headroom: stronger coder (JPEG XL, max effort) and inter-slice prediction
     vs per-slice JPEG-LS, the baseline already in the DICOM standard.
  2. Where the bits go: the share of lossless bits spent on fine-scale noise, estimated by
     coding a lightly denoised copy losslessly.
  3. Single-file lossy sizes (nothing else kept): bounded-error coding at several bounds.

usage: headroom.py <data_root> <out.json> [step]     (step: use every Nth slice, default 4)
Code-stream bytes only. Ratios are against 16-bit uncompressed pixels.
"""
import sys, os, json, collections, warnings
import numpy as np, pydicom, imagecodecs as ic
from scipy import ndimage
warnings.filterwarnings("ignore")
SER = {"CT chest 2.0mm": "ohif/dcm/scoord3d-and-scoord/scoord-bounding-box", "CT whole-body 5.0mm": "ohif/dcm/Juno", "CT 3.75mm (PET/CT)": "ohif/dcm/acrin"}
jxl = lambda a, e=9: len(ic.jpegxl_encode(a, lossless=True, effort=e))


def load(root, d):
    g = collections.defaultdict(list)
    for dp, _, fs in os.walk(os.path.join(root, d)):
        for f in fs:
            try: ds = pydicom.dcmread(os.path.join(dp, f))
            except Exception: continue
            if ds.get("Modality") == "CT" and "PixelData" in ds and "ImagePositionPatient" in ds: g[ds.SeriesInstanceUID].append(ds)
    sl = max(g.values(), key=len); sl.sort(key=lambda s: float(s.ImagePositionPatient[2]))
    return [s.pixel_array.astype(np.int64) for s in sl]


def main():
    root, out = sys.argv[1], sys.argv[2]; step = int(sys.argv[3]) if len(sys.argv) > 3 else 4; res = {}
    for name, d in SER.items():
        vol = load(root, d); t = collections.Counter()
        for k in range(1, len(vol), step):
            a = vol[k]; off = int(a.min()); u = (a - off).astype(np.uint16); t["raw"] += a.size * 2; t["n"] += 1
            t["jpegls"] += len(ic.jpegls_encode(u)); x = jxl(u); t["jxl"] += x
            # inter-slice prediction: code the difference from the previous slice when that is smaller
            dlt = a - vol[k - 1]; t["jxl_interslice"] += min(x, jxl((dlt - int(dlt.min())).astype(np.uint16)))
            # where the bits go: a lightly smoothed copy, coded losslessly
            for label, sm in (("gauss0.7", ndimage.gaussian_filter(u.astype(np.float32), 0.7)), ("median3", ndimage.median_filter(u, 3))):
                t["smooth_" + label] += jxl(np.clip(np.round(sm), 0, 65535).astype(np.uint16))
            # single-file lossy, residual discarded
            for n in (2, 4, 8, 16):
                t[f"near{n}_jpegls"] += len(ic.jpegls_encode(u, level=n))
                q = 2 * n + 1; t[f"near{n}_quant_jxl"] += jxl(((u.astype(np.int64) + n) // q).astype(np.uint16))
        r = {"slices": t["n"], "ratio": {k: round(t["raw"] / v, 2) for k, v in t.items() if k not in ("raw", "n")},
             "bits_per_pixel": {k: round(16 * v / t["raw"], 2) for k, v in t.items() if k not in ("raw", "n")}}
        r["jxl_gain_over_jpegls_pct"] = round(100 * (1 - t["jxl"] / t["jpegls"]), 1)
        r["interslice_gain_over_jpegls_pct"] = round(100 * (1 - t["jxl_interslice"] / t["jpegls"]), 1)
        r["share_of_lossless_bits_removed_by_light_smoothing_pct"] = {k[7:]: round(100 * (1 - t[k] / t["jxl"]), 1) for k in t if k.startswith("smooth_")}
        res[name] = r; print(name, json.dumps(r), flush=True)
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
