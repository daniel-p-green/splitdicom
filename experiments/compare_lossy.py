"""At equal total size, is a bounded-error file better than the lossy codecs already in use?

For each CT sample series this builds decoded volumes for:
  q4, q8        bounded-error coding (uniform quantizer + lossless JPEG XL), every pixel within +/-4 or +/-8
  j2k_N         lossy JPEG 2000 (irreversible 9/7 wavelet, OpenJPEG)
  htj2k_N       lossy HTJ2K (OpenJPH)
  jxl_N         lossy JPEG XL (VarDCT)
where each lossy codec uses one quality setting per series, tuned so its total bytes match the
bounded-error total for N = 4 and N = 8. It records bytes and pixel-error statistics, and writes
NIfTI volumes for experiments/compare.py (segmentation agreement).

usage: compare_lossy.py <data_root> <out_dir> <out.json>
"""
import sys, os, json, collections, warnings
import numpy as np, pydicom, nibabel as nib, imagecodecs as ic
warnings.filterwarnings("ignore")
from build_volumes import SER, CROP
BOUNDS = [4, 8]
SUB = 6            # tune quality on every SUB-th slice, then encode everything


def codecs(bits):
    return {
        "j2k":   (lambda u, p: ic.jpeg2k_encode(u, level=p, codecformat="J2K", bitspersample=bits, reversible=False), ic.jpeg2k_decode, 20.0, 120.0, +1),
        "htj2k": (lambda u, p: ic.htj2k_encode(u, level=p, bitspersample=bits, reversible=False), ic.htj2k_decode, 1e-6, 1.0, -1),
        "jxl":   (lambda u, p: ic.jpegxl_encode(u, distance=p, effort=7, bitspersample=bits, lossless=False), ic.jpegxl_decode, 0.01, 25.0, -1),
    }


def tune(enc, lo, hi, sign, slices, target):
    """Log-space bisection for the setting whose bytes on `slices` are just above `target`."""
    best = None
    for _ in range(14):
        mid = float(np.sqrt(lo * hi)); size = sum(len(enc(u, mid)) for u in slices)
        if size >= target: best = (mid, size)
        bigger = size < target                      # need more bytes
        if (sign > 0) == bigger: lo = mid
        else: hi = mid
    return best[0] if best else (hi if sign > 0 else lo)


def main():
    root, out, outjson = sys.argv[1:4]; os.makedirs(out, exist_ok=True); res = {}
    for name, d in SER.items():
        g = collections.defaultdict(list)
        for dp, _, fs in os.walk(os.path.join(root, d)):
            for f in fs:
                try: ds = pydicom.dcmread(os.path.join(dp, f))
                except Exception: continue
                if ds.get("Modality") == "CT" and "PixelData" in ds and "ImagePositionPatient" in ds: g[ds.SeriesInstanceUID].append(ds)
        sl = max(g.values(), key=len)
        iop = np.array(sl[0].ImageOrientationPatient, float); rc, cc = iop[:3], iop[3:]; nrm = np.cross(rc, cc)
        sl.sort(key=lambda s: float(np.dot(nrm, np.array(s.ImagePositionPatient, float))))
        ps = [float(x) for x in sl[0].PixelSpacing]; p0 = np.array(sl[0].ImagePositionPatient, float); p1 = np.array(sl[-1].ImagePositionPatient, float)
        aff = np.eye(4); aff[:3, 0] = rc * ps[1]; aff[:3, 1] = cc * ps[0]; aff[:3, 2] = (p1 - p0) / (len(sl) - 1); aff[:3, 3] = p0; aff[:2, :] *= -1
        slope, icpt = float(sl[0].get("RescaleSlope", 1)), float(sl[0].get("RescaleIntercept", 0))
        arrs = [s.pixel_array.astype(np.int64) for s in sl]
        off = min(0, min(int(a.min()) for a in arrs)); us = [(a - off).astype(np.uint16) for a in arrs]
        bits = max(8, max(int(u.max()) for u in us).bit_length()); raw = sum(u.size * 2 for u in us)
        outname, z0, z1 = CROP.get(name, (name, 0, len(us)))
        r = {"slices": len(us), "bits": bits, "lossless_jxl_ratio": None, "variants": {}}

        def record(tag, recon, nbytes, setting=None):
            e = np.concatenate([(x.astype(np.int64) - u.astype(np.int64)).ravel() for x, u in zip(recon, us)])
            ae = np.abs(e)
            r["variants"][tag] = {"bytes": int(nbytes), "ratio": round(raw / nbytes, 2), "setting": setting, "max_abs_err": int(ae.max()),
                                  "p999_abs_err": float(np.percentile(ae, 99.9)), "rmse": round(float(np.sqrt((e.astype(np.float64) ** 2).mean())), 3),
                                  "pixels_over_8_pct": round(100 * float((ae > 8).mean()), 4), "pixels_over_16_pct": round(100 * float((ae > 16).mean()), 5)}
            vol = np.stack([np.round((x.astype(np.int64) + off) * slope + icpt).astype(np.int16) for x in recon]).transpose(2, 1, 0)
            nib.save(nib.Nifti1Image(vol, aff).slicer[:, :, z0:z1], os.path.join(out, f"{outname}_{tag}.nii.gz"))
            print(name, tag, r["variants"][tag], flush=True)

        targets = {}
        for n in BOUNDS:
            q = 2 * n + 1; qi = [((u.astype(np.int64) + n) // q) for u in us]
            sizes = [len(ic.jpegxl_encode(x.astype(np.uint16), lossless=True, effort=7)) for x in qi]
            targets[n] = sizes; record(f"q{n}", [x * q for x in qi], sum(sizes))
        cs = codecs(bits); sub = list(range(0, len(us), SUB))
        for cname, (enc, dec, lo, hi, sign) in cs.items():
            for n in BOUNDS:
                p = tune(enc, lo, hi, sign, [us[i] for i in sub], sum(targets[n][i] for i in sub))
                streams = [enc(u, p) for u in us]
                record(f"{cname}_{n}", [np.asarray(dec(s)) for s in streams], sum(len(s) for s in streams), p)
        res[outname] = r
        json.dump(res, open(outjson, "w"), indent=1)


if __name__ == "__main__":
    main()
