"""Compare segmentations and HU measurements: variants vs the original (notes/03).

usage: compare.py <vols_dir> <seg_dir> <out.json> <series,series,...> <variant,variant,...>
Every listed series and variant is required: a missing file is an error, not a skip.
"""
import sys, os, json, numpy as np, nibabel as nib
from totalsegmentator.map_to_binary import class_map
MINVOX = 1000


def load(path):
    if not os.path.exists(path): sys.exit(f"error: required file is missing: {path}")
    return np.asarray(nib.load(path).dataobj)


def main():
    vols, segdir, outpath = sys.argv[1:4]; series = sys.argv[4].split(","); variants = sys.argv[5].split(",")
    names = class_map["total"]; out = {}
    for n in series:
        so = load(f"{segdir}/{n}_orig.nii.gz"); vo = load(f"{vols}/{n}_orig.nii.gz").astype(np.float64)
        cnt = np.bincount(so.ravel(), minlength=120); labs = [l for l in range(1, len(cnt)) if cnt[l] >= MINVOX]
        out[n] = {"labels": len(labs)}
        for k in variants:
            sk = load(f"{segdir}/{n}_{k}.nii.gz"); vk = load(f"{vols}/{n}_{k}.nii.gz").astype(np.float64)
            ck = np.bincount(sk.ravel(), minlength=120); d = {}; vol = {}; hu = {}
            for l in labs:
                mo = so == l; inter = np.count_nonzero(mo & (sk == l)); d[l] = 2 * inter / (cnt[l] + ck[l]); vol[l] = 100 * (ck[l] - cnt[l]) / cnt[l]
                hu[l] = abs(vk[mo].mean() - vo[mo].mean())
            dv = np.array(list(d.values())); worst = min(d, key=d.get)
            out[n][k] = {"dice_median": round(float(np.median(dv)), 4), "dice_p5": round(float(np.percentile(dv, 5)), 4), "dice_min": round(float(dv.min()), 4),
                         "worst": names.get(worst, str(worst)), "below_0.95": int((dv < 0.95).sum()), "vol_change_max_pct": round(float(np.max(np.abs(list(vol.values())))), 2),
                         "vol_change_median_pct": round(float(np.median(np.abs(list(vol.values())))), 3),
                         "voxels_changed_pct": round(100 * float(np.count_nonzero(so != sk)) / max(1, np.count_nonzero((so > 0) | (sk > 0))), 3),
                         "organ_meanHU_shift_max": round(float(max(hu.values())), 3), "organ_meanHU_shift_median": round(float(np.median(list(hu.values()))), 3),
                         "rmse_HU": round(float(np.sqrt(((vk - vo) ** 2).mean())), 2)}
            print(n, k, out[n][k], flush=True)
    json.dump(out, open(outpath, "w"), indent=1)


if __name__ == "__main__":
    main()
