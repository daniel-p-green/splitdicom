"""Summarise nodule_test.py output. usage: nodule_analyze.py <nodules.jsonl> <summary.json>"""
import sys, json, collections
import numpy as np

T = 0.5                       # detection score threshold
VARS = ["rand1", "4", "8", "16"]


BOX_OFFSET = 0.36             # see ctr(): calibrated on synthetic spheres


def tol(d): return max(d / 2, 3.0)   # match radius in mm: LUNA16 uses the nodule radius; floor of 3 mm for tiny nodules


def ctr(x):
    """Corrected detection centre (mm, RAS, relative to the annotated centre).

    With MONAI 1.6.1 this bundle returns boxes whose midpoint sits about 0.36 x box size
    past the true centre on every axis (measured on synthetic spheres, see notes/04).
    The offset is identical in every variant, so paired comparisons are unaffected;
    it only matters for matching detections to annotations.
    """
    return np.array(x[:3]) - BOX_OFFSET * x[4]


def matched(dets, d):
    m = [x for x in dets if np.linalg.norm(ctr(x)) <= tol(d)]
    return max(m, key=lambda x: x[3]) if m else None


def false_pos(dets, d, others):
    out = []
    for x in dets:
        if x[3] < T or np.linalg.norm(ctr(x)) <= tol(d): continue
        c = ctr(x)
        # others are LPS offsets; detections are RAS offsets
        if any(np.linalg.norm([c[0] + o[0], c[1] + o[1], c[2] - o[2]]) <= tol(o[3]) for o in others): continue
        out.append(x)
    return out


def main():
    scans = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
    used = [s for s in scans if "nodules" in s]; skipped = [s for s in scans if "skipped" in s]
    rows = []
    for s in used:
        for n in s["nodules"]:
            d = n["diameter_mm"]; v = n["variants"]; o = matched(v["orig"]["dets"], d)
            r = {"d": d, "orig_score": o[3] if o else 0.0, "orig_vol": v["orig"]["vol_mm3"], "fp_orig": len(false_pos(v["orig"]["dets"], d, n["others_lps_mm"]))}
            for k in VARS:
                m = matched(v[k]["dets"], d)
                r[k] = {"score": m[3] if m else 0.0,
                        "shift_mm": float(np.linalg.norm(np.array(m[:3]) - np.array(o[:3]))) if m and o else None,
                        "box_pct": 100 * (m[4] / o[4] - 1) if m and o else None,
                        "vol_pct": 100 * (v[k]["vol_mm3"] / v["orig"]["vol_mm3"] - 1) if v["orig"]["vol_mm3"] > 0 else None,
                        "hu_shift": (v[k]["mean_hu_in_orig_mask"] - v["orig"]["mean_hu_in_orig_mask"]) if v["orig"]["mean_hu_in_orig_mask"] is not None and v[k]["mean_hu_in_orig_mask"] is not None else None,
                        "fp": len(false_pos(v[k]["dets"], d, n["others_lps_mm"]))}
            rows.append(r)
    N = len(rows); det0 = sum(r["orig_score"] >= T for r in rows)
    summ = {"scans_used": len(used), "scans_skipped": len(skipped), "nodules": N, "nodules_under_6mm": sum(r["d"] < 6 for r in rows),
            "diameter_min_median_max": [round(float(x), 1) for x in (min(r["d"] for r in rows), np.median([r["d"] for r in rows]), max(r["d"] for r in rows))],
            "score_threshold": T, "detected_on_original": det0, "false_positives_on_original": sum(r["fp_orig"] for r in rows),
            "slice_spacing_mm": dict(collections.Counter(str(s["meta"]["slice_spacing"]) for s in used)),
            "manufacturers": dict(collections.Counter(s["meta"]["manufacturer"] for s in used)), "variants": {}}
    for k in VARS:
        sc = np.array([r[k]["score"] - r["orig_score"] for r in rows])
        vol = np.array([r[k]["vol_pct"] for r in rows if r[k]["vol_pct"] is not None])
        small = np.array([r[k]["vol_pct"] for r in rows if r[k]["vol_pct"] is not None and r["d"] < 6])
        box = np.array([r[k]["box_pct"] for r in rows if r[k]["box_pct"] is not None])
        sh = np.array([r[k]["shift_mm"] for r in rows if r[k]["shift_mm"] is not None])
        hu = np.array([r[k]["hu_shift"] for r in rows if r[k]["hu_shift"] is not None])
        ratios = [s["base_ratio"][k] for s in used if k in s.get("base_ratio", {})]
        flips = {str(t): {"lost": int(sum(r["orig_score"] >= t and r[k]["score"] < t for r in rows)), "gained": int(sum(r["orig_score"] < t and r[k]["score"] >= t for r in rows))} for t in (0.3, 0.5, 0.7)}
        summ["variants"][k] = {
            "base_ratio_min_median_max": [round(float(x), 2) for x in (min(ratios), np.median(ratios), max(ratios))] if ratios else None,
            "detected": int(sum(r[k]["score"] >= T for r in rows)), "flips_by_threshold": flips,
            "score_change_abs_median": round(float(np.median(np.abs(sc))), 5), "score_change_abs_p95": round(float(np.percentile(np.abs(sc), 95)), 4), "score_change_abs_max": round(float(np.abs(sc).max()), 4),
            "det_shift_mm_max": round(float(sh.max()), 3), "det_box_size_change_pct_abs_max": round(float(np.abs(box).max()), 2),
            "false_positives": int(sum(r[k]["fp"] for r in rows)),
            "vol_change_pct_abs_median": round(float(np.median(np.abs(vol))), 3), "vol_change_pct_abs_p95": round(float(np.percentile(np.abs(vol), 95)), 2), "vol_change_pct_abs_max": round(float(np.abs(vol).max()), 2),
            "vol_change_over_5pct": int((np.abs(vol) > 5).sum()), "vol_change_over_10pct": int((np.abs(vol) > 10).sum()), "vol_measured": int(len(vol)),
            "small_nodule_vol_change_pct_abs_median": round(float(np.median(np.abs(small))), 3), "small_nodule_vol_change_pct_abs_max": round(float(np.abs(small).max()), 2),
            "mean_hu_shift_abs_median": round(float(np.median(np.abs(hu))), 3), "mean_hu_shift_abs_max": round(float(np.abs(hu).max()), 2)}
        print(k, json.dumps(summ["variants"][k]))
    print({k: v for k, v in summ.items() if k != "variants"})
    json.dump(summ, open(sys.argv[2], "w"), indent=1)


if __name__ == "__main__":
    main()
