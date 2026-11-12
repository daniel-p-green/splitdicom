"""Summarise nodule_test.py output, overall and by radiologist texture group when present.
usage: nodule_analyze.py <nodules.jsonl> <summary.json>"""
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


def pct(a, b): return 100 * (a / b - 1) if b and b > 0 and a is not None else None


def group_of(n):
    t = n.get("info", {}).get("texture_median")
    if t is None: return "all"
    return "non-solid" if float(t) <= 2 else "part-solid"


def stats(vals):
    v = np.abs(np.array([x for x in vals if x is not None], float))
    if len(v) == 0: return None
    return {"n": int(len(v)), "abs_median": round(float(np.median(v)), 3), "abs_p95": round(float(np.percentile(v, 95)), 2), "abs_max": round(float(v.max()), 2),
            "over_5": int((v > 5).sum()), "over_10": int((v > 10).sum())}


def summarise(rows, used):
    N = len(rows)
    out = {"nodules": N, "nodules_under_6mm": sum(r["d"] < 6 for r in rows),
           "diameter_min_median_max": [round(float(x), 1) for x in (min(r["d"] for r in rows), np.median([r["d"] for r in rows]), max(r["d"] for r in rows))],
           "detected_on_original": sum(r["orig_score"] >= T for r in rows), "orig_score_between_0.2_and_0.8": sum(0.2 <= r["orig_score"] <= 0.8 for r in rows),
           "false_positives_on_original": sum(r["fp_orig"] for r in rows), "adaptive_measurable": sum(r["orig_avol"] is not None and r["orig_avol"] > 0 for r in rows), "variants": {}}
    for k in VARS:
        sc = np.array([r[k]["score"] - r["orig_score"] for r in rows])
        sh = [r[k]["shift_mm"] for r in rows if r[k]["shift_mm"] is not None]
        flips = {str(t): {"lost": int(sum(r["orig_score"] >= t and r[k]["score"] < t for r in rows)), "gained": int(sum(r["orig_score"] < t and r[k]["score"] >= t for r in rows))} for t in (0.3, 0.5, 0.7)}
        out["variants"][k] = {"detected": int(sum(r[k]["score"] >= T for r in rows)), "flips_by_threshold": flips,
            "score_change_abs_median": round(float(np.median(np.abs(sc))), 5), "score_change_abs_p95": round(float(np.percentile(np.abs(sc), 95)), 4), "score_change_abs_max": round(float(np.abs(sc).max()), 4),
            "det_shift_mm_max": round(float(max(sh)), 3) if sh else None, "det_box_size_change_pct": stats([r[k]["box_pct"] for r in rows]),
            "false_positives": int(sum(r[k]["fp"] for r in rows)),
            "fixed_threshold_vol_change_pct": stats([r[k]["vol_pct"] for r in rows]), "adaptive_threshold_vol_change_pct": stats([r[k]["avol_pct"] for r in rows]),
            "adaptive_vol_change_pct_regions_50mm3_plus": stats([r[k]["avol_pct"] for r in rows if r["orig_avol"] and r["orig_avol"] >= 50]),
            "mean_hu_shift": stats([r[k]["hu_shift"] for r in rows]), "adaptive_mean_hu_shift": stats([r[k]["ahu_shift"] for r in rows])}
    return out


def main():
    scans = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
    used = [s for s in scans if "nodules" in s]; skipped = [s for s in scans if "skipped" in s]
    groups = collections.defaultdict(list)
    for s in used:
        for n in s["nodules"]:
            d = n["diameter_mm"]; v = n["variants"]; o = matched(v["orig"]["dets"], d)
            diff = lambda a, b: (a - b) if a is not None and b is not None else None
            r = {"d": d, "orig_score": o[3] if o else 0.0, "fp_orig": len(false_pos(v["orig"]["dets"], d, n["others_lps_mm"])), "orig_avol": v["orig"].get("adaptive_vol_mm3")}
            for k in VARS:
                m = matched(v[k]["dets"], d)
                r[k] = {"score": m[3] if m else 0.0,
                        "shift_mm": float(np.linalg.norm(np.array(m[:3]) - np.array(o[:3]))) if m and o else None,
                        "box_pct": 100 * (m[4] / o[4] - 1) if m and o else None,
                        "vol_pct": pct(v[k]["vol_mm3"], v["orig"]["vol_mm3"]), "avol_pct": pct(v[k].get("adaptive_vol_mm3"), v["orig"].get("adaptive_vol_mm3")),
                        "hu_shift": diff(v[k]["mean_hu_in_orig_mask"], v["orig"]["mean_hu_in_orig_mask"]),
                        "ahu_shift": diff(v[k].get("adaptive_mean_hu_in_orig_mask"), v["orig"].get("adaptive_mean_hu_in_orig_mask")),
                        "fp": len(false_pos(v[k]["dets"], d, n["others_lps_mm"]))}
            groups[group_of(n)].append(r); 
            if group_of(n) != "all": groups["all"].append(r)
    ratios = {k: [s["base_ratio"][k] for s in used if k in s.get("base_ratio", {})] for k in VARS}
    summ = {"scans_used": len(used), "scans_skipped": len(skipped), "score_threshold": T,
            "slice_spacing_mm": dict(collections.Counter(str(s["meta"]["slice_spacing"]) for s in used)),
            "manufacturers": dict(collections.Counter(s["meta"]["manufacturer"] for s in used)),
            "base_ratio_min_median_max": {k: [round(float(x), 2) for x in (min(v), np.median(v), max(v))] for k, v in ratios.items() if v},
            "groups": {g: summarise(rows, used) for g, rows in groups.items()}}
    json.dump(summ, open(sys.argv[2], "w"), indent=1)
    print({k: v for k, v in summ.items() if k != "groups"})
    for g, s in summ["groups"].items():
        print("==", g, {k: v for k, v in s.items() if k != "variants"})
        for k, v in s["variants"].items(): print("  ", k, json.dumps(v))


if __name__ == "__main__":
    main()
