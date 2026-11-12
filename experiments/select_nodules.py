"""Select LIDC nodules by radiologist texture rating from the pylidc annotation database.

Clusters the per-reader annotations into nodules, then keeps nodules marked by at least two
readers. Output: one row per nodule with voxel-space centroid (first, second contour coordinate),
z position in mm, equivalent diameter, median texture, mean subtlety.
usage: select_nodules.py <pylidc.sqlite> <out.json>
"""
import sys, json, sqlite3, collections, hashlib
import numpy as np
c = sqlite3.connect(sys.argv[1])
scans = {r[0]: {"uid": r[1], "thick": r[2], "ps": r[3]} for r in c.execute("select id,series_instance_uid,slice_thickness,pixel_spacing from scans")}
zvals = collections.defaultdict(list)
for sid, v in c.execute("select scan_id,val from zvals"): zvals[sid].append(v)
anns = collections.defaultdict(list)
cont = collections.defaultdict(list)
for aid, inc, z, coords in c.execute("select annotation_id,inclusion,image_z_position,coords from contours"):
    pts = np.array([[float(v) for v in ln.split(",")] for ln in coords.strip().split("\n")])
    cont[aid].append((inc, z, pts))
for aid, sid, tex, sub in c.execute("select id,scan_id,texture,subtlety from annotations"):
    inc = [(z, p) for i, z, p in cont[aid] if i]
    if not inc: continue
    s = scans[sid]; zs = sorted(set(zvals[sid])); dz = float(np.median(np.diff(zs))) if len(zs) > 1 else s["thick"]
    area = 0.0; ws = []; cs = []
    for z, p in inc:
        x, y = p[:, 0], p[:, 1]
        a = 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))) * s["ps"] ** 2
        a = max(a, len(p) * s["ps"] ** 2 * 0.5) if len(p) < 4 else a
        area += a; ws.append(max(a, 1e-3)); cs.append([p[:, 0].mean(), p[:, 1].mean(), z])
    cs = np.array(cs); ws = np.array(ws); cen = (cs * ws[:, None]).sum(0) / ws.sum()
    vol = area * dz
    anns[sid].append({"c": cen, "d": float((6 * vol / np.pi) ** (1 / 3)), "tex": tex, "sub": sub})
out = []
for sid, al in anns.items():
    s = scans[sid]; used = [False] * len(al)
    for i, a in enumerate(al):
        if used[i]: continue
        grp = [a]; used[i] = True
        for j in range(i + 1, len(al)):
            b = al[j]
            if used[j]: continue
            dist = np.sqrt(((a["c"][0] - b["c"][0]) * s["ps"]) ** 2 + ((a["c"][1] - b["c"][1]) * s["ps"]) ** 2 + (a["c"][2] - b["c"][2]) ** 2)
            if dist <= max(a["d"], b["d"], 4.0) * 0.75: grp.append(b); used[j] = True
        cen = np.mean([g["c"] for g in grp], 0)
        out.append({"uid": s["uid"], "thick": s["thick"], "ps": s["ps"], "c0": float(cen[0]), "c1": float(cen[1]), "z": float(cen[2]),
                    "diameter_mm": float(np.mean([g["d"] for g in grp])), "readers": len(grp), "texture_median": float(np.median([g["tex"] for g in grp])),
                    "textures": [g["tex"] for g in grp], "subtlety_mean": float(np.mean([g["sub"] for g in grp]))})
json.dump(out, open(sys.argv[2], "w"))
ok = [n for n in out if n["readers"] >= 2 and n["thick"] <= 2.5 and n["diameter_mm"] >= 3.0]
print("nodules", len(out), "| >=2 readers, <=2.5mm slices, >=3mm:", len(ok))
for lab, f in [("non-solid (median texture <=2)", lambda n: n["texture_median"] <= 2), ("part-solid (2<median<4)", lambda n: 2 < n["texture_median"] < 4), ("solid (>=4)", lambda n: n["texture_median"] >= 4)]:
    g = [n for n in ok if f(n)]; print(lab, len(g), "nodules in", len({n["uid"] for n in g}), "scans; diam median", round(float(np.median([n["diameter_mm"] for n in g])), 1))
print(c.execute("select coords from contours limit 1").fetchone()[0][:40].replace("\n", " | "))
