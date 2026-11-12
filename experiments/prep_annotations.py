"""Convert picked LIDC nodules (voxel centroids) to a LUNA16-style CSV in patient coordinates.
usage: prep_annotations.py <picked.json> <dicom_zip_dir> <out.csv>"""
import sys, json, csv, io, os, zipfile, pydicom
pick = json.load(open(sys.argv[1])); geo = {}
w = csv.writer(open(sys.argv[3], "w", newline=""))
w.writerow(["seriesuid", "coordX", "coordY", "coordZ", "diameter_mm", "texture_median", "readers", "subtlety_mean"])
for n in pick:
    u = n["uid"]; p = os.path.join(sys.argv[2], u + ".zip")
    if not os.path.exists(p): continue
    if u not in geo:
        z = zipfile.ZipFile(p)
        ds = pydicom.dcmread(io.BytesIO(z.read([m for m in z.namelist() if m.lower().endswith(".dcm")][0])), stop_before_pixels=True)
        geo[u] = (float(ds.ImagePositionPatient[0]), float(ds.ImagePositionPatient[1]), float(ds.PixelSpacing[0]), float(ds.PixelSpacing[1]))
    x0, y0, pr, pc = geo[u]
    # first contour coordinate is the column (x), second the row (y): validated against LUNA16 positions
    w.writerow([u, x0 + n["c0"] * pc, y0 + n["c1"] * pr, n["z"], n["diameter_mm"], n["texture_median"], n["readers"], round(n["subtlety_mean"], 2)])
