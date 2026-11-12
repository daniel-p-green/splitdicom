"""End-to-end check of splitdicom.py on real sample series: encode, decode, verify every image.

usage: e2e.py <data_root> [--near N] [--out results.json]
Exits non-zero if any selected image fails or if nothing was selected. Unsupported inputs
that the codec rejects up front are counted separately as skipped, with reasons.
"""
import sys, os, json, argparse, tempfile, time, warnings, collections
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, ".."))
import numpy as np, pydicom, imagecodecs as ic, splitdicom as sd
from bench import SERIES


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("root"); ap.add_argument("--near", type=int, default=2); ap.add_argument("--out")
    args = ap.parse_args()
    if not os.path.isdir(args.root): print(f"error: data root {args.root} does not exist", file=sys.stderr); return 2
    tmp = tempfile.mkdtemp(); c = collections.Counter(); skip_reasons = collections.Counter(); fail = []; per = {}
    t0 = time.time()
    for label, (d, mod, desc) in SERIES.items():
        s = collections.Counter()
        for dp, _, fs in os.walk(os.path.join(args.root, d)):
            for f in sorted(fs):
                p = os.path.join(dp, f)
                try: ds = pydicom.dcmread(p, stop_before_pixels=True)
                except Exception: continue
                c["discovered"] += 1
                if ds.get("Modality") != mod or "Rows" not in ds: continue
                if desc and not str(ds.get("SeriesDescription", "")).startswith(desc): continue
                c["selected"] += 1
                b, r, o = (os.path.join(tmp, x) for x in ("b.dcm", "r.sdr", "o.dcm"))
                try:
                    m = sd.encode(p, b, r, args.near, overwrite=True)
                    sd.decode(b, r, o, overwrite=True)
                    diffs = sd.verify(p, o)
                    if diffs: raise AssertionError(f"verify differences: {diffs[:5]}")
                    # independent check: the base opens through pydicom's own decoder and honours the bound
                    bo = pydicom.dcmread(b); ba = bo.pixel_array.astype(np.int64); oa = pydicom.dcmread(p).pixel_array.astype(np.int64)
                    off = int(oa.min()) if oa.min() < 0 else 0
                    if np.abs((ba + off) - oa).max() > args.near: raise AssertionError("base exceeds the error bound in an independent decoder")
                    c["passed"] += 1
                    for k in ("pixel_bytes", "source_file_bytes", "base_file_bytes", "residual_file_bytes", "base_codestream_bytes"): s[k] += m[k]
                    # Like-for-like lossless baseline: the same DICOM envelope around a lossless JPEG-LS code stream.
                    u = (oa - off).astype(np.uint8 if bo.BitsAllocated == 8 else np.uint16)
                    s["lossless_jpegls_file_bytes_est"] += m["base_file_bytes"] - m["base_codestream_bytes"] + len(ic.jpegls_encode(u, level=0))
                    s["images"] += 1
                except sd.SplitDicomError as e:
                    if str(e).startswith("unsupported input"): c["skipped_unsupported"] += 1; skip_reasons[str(e)] += 1
                    else: c["failed"] += 1; fail.append([label, f, str(e)[:200]])
                except Exception as e:
                    c["failed"] += 1; fail.append([label, f, f"{type(e).__name__}: {str(e)[:200]}"])
        if s["images"]: per[label] = dict(s)
    out = {"near": args.near, "counts": dict(c), "skip_reasons": dict(skip_reasons), "failures": fail[:50], "per_series": per, "seconds": round(time.time() - t0, 1),
           "versions": {"pydicom": pydicom.__version__, "numpy": np.__version__}}
    print(json.dumps({k: v for k, v in out.items() if k != "per_series"}, indent=1))
    if args.out: json.dump(out, open(args.out, "w"), indent=1)
    if c["selected"] == 0: print("error: no images were selected", file=sys.stderr); return 2
    return 1 if c["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
