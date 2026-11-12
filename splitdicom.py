#!/usr/bin/env python3
"""splitdicom: two-layer, bit-exact DICOM storage prototype.

  encode  in.dcm  base.dcm  residual.sdr  [--near N]
      base.dcm      valid DICOM, JPEG-LS near-lossless (1.2.840.10008.1.2.4.81).
                    Every pixel is within +/-N stored units of the original.
      residual.sdr  exact per-pixel difference + the tags needed to restore the original.

  decode  base.dcm  residual.sdr  out.dcm
      Rebuilds an uncompressed DICOM whose pixels are bit-identical to the original
      (checked against a SHA-256 recorded at encode time).

  verify  original.dcm  restored.dcm
      Compares pixels and every non-pixel data element.

Scope: single-frame, single-sample (grayscale) images, 8 to 16 bits. Research prototype.
Not validated for clinical use.
"""
import argparse, hashlib, json, lzma, struct, sys
import numpy as np
import pydicom
import imagecodecs as ic
from pydicom.encaps import encapsulate, generate_frames
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

JPEGLS_NEAR = "1.2.840.10008.1.2.4.81"
MAGIC = b"SDR1"
# Tags the base layer changes; originals are kept in the sidecar and restored on decode.
CHANGED = ["SOPInstanceUID", "PixelRepresentation", "BitsStored", "HighBit", "RescaleIntercept",
           "RescaleSlope", "LossyImageCompression", "LossyImageCompressionMethod",
           "LossyImageCompressionRatio", "ImageType", "PixelPaddingValue"]


def _sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def encode(src, base_path, res_path, near):
    ds = pydicom.dcmread(src)
    a = ds.pixel_array
    if a.ndim != 2 or int(ds.get("SamplesPerPixel", 1)) != 1:
        sys.exit("prototype handles single-frame grayscale only")
    offset = int(a.min()) if a.min() < 0 else 0          # signed data is shifted to unsigned
    u = (a.astype(np.int64) - offset).astype(np.uint16)
    codestream = ic.jpegls_encode(u, level=near)
    b = ic.jpegls_decode(codestream).astype(np.int64)
    r = u.astype(np.int64) - b
    assert int(np.abs(r).max()) <= near, "error bound violated"

    meta = {"near": near, "offset": offset, "dtype": str(a.dtype), "shape": list(a.shape),
            "sha256": _sha(a), "transfer_syntax": str(ds.file_meta.TransferSyntaxUID),
            "orig": {}, "absent": []}
    for kw in CHANGED:
        if kw in ds:
            meta["orig"][kw] = ds.data_element(kw).to_json_dict(None, None)
        else:
            meta["absent"].append(kw)
    payload = lzma.compress((r + near).astype(np.uint8).tobytes(), preset=6)
    head = json.dumps(meta).encode()
    with open(res_path, "wb") as f:
        f.write(MAGIC + struct.pack("<I", len(head)) + head + payload)

    # Base layer: a standards-conformant lossy derived image with its own identity.
    slope = float(ds.get("RescaleSlope", 1) or 1)
    intercept = float(ds.get("RescaleIntercept", 0) or 0)
    ds.SOPInstanceUID = generate_uid()
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.file_meta.TransferSyntaxUID = JPEGLS_NEAR
    ds.PixelRepresentation = 0
    ds.BitsStored, ds.HighBit = 16, 15
    if offset:
        ds.RescaleSlope, ds.RescaleIntercept = slope, intercept + slope * offset
        if "PixelPaddingValue" in ds:
            del ds.PixelPaddingValue
    ds.LossyImageCompression = "01"
    ds.LossyImageCompressionMethod = "ISO_14495_1"
    ds.LossyImageCompressionRatio = round(a.size * 2 / len(codestream), 2)
    it = list(ds.get("ImageType", ["ORIGINAL", "PRIMARY"]))
    ds.ImageType = ["DERIVED"] + it[1:]
    ds.PixelData = encapsulate([codestream])
    ds["PixelData"].VR = "OB"
    ds.save_as(base_path, enforce_file_format=True)
    return a.size * 2, len(codestream), len(payload) + len(head) + 8


def decode(base_path, res_path, out_path):
    blob = open(res_path, "rb").read()
    if blob[:4] != MAGIC:
        sys.exit("not a splitdicom residual file")
    n = struct.unpack("<I", blob[4:8])[0]
    meta = json.loads(blob[8:8 + n])
    near, shape = meta["near"], tuple(meta["shape"])
    ds = pydicom.dcmread(base_path)
    frame = next(generate_frames(ds.PixelData, number_of_frames=1))
    b = ic.jpegls_decode(frame).astype(np.int64)
    r = np.frombuffer(lzma.decompress(blob[8 + n:]), np.uint8).reshape(shape).astype(np.int64) - near
    a = (b + r + meta["offset"]).astype(np.dtype(meta["dtype"]))
    if _sha(a) != meta["sha256"]:
        sys.exit("FAILED: restored pixels do not match the original hash")
    for kw in meta["absent"]:
        if kw in ds:
            delattr(ds, kw)
    for kw, js in meta["orig"].items():
        tag = pydicom.datadict.tag_for_keyword(kw)
        ds[tag] = pydicom.DataElement.from_json(pydicom.Dataset, tag, js["vr"], js.get("Value", [None]), "Value")
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    ds.PixelData = a.tobytes()
    ds["PixelData"].VR = "OW" if ds.BitsAllocated > 8 else "OB"
    ds.save_as(out_path, enforce_file_format=True)


def verify(orig_path, restored_path):
    o, r = pydicom.dcmread(orig_path), pydicom.dcmread(restored_path)
    if not np.array_equal(o.pixel_array, r.pixel_array):
        return ["pixel data differs"]
    diffs = []
    for tag in sorted(set(o.keys()) | set(r.keys())):
        if tag == 0x7FE00010 or tag.element == 0:   # skip pixels (checked above) and retired group-length elements,
            continue                                 # which are byte counts that change with the encoding
        if tag not in o or tag not in r or o[tag].value != r[tag].value:
            diffs.append(str(tag))
    return diffs


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest="cmd", required=True)
    e = s.add_parser("encode"); e.add_argument("src"); e.add_argument("base"); e.add_argument("residual")
    e.add_argument("--near", type=int, default=2)
    d = s.add_parser("decode"); d.add_argument("base"); d.add_argument("residual"); d.add_argument("out")
    v = s.add_parser("verify"); v.add_argument("original"); v.add_argument("restored")
    args = p.parse_args()
    if args.cmd == "encode":
        raw, base, res = encode(args.src, args.base, args.residual, args.near)
        print(f"raw pixels {raw:,} B | base {base:,} B ({raw/base:.2f}x) | residual {res:,} B ({raw/res:.2f}x)")
    elif args.cmd == "decode":
        decode(args.base, args.residual, args.out)
        print("restored; pixel hash matches original")
    else:
        diffs = verify(args.original, args.restored)
        print("IDENTICAL: pixels and all non-pixel elements match" if not diffs else f"DIFFERENCES: {diffs}")
        sys.exit(1 if diffs else 0)
