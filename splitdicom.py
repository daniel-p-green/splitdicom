#!/usr/bin/env python3
"""splitdicom: two-layer DICOM storage prototype.

  encode  in.dcm  base.dcm  residual.sdr  [--near N] [--overwrite]
      base.dcm      valid DICOM, JPEG-LS near-lossless (1.2.840.10008.1.2.4.81).
                    Every pixel is within +/-N stored units of the original
                    (+/- N x |RescaleSlope| in rescaled units such as HU).
      residual.sdr  exact per-pixel difference, the original values of every element
                    the base changes, and hashes binding it to this base and source.

  decode  base.dcm  residual.sdr  out.dcm  [--overwrite]
      Rebuilds an uncompressed Explicit VR Little Endian DICOM.

  verify  original.dcm  restored.dcm
      Compares pixels and every non-pixel data element.

What "restored" means here: exact decoded pixel samples and verified data-set metadata
for the supported inputs listed in README.md. The original file bytes and transfer
syntax are not recreated. Unsupported inputs are rejected before anything is written.

Research prototype. Not validated for clinical use.
"""
import argparse, base64, hashlib, io, json, lzma, os, struct, sys, tempfile
import numpy as np
import pydicom
import imagecodecs as ic
from pydicom.encaps import encapsulate, generate_frames
from pydicom.uid import ExplicitVRLittleEndian, generate_uid

JPEGLS_NEAR = "1.2.840.10008.1.2.4.81"
MAGIC = b"SDR2"
FORMAT_VERSION = 2
MAX_NEAR = 127                  # residual is stored as one unsigned byte per pixel: 0..2*near
MAX_HEADER = 4 * 1024 * 1024    # parser limits for untrusted sidecars
MAX_PIXELS = 1 << 28
PIXEL_DATA = 0x7FE00010
BIG_ENDIAN = {"1.2.840.10008.1.2.2"}


class SplitDicomError(Exception):
    """Raised for unsupported input, unsafe paths, or failed integrity checks."""


# ---------------------------------------------------------------- hashing

def _sha(b):
    return hashlib.sha256(b).hexdigest()


def pixel_hash(a):
    return _sha(np.ascontiguousarray(a).tobytes())


def _canon_el(el):
    if el.VR == "SQ":
        val = [_canon(item) for item in el.value]
    elif isinstance(el.value, (bytes, bytearray)):
        val = "bytes:" + _sha(bytes(el.value))
    elif isinstance(el.value, (list, tuple, pydicom.multival.MultiValue)):
        val = [str(v) for v in el.value]
    else:
        val = str(el.value)
    return [int(el.tag), el.VR, val]


def _canon(ds):
    # pixels are hashed separately; group lengths are byte counts that change with the encoding
    return [_canon_el(el) for el in ds if el.tag != PIXEL_DATA and el.tag.element != 0]


def metadata_hash(ds):
    """Digest of every data element except pixel data, including nested sequences."""
    return _sha(json.dumps(_canon(ds), sort_keys=True).encode())


def _pack_elements(ds, tags):
    """Serialise original elements exactly (binary DICOM, base64) so string formatting survives."""
    fm = pydicom.dataset.FileMetaDataset(); fm.TransferSyntaxUID = ExplicitVRLittleEndian
    fm.MediaStorageSOPClassUID = "1.2.840.10008.5.1.4.1.1.7"; fm.MediaStorageSOPInstanceUID = "1.2.3"
    c = pydicom.dataset.FileDataset(None, {}, file_meta=fm, preamble=b"\0" * 128)
    if "SpecificCharacterSet" in ds: c.SpecificCharacterSet = ds.SpecificCharacterSet
    for t in tags: c[t] = ds[t]
    buf = io.BytesIO(); c.save_as(buf, enforce_file_format=True)
    return base64.b64encode(buf.getvalue()).decode()


def _unpack_elements(b64):
    try: return pydicom.dcmread(io.BytesIO(base64.b64decode(b64)))
    except Exception: raise SplitDicomError("residual header holds unreadable original elements")


# ---------------------------------------------------------------- validation

def check_supported(ds, arr=None):
    """Return a list of reasons this data set is outside the supported subset (empty if supported)."""
    why = []
    ts = str(getattr(ds.file_meta, "TransferSyntaxUID", ""))
    if ts in BIG_ENDIAN: why.append("big-endian transfer syntax")
    if PIXEL_DATA not in ds: why.append("no integer Pixel Data element")
    if int(ds.get("SamplesPerPixel", 1)) != 1: why.append("more than one sample per pixel")
    if str(ds.get("PhotometricInterpretation", "")) not in ("MONOCHROME1", "MONOCHROME2"): why.append("photometric interpretation is not MONOCHROME1/2")
    if int(ds.get("NumberOfFrames", 1) or 1) != 1: why.append("multi-frame image")
    ba, bs, hb, pr = (ds.get(k) for k in ("BitsAllocated", "BitsStored", "HighBit", "PixelRepresentation"))
    if ba not in (8, 16): why.append(f"BitsAllocated={ba} (only 8 and 16 are supported)")
    if bs is None or ba is None or not (1 <= bs <= ba): why.append(f"BitsStored={bs} inconsistent with BitsAllocated={ba}")
    elif hb != bs - 1: why.append(f"HighBit={hb} is not BitsStored-1")
    if pr not in (0, 1): why.append(f"PixelRepresentation={pr}")
    for g in range(0x6000, 0x6020, 2):
        if any(t.group == g for t in ds.keys()) and (g, 0x3000) not in ds:
            why.append("overlay stored in unused pixel bits (retired encoding)")
            break
    if "ModalityLUTSequence" in ds and pr == 1: why.append("signed pixels with a Modality LUT Sequence")
    if arr is not None and not why:
        if arr.ndim != 2: why.append("pixel array is not two-dimensional")
        elif pr == 0 and int(arr.max()) >= (1 << bs): why.append("pixel values use bits above BitsStored")
        elif pr == 1 and (int(arr.min()) < -(1 << (bs - 1)) or int(arr.max()) >= (1 << (bs - 1))): why.append("pixel values outside the signed BitsStored range")
        elif bs < ba and not ds.file_meta.TransferSyntaxUID.is_compressed:
            # Decoders mask the bits above BitsStored, so anything stored there would be lost silently.
            n = arr.size * (ba // 8)
            raw = np.frombuffer(bytes(ds.PixelData[:n]), dtype="<u1" if ba == 8 else "<u2").astype(np.int64)
            high = raw >> bs
            if pr == 0: bad = high.any()
            else: bad = not ((high == 0) | (high == ((raw >> (bs - 1)) & 1) * ((1 << (ba - bs)) - 1))).all()
            if bad: why.append("data in the unused bits above BitsStored")
    return why


def _check_paths(inputs, outputs, overwrite):
    real = lambda p: os.path.realpath(os.path.abspath(p))
    ins, outs = [real(p) for p in inputs], [real(p) for p in outputs]
    if len(set(outs)) != len(outs): raise SplitDicomError("output paths must be different files")
    for o, shown in zip(outs, outputs):
        if o in ins or any(os.path.exists(o) and os.path.exists(i) and os.path.samefile(o, i) for i in ins):
            raise SplitDicomError(f"refusing to write {shown}: it is one of the input files")
        if os.path.exists(o) and not overwrite:
            raise SplitDicomError(f"refusing to overwrite existing file {shown} (pass --overwrite)")
        if not os.path.isdir(os.path.dirname(o)):
            raise SplitDicomError(f"output directory does not exist for {shown}")


def _stage(path):
    fd, tmp = tempfile.mkstemp(prefix=".splitdicom-", suffix=".tmp", dir=os.path.dirname(os.path.realpath(os.path.abspath(path))))
    os.close(fd)
    return tmp


def _jpegls_precision(codestream):
    i = codestream.find(b"\xff\xf7")
    return codestream[i + 4] if i >= 0 else None


# ---------------------------------------------------------------- core

def _restore(base_bytes, sidecar_bytes):
    """Rebuild the original data set in memory from a base file and its sidecar. Raises on any mismatch."""
    if sidecar_bytes[:4] != MAGIC:
        raise SplitDicomError("not a splitdicom v2 residual file")
    if len(sidecar_bytes) < 8: raise SplitDicomError("truncated residual file")
    n = struct.unpack("<I", sidecar_bytes[4:8])[0]
    if n > MAX_HEADER or 8 + n > len(sidecar_bytes): raise SplitDicomError("residual header size is out of bounds")
    try: meta = json.loads(sidecar_bytes[8:8 + n])
    except ValueError: raise SplitDicomError("residual header is not valid JSON")
    if meta.get("version") != FORMAT_VERSION: raise SplitDicomError(f"unsupported residual format version {meta.get('version')}")
    near, shape = int(meta["near"]), tuple(int(x) for x in meta["shape"])
    if not 0 <= near <= MAX_NEAR or len(shape) != 2 or min(shape) < 1 or shape[0] * shape[1] > MAX_PIXELS:
        raise SplitDicomError("residual header values are out of bounds")
    if meta["dtype"] not in ("uint8", "int8", "uint16", "int16"): raise SplitDicomError("unsupported pixel type in residual header")
    if _sha(base_bytes) != meta["base_sha256"]:
        raise SplitDicomError("base file does not match this residual (wrong pair, or the base was modified)")
    ds = pydicom.dcmread(io.BytesIO(base_bytes))
    if (int(ds.Rows), int(ds.Columns)) != shape: raise SplitDicomError("base dimensions do not match the residual")
    frame = next(generate_frames(ds.PixelData, number_of_frames=1))
    b = ic.jpegls_decode(frame).astype(np.int64)
    if b.shape != shape: raise SplitDicomError("decoded base dimensions do not match the residual")
    npx = shape[0] * shape[1]
    try: raw = lzma.LZMADecompressor().decompress(sidecar_bytes[8 + n:], max_length=npx + 1)
    except lzma.LZMAError: raise SplitDicomError("residual payload is corrupt")
    if len(raw) != npx: raise SplitDicomError("residual payload has the wrong size")
    r = np.frombuffer(raw, np.uint8).reshape(shape).astype(np.int64) - near
    a = (b + r + int(meta["offset"])).astype(np.dtype(meta["dtype"]))
    if pixel_hash(a) != meta["pixel_sha256"]:
        raise SplitDicomError("restored pixels do not match the original hash")
    for tag in meta["absent"]:
        if int(tag) in ds: del ds[int(tag)]
    packed = _unpack_elements(meta["orig_elements"])
    for tag in meta["orig_tags"]:
        if int(tag) not in packed: raise SplitDicomError("residual header is missing an original element")
        ds[int(tag)] = packed[int(tag)]
    ds.PixelData = a.tobytes()
    ds["PixelData"].VR = "OW" if ds.BitsAllocated > 8 else "OB"
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.file_meta.TransferSyntaxUID = ExplicitVRLittleEndian
    if metadata_hash(ds) != meta["metadata_sha256"]:
        raise SplitDicomError("restored metadata does not match the original (identity, geometry, or other elements differ)")
    return ds, a, meta


def encode(src, base_path, res_path, near=2, overwrite=False):
    """Write base + residual for one image. Returns a dict of byte counts. Raises SplitDicomError."""
    if not isinstance(near, int) or not 0 <= near <= MAX_NEAR:
        raise SplitDicomError(f"--near must be an integer from 0 to {MAX_NEAR}")
    _check_paths([src], [base_path, res_path], overwrite)
    ds = pydicom.dcmread(src)
    why = check_supported(ds)
    if why: raise SplitDicomError("unsupported input: " + "; ".join(why))
    try: a = ds.pixel_array
    except Exception as e: raise SplitDicomError(f"cannot decode source pixels: {e}")
    why = check_supported(ds, a)
    if why: raise SplitDicomError("unsupported input: " + "; ".join(why))
    orig = pydicom.dcmread(src, stop_before_pixels=True)
    orig_ts = str(ds.file_meta.TransferSyntaxUID)

    bits = int(ds.BitsAllocated)
    offset = int(a.min()) if a.min() < 0 else 0          # signed data is shifted to unsigned
    u = (a.astype(np.int64) - offset).astype(np.uint8 if bits == 8 else np.uint16)
    codestream = ic.jpegls_encode(u, level=near)
    b = ic.jpegls_decode(codestream).astype(np.int64)
    r = u.astype(np.int64) - b
    if int(np.abs(r).max()) > near: raise SplitDicomError("codec violated the error bound")
    if _jpegls_precision(codestream) != bits: raise SplitDicomError("unexpected JPEG-LS sample precision")

    # Base layer: a lossy derived image with its own identity. Header agrees with the code stream.
    slope = float(ds.get("RescaleSlope", 1) or 1)
    intercept = float(ds.get("RescaleIntercept", 0) or 0)
    ds.SOPInstanceUID = generate_uid()
    ds.file_meta.MediaStorageSOPInstanceUID = ds.SOPInstanceUID
    ds.file_meta.TransferSyntaxUID = JPEGLS_NEAR
    ds.PixelRepresentation = 0
    ds.BitsStored, ds.HighBit = bits, bits - 1
    if offset:
        if "RescaleSlope" not in ds: ds.RescaleSlope = 1
        ds.RescaleIntercept = intercept + slope * offset
    for kw in ("PixelPaddingValue", "PixelPaddingRangeLimit", "SmallestImagePixelValue", "LargestImagePixelValue"):
        if kw in ds: delattr(ds, kw)                       # stored-value statements that no longer hold for the base
    ds.LossyImageCompression = "01"
    ds.LossyImageCompressionMethod = "ISO_14495_1"
    ds.LossyImageCompressionRatio = round(a.nbytes / len(codestream), 2)
    ds.ImageType = ["DERIVED"] + list(ds.get("ImageType", ["ORIGINAL", "PRIMARY"]))[1:]
    ds.PixelData = encapsulate([codestream])
    ds["PixelData"].VR = "OB"

    base_tmp, res_tmp = _stage(base_path), _stage(res_path)
    try:
        ds.save_as(base_tmp, enforce_file_format=True)
        base_bytes = open(base_tmp, "rb").read()
        written = pydicom.dcmread(io.BytesIO(base_bytes), stop_before_pixels=True)
        changed = []; absent = []
        for tag in sorted(set(orig.keys()) | set(written.keys())):
            if tag == PIXEL_DATA or tag.element == 0: continue
            if tag not in orig: absent.append(int(tag))
            elif tag not in written or _canon_el(orig[tag]) != _canon_el(written[tag]): changed.append(int(tag))
        meta = {"version": FORMAT_VERSION, "near": near, "offset": offset, "dtype": str(a.dtype), "shape": list(a.shape),
                "pixel_sha256": pixel_hash(a), "base_sha256": _sha(base_bytes), "metadata_sha256": metadata_hash(orig),
                "source_sop_instance_uid": str(orig.get("SOPInstanceUID", "")), "source_transfer_syntax": orig_ts,
                "orig_tags": changed, "orig_elements": _pack_elements(orig, changed), "absent": absent}
        head = json.dumps(meta).encode()
        payload = lzma.compress((r + near).astype(np.uint8).tobytes(), preset=6)
        with open(res_tmp, "wb") as f:
            f.write(MAGIC + struct.pack("<I", len(head)) + head + payload)
        # Prove the staged pair restores before publishing either file.
        res_bytes = open(res_tmp, "rb").read()
        _restore(base_bytes, res_bytes)
        os.replace(res_tmp, res_path)
        try: os.replace(base_tmp, base_path)
        except OSError:
            os.unlink(res_path); raise
    finally:
        for t in (base_tmp, res_tmp):
            if os.path.exists(t): os.unlink(t)
    return {"pixel_bytes": int(a.nbytes), "source_file_bytes": os.path.getsize(src), "base_file_bytes": len(base_bytes),
            "residual_file_bytes": len(res_bytes), "base_codestream_bytes": len(codestream), "bound_rescaled_units": near * abs(slope)}


def decode(base_path, res_path, out_path, overwrite=False):
    """Restore an image. The output is reopened and re-verified before it is published."""
    _check_paths([base_path, res_path], [out_path], overwrite)
    ds, a, meta = _restore(open(base_path, "rb").read(), open(res_path, "rb").read())
    tmp = _stage(out_path)
    try:
        ds.save_as(tmp, enforce_file_format=True)
        back = pydicom.dcmread(tmp)
        if pixel_hash(back.pixel_array) != meta["pixel_sha256"] or metadata_hash(back) != meta["metadata_sha256"]:
            raise SplitDicomError("the written file does not read back as the original")
        os.replace(tmp, out_path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)


def verify(orig_path, restored_path):
    """Return a list of differences between two files' pixels and non-pixel elements (empty if none)."""
    o, r = pydicom.dcmread(orig_path), pydicom.dcmread(restored_path)
    if not np.array_equal(o.pixel_array, r.pixel_array):
        return ["pixel data differs"]
    diffs = []
    for tag in sorted(set(o.keys()) | set(r.keys())):
        if tag == PIXEL_DATA or tag.element == 0:   # skip pixels (checked above) and retired group-length elements,
            continue                                 # which are byte counts that change with the encoding
        if tag not in o or tag not in r or o[tag].value != r[tag].value:
            diffs.append(str(tag))
    return diffs


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_subparsers(dest="cmd", required=True)
    e = s.add_parser("encode"); e.add_argument("src"); e.add_argument("base"); e.add_argument("residual")
    e.add_argument("--near", type=int, default=2, help=f"error bound in stored units, 0 to {MAX_NEAR}")
    e.add_argument("--overwrite", action="store_true")
    d = s.add_parser("decode"); d.add_argument("base"); d.add_argument("residual"); d.add_argument("out")
    d.add_argument("--overwrite", action="store_true")
    v = s.add_parser("verify"); v.add_argument("original"); v.add_argument("restored")
    args = p.parse_args(argv)
    try:
        if args.cmd == "encode":
            m = encode(args.src, args.base, args.residual, args.near, args.overwrite)
            px = m["pixel_bytes"]
            print(f"pixels {px:,} B | base file {m['base_file_bytes']:,} B ({px / m['base_file_bytes']:.2f}x) | "
                  f"residual file {m['residual_file_bytes']:,} B ({px / m['residual_file_bytes']:.2f}x) | "
                  f"bound +/-{args.near} stored units = +/-{m['bound_rescaled_units']:g} rescaled units | pair verified")
        elif args.cmd == "decode":
            decode(args.base, args.residual, args.out, args.overwrite)
            print("restored; pixels and metadata match the original")
        else:
            diffs = verify(args.original, args.restored)
            print("IDENTICAL: pixels and all non-pixel elements match" if not diffs else f"DIFFERENCES: {diffs}")
            return 1 if diffs else 0
    except SplitDicomError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
