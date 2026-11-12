"""Regression tests for splitdicom.py. Synthetic fixtures only: no downloads, no patient data.

Run:  python -m pytest tests -q
"""
import hashlib, os, struct, sys
import numpy as np, pydicom, pytest
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import (CTImageStorage, ExplicitVRBigEndian, ExplicitVRLittleEndian, ImplicitVRLittleEndian,
                         SecondaryCaptureImageStorage, generate_uid)

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import splitdicom as sd

RNG = np.random.default_rng(7)


def make(path, arr, bits_alloc=16, bits_stored=None, signed=False, pid="P1", sop=CTImageStorage, ts=ExplicitVRLittleEndian, photometric="MONOCHROME2", edit=None, **extra):
    fm = FileMetaDataset(); fm.MediaStorageSOPClassUID = sop; fm.MediaStorageSOPInstanceUID = generate_uid(); fm.TransferSyntaxUID = ts
    be = ts == ExplicitVRBigEndian
    ds = FileDataset(str(path), {}, file_meta=fm, preamble=b"\0" * 128, is_implicit_VR=ts == ImplicitVRLittleEndian, is_little_endian=not be)
    ds.SOPClassUID = sop; ds.SOPInstanceUID = fm.MediaStorageSOPInstanceUID
    ds.PatientID = pid; ds.PatientName = pid; ds.Modality = "CT"
    ds.StudyInstanceUID = generate_uid(); ds.SeriesInstanceUID = generate_uid()
    ds.ImageType = ["ORIGINAL", "PRIMARY", "AXIAL"]
    ds.Rows, ds.Columns = arr.shape[-2], arr.shape[-1]; ds.SamplesPerPixel = 1; ds.PhotometricInterpretation = photometric
    ds.BitsAllocated = bits_alloc; ds.BitsStored = bits_stored or bits_alloc; ds.HighBit = ds.BitsStored - 1
    ds.PixelRepresentation = 1 if signed else 0
    ds.PixelSpacing = [0.7, 0.7]; ds.ImagePositionPatient = [1.5, -2.5, 30.0]; ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    ds.add_new(0x00091001, "LO", "private value")
    for k, v in extra.items(): setattr(ds, k, v)
    ds.PixelData = (arr.byteswap() if be else arr).tobytes()
    if edit: edit(ds)
    ds.save_as(str(path), enforce_file_format=True)
    return str(path)


def sha(path): return hashlib.sha256(open(path, "rb").read()).hexdigest()


def ct12(shape=(48, 64)): return RNG.integers(0, 4096, shape).astype(np.uint16)


def roundtrip(tmp_path, src, near):
    base, res, out = (str(tmp_path / n) for n in ("base.dcm", "res.sdr", "out.dcm"))
    info = sd.encode(src, base, res, near)
    sd.decode(base, res, out)
    assert sd.verify(src, out) == []
    return base, res, out, info


# ---------------------------------------------------------------- round trips

@pytest.mark.parametrize("near", [0, 1, 2, 16, 127])
def test_roundtrip_unsigned_12bit(tmp_path, near):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12, RescaleSlope="1", RescaleIntercept="-1024")
    base, res, out, info = roundtrip(tmp_path, src, near)
    b = pydicom.dcmread(base)
    assert b.file_meta.TransferSyntaxUID == sd.JPEGLS_NEAR and b.SOPInstanceUID != pydicom.dcmread(src).SOPInstanceUID
    assert (b.BitsAllocated, b.BitsStored, b.HighBit, b.PixelRepresentation) == (16, 16, 15, 0)
    restored = pydicom.dcmread(out)
    assert (restored.BitsStored, restored.HighBit) == (12, 11)
    assert info["pixel_bytes"] == 48 * 64 * 2 and info["base_file_bytes"] == os.path.getsize(base) and info["residual_file_bytes"] == os.path.getsize(res)


@pytest.mark.parametrize("near", [0, 2, 8])
def test_roundtrip_signed_with_padding(tmp_path, near):
    a = RNG.integers(-1024, 3000, (40, 40)).astype(np.int16); a[:4] = -2000
    def edit(ds): ds.add_new(0x00280121, "SS", -1900)        # PixelPaddingRangeLimit
    src = make(tmp_path / "in.dcm", a, signed=True, PixelPaddingValue=-2000, RescaleSlope="1", RescaleIntercept="0", edit=edit)
    base, res, out, _ = roundtrip(tmp_path, src, near)
    b = pydicom.dcmread(base)
    assert b.PixelRepresentation == 0 and "PixelPaddingValue" not in b and 0x00280121 not in b
    assert float(b.RescaleIntercept) == -2000.0
    r = pydicom.dcmread(out)
    assert r.PixelRepresentation == 1 and r.PixelPaddingValue == -2000 and r[0x00280121].value == -1900 and float(r.RescaleIntercept) == 0.0


def test_roundtrip_signed_non_negative(tmp_path):
    src = make(tmp_path / "in.dcm", RNG.integers(0, 3000, (32, 32)).astype(np.int16), signed=True)
    roundtrip(tmp_path, src, 2)


@pytest.mark.parametrize("near", [0, 2, 20])
def test_roundtrip_8bit(tmp_path, near):
    a = RNG.integers(0, 256, (64, 64)).astype(np.uint8)
    src = make(tmp_path / "in.dcm", a, bits_alloc=8, sop=SecondaryCaptureImageStorage)
    base, res, out, info = roundtrip(tmp_path, src, near)
    b = pydicom.dcmread(base)
    assert (b.BitsAllocated, b.BitsStored, b.HighBit) == (8, 8, 7)
    assert info["pixel_bytes"] == a.nbytes == 4096


def test_roundtrip_implicit_vr_source(tmp_path):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12, ts=ImplicitVRLittleEndian)
    roundtrip(tmp_path, src, 2)


def test_base_decodes_independently_within_bound(tmp_path):
    """The base must open through pydicom's own decoder path, with headers that agree with the code stream."""
    a = RNG.integers(-500, 1500, (40, 56)).astype(np.int16)
    src = make(tmp_path / "in.dcm", a, signed=True, RescaleSlope="1", RescaleIntercept="-24")
    base, res, out, _ = roundtrip(tmp_path, src, 3)
    b = pydicom.dcmread(base)
    try: px = b.pixel_array
    except Exception as e: pytest.skip(f"no independent JPEG-LS decoder installed: {e}")
    hu_base = px.astype(np.float64) * float(b.RescaleSlope) + float(b.RescaleIntercept)
    hu_orig = a.astype(np.float64) * 1 + (-24)
    assert np.abs(hu_base - hu_orig).max() <= 3


def test_bound_is_reported_in_rescaled_units(tmp_path):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12, RescaleSlope="2.5", RescaleIntercept="0")
    info = sd.encode(src, str(tmp_path / "b.dcm"), str(tmp_path / "r.sdr"), 4)
    assert info["bound_rescaled_units"] == 10.0


# ---------------------------------------------------------------- integrity (review finding A)

def test_modified_base_is_rejected(tmp_path):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12)
    base, res, out, _ = roundtrip(tmp_path, src, 2)
    b = pydicom.dcmread(base); b.PatientID = "SOMEONE_ELSE"; b.PixelSpacing = [9.0, 9.0]; b.save_as(base)
    with pytest.raises(sd.SplitDicomError, match="does not match this residual"):
        sd.decode(base, res, str(tmp_path / "out2.dcm"))
    assert not os.path.exists(tmp_path / "out2.dcm")


def test_wrong_pair_with_identical_pixels_is_rejected(tmp_path):
    a = ct12()
    x = make(tmp_path / "x.dcm", a, bits_stored=12, pid="ALICE"); y = make(tmp_path / "y.dcm", a, bits_stored=12, pid="BOB")
    sd.encode(x, str(tmp_path / "x.base"), str(tmp_path / "x.sdr"), 2); sd.encode(y, str(tmp_path / "y.base"), str(tmp_path / "y.sdr"), 2)
    with pytest.raises(sd.SplitDicomError):
        sd.decode(str(tmp_path / "x.base"), str(tmp_path / "y.sdr"), str(tmp_path / "xy.dcm"))
    assert not os.path.exists(tmp_path / "xy.dcm")


def test_metadata_hash_covers_identity_geometry_and_sequences():
    def ds(pid="A", spacing=(1.0, 1.0), code="X"):
        d = Dataset(); d.PatientID = pid; d.PixelSpacing = list(spacing); item = Dataset(); item.CodeValue = code; d.ProcedureCodeSequence = Sequence([item]); return d
    h = sd.metadata_hash(ds())
    assert h == sd.metadata_hash(ds())
    assert len({h, sd.metadata_hash(ds(pid="B")), sd.metadata_hash(ds(spacing=(1.0, 2.0))), sd.metadata_hash(ds(code="Y"))}) == 4


@pytest.mark.parametrize("damage", ["truncate", "flip", "header_len", "magic"])
def test_damaged_residual_is_rejected(tmp_path, damage):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12)
    base, res, out, _ = roundtrip(tmp_path, src, 2)
    raw = bytearray(open(res, "rb").read())
    if damage == "truncate": raw = raw[:len(raw) - 40]
    elif damage == "flip": raw[-20] ^= 0xFF
    elif damage == "header_len": raw[4:8] = struct.pack("<I", 0x7FFFFFFF)
    else: raw[:4] = b"SDR1"
    open(res, "wb").write(bytes(raw))
    with pytest.raises(sd.SplitDicomError):
        sd.decode(base, res, str(tmp_path / "out2.dcm"))
    assert not os.path.exists(tmp_path / "out2.dcm")


# ---------------------------------------------------------------- path safety (review finding B)

def test_same_path_for_everything_is_refused_and_source_survives(tmp_path):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12); before = sha(src)
    for base, res in [(src, src), (src, str(tmp_path / "r.sdr")), (str(tmp_path / "b.dcm"), src), (str(tmp_path / "same"), str(tmp_path / "same"))]:
        with pytest.raises(sd.SplitDicomError):
            sd.encode(src, base, res, 2)
    assert sha(src) == before and sorted(os.listdir(tmp_path)) == ["in.dcm"]


def test_symlink_alias_of_source_is_refused(tmp_path):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12); before = sha(src)
    os.symlink(src, tmp_path / "alias.dcm")
    with pytest.raises(sd.SplitDicomError):
        sd.encode(src, str(tmp_path / "alias.dcm"), str(tmp_path / "r.sdr"), 2, overwrite=True)
    assert sha(src) == before


def test_existing_outputs_need_overwrite(tmp_path):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12)
    base, res = str(tmp_path / "b.dcm"), str(tmp_path / "r.sdr")
    sd.encode(src, base, res, 2); first = sha(base)
    with pytest.raises(sd.SplitDicomError, match="overwrite"):
        sd.encode(src, base, res, 2)
    assert sha(base) == first
    sd.encode(src, base, res, 2, overwrite=True)
    sd.decode(base, res, str(tmp_path / "o.dcm"))
    with pytest.raises(sd.SplitDicomError, match="overwrite"):
        sd.decode(base, res, str(tmp_path / "o.dcm"))
    with pytest.raises(sd.SplitDicomError):
        sd.decode(base, res, base, overwrite=True)


def test_bad_destination_leaves_nothing_behind(tmp_path):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12)
    with pytest.raises(sd.SplitDicomError):
        sd.encode(src, str(tmp_path / "nodir" / "b.dcm"), str(tmp_path / "r.sdr"), 2)
    assert sorted(os.listdir(tmp_path)) == ["in.dcm"]


# ---------------------------------------------------------------- bounds and scope (review findings D, E)

@pytest.mark.parametrize("near", [-1, 128, 255, 100000])
def test_out_of_range_bound_is_refused_before_writing(tmp_path, near):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12)
    with pytest.raises(sd.SplitDicomError, match="near"):
        sd.encode(src, str(tmp_path / "b.dcm"), str(tmp_path / "r.sdr"), near)
    assert sorted(os.listdir(tmp_path)) == ["in.dcm"]


def _unsupported_cases():
    def lut(ds):
        item = Dataset(); item.add_new(0x00283002, "US", [16, 0, 16]); item.ModalityLUTType = "US"; item.add_new(0x00283006, "OW", np.arange(16, dtype=np.uint16).tobytes()); ds.ModalityLUTSequence = Sequence([item])
    def overlay(ds):
        ds.add_new(0x60000010, "US", 8); ds.add_new(0x60000011, "US", 8); ds.add_new(0x60000100, "US", 16); ds.add_new(0x60000102, "US", 12)
    return {
        "32-bit": dict(arr=RNG.integers(0, 100000, (16, 16)).astype(np.uint32), bits_alloc=32),
        "big-endian": dict(arr=ct12((16, 16)), ts=ExplicitVRBigEndian),
        "rgb-declared": dict(arr=ct12((16, 16)), photometric="PALETTE COLOR"),
        "multi-frame": dict(arr=np.stack([ct12((16, 16))] * 3), NumberOfFrames=3),
        "signed+modality-lut": dict(arr=RNG.integers(-50, 50, (16, 16)).astype(np.int16), signed=True, edit=lut),
        "embedded-overlay": dict(arr=ct12((16, 16)), bits_stored=12, edit=overlay),
        "bits-above-stored": dict(arr=(ct12((16, 16)) | 0x8000).astype(np.uint16), bits_stored=12),
    }


@pytest.mark.parametrize("name", list(_unsupported_cases()))
def test_unsupported_input_is_rejected_before_writing(tmp_path, name):
    kw = _unsupported_cases()[name]; arr = kw.pop("arr")
    src = make(tmp_path / "in.dcm", arr, **kw); before = sha(src)
    with pytest.raises(sd.SplitDicomError, match="unsupported input"):
        sd.encode(src, str(tmp_path / "b.dcm"), str(tmp_path / "r.sdr"), 2)
    assert sorted(os.listdir(tmp_path)) == ["in.dcm"] and sha(src) == before


# ---------------------------------------------------------------- CLI

def test_cli_exit_codes(tmp_path, capsys):
    src = make(tmp_path / "in.dcm", ct12(), bits_stored=12)
    b, r, o = (str(tmp_path / n) for n in ("b.dcm", "r.sdr", "o.dcm"))
    assert sd.main(["encode", src, b, r, "--near", "2"]) == 0
    assert sd.main(["encode", src, b, r, "--near", "2"]) == 2            # exists, no --overwrite
    assert sd.main(["encode", src, src, r, "--near", "2", "--overwrite"]) == 2
    assert sd.main(["decode", b, r, o]) == 0
    assert sd.main(["verify", src, o]) == 0
    other = make(tmp_path / "other.dcm", ct12(), bits_stored=12)
    assert sd.main(["verify", src, other]) == 1
