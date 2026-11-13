# splitdicom

Research notes and prototype code on one question: can DICOM medical images be stored for much less money without losing information that matters?

**Status: research prototype. Not validated for clinical use. No patient data is in this repository.**

## Short answer so far

- Storing pixels "as code" does not save space. Pixel data is almost all of a scan's bytes, and most of those bytes are measurement noise that no program can regenerate exactly.
- Plain lossless compression gives about 2.7x to 3.5x on the CT scans tested here.
- The cheapest lossless option is existing technology: lossless compression plus automatic storage tiering. **Whether the two-layer design below saves money over that is unproven** (see "Storage cost" below).
- The prototype splits each image into a base layer with a hard per-pixel error bound and an exact residual. Together they reproduce the original pixel samples and metadata.
- Measured on complete files from the shipping tool, a CT base file is about 64% to 69% of the size of a lossless file at a bound of ±2, and 41% to 46% at ±8. The base and residual together take 4% to 10% more space than one lossless file.
- In the samples tested, an open segmentation model produced near-identical results on base layers through ±8, and a nodule detector kept every detection it made on the original patches at ±4, ±8, and ±16. Threshold-based size measurement of very small nodules was not stable at any bound.

- The only large cut in total bytes comes from keeping one lossy file and discarding the exact original (57% to 66% fewer bytes than lossless at ±8). At the same file size, bounded-error coding has far lower pixel error than lossy JPEG 2000, HTJ2K, or JPEG XL, but it showed no advantage on the downstream tasks tested, and a plain quantizer was worse for threshold-based nodule measurement (note 08).

Details and limits are in [`notes/`](notes/). Note 06 records an independent code review and the fixes that followed.

## What the codec guarantees

**Exact decoded pixel samples and verified data-set metadata, for the supported inputs below.** The original file bytes and transfer syntax are not recreated: `decode` always writes uncompressed Explicit VR Little Endian.

A residual file is bound to one base file and one source image by SHA-256 hashes of the base file, the original pixels, and the original metadata. Decoding fails, and writes nothing, if the pair does not match, if the base was modified, or if the restored file does not read back as the original.

A bound of ±N is in stored pixel units. In rescaled units such as HU it is ±N x |RescaleSlope|. The bound limits per-pixel intensity error. It does not guarantee that thresholded regions, faint edges, or downstream measurements are unchanged (note 05 has a counterexample).

| Input | Status |
|---|---|
| Single-frame, one sample per pixel, MONOCHROME1 or MONOCHROME2 | Supported |
| 8 or 16 bits allocated, any bits stored, signed or unsigned | Supported |
| Little-endian transfer syntaxes that pydicom can decode | Supported |
| Linear rescale (slope and intercept), pixel padding value and range | Supported |
| Multi-frame, colour, palette colour | Rejected |
| 32-bit integer, float, or double pixel data | Rejected |
| Explicit VR Big Endian | Rejected |
| Signed pixels with a Modality LUT Sequence | Rejected |
| Overlays or other data in the unused bits above BitsStored | Rejected |
| Error bounds outside 0 to 127 | Rejected |

Unsupported inputs are rejected before anything is written.

## What is here

| Path | What it is |
|---|---|
| `splitdicom.py` | The two-layer codec: `encode`, `decode`, `verify` |
| `tests/` | Regression tests on synthetic images (no downloads) |
| `experiments/` | Benchmark and evaluation scripts |
| `results/` | Raw output from each experiment, and the package versions used |
| `notes/` | Findings, written up in order |

## Using the codec

```
pip install -r requirements.txt

python splitdicom.py encode in.dcm base.dcm residual.sdr --near 2
python splitdicom.py decode base.dcm residual.sdr restored.dcm
python splitdicom.py verify in.dcm restored.dcm
```

- `base.dcm` is a standard DICOM file (JPEG-LS near-lossless) with its own SOP Instance UID.
- `residual.sdr` holds the exact difference, the original value of every element the base changes, and the binding hashes.
- Existing files are never overwritten unless you pass `--overwrite`, and an output can never be one of the inputs.
- Both outputs are staged, the pair is restored once in memory, and only then are the files moved into place. The two moves are not one atomic step. If the process dies between them, the residual's hash of the base still detects an incomplete pair.

## Tests

```
pip install -r requirements-dev.txt
python -m pytest tests -q
```

## Reproducing the experiments

```
pip install -r requirements-research.txt
python experiments/e2e.py <data_root> --near 2          # every sample image: encode, decode, verify; non-zero exit on failure
experiments/run_segmentation.sh <data_root> <work_dir>  # note 03
python experiments/nodule_test.py ...                   # notes 04 and 05; see the script header
```

Exact package versions for the published runs are in `results/environment.txt`.

## Test data

The experiments use public data that is not redistributed here:

- [OHIF/viewer-testdata](https://github.com/OHIF/viewer-testdata) (CT, MR, radiograph), cloned as `<data_root>/ohif`
- [datalad/example-dicom-structural](https://github.com/datalad/example-dicom-structural) (MR), cloned as `<data_root>/mr`
- LIDC-IDRI chest CT from The Cancer Imaging Archive (CC BY 3.0), with LUNA16 and pylidc annotations, for notes 04 and 05

## Storage cost

An earlier version of these notes priced the design as "base in fast storage, residual in archive storage" using AWS HealthImaging rates. That comparison does not hold for that service:

- Both of its storage tiers give millisecond access, and data moves to the cheaper tier automatically after 30 days without access. Fast access is not a reason to keep a separate approximate copy there.
- It transcodes JPEG-LS near-lossless files to lossless HTJ2K on import, so the base-file sizes measured here are not what it would store or bill.

The claim still worth testing is narrower: frequently read approximate images might cost less to serve if exact recovery is rare and residuals are cheap to keep. That needs a concrete deployment and real access patterns. Note 06 gives the comparison to run.

## Limits

- A handful of public sample series plus 86 public chest CTs. Too few to generalize.
- The nodule tests look at patches around known nodules, with a detector trained on overlapping data. Whole-scan behaviour and held-out performance are unmeasured.
- No radiologist has reviewed any output.
- This repository has no license yet, so default copyright applies.
