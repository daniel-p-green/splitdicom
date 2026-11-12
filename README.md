# splitdicom

Research notes and prototype code on one question: can DICOM medical images be stored for much less money without losing information that matters?

**Status: research prototype. Not validated for clinical use. No patient data is in this repository.**

## Short answer so far

- Storing pixels "as code" does not save space. Pixel data is almost all of a scan's bytes, and most of those bytes are measurement noise that no program can regenerate exactly.
- Plain lossless compression gives about 2.8x to 3.5x on the CT scans tested here.
- The cheapest lossless option is existing technology: lossless compression plus archive-tier storage.
- The prototype in this repo splits each image into two layers: a small base layer with a hard per-pixel error bound, and an exact residual. Together they reproduce the original bit for bit. The base layer can sit in fast storage and the residual in cheap storage.
- Measured base-layer size on CT: about 4x to 5x at a bound of ±2, and about 6x to 8x at ±8.
- An open segmentation model produced near-identical results on base layers through ±8. Small structures started to drift at ±16.
- A lung-nodule detector kept every nodule it found on the originals (79 of 82 annotated) at ±4, ±8, and ±16. Threshold volumetry moved by under 2% for 95% of nodules through ±8. At ±16 one faint nodule's measured volume more than doubled.
- On 95 faint nodules (61 ground-glass, 34 part-solid) the detector again lost none at any bound. Measured volume stayed within 5% through ±8 for 94 of 95. One very small nodule's measured volume roughly doubled at ±4 and ±8.

Details and limits are in [`notes/`](notes/).

## What is here

| Path | What it is |
|---|---|
| `splitdicom.py` | The two-layer codec: `encode`, `decode`, `verify` |
| `experiments/` | Benchmark and test scripts |
| `results/` | Raw JSON output from each experiment |
| `notes/` | Findings, written up in order |

## Using the codec

```
pip install -r requirements.txt

python splitdicom.py encode in.dcm base.dcm residual.sdr --near 2
python splitdicom.py decode base.dcm residual.sdr restored.dcm
python splitdicom.py verify in.dcm restored.dcm
```

- `base.dcm` is a standard DICOM file (JPEG-LS near-lossless). Every pixel is within ±N stored units of the original.
- `residual.sdr` holds the exact difference plus the tags needed to restore the original.
- `decode` rebuilds an uncompressed DICOM and checks the pixels against a SHA-256 hash recorded at encode time.

Scope: single-frame grayscale images, 8 to 16 bits.

## Test data

The experiments use public sample data that is not redistributed here:

- [OHIF/viewer-testdata](https://github.com/OHIF/viewer-testdata) (CT, MR, radiograph)
- [datalad/example-dicom-structural](https://github.com/datalad/example-dicom-structural) (MR)
- LIDC-IDRI chest CT from The Cancer Imaging Archive (CC BY 3.0), with LUNA16 annotations, for the nodule test in note 04

Clone both into one folder as `ohif/` and `mr/`, then pass that folder to the scripts.

## Limits

- A handful of public sample series plus 96 public chest CTs. Too few to generalize.
- No radiologist has reviewed any output.
- Cost figures in the notes use list prices quoted during the investigation and were not re-verified.
- This repository has no license yet, so default copyright applies.
