# 04. Do lung nodules survive the base layer?

Date: 2026-11-12

## Question

Note 03 showed organ outlines hold up through ±8. Organ outlining is forgiving. This note tests small findings: annotated lung nodules, down to 3.3 mm.

## Method

- **Scans:** 40 chest CT series from LIDC-IDRI, chosen by sorting the LUNA16 series IDs by SHA-1 hash and taking the first 40 that have annotations. One was skipped for non-standard orientation, leaving 39 scans and 82 annotated nodules (36 under 6 mm; diameters 3.3 to 24.8 mm, median 6.2). Slice spacing ranged from 0.5 to 2.5 mm. Three scanner vendors.
- **Variants:** original, a control (random -1, 0, or +1 added to every pixel), and base layers at ±4, ±8, ±16 (JPEG-LS near-lossless on the full native slices, same codec as `splitdicom.py`).
- **Detection:** the MONAI `lung_nodule_ct_detection` model (3D RetinaNet), run on one 128 x 128 x 64 patch (90 x 90 x 80 mm) centred on each annotated nodule. A nodule counts as detected when a detection scoring 0.5 or higher lies within the nodule radius (3 mm minimum) of the annotated centre.
- **Size:** fixed-threshold volumetry at -400 HU on the native grid: the connected region nearest the annotated centre, inside a sphere of 0.75 x diameter + 1 mm.

Scripts: `experiments/nodule_test.py`, `experiments/nodule_analyze.py`. Output: `results/nodules_raw.jsonl`, `results/nodule_summary.json`.

## Results: detection (82 nodules)

| Variant | Base layer size (median, range) | Detected | Lost | Gained | Largest score change | False positives in patches |
|---|---|---|---|---|---|---|
| Original | n/a | 79 | n/a | n/a | n/a | 19 |
| Control (random ±1) | n/a | 79 | 0 | 0 | 0.002 | 19 |
| ±4 | 5.3x (3.2x to 8.3x) | 80 | 0 | 1 | 0.015 | 19 |
| ±8 | 6.9x (3.7x to 11.4x) | 80 | 0 | 1 | 0.013 | 19 |
| ±16 | 9.4x (4.5x to 16.0x) | 79 | 0 | 0 | 0.021 | 18 |

- No detected nodule was lost at any bound, at thresholds of 0.3, 0.5, or 0.7.
- The one "gained" nodule is a 5.0 mm nodule whose score moved from 0.495 to 0.510. That is a borderline case crossing the line, not a new finding.
- At ±16 one false positive dropped out: its score fell from 0.54 to below 0.5.
- Detection position moved by at most 0.04 mm through ±8 and 0.22 mm at ±16. Detected box size changed by under 0.4% through ±8.

## Results: size (81 nodules with a measurable region)

| Variant | Median volume change | 95th percentile | Largest | Over 5% | Over 10% |
|---|---|---|---|---|---|
| Control (random ±1) | 0.00% | 0.55% | 6.3% | 1 | 0 |
| ±4 | 0.00% | 1.45% | 2.2% | 0 | 0 |
| ±8 | 0.17% | 1.85% | 8.0% | 2 | 0 |
| ±16 | 0.35% | 2.63% | 132% | 2 | 1 |

- Changes over 5% were all on regions of a few voxels (6 to 18 mm3), where one voxel is several percent. The control produced one such change too.
- The 132% case at ±16 is a faint nodule whose density sits near the -400 HU threshold. Its measured region went from 6.5 to 15.1 mm3. Low-contrast findings are where a large bound first breaks threshold-based measurement.
- Mean density inside the original nodule mask shifted by a median of 0.2 HU at ±8 (largest 1.2 HU) and 0.5 HU at ±16 (largest 3.7 HU).

## Reading

- Through ±8, this detector's output on annotated nodules is unchanged for practical purposes, and threshold volumetry moves by under 2% for 95% of nodules.
- ±16 still detects everything here but distorts the measurement of one faint nodule. That fits note 03, where ±16 was the first bound to degrade small structures.
- Base layers on these scans were smaller than on the note 02 samples: a median of 6.9x at ±8.

## Limits

- **The model was trained on LUNA16, which includes these scans.** It may be unusually confident on them. A model seeing unfamiliar scans could be more sensitive to the bound.
- **Patch-based, not whole-scan.** Each run looks at an 80 to 90 mm block around a known nodule. Whole-scan false-positive behaviour was not measured.
- **Box offset.** With MONAI 1.6.1 this model returns boxes whose midpoint sits about 0.36 x box size past the true centre on each axis (confirmed on synthetic spheres). The cause was not traced. The offset is identical in every variant, so the comparisons above are unaffected. It is corrected only when matching detections to annotations.
- Solid nodules dominate LUNA16. Ground-glass and part-solid nodules, the low-contrast cases, are under-represented, and the one clear failure here was of that kind.
- 82 nodules cannot rule out rare failures. Zero lost out of 79 is consistent with a true loss rate of up to roughly 4%.
- One detector and one measurement method. No radiologist reviewed anything.

## Data and model credits

- LIDC-IDRI: Armato et al., The Cancer Imaging Archive, CC BY 3.0. Downloaded per series from the TCIA API. Not redistributed here.
- LUNA16 annotations: Setio et al., [Zenodo record 3723295](https://zenodo.org/records/3723295).
- Detector: [MONAI/lung_nodule_ct_detection](https://huggingface.co/MONAI/lung_nodule_ct_detection).
