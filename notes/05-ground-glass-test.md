# 05. Faint nodules: non-solid and part-solid

Date: 2026-11-12

## Question

Note 04 found no lost nodules through ±16, but the one measurement failure was a faint nodule, and that dataset is mostly solid nodules. This note tests the faint ones directly.

## Method

- **Nodule selection:** radiologist texture ratings from the LIDC-IDRI annotations (via the pylidc database). Per-reader annotations were clustered into nodules. Kept: nodules marked by at least two readers, 3 mm or larger, on scans with slices of 2.5 mm or thinner.
  - **Non-solid** (ground-glass): median texture rating 1 or 2.
  - **Part-solid:** median texture rating above 2 and below 4.
- **Scans:** the first 40 scans (by SHA-1 of series ID) containing a non-solid nodule, plus the first 20 others containing a part-solid nodule. Three were skipped for non-standard orientation, leaving 57 scans and 95 nodules: 61 non-solid, 34 part-solid. Four scanner vendors. Ten of the 57 scans also appear in note 04, so the two notes cover 86 distinct scans.
- **Variants, detector, and matching:** same as note 04.
- **Size:** fixed -400 HU volumetry cannot see ground-glass, which is less dense than that. This note adds an adaptive threshold set half-way between the nodule core and the surrounding lung, measured once on the original and reused for every variant.

Scripts: `experiments/select_nodules.py`, `experiments/prep_annotations.py`, `experiments/nodule_test.py`, `experiments/nodule_analyze.py`. Output: `results/ggo_raw.jsonl`, `results/ggo_summary.json`.

Base-layer code-stream ratio for the slices read around these nodules, headers excluded (median, range): 4.9x (3.3x to 6.7x) at ±4, 6.3x (3.9x to 8.9x) at ±8, 8.8x (4.7x to 13.0x) at ±16.

## Results: non-solid nodules (61)

| Variant | Detected | Lost | Largest score change | False positives | Volume change, 95th pct | Largest volume change |
|---|---|---|---|---|---|---|
| Original | 49 | n/a | n/a | 15 | n/a | n/a |
| Control (random ±1) | 49 | 0 | 0.004 | 15 | 0.7% | 3.7% |
| ±4 | 49 | 0 | 0.029 | 15 | 1.4% | 3.7% |
| ±8 | 49 | 0 | 0.022 | 16 | 2.5% | 4.3% |
| ±16 | 49 | 0 | 0.077 | 17 | 5.3% | 15.7% |

## Results: part-solid nodules (34)

| Variant | Detected | Lost | Largest score change | False positives | Volume change, 95th pct | Largest volume change |
|---|---|---|---|---|---|---|
| Original | 27 | n/a | n/a | 19 | n/a | n/a |
| Control (random ±1) | 27 | 0 | 0.003 | 19 | 0.8% | 1.1% |
| ±4 | 27 | 0 | 0.015 | 19 | 0.9% | 91% |
| ±8 | 27 | 0 | 0.026 | 19 | 2.1% | 95% |
| ±16 | 27 | 0 | 0.014 | 18 | 3.7% | 4.7% |

Volume figures use the adaptive threshold. Detection uses a score threshold of 0.5.

## Reading

Scope of every statement below: in this sample, on patches centred on known nodules, under the recorded matching rule.

- **Detection held.** The detector is much less sure of faint nodules (49 of 61 non-solid found on the originals, against 79 of 82 in note 04, and 16 of 95 scores in the uncertain 0.2 to 0.8 range), which makes this a harder test. No nodule was lost or gained at a 0.5 threshold at any bound.
- **Borderline scores do cross other thresholds.** At a 0.7 threshold one non-solid nodule dropped below the line at ±8. At 0.3, one part-solid nodule dropped and one rose at ±8, and one rose in the control as well. Every score change was 0.03 or less through ±8.
- **False positives crept up on non-solid patches:** 15 on the originals, 16 at ±8, 17 at ±16.
- **Adaptive volumetry held through ±8 with one exception.** For 94 of 95 nodules the change was under 5%.
- **The exception matters.** A 3.4 mm part-solid nodule on 2.5 mm slices measured 21 mm3 on the original and about 40 mm3 at both ±4 and ±8. The measured region is under 20 voxels, and the base layer let it connect to a neighbouring structure. It measured 22 mm3 at ±16, so this is instability, not a trend. A threshold-based measurement of a very small region can jump at any bound tested.
- **±16 is visibly worse:** 4 of 61 non-solid nodules changed by more than 5%, the largest by 15.7%.
- **Fixed -400 HU volumetry is unusable on ground-glass** with or without compression. It sees only a few stray voxels, and those flip at ±4. It is reported in the raw output and excluded here.

## What this changes

- The evidence still supports ±4 to ±8 for detection on this model.
- For measurement, a rule keyed to the bound alone is not enough. Very small regions (roughly under 50 mm3) on thick slices can change at ±4. Among regions of 50 mm3 or more, the largest change through ±8 was 4.3%.
- Any serious use needs the exact residual available on demand for measurement tasks. That is what the two-layer design provides.

## Limits

- Same as note 04: the detector was trained on LUNA16, which overlaps these scans. Each run covers a block around a known nodule, not a whole scan. One detector, one measurement method, no radiologist review.
- Zero lost out of 76 detected allows a true loss rate of up to roughly 4%. For the 49 non-solid nodules alone, up to roughly 6%.
- Texture ratings are subjective, and only two readers were required.
- The adaptive threshold is a simple method chosen for this test. Clinical volumetry software differs.

## Data credits

- LIDC-IDRI: Armato et al., The Cancer Imaging Archive, CC BY 3.0. Not redistributed here.
- Annotation database: [pylidc](https://pylidc.github.io/).
- Detector: [MONAI/lung_nodule_ct_detection](https://huggingface.co/MONAI/lung_nodule_ct_detection).
