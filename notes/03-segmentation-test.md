# 03. Does the error bound change a downstream model? Organ segmentation

Date: 2026-11-11 to 2026-11-12

## Question

If viewers and AI read the base layer in place of the original, do results change?

## Method

- Model: TotalSegmentator 2.18.0 (open CT model that outlines about 80 to 90 anatomical structures per scan), CPU.
- Scans: the three CT series from note 02. The whole-body scan was cropped to 88 torso slices because the full volume ran out of memory.
- Variants: original, and base layers at ±2, ±4, ±8, ±16 (JPEG-LS near-lossless, same codec as `splitdicom.py`).
- Control: the original with random noise of -1, 0, or +1 added to every voxel. This shows how much the model's output moves from a trivial change.
- Metrics, each variant vs the original's segmentation: Dice overlap per structure (structures with at least 1,000 voxels), and the shift in mean density (HU) inside each original structure mask.

Scripts: `experiments/build_volumes.py`, `experiments/compare.py`. Raw output: `results/compare_fast.json`, `results/compare_full.json`.

## Results: fast model (3 mm), three scans, ranges across scans

| Variant | Base layer size | Median Dice | Worst-structure Dice | Largest organ mean shift |
|---|---|---|---|---|
| Control (random ±1) | n/a | 0.999 to 1.000 | 0.966 to 0.996 | 0.03 to 0.05 HU |
| ±2 | 4.0x to 5.2x | 0.999 to 1.000 | 0.987 to 0.992 | 0.05 to 0.07 HU |
| ±4 | 4.8x to 6.4x | 0.998 to 0.999 | 0.957 to 0.986 | 0.09 to 0.13 HU |
| ±8 | 6.0x to 8.4x | 0.997 to 0.999 | 0.965 to 0.975 | 0.12 to 0.19 HU |
| ±16 | 8.0x to 11.8x | 0.995 to 0.998 | 0.880 to 0.960 | 0.30 to 0.66 HU |

## Results: full-resolution model (1.5 mm), chest scan only

| Variant | Median Dice | Worst-structure Dice |
|---|---|---|
| Control (random ±1) | 1.000 | 0.995 |
| ±8 | 0.999 | 0.984 |
| ±16 | 0.998 | 0.972 |

## Reading

- ±2 is at the model's own jitter level. The random ±1 control moved the output about as much.
- Through ±8, every structure stayed at 0.957 or above and organ mean density moved by under 0.2 HU.
- At ±16, one small structure fell below 0.95 in two of three scans (adrenal gland, carotid artery).
- The worst cases at every bound were small structures: adrenal glands and small vessels.

## Limits

- Organ outlining is a forgiving task. It says little about small or low-contrast findings such as lung nodules, which are the real clinical risk. Note 04 addresses that.
- Three scans.
- No radiologist reviewed anything.
