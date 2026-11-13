# 08. At the same file size, is bounded-error coding better than existing lossy codecs?

Date: 2026-11-12

## Why this matters

Note 07 showed that the only large cut in total bytes comes from keeping one lossy file. Lossy codecs at these ratios already exist and are already used for medical images. If bounded-error coding is no better at the same size, this project has no technical advantage, only a policy choice that is already available.

## Method

Each codec is given the same number of bytes, then judged on the same tests as notes 03 to 05.

- **Bounded-error (q4, q8):** uniform quantizer with step 2N + 1, then lossless JPEG XL. Every pixel within ±4 or ±8. This is the smallest bounded-error coder measured in note 07.
- **Lossy JPEG 2000:** irreversible 9/7 wavelet (OpenJPEG).
- **Lossy HTJ2K** (OpenJPH).
- **Lossy JPEG XL** (VarDCT). It is tuned for how photographs look on a display, which is not what CT values are.
- **Size matching:** each lossy codec gets one quality setting per series or scan, tuned on a sample of slices to match the bounded-error bytes. Achieved totals are within 0.2% of the target.

Scripts: `experiments/compare_lossy.py`, `experiments/nodule_test_lossy.py`. Output: `results/lossy_compare_pixel.json`, `results/lossy_compare_segmentation.json`, `results/lossy_compare_nodules_raw.jsonl`, `results/lossy_compare_nodules_summary.json`.

## Result 1: pixel error (three CT series, all 452 slices)

At the size of the ±8 bounded-error file:

| Codec | Largest pixel error | RMSE | Pixels off by more than 8 |
|---|---|---|---|
| Bounded-error | 8 | 4.3 | 0% |
| Lossy JPEG 2000 | 39 to 48 | 4.9 to 6.1 | 10% to 17% |
| Lossy HTJ2K | 31 to 43 | 5.3 to 6.9 | 12% to 22% |
| Lossy JPEG XL | 67 to 153 | 5.3 to 7.9 | 12% to 17% |

At the ±4 size the pattern is the same: largest error 4 against 20 to 92, RMSE 2.3 against 2.8 to 5.0.

Bounded-error coding wins clearly on both worst-case and average pixel error.

## Result 2: organ segmentation (three scans, fast model)

| Codec, at the ±8 size | Voxels whose label changed (mean of 3 scans) | Structures below 0.95 Dice | Worst-structure Dice |
|---|---|---|---|
| Bounded-error | 0.34% | 0 | 0.970 to 0.975 |
| Lossy JPEG 2000 | 0.38% | 2 | 0.942 to 0.974 |
| Lossy HTJ2K | 0.44% | 2 | 0.946 to 0.983 |
| Lossy JPEG XL | 0.68% | 2 | 0.943 to 0.960 |

A small edge for bounded-error coding. On one of the three scans JPEG 2000 changed fewer voxels. Three scans cannot separate these.

## Result 3: faint nodules (95 nodules, 57 scans, at the ±8 size)

Detection at a 0.5 score threshold: all three codecs detected the same 76 nodules as the original. None lost, none gained.

Volume by adaptive threshold, change against the original:

| Codec | Median | 95th percentile | Over 5% | Over 10% | Largest |
|---|---|---|---|---|---|
| Bounded-error (quantizer) | 2.0% | 7.5% | 16 | 2 | 100% |
| Lossy JPEG 2000 | 0.4% | 4.3% | 4 | 2 | 100% |
| Lossy JPEG XL | 0.6% | 3.4% | 2 | 0 | 5.3% |

Pixel error on these same scans still favoured the quantizer: largest error 8 against 30 to 52 (JPEG 2000) and 41 to 207 (JPEG XL), median RMSE 4.4 against 6.0 and 6.5.

**The quantizer is worse for measurement despite having lower pixel error.** The likely reason is the shape of its error. A uniform quantizer snaps every pixel to a multiple of 17. Wherever a measurement threshold falls between two levels, every pixel in that 17-unit band moves to the same side together, which shifts the measured edge. Transform codecs spread their error like noise, which mostly averages out across a region. This explanation fits the data and was not tested directly.

Two more observations support it:

- JPEG-LS near-lossless at the same ±8 bound, which does not snap to a fixed grid, did much better in note 05 on the same nodules: 95th percentile 2.5% and largest 4.3% for the non-solid group, against 7.5% and 9.4% for the quantizer here. Its files are 9% to 23% larger at that bound, so that is not an equal-size comparison.
- Mean density inside each nodule moved least with the quantizer (median 0.17 HU against 0.28 and 0.47). Averages are fine. Thresholds are not.

The 100% case is one very small nodule whose measured region disappeared under both the quantizer and JPEG 2000.

## Conclusion

- **A hard bound on pixel error is real and existing lossy codecs do not offer it.** At the same size their worst pixel is 4 to 20 times further off.
- **That advantage did not carry through to the tasks tested.** Detection was identical, segmentation was nearly identical, and threshold-based measurement was worse with the plain quantizer.
- **Smaller pixel error is not the same as better downstream results.** The structure of the error matters more than its size. The claim in note 07 that the quantizer coder is simply better than JPEG-LS near-lossless because it is smaller does not hold once measurement is considered.
- On this evidence, bounded-error coding is not a clear improvement over lossy JPEG 2000 for cutting storage. Lossy JPEG 2000 at 6x to 10x is available today.

## What could still change this

- **Bounded error without the grid.** Subtractive dither (add a reproducible pseudo-random offset before quantizing, subtract it after decoding) keeps the same hard bound and makes the error independent of the image, which should remove the threshold bias. It will cost bytes in flat regions. Whether it keeps the pixel-error advantage and matches JPEG 2000 on measurement at equal size is testable with this harness.
- **JPEG-LS near-lossless against JPEG 2000 at equal size.** Not run. It measured well in note 05 but was never size-matched.

## Limits

- Three scans for pixel error and segmentation, 57 scans for nodules, patches around known nodules, one detector trained on overlapping data, one measurement method, no radiologist review.
- Lossy codec settings were tuned to match size, not tuned for CT. Better-configured encoders may do better than shown.
- The nodule run was completed in two parts after an out-of-memory failure, so the raw file contains two run records. The solid-nodule set from note 04 was not run for this comparison.
