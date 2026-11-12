# 07. Where a real reduction in total bytes can come from

Date: 2026-11-12

## Why this note exists

The two-layer design keeps a base and a residual, and together they are 4% to 10% larger than one lossless file (note 02). It moves bytes between price tiers. It never reduces them. The goal is lower storage cost, so this note measures what could shrink the total.

Script: `experiments/headroom.py`. Output: `results/headroom.json`. Three CT series, every fourth slice (114 slices), code-stream bytes, ratios against 16-bit uncompressed pixels.

## 1. Lossless is close to its limit

| Series | JPEG-LS (in the DICOM standard) | JPEG XL, maximum effort | JPEG XL with inter-slice prediction |
|---|---|---|---|
| CT chest 2.0 mm | 2.74x | 2.81x | 2.81x |
| CT whole-body 5.0 mm | 3.05x | 3.14x | 3.14x |
| CT 3.75 mm (PET/CT) | 3.29x | 3.53x | 3.53x |

- A stronger coder gains 2% to 7%.
- Predicting each slice from the previous one gained nothing on any slice. The noise in neighbouring slices is independent, so the difference between slices is no cheaper to code than the slice itself.
- Published learned lossless coders report roughly 9% to 15% over standard codecs (note 01). That was not tested here.

Conclusion: lossless methods offer at most a further 5% to 15%. Nothing lossless turns 3x into 10x.

## 2. A single bounded-error file, residual discarded

This is the only option measured so far that cuts total bytes substantially. The original cannot be recovered afterwards.

Coder: uniform quantizer (step 2N + 1) followed by lossless JPEG XL. Every pixel stays within ±N stored units.

| Bound | Ratio vs uncompressed | Bytes saved vs best lossless |
|---|---|---|
| Lossless (JPEG XL) | 2.81x to 3.53x | 0% |
| ±2 | 4.19x to 5.88x | 33% to 40% |
| ±4 | 5.08x to 7.51x | 45% to 53% |
| ±8 | 6.51x to 10.29x | 57% to 66% |
| ±16 | 8.95x to 14.76x | 69% to 76% |

The quantizer plus JPEG XL coder is 5% to 19% smaller than JPEG-LS near-lossless at the same bound, and the gap widens as the bound grows.

## 3. What that is worth

Per 1 TB of uncompressed CT per month, middle series, at $0.105 and $0.006 per GB-month, storage charges only:

| Stored as | Expensive tier | Cheap tier |
|---|---|---|
| Lossless | about $34 | about $1.96 |
| ±4, single file | about $17 | about $0.99 |
| ±8, single file | about $13 | about $0.73 |

The percentage saving is the same in both tiers. The money is very different. In the cheap tier, ±8 saves about $1.20 per TB a month, roughly $15,000 a year per petabyte of uncompressed data. In the expensive tier it saves about $21 per TB a month, roughly $260,000 a year per petabyte. The value of any further compression depends almost entirely on the price paid per gigabyte.

## 4. What stands between this and a usable result

- **It is lossy.** Evidence so far (notes 03 to 05): two AI tasks were stable through ±8 on these samples, and threshold-based measurement of very small nodules was not stable at any bound. No radiologist has reviewed anything.
- **Rules on lossy storage of diagnostic images vary by country and modality.** They were not checked for this note.
- **It has not been compared with the lossy codecs already in use.** Lossy JPEG 2000 at similar ratios has been deployed for years. If a bounded-error file is no better at the same size, there is no technical advantage here, only a policy choice anyone can already make. That comparison is the next experiment.

## Next experiment

At equal file size, compare bounded-error coding with lossy JPEG 2000, lossy HTJ2K, and lossy JPEG XL on the existing harness: maximum pixel error, organ segmentation agreement, nodule detection, and nodule measurement. The question is whether a guaranteed error bound preserves downstream results better than codecs that only control average error.
