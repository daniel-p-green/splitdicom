# 02. Split-fidelity storage: design and measured ratios

Date: 2026-11-11. Corrected 2026-11-12 after an independent review (see note 06). Corrections are marked.

## Design

Two layers per image, nothing discarded.

- **Base layer.** Encoded with a hard per-pixel error bound: no pixel differs from the original by more than ±N stored units. Uses JPEG-LS near-lossless, which DICOM already defines. The bound is a guarantee on per-pixel intensity error. **Correction:** an earlier version said the base layer therefore "cannot invent structure". That overstated it. A bounded intensity error can still change a thresholded region or a faint edge, and note 05 records a case where it did.
- **Residual layer.** The exact difference between base and original, compressed losslessly.
- Base plus residual reproduces the original pixels bit for bit, checked by SHA-256.

The intended use is tiering: base layer in fast storage for viewing and AI, residual in archive storage for exact recovery.

## Round-trip check

`experiments/e2e.py` runs `splitdicom.py` encode, decode, and verify on all 988 sample images (3 CT series, 3 MR series, 1 chest radiograph), and opens every base file through pydicom's own JPEG-LS decoder to confirm the bound. With the corrected codec (format version 2), all 988 passed at ±2, ±4, and ±8. Run records: `results/e2e_near2.json`, `results/e2e_near4.json`, `results/e2e_near8.json`.

**Correction:** the first version of the checker could not fail. It printed a count and always exited with success. It now exits non-zero on any failure or when no images are selected.

Hashes and comparisons ignore retired group-length elements. Those are byte counts that change with the encoding.

## Lossless baselines (ratio vs 16-bit uncompressed pixels)

| Series | Slices | JPEG-LS | JPEG 2000 | HTJ2K | JPEG XL |
|---|---|---|---|---|---|
| CT chest 2.0 mm | 143 | 2.74 | 2.58 | 2.47 | 2.79 |
| CT whole-body 5.0 mm | 176 | 3.07 | 2.97 | 2.83 | 3.13 |
| CT 3.75 mm (PET/CT) | 136 | 3.28 | 3.12 | 2.90 | 3.46 |
| MR MPRAGE 1.33 mm | 128 | 2.75 | 2.77 | 2.66 | 2.90 |
| MR T1w 0.67 mm | 384 | 7.96 | 8.06 | 7.71 | 9.02 |
| MR T1 3 mm, 512 px | 20 | 3.36 | 3.41 | 3.22 | 3.57 |
| DX chest radiograph | 1 | 2.06 | 2.03 | 1.92 | 2.03 |

The MR T1w 0.67 mm series is an outlier, likely because of a blank background.

## Complete files from the shipping tool (added 2026-11-12)

The tables further down are exploratory: they count code-stream bytes only, and one of them picks the best of three residual coders per image. These figures are what `splitdicom.py` writes to disk, headers included, compared with a like-for-like lossless file (the same DICOM header around a lossless JPEG-LS code stream).

CT, three series, 455 images:

| Bound | Base file vs 16-bit pixels | Base file vs lossless file | Residual file vs lossless file | Base + residual vs lossless file |
|---|---|---|---|---|
| Lossless JPEG-LS file | 2.71x to 3.20x | 100% | n/a | n/a |
| ±2 | 3.92x to 4.98x | 64% to 69% | 35% to 42% | +4.6% to +6.1% |
| ±4 | 4.68x to 6.09x | 52% to 58% | 46% to 55% | +4.4% to +7.5% |
| ±8 | 5.84x to 7.83x | 41% to 46% | 58% to 69% | +4.6% to +9.6% |

Headers matter for small images. In the MR MPRAGE series each file carries about 67 KB of metadata, which is 70% to 82% of the base file. There the base is only 16% to 28% smaller than a lossless file. The split helps least where pixels are a small share of the file.

## Two-layer results on CT (exploratory, code stream only)

Base coder: uniform quantizer plus JPEG XL lossless. Every fourth slice. Source: `results/sweep.json`.

| Bound | Base layer | Residual | Total size vs lossless JPEG XL |
|---|---|---|---|
| ±2 | 4.1x to 5.7x | 8.3x to 8.4x | +1.1% to +2.4% |
| ±4 | 5.0x to 7.2x | 6.0x to 6.1x | +2.0% to +4.9% |
| ±8 | 6.4x to 9.7x | 4.7x to 4.8x | +2.5% to +8.2% |
| ±16 | 8.6x to 13.4x | 3.9x to 4.0x | +3.9% to +13.7% |

With the JPEG-LS near-lossless base coder that `splitdicom.py` ships (all slices, `results/results.json`), the base layer at ±2 was 3.98x to 5.18x and the residual about 8.5x.

## What this means

- An early estimate assumed an 8x base layer at a tight bound. That was wrong. At ±2 the base is about 4x to 5x.
- Reaching about 8x on CT takes a bound near ±8.
- A rough noise estimate for these scans was 7 to 19 units, so ±8 is about half to all of the noise level. That needs validation (see note 03).
- The residual is the smaller layer at tight bounds and becomes the larger one near ±8.

## Illustrative monthly cost per 1 TB of uncompressed CT

**Correction: this table does not describe AWS HealthImaging, and the saving it implies is unproven.** Two facts checked on 2026-11-12 against AWS's pricing page and developer guide:

- Both HealthImaging tiers give millisecond access, and data moves to the cheaper tier automatically after 30 days without access. Speed is not a reason to keep a separate base layer there.
- HealthImaging transcodes JPEG-LS near-lossless files to lossless HTJ2K on import. The base-file sizes measured here are not what it would store or bill.

Read the table only as arithmetic for a generic store with two prices where the base stays at the high price and the residual at the low price. Note 06 gives the comparison that would test whether the split beats ordinary lossless storage with tiering.

Uses AWS HealthImaging list prices quoted during the investigation ($0.105 per GB-month frequent access, $0.006 archive). Not re-verified. Middle CT series.

| Approach | Cost |
|---|---|
| Uncompressed, frequent access | about $108 |
| Lossless, frequent access | about $34 |
| Split at ±2 (base frequent, residual archive) | about $24 |
| Split at ±8 | about $15 |
| Lossless, all in archive | about $2 |

Plain lossless plus archive tiering beats everything on cost. The split helps only for data that must stay in fast storage.

## Limits

- Seven public sample series.
- Single-frame grayscale only.
- For signed pixel data the base file shifts values to unsigned and adjusts the rescale intercept. Originals are restored on decode.
- Base-plus-residual coding is a known technique in general imaging. Whether anyone ships it as a tiered DICOM archive was not confirmed.
