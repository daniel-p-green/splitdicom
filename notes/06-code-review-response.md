# 06. Independent code review: what it found and what changed

Date: 2026-11-12

An independent review of this repository at commit `fbad0e3` reported correctness bugs in the codec, gaps in reproducibility, and claims that outran the evidence. Every codec defect it reported was reproduced here with synthetic images before any change was made. This note records each finding, the fix, and what is still open.

## Codec defects (all reproduced, all fixed, all covered by tests)

| Finding | Reproduced | Fix |
|---|---|---|
| The pixel hash passed while patient identity or geometry was wrong. A modified base, or a base paired with another image's residual, decoded "successfully". | Yes | The residual now records SHA-256 hashes of the base file, the original pixels, and the original metadata (every element, nested sequences included). Decode checks all three, writes to a temporary file, reads it back, and re-checks before publishing. |
| Passing the same path for source, base, and residual overwrote the source and returned success. A bad base path left an orphaned residual. | Yes | Outputs may not be inputs or each other (symlinks resolved). Existing files are refused without `--overwrite`. Both outputs are staged, the pair is restored once in memory, then moved into place. |
| 8-bit input produced a base with BitsStored 16 and BitsAllocated 8, which pydicom refuses to read. Byte counts assumed two bytes per pixel. | Yes | 8-bit images are coded as 8-bit. The base header always agrees with the JPEG-LS stream precision, and encode checks this. Byte counts are real. |
| Bounds above 127 encoded without error and could not be restored. | Yes | Bounds outside 0 to 127 are refused before anything is written. |
| 32-bit and big-endian inputs failed late. Signed images with a Modality LUT got a second, conflicting transform. Pixel padding range was left behind in an unsigned base. Overlay bits stored above BitsStored were dropped silently. | Yes (32-bit, padding range). Accepted on inspection (the rest). | A preflight check rejects everything outside a stated supported subset before any output. Padding value, padding range, and smallest/largest pixel value are removed from the base and restored on decode. |

Two further problems surfaced while fixing these:

- Original element values were stored as JSON, which rewrote a decimal string such as `-24` as `-24.0`. The old verifier compared numerically and missed it. Originals are now stored as binary DICOM, and the metadata hash compares text exactly.
- pydicom masks the bits above BitsStored when decoding, so stray data there cannot be seen in the pixel array. The check reads the raw bytes.

`tests/test_splitdicom.py` holds 38 tests on synthetic images. With the new codec, all 988 real sample images pass encode, decode, verify, and an independent decode of the base at ±2, ±4, and ±8.

The residual format changed (version 2). Version 1 files are not readable by this version.

## What "lossless" means here

The review separated three promises. The project now states which one it makes:

| Promise | Status |
|---|---|
| Exact decoded pixel samples | Guaranteed for supported inputs, checked by hash |
| Identity, geometry, calibration, and all other data-set elements | Guaranteed for supported inputs, checked by hash |
| The original file's bytes and transfer syntax | Not provided |

## Reproducibility

| Finding | Fix |
|---|---|
| The end-to-end checker always exited with success, even with no data. | It now reports discovered, selected, passed, skipped, and failed counts, and exits non-zero on any failure or an empty selection. |
| The segmentation scripts did not form a complete path to note 03: the control volumes, the torso crop, and the model commands were run by hand. | `experiments/build_volumes.py` now writes the control and the crop, `experiments/run_segmentation.sh` holds every command, and `experiments/compare.py` fails on a missing file. The rebuilt volumes are identical to the ones used for the published run (18 of 18), and the full-resolution comparison reproduces `results/compare_full.json` exactly. |
| The nodule runner resumed any existing output file regardless of configuration. | Each output starts with a run record (hashes of code, weights, and annotations, settings, package versions). Resuming under a different fingerprint is refused. Missing scans are recorded as skipped. |
| Dependencies were unpinned and the core codec pulled in the research stack. | Requirements are split into core, dev, and research. `results/environment.txt` records the versions used. |

Not done: the published `nodules_raw.jsonl` and `ggo_raw.jsonl` predate the run record and were not regenerated. The fast-model segmentation comparison was not rerun end to end (about 10 minutes of compute plus 18 model runs).

## Claims corrected

- **"The base layer cannot invent structure."** Withdrawn. A bounded intensity error does not guarantee that thresholded regions or faint edges survive. Note 05 already contained a counterexample.
- **Compression figures.** The earlier benchmark picked the best of three residual coders per image and counted code streams only, while the tool always uses one coder and writes headers. Note 02 now leads with complete-file figures from the shipping tool. The combined base and residual cost 4% to 10% more than a lossless file on CT, not the 1% to 8% reported before.
- **Nodule-test compression ratios** covered only the slices read around each nodule. Notes 04 and 05 now say so.
- **Storage cost.** The AWS HealthImaging comparison was wrong in two ways, both confirmed against AWS's own pages: both tiers give millisecond access with automatic tiering after 30 days, and near-lossless JPEG-LS is transcoded to lossless HTJ2K on import. The table in note 02 is now labelled as generic arithmetic only.
- **Detection counts** depend on an unexplained box-offset correction. Note 04 now reports the counts with and without it. The paired result (no detection lost) holds either way.

## Still open

1. **The box offset.** The detector's boxes sit about 0.36 x box size past the true centre. The review found no missing setting against the official configuration. The next step is to run the official pipeline and this one on synthetic objects with known positions. All six box coordinates are now saved for that.
2. **Held-out, whole-scan evaluation.** The detector was trained on data that overlaps the test scans, and only patches around known nodules were run. Needed: held-out cases, full-scan inference, one-to-one matching, false positives per scan, uncertainty that accounts for several nodules per scan, and tolerances fixed in advance.
3. **Whether the split saves money.** Unproven. The comparison to run, for a chosen deployment:

   - Split: `H x B + C x R + overhead`
   - Ordinary lossless: `L x (p x H + (1 - p) x C) + overhead`

   `B`, `R`, and `L` are the stored bytes of the base, the residual, and a conventional lossless image. `H` and `C` are the high and low storage prices. `p` is the share of time conventional data spends at the high price. Overhead covers requests, transfer, compute, and operations. On the 455 CT files measured here at ±8, `B` is 0.43 x `L` and `R` is 0.63 x `L`. With prices of $0.105 and $0.006 per GB-month and overhead ignored, the split breaks even at `p` of about 0.43. It wins only for data that would otherwise spend more than about 43% of its time at the high price, and overhead moves that threshold up.
4. **Atomic pair publication.** The two output files are moved one after the other. A manifest or committed directory would make completion explicit.
5. **Packaging.** No `pyproject.toml` or installed command yet. No license yet.
6. **A study-level workflow** (inventory a study, flag unsupported images, encode, verify, restore) is the natural next feature and would exercise the intended storage lifecycle.

## Sources for the storage facts

- [AWS HealthImaging pricing](https://aws.amazon.com/healthimaging/pricing/)
- [AWS HealthImaging supported transfer syntaxes](https://docs.aws.amazon.com/healthimaging/latest/devguide/supported-transfer-syntaxes.html)
