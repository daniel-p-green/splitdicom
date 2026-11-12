#!/usr/bin/env bash
# Reproduces notes/03 end to end. Needs about 8 GB RAM; CPU is enough.
#   usage: experiments/run_segmentation.sh <data_root> <work_dir>
# <data_root> holds ohif/ (OHIF/viewer-testdata) as described in README.md.
# TotalSegmentator downloads its own weights on first use. Published run: TotalSegmentator 2.18.0, torch 2.15.0 (CPU).
set -euo pipefail
ROOT=$1; WORK=$2; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK/vols" "$WORK/seg" "$WORK/segfull"
python3 "$HERE/build_volumes.py" "$ROOT" "$WORK/vols"
seg() {  # one run at a time: a killed run can leave worker processes holding memory
  TotalSegmentator -i "$1" -o "$2" --ml --device cpu --nr_thr_resamp 1 --nr_thr_saving 1 ${3:-}
}
for n in chest2mm wbtorso petct375; do for k in orig rand1 2 4 8 16; do
  [ -f "$WORK/seg/${n}_${k}.nii.gz" ] || seg "$WORK/vols/${n}_${k}.nii.gz" "$WORK/seg/${n}_${k}.nii.gz" --fast
done; done
for k in orig rand1 8 16; do   # full-resolution model, chest scan only (about 7 minutes per run on 2 cores)
  [ -f "$WORK/segfull/chest2mm_${k}.nii.gz" ] || seg "$WORK/vols/chest2mm_${k}.nii.gz" "$WORK/segfull/chest2mm_${k}.nii.gz"
done
python3 "$HERE/compare.py" "$WORK/vols" "$WORK/seg" "$WORK/compare_fast.json" chest2mm,wbtorso,petct375 rand1,2,4,8,16
python3 "$HERE/compare.py" "$WORK/vols" "$WORK/segfull" "$WORK/compare_full.json" chest2mm rand1,8,16
