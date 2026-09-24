#!/usr/bin/env bash
# Usage: bash scripts/make_pseudo_inputs.sh <bam_dir> <out_dir>
set -euo pipefail
bam_dir=$1; out_dir=$2
mark=${MARK:-H3K4me3}; threads=${THREADS:-8}
mapfile -t cells < cell_types.txt
(( ${#cells[@]} >= 2 )) || { echo "Need >= 2 cell types" >&2; exit 1; }
mkdir -p "$out_dir"
for ct in "${cells[@]}"; do
  samtools quickcheck -v "$bam_dir/${mark}_${ct}.bam"
done
for target in "${cells[@]}"; do
  sources=()
  for ct in "${cells[@]}"; do
    [[ "$ct" == "$target" ]] && continue
    sources+=("$bam_dir/${mark}_${ct}.bam")
  done
  samtools merge -@ "$threads" -u - "${sources[@]}" |
    samtools sort -@ "$threads" -o "$out_dir/input_${target}.bam" -
  samtools index "$out_dir/input_${target}.bam"
done
