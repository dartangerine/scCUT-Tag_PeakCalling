#!/usr/bin/env bash
# Usage: bash scripts/call_macs2.sh <bam_dir> <out_dir> [control_dir]
set -euo pipefail
bam_dir=$1; out_dir=$2; control_dir=${3:-}
mark=${MARK:-H3K4me3}; mode=${MODE:-narrow}; gsize=${GENOME_SIZE:-2.7e9}
case "$mode" in narrow|broad) ;; *) echo "MODE must be narrow|broad" >&2; exit 1;; esac
mkdir -p "$out_dir"
while IFS= read -r ct; do
  sample="${mark}_${ct}"
  samtools quickcheck -v "$bam_dir/$sample.bam"
  args=(-t "$bam_dir/$sample.bam" -f BAMPE -g "$gsize" -q 0.05 \
        --keep-dup all -n "$sample" --outdir "$out_dir")
  [[ "$mode" == broad ]] && args+=(--broad --broad-cutoff 0.05)
  [[ -n "$control_dir" ]] && args+=(-c "$control_dir/input_${ct}.bam")
  macs2 callpeak "${args[@]}" > "$out_dir/$sample.log" 2>&1
done < cell_types.txt
