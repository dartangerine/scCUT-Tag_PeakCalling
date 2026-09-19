#!/usr/bin/env bash
# Usage: bash scripts/call_sicer2.sh <bam_dir> <out_dir> [control_dir]
set -euo pipefail
bam_dir=$1; out_dir=$2; control_dir=${3:-}
mark=${MARK:-H3K4me3}; mode=${MODE:-narrow}
genome=${GENOME:-hg38}; threads=${THREADS:-8}
case "$mode" in narrow) w=50; g=100;; broad) w=100; g=200;; *) exit 1;; esac
while IFS= read -r ct; do
  sample="${mark}_${ct}"; run="$out_dir/$sample"
  mkdir -p "$run"
  bedtools bamtobed -i "$bam_dir/$sample.bam" > "$run/$sample.bed"
  args=(-t "$run/$sample.bed" -s "$genome" -w "$w" -g "$g" -f 150 \
        -egf 0.80 -rt 1 -e 1000 -cpu "$threads" -o "$run")
  stem="$run/$sample-W${w}-G${g}"
  if [[ -n "$control_dir" ]]; then
    bedtools bamtobed -i "$control_dir/input_${ct}.bam" > "$run/control.bed"
    args+=(-c "$run/control.bed" -fdr 0.05)
    peaks="$stem-FDR0.05-island.bed"
  else
    peaks="$stem.scoreisland"
  fi
  sicer "${args[@]}" > "$run/sicer.log" 2>&1
  [[ -f "$peaks" ]] || { echo "Missing output: $peaks" >&2; exit 1; }
  awk 'BEGIN{OFS="\t"} {print $1, $2, $3+1}' "$peaks" > "$run/peaks.bed"
done < cell_types.txt
