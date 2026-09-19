#!/usr/bin/env bash
# Usage: bash scripts/split_bam.sh [input.bam] [out_dir]
set -euo pipefail
bam=${1:-data/H3K4me3.bam}
out=${2:-results/aligned/bam}
mark=${MARK:-H3K4me3}
threads=${THREADS:-8}
mkdir -p "$out"
samtools quickcheck -v "$bam"
printf 'cell_type\tretained_pairs\n' > "$out/counts.tsv"
while IFS= read -r ct; do
  list="barcodes/${ct}_barcodes.txt"; outbam="$out/${mark}_${ct}.bam"
  [[ -s "$list" ]]
  samtools view -u -f 2 -F 3844 -D "CB:$list" "$bam" |
    samtools sort -@ "$threads" -o "$outbam" -
  samtools quickcheck -v "$outbam"; samtools index "$outbam"
  pairs=$(samtools view -c -f 64 "$outbam")
  (( pairs > 0 )) || { echo "No retained pairs for $ct" >&2; exit 1; }
  printf '%s\t%s\n' "$ct" "$pairs" >> "$out/counts.tsv"
done < cell_types.txt
