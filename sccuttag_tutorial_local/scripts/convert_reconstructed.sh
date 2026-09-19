#!/usr/bin/env bash
set -euo pipefail
out=results/reconstructed/bam; mkdir -p "$out"
while IFS= read -r ct; do
  sample="H3K4me3_${ct}"
  sam="results/reconstructed/sam/$sample.sam"
  [[ -s "$sam" ]]
  samtools sort -@ 8 -o "$out/$sample.bam" "$sam"
  samtools quickcheck -v "$out/$sample.bam"; samtools index "$out/$sample.bam"
  (( $(samtools view -c -f 64 "$out/$sample.bam") > 0 ))
done < cell_types.txt
