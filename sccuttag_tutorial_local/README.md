# scCUT&Tag peak-calling tutorial

This is the project for the scCUT&Tag peak-calling section. Once the environment is properly configured, you can quickly complete the entire research example by following the workflow below.


```bash
Rscript scripts/export_barcodes.R
cat cell_counts.tsv
Rscript scripts/reconstruct_fragments.R
bash scripts/convert_reconstructed.sh

export MARK=H3K4me3 MODE=narrow GENOME=hg38 GENOME_SIZE=2.7e9 THREADS=8
bash scripts/call_macs2.sh results/reconstructed/bam results/reconstructed/MACS2/narrow/no_control
bash scripts/call_sicer2.sh results/reconstructed/bam results/reconstructed/SICER2/narrow/no_control
bash scripts/make_pseudo_inputs.sh results/reconstructed/bam results/reconstructed/pseudo_input
bash scripts/call_macs2.sh results/reconstructed/bam results/reconstructed/MACS2/narrow/pseudo_input results/reconstructed/pseudo_input
bash scripts/call_sicer2.sh results/reconstructed/bam results/reconstructed/SICER2/narrow/pseudo_input results/reconstructed/pseudo_input
REPRESENTATIONS=reconstructed Rscript scripts/summarize_peaks.R
```

Inspect `peak_summary.csv` and caller logs under `results/reconstructed/`.
