library(Seurat)
object <- readRDS("data/H3K4me3.rds")
label <- "predicted.celltype.l1"
meta <- object[[]]
stopifnot(label %in% colnames(meta))
barcodes <- if ("bam_barcode" %in% colnames(meta)) as.character(meta$bam_barcode) else rownames(meta)
labels <- as.character(meta[[label]])
labels[labels == "CD4 T"] <- "CD4T"
labels[labels == "CD8 T"] <- "CD8T"
labels[labels == "other T"] <- "otherT"
keep <- !is.na(labels) & nzchar(labels) & labels != "Unknown"
cells <- data.frame(barcode = barcodes[keep], cell_type = labels[keep])
stopifnot(nrow(cells) > 0, !anyNA(cells$barcode), !anyDuplicated(cells$barcode))
stopifnot(all(grepl("^[A-Za-z0-9_-]+$", cells$cell_type)))
dir.create("barcodes", showWarnings = FALSE)
write.table(cells, "cells.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
types <- sort(unique(cells$cell_type))
writeLines(types, "cell_types.txt")
for (ct in types) {
  writeLines(cells$barcode[cells$cell_type == ct],
             file.path("barcodes", paste0(ct, "_barcodes.txt")))
}
write.table(as.data.frame(table(cells$cell_type)),
            "cell_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)