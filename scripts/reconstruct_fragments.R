cells <- read.delim("cells.tsv", colClasses = "character")
sizes <- read.delim("data/hg38.chrom.sizes", header = FALSE,
                    col.names = c("chrom", "length"))
stopifnot(!anyDuplicated(cells$barcode), !anyDuplicated(sizes$chrom),
          all(sizes$length > 0))
cell_type <- setNames(cells$cell_type, cells$barcode)
chrom_limit <- setNames(sizes$length, sizes$chrom)
types <- sort(unique(cells$cell_type))
out_dir <- "results/reconstructed/sam"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
header <- c("@HD\tVN:1.6\tSO:unsorted",
            paste0("@SQ\tSN:", sizes$chrom, "\tLN:", sizes$length))
write_pairs <- function(fr, first_id, con) {
  width <- fr$end - fr$start; rl <- pmin(50L, width)
  left <- fr$start + 1L; right <- fr$end - rl + 1L
  ids <- paste0("fragment_", first_id + seq_len(nrow(fr)))
  cigar <- paste0(rl, "M"); cb <- paste0("CB:Z:", fr$barcode)
  r1 <- paste(ids, 99L, fr$chrom, left, 255L, cigar, "=", right, width, "*", "*", cb, sep = "\t")
  r2 <- paste(ids, 147L, fr$chrom, right, 255L, cigar, "=", left, -width, "*", "*", cb, sep = "\t")
  writeLines(as.vector(rbind(r1, r2)), con)
}
paths <- file.path(out_dir, paste0("H3K4me3_", types, ".sam"))
outs <- setNames(lapply(paths, file, open = "wt"), types)
on.exit(invisible(lapply(outs, close)), add = TRUE)
for (con in outs) writeLines(header, con)
input <- gzfile("data/H3K4me3_fragments.tsv.gz", "rt")
on.exit(close(input), add = TRUE)
written <- setNames(numeric(length(types)), types); next_id <- 0
repeat {
  lines <- readLines(input, n = 100000L, warn = FALSE)
  if (!length(lines)) break
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (!length(lines)) next
  x <- read.table(text = lines, sep = "\t", header = FALSE,
                  colClasses = "character", quote = "", comment.char = "")
  x <- x[x$V4 %in% cells$barcode, 1:4]
  if (!nrow(x)) next
  names(x) <- c("chrom", "start", "end", "barcode")
  x$start <- as.integer(x$start); x$end <- as.integer(x$end)
  limit <- unname(chrom_limit[x$chrom])
  stopifnot(!anyNA(x), !anyNA(limit), all(x$start >= 0),
            all(x$end > x$start), all(x$end <= limit))
  for (ct in types) {
    part <- x[cell_type[x$barcode] == ct, , drop = FALSE]
    if (!nrow(part)) next
    write_pairs(part, next_id, outs[[ct]])
    next_id <- next_id + nrow(part); written[ct] <- written[ct] + nrow(part)
  }
}
stopifnot(all(written > 0))
write.table(data.frame(cell_type = types, fragments = written),
            file.path(out_dir, "counts.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
