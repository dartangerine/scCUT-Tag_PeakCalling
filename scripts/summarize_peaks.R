reps <- strsplit(Sys.getenv("REPRESENTATIONS", "aligned"), ",", fixed = TRUE)[[1]]
types <- readLines("cell_types.txt")
manifest <- expand.grid(cell_type = types, caller = c("MACS2", "SICER2"),
                        control = c("no_control", "pseudo_input"),
                        representation = reps, stringsAsFactors = FALSE)
manifest$path <- apply(manifest, 1, function(r) {
  base <- file.path("results", r["representation"], r["caller"],
                    "narrow", r["control"])
  s <- paste0("H3K4me3_", r["cell_type"])
  if (r["caller"] == "MACS2") file.path(base, paste0(s, "_peaks.narrowPeak"))
  else file.path(base, s, "peaks.bed")
})
summ <- function(path) {
  if (!file.exists(path)) stop("Missing peak file: ", path)
  if (file.size(path) == 0) return(c(n_peaks = 0, mean_width_bp = NA_real_))
  x <- read.table(path, sep = "\t", comment.char = "#", quote = "")
  stopifnot(ncol(x) >= 3, all(x$V3 > x$V2))
  c(n_peaks = nrow(x), mean_width_bp = mean(x$V3 - x$V2))
}
stats <- t(vapply(manifest$path, summ, numeric(2)))
write.csv(cbind(manifest, stats), "peak_summary.csv", row.names = FALSE)
