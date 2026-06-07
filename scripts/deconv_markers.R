#!/usr/bin/env Rscript
# deconv_markers.R -- base-R marker z-score deconvolution (Arm 5, R side).
#
#   Rscript --vanilla deconv_markers.R <tpm_in.tsv> <scores_out.tsv>
#
# Input : a gene_name x sample TPM matrix (TSV; first column = gene_name).
# Output: a sample x cell_type score matrix (TSV; first column = sample).
#
# Method (deliberately DIFFERENT from the Python ssGSEA mean-rank, so the
# cross-method agreement check is meaningful): per gene, z-score log2(TPM+1)
# across samples; per cell type, a sample's score is the mean z of its present
# marker genes. Base R only -- no CRAN/Bioconductor packages -- so it runs
# wherever R is installed.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: deconv_markers.R <tpm_in.tsv> <scores_out.tsv>")
in_tsv <- args[[1]]; out_tsv <- args[[2]]

# Marker panels -- must match hnscc_time.genomics.IMMUNE_SIGNATURES (same cell
# types / genes) so the two methods are comparable.
markers <- list(
  CD3   = c("CD3D", "CD3E", "CD3G"),
  CD8   = c("CD8A", "CD8B", "GZMB", "PRF1"),
  FoxP3 = c("FOXP3", "IL2RA", "CTLA4"),
  PanCK = c("KRT5", "KRT6A", "KRT14", "KRT17")
)

# Read explicitly: column 1 = gene_name, remaining columns = samples. (Avoids
# row.names=1 header-raggedness ambiguity across writers.)
raw <- read.delim(in_tsv, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
genes <- as.character(raw[[1]])
m <- as.matrix(raw[, -1, drop = FALSE])
storage.mode(m) <- "double"
rownames(m) <- genes

# log2(TPM+1) then per-gene (row) z-score across samples.
logm <- log2(m + 1)
mu <- rowMeans(logm)
sdv <- apply(logm, 1L, sd)
z <- sweep(logm, 1L, mu, "-")
z <- sweep(z, 1L, ifelse(sdv == 0, 1, sdv), "/")
z[sdv == 0, ] <- 0  # genes with no variance carry no signal

samples <- colnames(m)
out <- matrix(0.0, nrow = length(samples), ncol = length(markers),
              dimnames = list(samples, names(markers)))
for (cell in names(markers)) {
  present <- intersect(markers[[cell]], rownames(z))
  if (length(present) > 0L) {
    out[, cell] <- colMeans(z[present, , drop = FALSE])
  }
}

df <- data.frame(sample = samples, out, check.names = FALSE, row.names = NULL)
write.table(df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
