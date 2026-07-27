#!/usr/bin/env Rscript
# Negative-binomial test for chromatin retention, per (subfamily, orientation).
#
# 06_enrich_chrrna.py stops at effect sizes on purpose. This supplies the p-value.
# 2 vs 2 replicates over a few hundred count features is exactly what DESeq2's
# dispersion shrinkage is for; a Poisson or chi-square test on pooled counts would
# treat between-replicate variance as zero and manufacture significance.
#
# Sense and antisense are tested in ONE model with orientation as a covariate is
# NOT what happens here -- they are fitted separately, because the antisense
# compartment has a different baseline (see the mintseq tree, where the antisense
# background median sat 1.4 log2 units below sense). Pooling them would let the
# sense baseline set the dispersion prior for antisense.
#
# DESeq2 is not in deps/ by default:  cd deps && pixi add bioconductor-deseq2
#
# Usage: Rscript 07_deseq2.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
res_dir <- if (length(args) >= 1) args[1] else "results/chrrna"
infile <- file.path(res_dir, "subfamily_counts.tsv")

if (!file.exists(infile)) {
  stop("missing ", infile, " -- run 06_enrich_chrrna.py first")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 not installed. Run:  cd deps && pixi add bioconductor-deseq2")
}
suppressPackageStartupMessages(library(DESeq2))

tab <- read.delim(infile, check.names = FALSE)
count_cols <- c("K562_chr_Takara_rep1", "K562_chr_Takara_rep2",
                "K562_cyto_Takara_rep1", "K562_cyto_Takara_rep2")
missing <- setdiff(count_cols, colnames(tab))
if (length(missing)) stop("missing columns: ", paste(missing, collapse = ", "))

coldata <- data.frame(
  fraction = factor(c("chromatin", "chromatin", "cytoplasm", "cytoplasm"),
                    levels = c("cytoplasm", "chromatin")),
  row.names = count_cols
)

out <- list()
for (orient in c("sense", "antisense")) {
  sub <- tab[tab$orientation == orient, , drop = FALSE]
  mat <- as.matrix(sub[, count_cols])
  rownames(mat) <- sub$subfamily
  mode(mat) <- "integer"

  # Drop categories with essentially no signal; they only inflate the BH
  # denominator. 10 reads summed across the 4 libraries is a floor, not a filter
  # on effect size.
  keep <- rowSums(mat) >= 10
  cat(sprintf("%s: %d of %d subfamilies pass the 10-read floor\n",
              orient, sum(keep), nrow(mat)))
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2) next

  dds <- DESeqDataSetFromMatrix(mat, coldata, ~ fraction)
  dds <- DESeq(dds, quiet = TRUE)
  r <- results(dds, contrast = c("fraction", "chromatin", "cytoplasm"))
  df <- as.data.frame(r)
  df$subfamily <- rownames(df)
  df$orientation <- orient
  df$family <- sub$family[match(df$subfamily, sub$subfamily)]
  out[[orient]] <- df
}

res <- do.call(rbind, out)
res <- res[order(-res$log2FoldChange), ]
cols <- c("subfamily", "orientation", "family", "baseMean",
          "log2FoldChange", "lfcSE", "pvalue", "padj")
res <- res[, cols]

outfile <- file.path(res_dir, "chromatin_retention_deseq2.tsv")
write.table(res, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote ", outfile, "\n", sep = "")

sig <- subset(res, !is.na(padj) & padj < 0.05 & log2FoldChange > 0)
cat(sprintf("\n%d categories chromatin-retained at padj < 0.05\n", nrow(sig)))
print(utils::head(sig, 15), row.names = FALSE)
