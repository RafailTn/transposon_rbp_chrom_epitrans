#!/usr/bin/env Rscript
# Chromatin retention per (subfamily, orientation, context), under both
# normalisations from 10_context.py.
#
# Same model as 07_deseq2.R (~ fraction, Takara chromatin vs cytoplasm, 10-read
# floor, sense and antisense fitted separately), with two changes:
#
#   * fitted separately per genomic context as well (utr3/exon/intron x same/opp,
#     intergenic, plus 'all' = contexts pooled), so a context's baseline cannot
#     set the dispersion prior for another -- the same reason 07 splits by
#     orientation;
#   * size factors are FIXED from size_factors.tsv rather than re-estimated per
#     subset. Re-estimating inside a subset would make every context
#     compositional on its own and erase exactly the between-context shift
#     this split exists to show. Both 'te' and 'gene' factors are reported.
#
# For intron_* contexts, delta_vs_intron = log2FC - log2(chromatin/cytoplasm of
# all host-intron fragments, under the same size factors): positive means the TE
# is retained MORE than the introns around it, i.e. beyond what host pre-mRNA
# would explain.
#
# Usage: Rscript 11_context_deseq2.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
# Project root = the parent of src/, wherever this is run from.
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
res_dir <- if (length(args) >= 1) args[1] else file.path(root, "results")
suppressPackageStartupMessages(library(DESeq2))

tab <- read.delim(file.path(res_dir, "context_counts.tsv"), check.names = FALSE)
sf <- read.delim(file.path(res_dir, "size_factors.tsv"), row.names = 1)
count_cols <- c("K562_chr_Takara_rep1", "K562_chr_Takara_rep2",
                "K562_cyto_Takara_rep1", "K562_cyto_Takara_rep2")
fraction <- factor(c("chromatin", "chromatin", "cytoplasm", "cytoplasm"),
                   levels = c("cytoplasm", "chromatin"))
coldata <- data.frame(fraction = fraction, row.names = count_cols)

pooled <- aggregate(tab[, count_cols], tab[, c("subfamily", "orientation", "family")], sum)
pooled$context <- "all"
tab <- rbind(tab[, c("subfamily", "orientation", "context", "family", count_cols)],
             pooled[, c("subfamily", "orientation", "context", "family", count_cols)])

out <- list()
for (norm in c("te", "gene")) {
  s <- sf[count_cols, paste0("sf_", norm)]
  ni <- sf[count_cols, "intron"] / s
  intron_lfc <- log2(mean(ni[1:2]) / mean(ni[3:4]))
  cat(sprintf("[%s] host-intron log2(chromatin/cytoplasm) = %.3f\n", norm, intron_lfc))
  for (orient in c("sense", "antisense")) {
    for (ctx in unique(tab$context)) {
      sub <- tab[tab$orientation == orient & tab$context == ctx, , drop = FALSE]
      mat <- as.matrix(sub[, count_cols])
      rownames(mat) <- sub$subfamily
      mode(mat) <- "integer"
      keep <- rowSums(mat) >= 10
      if (sum(keep) < 2) next
      mat <- mat[keep, , drop = FALSE]
      dds <- DESeqDataSetFromMatrix(mat, coldata, ~ fraction)
      sizeFactors(dds) <- s
      dds <- DESeq(dds, quiet = TRUE, fitType = if (nrow(mat) < 30) "mean" else "parametric")
      df <- as.data.frame(results(dds, contrast = c("fraction", "chromatin", "cytoplasm")))
      df$subfamily <- rownames(df)
      df$family <- sub$family[match(df$subfamily, sub$subfamily)]
      df$orientation <- orient
      df$context <- ctx
      df$norm <- norm
      df$delta_vs_intron <- if (startsWith(ctx, "intron")) df$log2FoldChange - intron_lfc else NA
      out[[length(out) + 1]] <- df
    }
  }
}

res <- do.call(rbind, out)
cols <- c("norm", "context", "orientation", "subfamily", "family", "baseMean",
          "log2FoldChange", "lfcSE", "pvalue", "padj", "delta_vs_intron")
res <- res[order(res$norm, res$context, res$orientation, -res$log2FoldChange), cols]
outfile <- file.path(res_dir, "chromatin_retention_context.tsv")
write.table(res, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote ", outfile, "\n\n", sep = "")

# Family-level summary: median log2FC and significant up/down, per cell.
sig <- !is.na(res$padj) & res$padj < 0.05
cells <- split(seq_len(nrow(res)), res[, c("norm", "context", "orientation", "family")], drop = TRUE)
summ <- do.call(rbind, lapply(cells, function(i) {
  r <- res[i[1], ]
  data.frame(norm = r$norm, family = r$family, orientation = r$orientation,
             context = r$context, n = length(i),
             up = sum(sig[i] & res$log2FoldChange[i] > 0),
             down = sum(sig[i] & res$log2FoldChange[i] < 0),
             med_lfc = round(median(res$log2FoldChange[i]), 2),
             med_delta_intron = round(median(res$delta_vs_intron[i]), 2))
}))
summ <- summ[order(summ$norm, summ$family, summ$orientation, summ$context), ]
write.table(summ, file.path(res_dir, "chromatin_retention_context_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)
