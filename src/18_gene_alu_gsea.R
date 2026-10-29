#!/usr/bin/env Rscript
# Host-gene function test for chromatin-retained intronic Alu RNA.
#
# Per gene (counts from 17_gene_alu.py, unique fragments, Takara libraries):
#
#   score = log2(chromatin/cytoplasm) of the gene's intronic Alu RNA
#         - log2(chromatin/cytoplasm) of the gene's own non-repeat intronic RNA
#
# fitted in ONE DESeq2 model over 8 columns (4 libraries x {alu, intron}):
#
#   ~ lib + type + aluChr        aluChr = 1 for the alu column of a chromatin library
#
# 'lib' absorbs everything library-specific (depth, rRNA, normalisation), so the
# aluChr coefficient compares Alu to intron WITHIN each library: it needs no
# size factors at all, and is immune to the TE-vs-gene normalisation problem that
# 11 had to solve. Positive = the Alu RNA is held on chromatin beyond the pre-mRNA
# it sits in; negative = it leaves chromatin more readily than its introns.
#
# Genes are ranked by the Wald statistic and tested with fgsea against MSigDB
# Hallmark, GO:BP and Reactome. The universe is the tested genes themselves
# (expressed genes carrying unique-read-supported intronic Alus), so the generic
# biases of where Alus sit -- long, GC-rich, highly expressed genes -- are already
# in the background. The statistic is also checked against gene length, intron
# bp, Alu number, Alu age mix, expression and Alu depth; GSEA is run on the raw
# statistic AND on its residual after regressing those out, and only themes that
# survive the adjustment should be trusted.
#
# Usage: Rscript 18_gene_alu_gsea.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
# Project root = the parent of src/, wherever this is run from.
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
res_dir <- if (length(args) >= 1) args[1] else file.path(root, "results")
suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
})
set.seed(1)

tab <- read.delim(file.path(res_dir, "gene_alu_counts.tsv"), check.names = FALSE)
libs <- c("K562_chr_Takara_rep1", "K562_chr_Takara_rep2",
          "K562_cyto_Takara_rep1", "K562_cyto_Takara_rep2")
sf <- read.delim(file.path(res_dir, "size_factors.tsv"), row.names = 1)[libs, "sf_gene"]

alu <- as.matrix(tab[, paste0("alu_", libs)])
intr <- as.matrix(tab[, paste0("intron_", libs)])
expr <- as.matrix(tab[, paste0("expr_", libs)])
# Floors, not effect filters: enough Alu and intron reads for a per-gene ratio,
# and a gene actually expressed in K562.
keep <- rowSums(alu) >= 20 & rowSums(intr) >= 50 & rowSums(expr) >= 40 &
        rowSums(alu[, 3:4]) + rowSums(intr[, 3:4]) > 0 &
        !grepl("_PAR_Y", tab$gene_id)    # chrY PAR copies share the X gene's Ensembl id
cat(sprintf("%d genes with intron_same Alus; %d pass the read floors\n", nrow(tab), sum(keep)))
tab <- tab[keep, ]
alu <- alu[keep, ]
intr <- intr[keep, ]
expr <- expr[keep, ]

mat <- cbind(alu, intr)
colnames(mat) <- c(paste0("alu.", libs), paste0("intron.", libs))
rownames(mat) <- sub("\\..*$", "", tab$gene_id)
mode(mat) <- "integer"
cd <- data.frame(
  lib = factor(rep(libs, 2)),
  type = factor(rep(c("alu", "intron"), each = 4), levels = c("intron", "alu")),
  row.names = colnames(mat))
cd$aluChr <- as.numeric(cd$type == "alu" & grepl("_chr_", cd$lib))

dds <- DESeqDataSetFromMatrix(mat, cd, ~ lib + type + aluChr)
sizeFactors(dds) <- rep(sf, 2)
dds <- DESeq(dds, quiet = TRUE)
r <- results(dds, name = "aluChr")
cat(sprintf("score (Alu minus intron retention): median %.2f, %d genes padj<0.05 (%d up, %d down)\n",
            median(r$log2FoldChange, na.rm = TRUE), sum(r$padj < 0.05, na.rm = TRUE),
            sum(r$padj < 0.05 & r$log2FoldChange > 0, na.rm = TRUE),
            sum(r$padj < 0.05 & r$log2FoldChange < 0, na.rm = TRUE)))

# Covariates and the adjusted statistic.
nexpr <- sweep(expr, 2, sf, "/")
cov <- data.frame(
  log_gene_len = log10(tab$gene_len),
  log_intron_bp = log10(tab$intron_bp + 1),
  log_n_alu = log10(tab$n_alu),
  frac_aluJ = tab$n_aluJ / tab$n_alu,
  frac_aluY = tab$n_aluY / tab$n_alu,
  log_expr_chr = log10(rowMeans(nexpr[, 1:2]) + 1),
  log_expr_cyto = log10(rowMeans(nexpr[, 3:4]) + 1),
  log_alu_depth = log10(rowSums(alu) + 1))
ok <- !is.na(r$stat)
stat <- r$stat[ok]
names(stat) <- rownames(r)[ok]
cv <- cov[ok, ]
cat("\nSpearman of the statistic with covariates:\n")
for (v in colnames(cv)) cat(sprintf("  %-14s %+.3f\n", v, cor(stat, cv[[v]], method = "spearman")))
fit <- lm(stat ~ splines::ns(log_gene_len, 3) + splines::ns(log_intron_bp, 3) +
            splines::ns(log_n_alu, 3) + frac_aluJ + frac_aluY +
            splines::ns(log_expr_chr, 3) + splines::ns(log_expr_cyto, 3) +
            splines::ns(log_alu_depth, 3), data = cv)
cat(sprintf("covariates explain R^2 = %.3f of the statistic\n", summary(fit)$r.squared))
stat_adj <- residuals(fit)
names(stat_adj) <- names(stat)

genes_out <- data.frame(gene_id = tab$gene_id, gene_name = tab$gene_name,
                        n_alu = tab$n_alu, baseMean = r$baseMean,
                        score_log2 = r$log2FoldChange, lfcSE = r$lfcSE,
                        stat = r$stat, padj = r$padj)
genes_out$stat_adj <- NA
genes_out$stat_adj[ok] <- stat_adj
write.table(genes_out[order(-genes_out$stat), ], file.path(res_dir, "gene_alu_retention.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Gene sets.
get_sets <- function(coll, sub = NULL) {
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub),
                error = function(e) msigdbr(species = "Homo sapiens", category = coll, subcategory = sub))
  split(m$ensembl_gene, m$gs_name)
}
sets <- list(H = get_sets("H"), GOBP = get_sets("C5", "GO:BP"),
             REACTOME = get_sets("C2", "CP:REACTOME"))
sym <- setNames(tab$gene_name, sub("\\..*$", "", tab$gene_id))

all_res <- list()
for (which in c("raw", "adjusted")) {
  s <- if (which == "raw") stat else stat_adj
  for (coll in names(sets)) {
    g <- fgsea(sets[[coll]], s, minSize = 15, maxSize = 500, eps = 0)
    if (!nrow(g)) next
    g$collection <- coll
    g$statistic <- which
    g$leadingEdge <- vapply(g$leadingEdge, function(x) paste(head(sym[x], 15), collapse = ","), "")
    all_res[[length(all_res) + 1]] <- as.data.frame(g)
  }
}
gs <- do.call(rbind, all_res)
gs <- gs[order(gs$statistic, gs$collection, gs$padj), ]
write.table(gs[, c("statistic", "collection", "pathway", "size", "NES", "pval", "padj", "leadingEdge")],
            file.path(res_dir, "gene_alu_gsea.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nuniverse: %d genes\n", length(stat)))
for (which in c("raw", "adjusted")) {
  for (coll in names(sets)) {
    x <- gs[gs$statistic == which & gs$collection == coll, ]
    cat(sprintf("\n== %s / %s: %d sets tested, %d padj<0.05, %d padj<0.25\n", which, coll,
                nrow(x), sum(x$padj < 0.05), sum(x$padj < 0.25)))
    x <- head(x[x$padj < 0.25, ], 12)
    if (nrow(x)) for (i in seq_len(nrow(x)))
      cat(sprintf("  %-60s n=%3d NES=%+5.2f padj=%.2g  %s\n", substr(x$pathway[i], 1, 60), x$size[i],
                  x$NES[i], x$padj[i], substr(x$leadingEdge[i], 1, 60)))
  }
}
