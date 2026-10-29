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
# Other families: pass the family as the second argument (Rscript
# 18_gene_alu_gsea.R results L1) to read gene_l1_counts.tsv from 17 run with
# FAMILY=L1; outputs are named gene_l1_*. The lineage-mix covariates are the
# fractions of every n_<lineage> column 17 wrote except the most common one (for
# Alu: AluJ and AluY, AluS the reference). If 17 also wrote
# gene_<fam>_group_counts.tsv, each lineage is scored on its own with the same
# model, into gene_<fam>_group_retention.tsv, for 19 to test with expression
# matching.
#
# Usage: Rscript 18_gene_alu_gsea.R [results_dir] [family]

args <- commandArgs(trailingOnly = TRUE)
# Project root = the parent of src/, wherever this is run from.
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
res_dir <- if (length(args) >= 1) args[1] else file.path(root, "results")
family <- if (length(args) >= 2) args[2] else "Alu"
te <- tolower(family)
suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
})
set.seed(1)

tab0 <- read.delim(file.path(res_dir, sprintf("gene_%s_counts.tsv", te)), check.names = FALSE)
libs <- c("K562_chr_Takara_rep1", "K562_chr_Takara_rep2",
          "K562_cyto_Takara_rep1", "K562_cyto_Takara_rep2")
sf <- read.delim(file.path(res_dir, "size_factors.tsv"), row.names = 1)[libs, "sf_gene"]

# Floors, not effect filters: enough TE and intron reads for a per-gene ratio,
# and a gene actually expressed in K562.
passes <- function(tcnt, icnt, ecnt, ids) {
  rowSums(tcnt) >= 20 & rowSums(icnt) >= 50 & rowSums(ecnt) >= 40 &
    rowSums(tcnt[, 3:4]) + rowSums(icnt[, 3:4]) > 0 &
    !grepl("_PAR_Y", ids)    # chrY PAR copies share the X gene's Ensembl id
}

# One DESeq2 fit of TE vs own-intron retention; returns the aluChr results.
fit_score <- function(tcnt, icnt, ids) {
  mat <- cbind(tcnt, icnt)
  colnames(mat) <- c(paste0(te, ".", libs), paste0("intron.", libs))
  rownames(mat) <- sub("\\..*$", "", ids)
  mode(mat) <- "integer"
  cd <- data.frame(
    lib = factor(rep(libs, 2)),
    type = factor(rep(c(te, "intron"), each = 4), levels = c("intron", te)),
    row.names = colnames(mat))
  cd$aluChr <- as.numeric(cd$type == te & grepl("_chr_", cd$lib))
  dds <- DESeqDataSetFromMatrix(mat, cd, ~ lib + type + aluChr)
  sizeFactors(dds) <- rep(sf, 2)
  dds <- DESeq(dds, quiet = TRUE)
  results(dds, name = "aluChr")
}

expr0 <- as.matrix(tab0[, paste0("expr_", libs)])
intr0 <- as.matrix(tab0[, paste0("intron_", libs)])
keep <- passes(as.matrix(tab0[, paste0(te, "_", libs)]), intr0, expr0, tab0$gene_id)
cat(sprintf("%d genes with intron_same %s; %d pass the read floors\n", nrow(tab0), family, sum(keep)))
tab <- tab0[keep, ]
alu <- as.matrix(tab[, paste0(te, "_", libs)])
intr <- intr0[keep, ]
expr <- expr0[keep, ]

r <- fit_score(alu, intr, tab$gene_id)
cat(sprintf("score (%s minus intron retention): median %.2f, %d genes padj<0.05 (%d up, %d down)\n",
            family, median(r$log2FoldChange, na.rm = TRUE), sum(r$padj < 0.05, na.rm = TRUE),
            sum(r$padj < 0.05 & r$log2FoldChange > 0, na.rm = TRUE),
            sum(r$padj < 0.05 & r$log2FoldChange < 0, na.rm = TRUE)))

# Covariates and the adjusted statistic.
n_te <- tab[[paste0("n_", te)]]
grp_cols <- setdiff(grep("^n_", names(tab), value = TRUE), paste0("n_", te))
ref_col <- grp_cols[which.max(colSums(tab[, grp_cols, drop = FALSE]))]
nexpr <- sweep(expr, 2, sf, "/")
cov <- data.frame(log_gene_len = log10(tab$gene_len),
                  log_intron_bp = log10(tab$intron_bp + 1))
cov[[paste0("log_n_", te)]] <- log10(n_te)
for (g in setdiff(grp_cols, ref_col)) cov[[sub("^n_", "frac_", g)]] <- tab[[g]] / n_te
cov$log_expr_chr <- log10(rowMeans(nexpr[, 1:2]) + 1)
cov$log_expr_cyto <- log10(rowMeans(nexpr[, 3:4]) + 1)
cov[[paste0("log_", te, "_depth")]] <- log10(rowSums(alu) + 1)
ok <- !is.na(r$stat)
stat <- r$stat[ok]
names(stat) <- rownames(r)[ok]
cv <- cov[ok, ]
cat("\nSpearman of the statistic with covariates:\n")
for (v in colnames(cv)) cat(sprintf("  %-14s %+.3f\n", v, cor(stat, cv[[v]], method = "spearman")))
# Splines on the continuous covariates, linear in the lineage fractions.
terms <- ifelse(grepl("^frac_", colnames(cv)), colnames(cv),
                sprintf("splines::ns(%s, 3)", colnames(cv)))
fit <- lm(as.formula(paste("stat ~", paste(terms, collapse = " + "))), data = cv)
cat(sprintf("covariates explain R^2 = %.3f of the statistic\n", summary(fit)$r.squared))
stat_adj <- residuals(fit)
names(stat_adj) <- names(stat)

genes_out <- data.frame(gene_id = tab$gene_id, gene_name = tab$gene_name,
                        n_te = n_te, baseMean = r$baseMean,
                        score_log2 = r$log2FoldChange, lfcSE = r$lfcSE,
                        stat = r$stat, padj = r$padj)
names(genes_out)[3] <- paste0("n_", te)
genes_out$stat_adj <- NA
genes_out$stat_adj[ok] <- stat_adj
write.table(genes_out[order(-genes_out$stat), ], file.path(res_dir, sprintf("gene_%s_retention.tsv", te)),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Per-lineage scores, same model and floors, each lineage's own reads only.
gfile <- file.path(res_dir, sprintf("gene_%s_group_counts.tsv", te))
if (file.exists(gfile)) {
  gcnt <- read.delim(gfile, check.names = FALSE)
  cat("\nper-lineage scores:\n")
  gres <- list()
  for (g in sort(unique(gcnt$group))) {
    x <- gcnt[gcnt$group == g, ]
    m <- match(x$gene_id, tab0$gene_id)
    tc <- as.matrix(x[, paste0(te, "_", libs)])
    k <- passes(tc, intr0[m, , drop = FALSE], expr0[m, , drop = FALSE], x$gene_id)
    if (sum(k) < 100) {
      cat(sprintf("  %-6s %5d genes pass the floors -- too few, skipped\n", g, sum(k)))
      next
    }
    rg <- fit_score(tc[k, ], intr0[m[k], ], x$gene_id[k])
    cat(sprintf("  %-6s %5d genes  median score %+.2f  %d padj<0.05 (%d up, %d down)\n", g, sum(k),
                median(rg$log2FoldChange, na.rm = TRUE), sum(rg$padj < 0.05, na.rm = TRUE),
                sum(rg$padj < 0.05 & rg$log2FoldChange > 0, na.rm = TRUE),
                sum(rg$padj < 0.05 & rg$log2FoldChange < 0, na.rm = TRUE)))
    gres[[g]] <- data.frame(group = g, gene_id = x$gene_id[k], gene_name = tab0$gene_name[m[k]],
                            baseMean = rg$baseMean, score_log2 = rg$log2FoldChange,
                            lfcSE = rg$lfcSE, stat = rg$stat, padj = rg$padj)
  }
  write.table(do.call(rbind, gres), file.path(res_dir, sprintf("gene_%s_group_retention.tsv", te)),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

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
            file.path(res_dir, sprintf("gene_%s_gsea.tsv", te)), sep = "\t", quote = FALSE, row.names = FALSE)

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
