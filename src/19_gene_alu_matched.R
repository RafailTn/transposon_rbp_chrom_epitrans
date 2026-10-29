#!/usr/bin/env Rscript
# Expression-matched version of the host-gene function test (18_gene_alu_gsea.R).
#
# In 18 the per-gene score (intronic Alu retention minus own-intron retention)
# correlated +0.62 with cytoplasmic expression, and the raw GSEA was dominated by
# the most highly expressed gene classes. Regressing covariates out (18's
# 'adjusted' run) can over- or under-correct at the extremes -- it turned the
# hugely expressed histone genes into strong negative outliers. Two
# non-parametric alternatives here, comparing genes only with genes of similar
# expression:
#
#   1. STRATIFIED GSEA. Genes are binned into deciles of cytoplasmic expression
#      (the strongest covariate); the Wald statistic is turned into a normal score
#      of its rank WITHIN the bin, and fgsea runs on that. A second run bins on
#      expression (chromatin x cytoplasm, 5 x 5) crossed with gene-length
#      tertiles. A gene scores high only by beating its expression peers.
#
#   2. MATCHED CONTROLS for named candidate sets (top raw hits from 18, plus
#      TGF-beta and OXPHOS): each member gets its K nearest non-members in
#      (log chromatin expr, log cytoplasm expr, log gene length), standardised;
#      the set's mean statistic is compared with the mean over matched controls,
#      and an empirical p-value comes from redrawing one control per member
#      10,000 times.
#
# Writes results/gene_alu_gsea_matched.tsv and results/gene_alu_matched_sets.tsv.
#
# Other families: second argument, as in 18 (Rscript 19_gene_alu_matched.R
# results L1). The named candidates above came from the Alu run, so for any other
# family every Hallmark set is checked with matched controls as well. When 18
# wrote gene_<fam>_group_retention.tsv, each lineage also gets the stratified
# GSEA (z_joint, strata built within that lineage's genes), into
# gene_<fam>_group_gsea_matched.tsv.
#
# Usage: Rscript 19_gene_alu_matched.R [results_dir] [family]

args <- commandArgs(trailingOnly = TRUE)
# Project root = the parent of src/, wherever this is run from.
root <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
res_dir <- if (length(args) >= 1) args[1] else file.path(root, "results")
family <- if (length(args) >= 2) args[2] else "Alu"
te <- tolower(family)
suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
})
set.seed(1)

libs <- c("K562_chr_Takara_rep1", "K562_chr_Takara_rep2",
          "K562_cyto_Takara_rep1", "K562_cyto_Takara_rep2")
sf <- read.delim(file.path(res_dir, "size_factors.tsv"), row.names = 1)[libs, "sf_gene"]
cnt <- read.delim(file.path(res_dir, sprintf("gene_%s_counts.tsv", te)), check.names = FALSE)
ret <- read.delim(file.path(res_dir, sprintf("gene_%s_retention.tsv", te)))
ret <- ret[!is.na(ret$stat), ]
d <- merge(ret, cnt, by = c("gene_id", "gene_name", paste0("n_", te)))
nexpr <- sweep(as.matrix(d[, paste0("expr_", libs)]), 2, sf, "/")
d$x_chr <- log10(rowMeans(nexpr[, 1:2]) + 1)
d$x_cyto <- log10(rowMeans(nexpr[, 3:4]) + 1)
d$x_len <- log10(d$gene_len)
d$ens <- sub("\\..*$", "", d$gene_id)
cat(sprintf("%d genes\n", nrow(d)))

within_score <- function(stat, strata) {
  z <- numeric(length(stat))
  for (s in unique(strata)) {
    i <- which(strata == s)
    z[i] <- qnorm((rank(stat[i]) - 0.5) / length(i))
  }
  # Equal ranks in different strata give identical scores; break those ties by
  # the raw statistic so fgsea's ordering is not arbitrary.
  z + 1e-6 * rank(stat) / length(stat)
}
q <- function(x, n) cut(x, unique(quantile(x, seq(0, 1, length.out = n + 1))), include.lowest = TRUE)
d$strat_expr <- q(d$x_cyto, 10)
d$strat_joint <- interaction(q(d$x_chr, 5), q(d$x_cyto, 5), q(d$x_len, 3), drop = TRUE)
d$z_expr <- within_score(d$stat, d$strat_expr)
d$z_joint <- within_score(d$stat, d$strat_joint)

cat("\nSpearman with covariates (raw stat -> within-decile -> within joint strata):\n")
for (v in c("x_chr", "x_cyto", "x_len")) {
  cat(sprintf("  %-7s %+.3f  %+.3f  %+.3f\n", v,
              cor(d$stat, d[[v]], method = "spearman"),
              cor(d$z_expr, d[[v]], method = "spearman"),
              cor(d$z_joint, d[[v]], method = "spearman")))
}
cat(sprintf("joint strata: %d, median %d genes each\n",
            nlevels(d$strat_joint), as.integer(median(table(d$strat_joint)))))

get_sets <- function(coll, sub = NULL) {
  m <- tryCatch(msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub),
                error = function(e) msigdbr(species = "Homo sapiens", category = coll, subcategory = sub))
  split(m$ensembl_gene, m$gs_name)
}
sets <- list(H = get_sets("H"), GOBP = get_sets("C5", "GO:BP"),
             REACTOME = get_sets("C2", "CP:REACTOME"))
sym <- setNames(d$gene_name, d$ens)

all_res <- list()
for (which in c("z_expr", "z_joint")) {
  s <- setNames(d[[which]], d$ens)
  for (coll in names(sets)) {
    g <- fgsea(sets[[coll]], s, minSize = 15, maxSize = 500, eps = 0)
    g$collection <- coll
    g$statistic <- which
    g$leadingEdge <- vapply(g$leadingEdge, function(x) paste(head(sym[x], 15), collapse = ","), "")
    all_res[[length(all_res) + 1]] <- as.data.frame(g)
  }
}
gs <- do.call(rbind, all_res)
gs <- gs[order(gs$statistic, gs$collection, gs$padj), ]
write.table(gs[, c("statistic", "collection", "pathway", "size", "NES", "pval", "padj", "leadingEdge")],
            file.path(res_dir, sprintf("gene_%s_gsea_matched.tsv", te)), sep = "\t", quote = FALSE, row.names = FALSE)
for (which in c("z_expr", "z_joint")) {
  for (coll in names(sets)) {
    x <- gs[gs$statistic == which & gs$collection == coll & !is.na(gs$padj), ]
    cat(sprintf("\n== %s / %s: %d sets, %d padj<0.05, %d padj<0.25\n", which, coll,
                nrow(x), sum(x$padj < 0.05), sum(x$padj < 0.25)))
    x <- head(x[x$padj < 0.25, ], 10)
    if (nrow(x)) for (i in seq_len(nrow(x)))
      cat(sprintf("  %-58s n=%3d NES=%+5.2f padj=%.2g  %s\n", substr(x$pathway[i], 1, 58), x$size[i],
                  x$NES[i], x$padj[i], substr(x$leadingEdge[i], 1, 55)))
  }
}

# Matched controls for named candidates.
cand <- c("HALLMARK_MYC_TARGETS_V1", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
          "HALLMARK_MTORC1_SIGNALING", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
          "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_TGF_BETA_SIGNALING",
          "HALLMARK_UV_RESPONSE_DN", "GOBP_TRANSLATION", "GOBP_RNA_SPLICING",
          "GOBP_RIBOSOME_BIOGENESIS", "REACTOME_MRNA_SPLICING")
pool <- c(sets$H, sets$GOBP, sets$REACTOME)
if (family != "Alu") cand <- union(cand, names(sets$H))
X <- scale(as.matrix(d[, c("x_chr", "x_cyto", "x_len")]))
K <- 10
B <- 10000
rows <- list()
cat(sprintf("\n== matched controls (K=%d nearest non-members on chr expr, cyto expr, length)\n", K))
cat(sprintf("  %-38s %4s %8s %8s %7s %9s\n", "set", "n", "mean_set", "mean_ctl", "diff", "p_emp"))
for (nm in cand) {
  mem <- which(d$ens %in% pool[[nm]])
  if (length(mem) < 10) next
  non <- setdiff(seq_len(nrow(d)), mem)
  ctl <- t(vapply(mem, function(i) {
    dd <- colSums((t(X[non, , drop = FALSE]) - X[i, ])^2)
    non[order(dd)[1:K]]
  }, integer(K)))
  ms <- mean(d$stat[mem])
  mc <- mean(d$stat[ctl])
  draws <- replicate(B, mean(d$stat[ctl[cbind(seq_along(mem), sample.int(K, length(mem), TRUE))]]))
  # Two-sided: how often does a random matched-control set sit as far from the
  # control mean as the real set does.
  p <- (sum(abs(draws - mc) >= abs(ms - mc)) + 1) / (B + 1)
  rows[[length(rows) + 1]] <- data.frame(set = nm, n = length(mem), mean_set = ms,
                                         mean_ctl = mc, diff = ms - mc, p_emp = p)
  cat(sprintf("  %-38s %4d %8.2f %8.2f %+7.2f %9.4f\n", substr(nm, 1, 38), length(mem), ms, mc, ms - mc, p))
}
ms_tab <- do.call(rbind, rows)
ms_tab$padj <- p.adjust(ms_tab$p_emp, "BH")
write.table(ms_tab, file.path(res_dir, sprintf("gene_%s_matched_sets.tsv", te)), sep = "\t", quote = FALSE, row.names = FALSE)

# Per-lineage stratified GSEA: each lineage's genes ranked only against
# lineage-mates of similar expression and length.
gfile <- file.path(res_dir, sprintf("gene_%s_group_retention.tsv", te))
if (file.exists(gfile)) {
  gret <- read.delim(gfile)
  gret <- gret[!is.na(gret$stat), ]
  cov_g <- cnt
  nx <- sweep(as.matrix(cov_g[, paste0("expr_", libs)]), 2, sf, "/")
  cov_g <- data.frame(gene_id = cov_g$gene_id, x_chr = log10(rowMeans(nx[, 1:2]) + 1),
                      x_cyto = log10(rowMeans(nx[, 3:4]) + 1), x_len = log10(cov_g$gene_len))
  gall <- list()
  for (g in unique(gret$group)) {
    dg <- merge(gret[gret$group == g, ], cov_g, by = "gene_id")
    dg$ens <- sub("\\..*$", "", dg$gene_id)
    dg$z_joint <- within_score(dg$stat, interaction(q(dg$x_chr, 5), q(dg$x_cyto, 5),
                                                     q(dg$x_len, 3), drop = TRUE))
    symg <- setNames(dg$gene_name, dg$ens)
    s <- setNames(dg$z_joint, dg$ens)
    for (coll in names(sets)) {
      x <- fgsea(sets[[coll]], s, minSize = 15, maxSize = 500, eps = 0)
      if (!nrow(x)) next
      x$group <- g
      x$n_genes <- nrow(dg)
      x$collection <- coll
      x$leadingEdge <- vapply(x$leadingEdge, function(y) paste(head(symg[y], 15), collapse = ","), "")
      gall[[length(gall) + 1]] <- as.data.frame(x)
    }
  }
  gg <- do.call(rbind, gall)
  gg <- gg[order(gg$group, gg$collection, gg$padj), ]
  write.table(gg[, c("group", "n_genes", "collection", "pathway", "size", "NES", "pval", "padj", "leadingEdge")],
              file.path(res_dir, sprintf("gene_%s_group_gsea_matched.tsv", te)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n== per-lineage z_joint GSEA (padj<0.25 shown)\n")
  for (g in unique(gg$group)) {
    x <- gg[gg$group == g & !is.na(gg$padj), ]
    cat(sprintf("\n  %s: %d genes, %d sets padj<0.05, %d padj<0.25\n", g, x$n_genes[1],
                sum(x$padj < 0.05), sum(x$padj < 0.25)))
    x <- head(x[order(x$padj), ][x$padj[order(x$padj)] < 0.25, ], 8)
    if (nrow(x)) for (i in seq_len(nrow(x)))
      cat(sprintf("    %-8s %-52s n=%3d NES=%+5.2f padj=%.2g  %s\n", x$collection[i], substr(x$pathway[i], 1, 52),
                  x$size[i], x$NES[i], x$padj[i], substr(x$leadingEdge[i], 1, 50)))
  }
}
