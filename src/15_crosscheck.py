#!/usr/bin/env python
"""Cross-checks for the retention results. Prints a report and writes
results/crosscheck_*.tsv.

1. Against the authors (GSE296212 per-sample tables in data/geo_tables/).
   Their tables are htseq-count output: gene symbols, ~10k named L1 loci
   (L1HS_1..332, L1PA*, L1M* -- a curated, mostly full-length subset, with no
   coordinates), and __no_feature / __alignment_not_unique rows. htseq drops
   NH>1 alignments, and their totals equal READ counts. So absolute numbers are
   not expected to match ours; what should match is
     a. gene expression (log CPM over shared genes), per sample;
     b. the ranking of full-length L1 subfamilies by abundance (their loci
        summed by name prefix vs our >= 5 kb copies summed per subfamily);
     c. full-length L1 retention, log2(chromatin/cytoplasm) on gene totals,
        computed the same way on both sides.
2. Random-placement arm (results/) vs EM arm (results/em/): per-category
   log2FC agreement; categories with |delta| > 1 are quantification-sensitive.
3. NEB chrRNA vs Takara chromatin: antisense share of TE fragments, and the
   per-subfamily sense/antisense ratio agreement (the open question from the
   v1 review).

Usage: 15_crosscheck.py [results_dir] [results_em_dir]
(per-copy counts are read from work/counts/, gene counts from work/gene_counts/)
"""
import glob
import gzip
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # project root; scripts live in src/
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results")
RES_EM = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "results", "em")
GEO = os.path.join(ROOT, "data", "geo_tables")
GENES = os.path.join(ROOT, "work", "gene_counts")
COUNTS = os.path.join(ROOT, "work", "counts")
GTF = "/mnt/raid1/rafail/aluaca/ref/gencode.v47.primary_assembly.annotation.gtf.gz"
FL_MIN = 5000


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh]


def enrich_counts(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            f = line.rstrip("\n").split("\t")
            out[f[0]] = float(f[-1])
    return out


def spearman(a, b):
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


def fnum(x):
    return float("nan") if x in ("NA", "") else float(x)


def geo_tables():
    out = {}
    for path in sorted(glob.glob(os.path.join(GEO, "GSM*.txt.gz"))):
        sample = re.sub(r"^GSM\d+_|\.txt\.gz$", "", os.path.basename(path))
        d = {}
        with gzip.open(path, "rt") as fh:
            next(fh)
            for line in fh:
                k, v = line.rstrip("\n").split("\t")
                d[k] = float(v)
        out[sample] = d
    return out


def gene_names():
    m = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if "\tgene\t" not in line:
                continue
            a = line.split("\t")[8]
            gid = re.search(r'gene_id "([^"]+)"', a).group(1)
            m[gid] = re.search(r'gene_name "([^"]+)"', a).group(1)
    return m


def check_geo(rows_out):
    geo = geo_tables()
    samples = list(geo)
    print("== 1. against the authors' GEO tables (%d samples)" % len(samples))

    # a. genes
    if os.path.isdir(GENES):
        names = gene_names()
        print("  a. gene log-CPM Pearson r (shared genes, both > 0):")
        for s in samples:
            ours = {}
            with open(os.path.join(GENES, "%s.genes.txt" % s)) as fh:
                for line in fh:
                    if line.startswith("#") or line.startswith("Geneid"):
                        continue
                    f = line.rstrip("\n").split("\t")
                    n = names.get(f[0])
                    if n:
                        ours[n] = ours.get(n, 0.0) + float(f[-1])
            theirs = {k: v for k, v in geo[s].items()
                      if not k.startswith("__") and not k.startswith("L1")}
            shared = [g for g in theirs if g in ours and ours[g] > 0 and theirs[g] > 0]
            to, tt = sum(ours.values()), sum(theirs.values())
            a = np.log2([ours[g] / to * 1e6 + 1 for g in shared])
            b = np.log2([theirs[g] / tt * 1e6 + 1 for g in shared])
            r = float(np.corrcoef(a, b)[0, 1])
            print("     %-26s r=%.3f  (n=%d)" % (s, r, len(shared)))
            rows_out.append(("geo_genes", s, "pearson_logcpm", "%.4f" % r))

    # b/c. L1. Their ~10k loci are, by count per subfamily, the full-length
    # copies: our L1 copies >= 5 kb reproduce their per-subfamily locus numbers
    # (L1HS 332 vs 356, L1PA2 1105 vs 1165, L1M1 65 vs 66). So sum OUR per-copy
    # sense counts over copies >= FL_MIN only, not the whole subfamily.
    key = {}
    with open(os.path.join(ROOT, "data", "te_copy_key.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[2] == "L1" and int(f[5]) - int(f[4]) + 1 >= FL_MIN:
                key[f[0]] = f[1]
    ours = {}
    for s in samples:
        agg = {}
        with open(os.path.join(COUNTS, "%s.sense.txt" % s)) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] in key:
                    agg[key[f[0]]] = agg.get(key[f[0]], 0.0) + float(f[-1])
        ours[s] = agg

    print("  b. full-length L1 (>= %d bp) subfamily abundance, Spearman vs theirs:" % FL_MIN)
    geo_sub = {}
    for s in samples:
        agg = {}
        for k, v in geo[s].items():
            if k.startswith("L1"):
                sub = re.sub(r"_\d+$", "", k)
                agg[sub] = agg.get(sub, 0.0) + v
        geo_sub[s] = agg
        shared = sorted(x for x in agg if x in ours[s])
        rho = spearman([agg[x] for x in shared], [ours[s][x] for x in shared])
        print("     %-26s rho=%.3f  (n=%d subfamilies)" % (s, rho, len(shared)))
        rows_out.append(("geo_L1_abundance", s, "spearman", "%.4f" % rho))

    # c. Retention direction on the same full-length L1 set, both sides
    # normalised to their own gene totals so the denominators are alike.
    def lfc(tab, gtot, sub):
        c = sum(tab[s].get(sub, 0.0) / gtot[s] for s in samples if "_chr_Takara" in s)
        y = sum(tab[s].get(sub, 0.0) / gtot[s] for s in samples if "_cyto_" in s)
        return c, y

    geo_gtot = {s: sum(v for k, v in geo[s].items()
                       if not k.startswith("__") and not k.startswith("L1")) for s in samples}
    our_gtot = {s: sum(enrich_counts(os.path.join(GENES, "%s.genes.txt" % s)).values())
                for s in samples} if os.path.isdir(GENES) else None
    if our_gtot is None:
        print("  c. skipped (no %s)" % GENES)
        return
    t_l, o_l, kept = [], [], []
    for x in sorted(set(geo_sub[samples[0]]) & set(ours[samples[0]])):
        tc, ty = lfc(geo_sub, geo_gtot, x)
        oc, oy = lfc(ours, our_gtot, x)
        if min(tc, ty, oc, oy) <= 0:
            continue
        t_l.append(math.log2(tc / ty))
        o_l.append(math.log2(oc / oy))
        kept.append(x)
    rho = spearman(t_l, o_l)
    agree = float(np.mean(np.sign(t_l) == np.sign(o_l)))
    print("  c. full-length L1 log2(chromatin/cytoplasm), gene-normalised: Spearman %.3f, "
          "sign agreement %.0f%%, median theirs %.2f vs ours %.2f (n=%d)"
          % (rho, 100 * agree, np.median(t_l), np.median(o_l), len(kept)))
    rows_out.append(("geo_L1_retention", "all", "spearman", "%.4f" % rho))
    rows_out.append(("geo_L1_retention", "all", "sign_agreement", "%.4f" % agree))


def check_arms(rows_out):
    p_em = os.path.join(RES_EM, "chromatin_retention_deseq2.tsv")
    if not os.path.exists(p_em):
        print("== 2. EM arm: %s missing, skipped" % p_em)
        return
    a = {(r["subfamily"], r["orientation"]): r
         for r in read_tsv(os.path.join(RES, "chromatin_retention_deseq2.tsv"))}
    b = {(r["subfamily"], r["orientation"]): r for r in read_tsv(p_em)}
    print("== 2. random-placement vs EM arm")
    out = open(os.path.join(RES, "crosscheck_arms.tsv"), "w")
    out.write("subfamily\torientation\tfamily\tlfc_random\tpadj_random\tlfc_em\tpadj_em\tdelta\n")
    for orient in ("sense", "antisense"):
        keys = [k for k in a if k in b and k[1] == orient]
        x = np.array([fnum(a[k]["log2FoldChange"]) for k in keys])
        y = np.array([fnum(b[k]["log2FoldChange"]) for k in keys])
        ok = ~(np.isnan(x) | np.isnan(y))
        r = float(np.corrcoef(x[ok], y[ok])[0, 1])
        big = [k for k, xi, yi, o in zip(keys, x, y, ok) if o and abs(xi - yi) > 1]
        print("  %-9s r=%.3f  n=%d  |delta|>1: %d" % (orient, r, ok.sum(), len(big)))
        for k in sorted(big, key=lambda k: -abs(fnum(a[k]["log2FoldChange"]) - fnum(b[k]["log2FoldChange"])))[:10]:
            print("     %-12s %-5s random %6.2f  EM %6.2f"
                  % (k[0], a[k]["family"], fnum(a[k]["log2FoldChange"]), fnum(b[k]["log2FoldChange"])))
        rows_out.append(("arms", orient, "pearson_lfc", "%.4f" % r))
        rows_out.append(("arms", orient, "n_delta_gt1", str(len(big))))
        for k in keys:
            out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n" % (
                k[0], k[1], a[k]["family"], a[k]["log2FoldChange"], a[k]["padj"],
                b[k]["log2FoldChange"], b[k]["padj"],
                fnum(b[k]["log2FoldChange"]) - fnum(a[k]["log2FoldChange"])))
    out.close()


def check_kits(rows_out):
    rows = read_tsv(os.path.join(RES, "subfamily_counts.tsv"))
    neb = ("K562_chrRNA_rep1", "K562_chrRNA_rep2")
    tak = ("K562_chr_Takara_rep1", "K562_chr_Takara_rep2")
    print("== 3. NEB chrRNA vs Takara chromatin")
    for s in neb + tak + ("K562_cyto_Takara_rep1", "K562_cyto_Takara_rep2"):
        tot = sum(float(r[s]) for r in rows)
        anti = sum(float(r[s]) for r in rows if r["orientation"] == "antisense")
        print("  %-26s antisense share %.2f" % (s, anti / tot))
        rows_out.append(("kits", s, "antisense_share", "%.4f" % (anti / tot)))
    d = {(r["subfamily"], r["orientation"]): r for r in rows}
    subs = sorted({r["subfamily"] for r in rows
                   if all((r["subfamily"], o) in d for o in ("sense", "antisense"))})

    def ratio(sub, grp):
        sn = sum(float(d[(sub, "sense")][s]) for s in grp)
        an = sum(float(d[(sub, "antisense")][s]) for s in grp)
        return math.log2((sn + 1) / (an + 1))
    rn = [ratio(x, neb) for x in subs]
    rt = [ratio(x, tak) for x in subs]
    r = float(np.corrcoef(rn, rt)[0, 1])
    print("  per-subfamily log2(sense/antisense): NEB median %.2f, Takara median %.2f, r=%.3f"
          % (np.median(rn), np.median(rt), r))
    rows_out.append(("kits", "all", "sense_antisense_ratio_r", "%.4f" % r))


def main():
    rows_out = []
    check_geo(rows_out)
    check_arms(rows_out)
    check_kits(rows_out)
    with open(os.path.join(RES, "crosscheck_summary.tsv"), "w") as fh:
        fh.write("check\tsample\tmetric\tvalue\n")
        for r in rows_out:
            fh.write("\t".join(r) + "\n")
    print("\nwrote %s" % os.path.join(RES, "crosscheck_summary.tsv"))


if __name__ == "__main__":
    main()
