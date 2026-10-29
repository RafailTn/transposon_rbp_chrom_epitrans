#!/usr/bin/env python
"""Per-gene counts for the host-gene function test: intronic Alu RNA vs the
gene's own intronic RNA, chromatin vs cytoplasm.

QUESTION. Among genes that carry sense-oriented intronic Alus (context
intron_same), is there a functional theme to the genes whose Alu sequence stays
on chromatin MORE than the rest of their introns do? 18_gene_alu_gsea.R scores
each gene and runs GSEA; this script makes the counts.

WHY UNIQUE FRAGMENTS ONLY. Per-copy counts from 05 are random partitions of
multi-mapping reads (03_align.sh), so a gene's Alu count there would include
reads from Alus anywhere in the genome. Here featureCounts runs WITHOUT -M, so
only NH==1 fragments count: a read is credited to an Alu only if it can only be
that copy. Old and middle-aged Alus (AluJ, most AluS) are diverged enough to
carry such reads; young AluY mostly are not, and drop out by construction.

Per gene, on the gene's strand (-s 2, --countReadPairs):
  ALU     the gene's intron_same Alu copies (copy midpoint inside exactly one
          protein_coding/lncRNA gene body on that strand)
  INTRON  the gene's introns MINUS every TE copy (any family, either strand),
          minus exons of all gene types and the rDNA mask -- the baseline is
          non-repeat intronic RNA from the same pre-mRNA
Intronic sequence shared by two same-strand genes is dropped from both.

Writes data/gene_alu/{alu,intron}.saf, counts in work/gene_alu_counts/, and
results/gene_alu_counts.tsv with covariates (gene length, intron bp, n Alu,
gene expression) for 18.

Usage: 17_gene_alu.py [bam_dir] [gencode.gtf.gz]
"""
import bisect
import gzip
import importlib.util
import os
import subprocess
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # project root; scripts live in src/
BAMDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "work", "bam")
GTF = sys.argv[2] if len(sys.argv) > 2 else (
    "/mnt/raid1/rafail/aluaca/ref/gencode.v47.primary_assembly.annotation.gtf.gz")
DATA = os.path.join(ROOT, "data", "gene_alu")
CNT = os.path.join(ROOT, "work", "gene_alu_counts")
RES = os.path.join(ROOT, "results")
THREADS = os.environ.get("THREADS", "16")

# Interval helpers and the gene-type / rDNA definitions from 08, not re-derived.
_spec = importlib.util.spec_from_file_location("ctx", os.path.join(HERE, "08_context.py"))
ctx = importlib.util.module_from_spec(_spec)
_argv, sys.argv = sys.argv, sys.argv[:1]
_spec.loader.exec_module(ctx)
sys.argv = _argv


def samples():
    with open(os.path.join(ROOT, "config", "runs.tsv")) as fh:
        next(fh)
        return [l.split("\t")[0] for l in fh if l.strip()]


def build():
    os.makedirs(DATA, exist_ok=True)
    genes = {}                                   # gene_id -> (chrom, strand, s, e, name)
    all_exons = defaultdict(list)
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[2] == "exon":
                all_exons[(f[0], f[6])].append([int(f[3]), int(f[4])])
            elif f[2] == "gene" and ctx.attr(f[8], "gene_type") in ctx.GENE_TYPES:
                if not f[0].startswith("chr") or f[0] == "chrM":
                    continue
                genes[ctx.attr(f[8], "gene_id")] = (f[0], f[6], int(f[3]), int(f[4]),
                                                    ctx.attr(f[8], "gene_name"))

    # Same-strand gene-body overlaps: sequence covered by >1 gene is dropped.
    by_cs = defaultdict(list)
    for gid, (c, st, s, e, _) in genes.items():
        by_cs[(c, st)].append((s, e, gid))
    shared = defaultdict(list)
    for k, v in by_cs.items():
        ev = sorted([(s, 1) for s, e, _ in v] + [(e + 1, -1) for s, e, _ in v])
        depth, start = 0, None
        for pos, d in ev:
            if depth >= 2 and start is not None and pos > start:
                shared[k].append([start, pos - 1])
            depth += d
            start = pos

    te_all = defaultdict(list)
    alu = []
    ctx_lab = {}
    with open(os.path.join(ROOT, "data", "te_context.tsv")) as fh:
        next(fh)
        for line in fh:
            g, c = line.rstrip("\n").split("\t")
            ctx_lab[g] = c
    with open(os.path.join(ROOT, "data", "te_copy_key.tsv")) as fh:
        next(fh)
        for line in fh:
            gid, sub, fam, chrom, s, e, strand = line.rstrip("\n").split("\t")
            s, e = int(s), int(e)
            te_all[chrom].append([s, e])
            if fam == "Alu" and ctx_lab[gid] == "intron_same":
                alu.append((gid, sub, chrom, s, e, strand))

    # Assign each Alu to the single same-strand gene whose body holds its
    # midpoint: one sweep per (chrom, strand) over genes sorted by start.
    alu_by = defaultdict(list)
    for gid, sub, chrom, s, e, strand in alu:
        alu_by[(chrom, strand)].append(((s + e) // 2, chrom, s, e, strand, sub))
    alu_gene = defaultdict(list)
    for k, items in alu_by.items():
        gl = sorted(by_cs.get(k, []))
        items.sort()
        i, active = 0, []
        for mid, chrom, s, e, strand, sub in items:
            while i < len(gl) and gl[i][0] <= mid:
                active.append(gl[i])
                i += 1
            active = [x for x in active if x[1] >= mid]
            if len(active) == 1:
                alu_gene[active[0][2]].append((chrom, s, e, strand, sub))

    with open(os.path.join(DATA, "alu.saf"), "w") as fh:
        fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for g in sorted(alu_gene):
            for chrom, s, e, strand, _ in alu_gene[g]:
                fh.write("%s\t%s\t%d\t%d\t%s\n" % (g, chrom, s, e, strand))

    # One merged mask per (chrom, strand); each gene subtracts only the slice
    # of it that overlaps its body.
    masks = {}
    for k in by_cs:
        chrom, strand = k
        masks[k] = ctx.merge(all_exons.get(k, []) + shared.get(k, []) + te_all.get(chrom, [])
                             + [[a, b] for mc, a, b in ctx.RDNA_MASK if mc == chrom])
    mstarts = {k: [x[0] for x in v] for k, v in masks.items()}
    cov = open(os.path.join(DATA, "covariates.tsv"), "w")
    cov.write("gene_id\tgene_name\tgene_len\tintron_bp\tn_alu\tn_aluJ\tn_aluS\tn_aluY\n")
    with open(os.path.join(DATA, "intron.saf"), "w") as fh:
        fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for g in sorted(alu_gene):
            chrom, strand, s, e, name = genes[g]
            k = (chrom, strand)
            lo = max(0, bisect.bisect_right(mstarts[k], s) - 1)
            hi = bisect.bisect_right(mstarts[k], e)
            ivs = ctx.subtract([[s, e]], masks[k][lo:hi])
            bp = sum(b - a + 1 for a, b in ivs)
            for a, b in ivs:
                fh.write("%s\t%s\t%d\t%d\t%s\n" % (g, chrom, a, b, strand))
            subs = [x[4] for x in alu_gene[g]]
            cov.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (
                g, name, e - s + 1, bp, len(subs),
                sum(x.startswith("AluJ") for x in subs),
                sum(x.startswith("AluS") for x in subs),
                sum(x.startswith("AluY") for x in subs)))
    cov.close()
    print("%d genes carry intron_same Alus assigned to a single gene (%d Alu copies)"
          % (len(alu_gene), sum(len(v) for v in alu_gene.values())))


def count():
    os.makedirs(CNT, exist_ok=True)
    for s in samples():
        bam = os.path.join(BAMDIR, "%s.Aligned.sortedByCoord.out.bam" % s)
        strand = os.environ.get("STRAND_" + s, os.environ.get("STRAND_DEFAULT", "2"))
        for kind in ("alu", "intron"):
            out = os.path.join(CNT, "%s.%s.txt" % (s, kind))
            if os.path.exists(out) and os.path.getsize(out):
                continue
            print("== %s %s" % (s, kind), flush=True)
            # No -M: unique (NH==1) fragments only -- see the docstring.
            subprocess.run(["featureCounts", "-F", "SAF", "-T", THREADS, "-p",
                            "--countReadPairs", "-s", strand,
                            "-a", os.path.join(DATA, kind + ".saf"), "-o", out, bam],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read_fc(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            f = line.rstrip("\n").split("\t")
            out[f[0]] = int(float(f[-1]))
    return out


def assemble():
    ss = samples()
    alu = {s: read_fc(os.path.join(CNT, "%s.alu.txt" % s)) for s in ss}
    intr = {s: read_fc(os.path.join(CNT, "%s.intron.txt" % s)) for s in ss}
    expr = {s: read_fc(os.path.join(ROOT, "work", "gene_counts", "%s.genes.txt" % s)) for s in ss}
    with open(os.path.join(DATA, "covariates.tsv")) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        cov = [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh]
    out = os.path.join(RES, "gene_alu_counts.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(head + ["alu_" + s for s in ss] + ["intron_" + s for s in ss]
                           + ["expr_" + s for s in ss]) + "\n")
        for r in cov:
            g = r["gene_id"]
            fh.write("\t".join([r[h] for h in head]
                               + [str(alu[s].get(g, 0)) for s in ss]
                               + [str(intr[s].get(g, 0)) for s in ss]
                               + [str(expr[s].get(g, 0)) for s in ss]) + "\n")
    tak = [s for s in ss if "Takara" in s]
    tot_a = {s: sum(alu[s].values()) for s in tak}
    tot_i = {s: sum(intr[s].values()) for s in tak}
    print("unique fragments, Takara:")
    for s in tak:
        print("  %-24s alu %10d  intron %12d" % (s, tot_a[s], tot_i[s]))
    print("wrote %s" % out)


if __name__ == "__main__":
    if not os.path.exists(os.path.join(DATA, "covariates.tsv")):
        build()
    count()
    assemble()
