#!/usr/bin/env python
"""Aggregate TE counts to (subfamily, orientation, context) and compute the two
normalisations 11_context_deseq2.R needs.

Context comes from 08_context.py (data/te_context.tsv). Writes to results_dir:

  context_counts.tsv   subfamily, orientation, context, family, one count column
                       per sample -- the 07-style matrix, split by context.
  size_factors.tsv     sample, sf_te, sf_gene, intron
      sf_te    total TE-assigned fragments / geometric mean (the 06 factor).
               Compositional: it forces the TE compartment to the same total in
               every library, so it measures retention RELATIVE to other TEs.
      sf_gene  DESeq2 median-of-ratios over GENCODE gene counts (genes with a
               count in every library). Puts the TE compartment against the
               transcriptome, so a global shift of TE reads between fractions
               survives instead of being normalised away.
      intron   total fragments on host introns (09_gene_count.sh), the baseline
               an intronic TE's retention is read against.

Usage: 10_context.py [counts_dir] [gene_counts_dir] [results_dir]
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # project root; scripts live in src/
COUNTS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "work", "counts")
GENES = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "work", "gene_counts")
RES = sys.argv[3] if len(sys.argv) > 3 else os.path.join(ROOT, "results")
CTX = os.path.join(ROOT, "data", "te_context.tsv")

# Reuse the sample list and parsers from 06 rather than re-deriving them.
_spec = importlib.util.spec_from_file_location(
    "enrich", os.path.join(HERE, "06_enrich_chrrna.py"))
enrich = importlib.util.module_from_spec(_spec)
_argv, sys.argv = sys.argv, sys.argv[:1]
_spec.loader.exec_module(enrich)
sys.argv = _argv
SAMPLES = [s for s, _, _ in enrich.SAMPLES]


def main():
    os.makedirs(RES, exist_ok=True)
    gid_sub, gid_fam = enrich.read_key()
    ctx = {}
    with open(CTX) as fh:
        next(fh)
        for line in fh:
            g, c = line.rstrip("\n").split("\t")
            ctx[g] = c

    agg = {}
    te_total = dict.fromkeys(SAMPLES, 0.0)
    for s in SAMPLES:
        for orient in enrich.ORIENTS:
            c = enrich.read_counts(os.path.join(COUNTS, "%s.%s.txt" % (s, orient)))
            for gid, v in c.items():
                te_total[s] += v
                if not v or gid not in gid_sub:
                    continue
                k = (gid_sub[gid], orient, ctx[gid], gid_fam[gid])
                agg.setdefault(k, dict.fromkeys(SAMPLES, 0.0))[s] += v

    with open(os.path.join(RES, "context_counts.tsv"), "w") as fh:
        fh.write("subfamily\torientation\tcontext\tfamily\t" + "\t".join(SAMPLES) + "\n")
        for k in sorted(agg):
            fh.write("\t".join(k) + "\t"
                     + "\t".join("%d" % round(agg[k][s]) for s in SAMPLES) + "\n")

    gene = np.array([[v for v in enrich.read_counts(
        os.path.join(GENES, "%s.genes.txt" % s)).values()] for s in SAMPLES])
    gene = gene[:, (gene > 0).all(axis=0)]
    logg = np.log(gene)
    sf_gene = np.exp(np.median(logg - logg.mean(axis=0), axis=1))
    gm = np.exp(np.mean(np.log(list(te_total.values()))))
    intron = {s: sum(enrich.read_counts(
        os.path.join(GENES, "%s.intron.txt" % s)).values()) for s in SAMPLES}

    with open(os.path.join(RES, "size_factors.tsv"), "w") as fh:
        fh.write("sample\tsf_te\tsf_gene\tintron\n")
        for i, s in enumerate(SAMPLES):
            fh.write("%s\t%.6f\t%.6f\t%d\n" % (s, te_total[s] / gm, sf_gene[i], intron[s]))

    print("%d (subfamily, orientation, context) categories; %d genes in sf_gene"
          % (len(agg), gene.shape[1]))
    print("%-26s %8s %8s %12s" % ("sample", "sf_te", "sf_gene", "intron"))
    for i, s in enumerate(SAMPLES):
        print("%-26s %8.3f %8.3f %12d" % (s, te_total[s] / gm, sf_gene[i], intron[s]))
    print("wrote context_counts.tsv and size_factors.tsv to %s" % RES)


if __name__ == "__main__":
    main()
