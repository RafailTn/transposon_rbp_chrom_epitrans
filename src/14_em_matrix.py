#!/usr/bin/env python
"""Turn TEcount tables (13_em_count.sh) into the subfamily_counts.tsv layout of
06_enrich_chrrna.py, so 07_deseq2.R runs on the EM arm unchanged:

    python src/14_em_matrix.py
    Rscript src/07_deseq2.R results/em

TEcount rows are genes (ENSG ids) and TEs named 'subfamily:family:class'. Only
the TE rows go in the matrix; gene rows are identical between the sense and
antisense runs and are written once to results/em/gene_counts_em.tsv.

Counts are EM-fractional; they are rounded, as DESeq2 needs integers.

Usage: 14_em_matrix.py [counts_em_dir] [results_dir]
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)          # project root; scripts live in src/
CNT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "work", "counts_em")
RES = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "results", "em")

_spec = importlib.util.spec_from_file_location(
    "enrich", os.path.join(HERE, "06_enrich_chrrna.py"))
enrich = importlib.util.module_from_spec(_spec)
_argv, sys.argv = sys.argv, sys.argv[:1]
_spec.loader.exec_module(enrich)
sys.argv = _argv
SAMPLES = [s for s, _, _ in enrich.SAMPLES]


def read_cnt(path):
    te, genes = {}, {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            name, v = line.rstrip("\n").split("\t")
            parts = name.split(":")
            if len(parts) == 3:
                te[(parts[0], parts[1])] = float(v)
            else:
                genes[name] = float(v)
    return te, genes


def main():
    os.makedirs(RES, exist_ok=True)
    mat, gmat = {}, {}
    for s in SAMPLES:
        for orient in enrich.ORIENTS:
            path = os.path.join(CNT, "%s.%s.cntTable" % (s, orient))
            if not os.path.exists(path):
                sys.exit("missing %s -- run 13_em_count.sh" % path)
            te, genes = read_cnt(path)
            for (sub, fam), v in te.items():
                mat.setdefault((sub, orient, fam), {})[s] = v
            if orient == "sense":
                for g, v in genes.items():
                    gmat.setdefault(g, {})[s] = v

    with open(os.path.join(RES, "subfamily_counts.tsv"), "w") as fh:
        fh.write("subfamily\torientation\tfamily\t" + "\t".join(SAMPLES) + "\n")
        for k in sorted(mat):
            fh.write("\t".join(k) + "\t"
                     + "\t".join("%d" % round(mat[k].get(s, 0.0)) for s in SAMPLES) + "\n")
    with open(os.path.join(RES, "gene_counts_em.tsv"), "w") as fh:
        fh.write("gene_id\t" + "\t".join(SAMPLES) + "\n")
        for g in sorted(gmat):
            fh.write(g + "\t" + "\t".join("%d" % round(gmat[g].get(s, 0.0)) for s in SAMPLES) + "\n")

    print("%d (subfamily, orientation) rows, %d genes -> %s" % (len(mat), len(gmat), RES))
    for s in SAMPLES:
        tot = sum(v.get(s, 0.0) for v in mat.values())
        print("  %-26s TE fragments %12.0f" % (s, tot))


if __name__ == "__main__":
    main()
