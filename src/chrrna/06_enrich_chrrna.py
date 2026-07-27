#!/usr/bin/env python
"""Aggregate featureCounts output to (subfamily, orientation) and score chromatin
retention of L1 / Alu / SVA in K562.

WHY SUBFAMILY AND NOT PER COPY. Everywhere else in this project the counting unit
is the individual TE copy, deliberately. Here it cannot be. STAR was run with
--outSAMmultNmax 1 --outMultimapperOrder Random (the authors' published setting,
kept for comparability), so a read that maps equally well to 500 AluY copies is
credited to one of them at random. Per-copy counts are therefore a random
partition of a subfamily total, not locus-level expression. Aggregating over
copies undoes the randomisation and restores an approximately unbiased number;
reading the per-copy column as expression does not. The copy level survives only
as a support filter (how many copies carry any read at all).

THE CONTRAST. chrRNA abundance on its own has no denominator -- "Alu is highly
represented in chromatin RNA" is close to a restatement of Alu copy number. The
Takara chromatin/cytoplasm pair supplies the matched denominator, the same role
TT-Seq played for MINT-Seq:

    retention = log2( chromatin_takara / cytoplasm_takara )

A subfamily high on this axis is retained on chromatin rather than exported,
which is the interesting property for a TE transcript. The NEB chrRNA samples are
carried through as an independent check that the Takara chromatin fraction agrees
with the primary chrRNA libraries.

STATISTICS. 2 vs 2 replicates over ~183 subfamilies x 2 orientations is a
negative-binomial differential-expression problem, and the right tool is DESeq2
(07_deseq2.R). This script deliberately stops short of inventing a p-value: it
emits the count matrix DESeq2 consumes, plus size-factor-normalised effect sizes
and a per-replicate-pair concordance flag, which are the parts that can be
computed exactly and checked. Do not quote the log2FC here as significant on its
own -- pair it with the DESeq2 q-value.

Usage: 06_enrich_chrrna.py <counts_dir> [results_dir]
"""
import os
import sys

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
COUNTS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJ, "data", "chrrna", "counts")
RES = sys.argv[2] if len(sys.argv) > 2 else os.path.join(PROJ, "results", "chrrna")
KEY = os.path.join(PROJ, "data", "chrrna", "te_copy_key.tsv")

SAMPLES = [
    ("K562_chrRNA_rep1", "chrRNA", 1),
    ("K562_chrRNA_rep2", "chrRNA", 2),
    ("K562_chr_Takara_rep1", "chromatin", 1),
    ("K562_chr_Takara_rep2", "chromatin", 2),
    ("K562_cyto_Takara_rep1", "cytoplasm", 1),
    ("K562_cyto_Takara_rep2", "cytoplasm", 2),
]
ORIENTS = ("sense", "antisense")


def read_key():
    gid_sub, gid_fam = {}, {}
    with open(KEY) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            gid_sub[f[0]] = f[1]
            gid_fam[f[0]] = f[2]
    return gid_sub, gid_fam


def read_counts(path):
    """featureCounts output: 1 comment line, then a header, then Geneid..count."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "Geneid":
                continue
            out[f[0]] = float(f[-1])
    if not out:
        sys.exit("no counts parsed from %s" % path)
    return out


def main():
    if not os.path.isdir(COUNTS):
        sys.exit("no counts dir at %s -- run 05_count.sh on the server first, "
                 "then copy the counts back" % COUNTS)
    os.makedirs(RES, exist_ok=True)
    gid_sub, gid_fam = read_key()

    # (subfamily, orientation) -> {sample: summed count}, and copies with support.
    agg = {}
    support = {}
    libsize = {}
    for sample, _, _ in SAMPLES:
        for orient in ORIENTS:
            path = os.path.join(COUNTS, "%s.%s.txt" % (sample, orient))
            if not os.path.exists(path):
                sys.exit("missing %s" % path)
            c = read_counts(path)
            libsize[(sample, orient)] = sum(c.values())
            for gid, v in c.items():
                sub = gid_sub.get(gid)
                if sub is None:
                    continue
                k = (sub, orient)
                agg.setdefault(k, {}).setdefault(sample, 0.0)
                agg[k][sample] += v
                if v > 0:
                    support.setdefault(k, set()).add(gid)

    keys = sorted(agg)
    print("aggregated %d (subfamily, orientation) categories" % len(keys))

    # Size factors: total TE-assigned reads per library, summed over orientations.
    # Deliberately not a genome-wide total -- the analysis lives entirely inside
    # the TE compartment, and a genome-wide factor would import gene-level
    # composition shifts between the chromatin and cytoplasmic fractions.
    total = {s: sum(libsize[(s, o)] for o in ORIENTS) for s, _, _ in SAMPLES}
    gm = np.exp(np.mean(np.log([v for v in total.values() if v > 0])))
    sf = {s: total[s] / gm for s in total}
    print("\nsize factors (total TE-assigned reads / geometric mean):")
    for s, _, _ in SAMPLES:
        print("  %-26s %12.0f  sf=%.3f" % (s, total[s], sf[s]))

    matrix_path = os.path.join(RES, "subfamily_counts.tsv")
    with open(matrix_path, "w") as fh:
        fh.write("subfamily\torientation\tfamily\t"
                 + "\t".join(s for s, _, _ in SAMPLES) + "\n")
        for sub, orient in keys:
            row = agg[(sub, orient)]
            gid = next(iter(support.get((sub, orient), {"TE0"})))
            fam = gid_fam.get(gid, "NA")
            fh.write("%s\t%s\t%s\t%s\n" % (
                sub, orient, fam,
                "\t".join("%d" % round(row.get(s, 0.0)) for s, _, _ in SAMPLES)))
    print("\nwrote %s  (feed this to 07_deseq2.R)" % matrix_path)

    # Descriptive effect sizes: normalised chromatin/cytoplasm ratio, per
    # replicate pair, so concordance is visible rather than averaged away.
    eps = 1.0
    rows = []
    for sub, orient in keys:
        row = agg[(sub, orient)]

        def n(s):
            return row.get(s, 0.0) / sf[s]

        c1, c2 = n("K562_chr_Takara_rep1"), n("K562_chr_Takara_rep2")
        y1, y2 = n("K562_cyto_Takara_rep1"), n("K562_cyto_Takara_rep2")
        l1 = np.log2((c1 + eps) / (y1 + eps))
        l2 = np.log2((c2 + eps) / (y2 + eps))
        chr_nb = 0.5 * (n("K562_chrRNA_rep1") + n("K562_chrRNA_rep2"))
        rows.append({
            "subfamily": sub, "orientation": orient,
            "family": gid_fam.get(next(iter(support.get((sub, orient), {"TE0"}))), "NA"),
            "n_copies_with_reads": len(support.get((sub, orient), ())),
            "chr_mean": 0.5 * (c1 + c2), "cyto_mean": 0.5 * (y1 + y2),
            "chrRNA_neb_mean": chr_nb,
            "log2_retention_rep1": l1, "log2_retention_rep2": l2,
            "log2_retention": 0.5 * (l1 + l2),
            "reps_agree": "yes" if (l1 > 0) == (l2 > 0) else "NO",
        })

    cols = ["subfamily", "orientation", "family", "n_copies_with_reads",
            "chr_mean", "cyto_mean", "chrRNA_neb_mean",
            "log2_retention_rep1", "log2_retention_rep2", "log2_retention",
            "reps_agree"]
    rows.sort(key=lambda r: -r["log2_retention"])
    out = os.path.join(RES, "chromatin_retention.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                ("%.4f" % r[c]) if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")
    print("wrote %s" % out)

    disagree = sum(1 for r in rows if r["reps_agree"] == "NO")
    print("\n%d/%d categories have replicate pairs disagreeing in sign"
          % (disagree, len(rows)))
    print("\ntop 15 chromatin-retained (replicate-concordant, >=20 copies with reads):")
    print("  %-12s %-10s %-5s %8s %9s %9s"
          % ("subfamily", "orient", "fam", "copies", "log2ret", "cyto"))
    shown = 0
    for r in rows:
        if r["reps_agree"] != "yes" or r["n_copies_with_reads"] < 20:
            continue
        print("  %-12s %-10s %-5s %8d %9.3f %9.1f"
              % (r["subfamily"], r["orientation"], r["family"],
                 r["n_copies_with_reads"], r["log2_retention"], r["cyto_mean"]))
        shown += 1
        if shown >= 15:
            break


if __name__ == "__main__":
    main()
