#!/usr/bin/env python
"""Partition TE copies by host-gene context, to separate real TE transcription
from host pre-mRNA passenger reads.

WHY THIS EXISTS. The chromatin/cytoplasm contrast in 06 has a confounder it
cannot see. Chromatin RNA is nascent and largely unspliced; cytoplasmic RNA is
mature and spliced. A TE copy sitting inside an intron, on the SAME strand as its
host gene, therefore collects host pre-mRNA reads into its own *sense* channel --
present in the chromatin fraction, absent from the cytoplasmic one -- and scores
as "chromatin-retained" with no TE transcription involved anywhere. L1 copies are
overwhelmingly intronic/intergenic, so the observed L1-sense retention is exactly
what that artefact would produce.

The mirror case matters too, and cuts the other way: a TE overlapping an EXON on
the host's strand rides in mature mRNA and is exported. Alu is the classic
exonized element, abundant in 3'UTRs, so Alu-sense export is exactly what
exonization alone would produce.

Both stories are testable by the same partition. Re-run 06 on the `intergenic`
class:

  * L1-sense retention survives at similar magnitude -> the signal is about L1.
  * It collapses -> the headline was host-gene pre-mRNA.

Copies whose host lies on the OPPOSITE strand are a useful control rather than a
nuisance: host pre-mRNA lands in their antisense channel, so their sense channel
is as clean as an intergenic copy's.

CLASSES

  intergenic       overlaps no annotated gene on either strand
  intronic_sense   inside a gene body on the TE's own strand, no exon overlap
  exonic_sense     overlaps an exon of a gene on the TE's own strand
  intronic_anti    inside a gene body on the opposite strand only, no exon overlap
  exonic_anti      overlaps an exon of a gene on the opposite strand only
  ambiguous        overlapped by genes on BOTH strands -- neither channel is clean
  no_annotation    on a sequence GENCODE carries no genes for -- EXCLUDED, not
                   silently folded into intergenic (see below)

COORDINATES. This is the one script in the project where nothing is converted.
te_copy_key.tsv carries the GTF coordinates through unchanged (make_saf.py copies
them straight across, since SAF is GTF-like), and GENCODE is a GTF too, so both
sides are 1-based inclusive. Overlap is `a_start <= b_end and a_end >= b_start`.
Do not paste the half-open test from build_te_index.py or the eclip scripts in
here -- it is off by one against these inputs.

"Genic" is defined against ALL GENCODE gene records, not just protein_coding.
That deliberately makes `intergenic` the strictest class: a copy called
intergenic is one no annotation of any biotype can explain, which is what the
test needs. It shrinks the intergenic set rather than contaminating it.

GENCODE is not shipped with this project. Fetch the hg38 primary annotation:

  wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/\\
gencode.v48.annotation.gtf.gz
  # (or browse .../latest_release/ for the current version number)

Chromosome names must match the TE key's UCSC style (chr1, chr2, ...), which the
standard GENCODE human release uses for the 25 real chromosomes. It names unplaced
scaffolds Ensembl-style (GL000008.2) where the TE key uses UCSC (chrUn_GL000008v1),
so those never match, and neither do _alt/_fix patch contigs, which the primary
assembly omits entirely.

Copies on any such sequence are classed `no_annotation` and are NOT part of
`intergenic`. The distinction is load-bearing: "no gene overlaps this copy" and
"nothing here was ever examined for genes" are different statements, and merging
them would quietly seed the control class -- the one the whole deconfounding test
rests on -- with unexamined copies. They keep a row in the output so the partition
stays auditable against the key, but no TE_KEEP class selects them.

Usage: 08_partition_copies.py <gencode.gtf[.gz]> [te_copy_key.tsv] [out.tsv]
"""
import gzip
import os
import sys
from collections import defaultdict

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
GENCODE = sys.argv[1] if len(sys.argv) > 1 else None
KEY = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    PROJ, "data", "chrrna", "te_copy_key.tsv")
OUT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(
    PROJ, "data", "chrrna", "te_copy_class.tsv")

# Brute-forced against a naive scan at the end. chr21 is small enough to hold
# raw unmerged intervals for, and carries all three families.
CHECK_CHROM = "chr21"
N_CHECK = 400

CLASSES = ["intergenic", "intronic_sense", "exonic_sense",
           "intronic_anti", "exonic_anti", "ambiguous", "no_annotation"]
(INTERGENIC, INTRONIC_S, EXONIC_S, INTRONIC_A, EXONIC_A, AMBIGUOUS,
 NO_ANNOT) = range(7)


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def merge(intervals):
    """Sort and merge 1-based inclusive intervals into non-overlapping runs.

    Merging is what makes the lookup a single indexed comparison rather than a
    window scan: for sorted NON-OVERLAPPING intervals, if any interval overlaps
    a query then the last one starting at or before the query's end must be the
    one that does. That identity fails on raw gene records, which overlap freely
    (nested genes, readthrough transcripts), hence the merge.
    """
    if not intervals:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    intervals.sort()
    starts, ends = [], []
    cs, ce = intervals[0]
    for s, e in intervals[1:]:
        if s <= ce + 1:          # overlapping or exactly abutting
            if e > ce:
                ce = e
        else:
            starts.append(cs)
            ends.append(ce)
            cs, ce = s, e
    starts.append(cs)
    ends.append(ce)
    return np.array(starts, np.int64), np.array(ends, np.int64)


def read_gencode(path):
    """-> merged gene and exon intervals, keyed by (chrom, strand)."""
    genes = defaultdict(list)
    exons = defaultdict(list)
    raw_genes = defaultdict(list)   # CHECK_CHROM only, unmerged, for the assert
    raw_exons = defaultdict(list)
    n_gene = n_exon = 0
    with opener(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 8:
                continue
            feat = f[2]
            if feat != "gene" and feat != "exon":
                continue
            chrom, strand = f[0], f[6]
            if strand != "+" and strand != "-":
                continue
            s, e = int(f[3]), int(f[4])
            if feat == "gene":
                genes[(chrom, strand)].append((s, e))
                n_gene += 1
                if chrom == CHECK_CHROM:
                    raw_genes[strand].append((s, e))
            else:
                exons[(chrom, strand)].append((s, e))
                n_exon += 1
                if chrom == CHECK_CHROM:
                    raw_exons[strand].append((s, e))
    if n_gene == 0:
        sys.exit("no 'gene' records in %s -- is this a GENCODE GTF?" % path)
    print("read %d gene and %d exon records from %s"
          % (n_gene, n_exon, os.path.basename(path)))
    gm = {k: merge(v) for k, v in genes.items()}
    em = {k: merge(v) for k, v in exons.items()}
    print("  merged to %d gene and %d exon intervals"
          % (sum(len(v[0]) for v in gm.values()),
             sum(len(v[0]) for v in em.values())))
    return gm, em, raw_genes, raw_exons


EMPTY = (np.empty(0, np.int64), np.empty(0, np.int64))


def overlaps(merged, te_s, te_e):
    """Vectorised 1-based-inclusive overlap against merged intervals."""
    ms, me = merged
    out = np.zeros(te_s.size, bool)
    if ms.size == 0:
        return out
    idx = np.searchsorted(ms, te_e, side="right") - 1
    ok = idx >= 0
    if not ok.any():
        return out
    out[ok] = me[idx[ok]] >= te_s[ok]
    return out


def classify(gm, em, chrom, te_s, te_e, te_plus):
    g_p = overlaps(gm.get((chrom, "+"), EMPTY), te_s, te_e)
    g_m = overlaps(gm.get((chrom, "-"), EMPTY), te_s, te_e)
    e_p = overlaps(em.get((chrom, "+"), EMPTY), te_s, te_e)
    e_m = overlaps(em.get((chrom, "-"), EMPTY), te_s, te_e)

    host_same = np.where(te_plus, g_p, g_m)
    host_opp = np.where(te_plus, g_m, g_p)
    exon_same = np.where(te_plus, e_p, e_m)
    exon_opp = np.where(te_plus, e_m, e_p)

    cls = np.full(te_s.size, INTERGENIC, np.int8)
    both = host_same & host_opp
    only_s = host_same & ~host_opp
    only_a = host_opp & ~host_same
    cls[only_s] = np.where(exon_same[only_s], EXONIC_S, INTRONIC_S)
    cls[only_a] = np.where(exon_opp[only_a], EXONIC_A, INTRONIC_A)
    cls[both] = AMBIGUOUS
    return cls


def naive(raw_g, raw_e, s, e, strand):
    """Independent O(n) reimplementation, for the correctness assert."""
    def hit(iv):
        return any(a <= e and b >= s for a, b in iv)
    opp = "-" if strand == "+" else "+"
    hs, ho = hit(raw_g[strand]), hit(raw_g[opp])
    if hs and ho:
        return AMBIGUOUS
    if hs:
        return EXONIC_S if hit(raw_e[strand]) else INTRONIC_S
    if ho:
        return EXONIC_A if hit(raw_e[opp]) else INTRONIC_A
    return INTERGENIC


def main():
    if GENCODE is None or not os.path.exists(GENCODE):
        sys.exit("usage: 08_partition_copies.py <gencode.gtf.gz> [key] [out]\n"
                 "GENCODE is not shipped with this project; see the module "
                 "docstring for the download URL.")
    if not os.path.exists(KEY):
        sys.exit("missing %s -- run make_saf.py first" % KEY)

    gm, em, raw_g, raw_e = read_gencode(GENCODE)

    # Read the key once, grouped by chromosome, keeping the original row order
    # so the output can be written back aligned to it.
    by_chrom = defaultdict(lambda: ([], [], [], []))   # idx, start, end, plus
    fam_of = []
    n = 0
    with open(KEY) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            i, st, en, pl = by_chrom[f[3]]
            i.append(n)
            st.append(int(f[4]))
            en.append(int(f[5]))
            pl.append(f[6] == "+")
            fam_of.append(f[2])
            n += 1
    print("read %d TE copies across %d sequences" % (n, len(by_chrom)))

    cls_all = np.full(n, INTERGENIC, np.int8)
    missing = []
    for chrom, (idx, st, en, pl) in by_chrom.items():
        ix = np.array(idx, np.int64)
        if (chrom, "+") not in gm and (chrom, "-") not in gm:
            # GENCODE never looked at this sequence. Not the same thing as
            # "no gene overlaps these copies", so it does not get to be
            # intergenic -- that class is the control for the whole test.
            missing.append((chrom, len(idx)))
            cls_all[ix] = NO_ANNOT
            continue
        cls_all[ix] = classify(gm, em, chrom,
                               np.array(st, np.int64), np.array(en, np.int64),
                               np.array(pl, bool))

    # ---- correctness check, same spirit as the other scripts in this tree ----
    if CHECK_CHROM in by_chrom and raw_g:
        idx, st, en, pl = by_chrom[CHECK_CHROM]
        rng = np.random.default_rng(0)
        pick = rng.choice(len(idx), size=min(N_CHECK, len(idx)), replace=False)
        bad = 0
        for j in pick:
            want = naive(raw_g, raw_e, st[j], en[j], "+" if pl[j] else "-")
            if want != cls_all[idx[j]]:
                bad += 1
        if bad:
            sys.exit("BUG: vectorised classification disagrees with the naive "
                     "scan on %d/%d %s copies" % (bad, len(pick), CHECK_CHROM))
        print("check: %d random %s copies match a naive interval scan"
              % (len(pick), CHECK_CHROM))
    else:
        print("WARNING: %s absent -- correctness check skipped" % CHECK_CHROM)

    # ---- write, aligned to the key's row order ----
    counts = defaultdict(int)
    with open(KEY) as fh, open(OUT, "w") as out:
        next(fh)
        out.write("gene_id\tclass\n")
        for row, line in enumerate(fh):
            gid = line.split("\t", 1)[0]
            c = CLASSES[cls_all[row]]
            out.write("%s\t%s\n" % (gid, c))
            counts[(fam_of[row], c)] += 1
    print("\nwrote %s" % OUT)

    fams = [f for f in ("Alu", "L1", "SVA")
            if sum(counts[(f, c)] for c in CLASSES)]
    # Percentages are of the ANNOTATED copies on each family, i.e. excluding
    # no_annotation, so they describe the population the test actually runs on.
    denom = {f: sum(counts[(f, c)] for c in CLASSES if c != "no_annotation")
             for f in fams}
    print("\n=== copies per class (%% of annotated copies in that family) ===")
    print("%-16s %s" % ("class", "  ".join("%18s" % f for f in fams)))
    for c in CLASSES:
        cells = []
        for f in fams:
            n_ = counts[(f, c)]
            cells.append("%18s" % ("%d" % n_ if c == "no_annotation"
                                   else "%d (%.1f%%)" % (n_, 100 * n_ / denom[f])))
        print("%-16s %s" % (c, "  ".join(cells)))
    print("%-16s %s" % ("  annotated", "  ".join("%18d" % denom[f] for f in fams)))

    if missing:
        tot = sum(k for _, k in missing)
        print("\n%d copies (%.1f%%) on %d sequences carrying no GENCODE genes were "
              "classed no_annotation" % (tot, 100.0 * tot / n, len(missing)))
        print("and are EXCLUDED -- no TE_KEEP class selects them. They are not "
              "intergenic:\nnothing there was examined, which is a different claim.")
        for c, k in sorted(missing, key=lambda x: -x[1])[:5]:
            print("      %-24s %d" % (c, k))

    print("\nNext: re-run 06 on one class at a time, e.g.\n"
          "  TE_CLASSES=%s TE_KEEP=intergenic \\\n"
          "    python src/chrrna/06_enrich_chrrna.py <counts_dir> "
          "results/chrrna_intergenic\n"
          "then 07_deseq2.R on that directory. Compare the L1-sense log2FC "
          "against the all-copies run." % OUT)


if __name__ == "__main__":
    main()
