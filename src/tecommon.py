"""Primitives shared by the enrichment stages: BED parsing and BH correction.

Imported by enrich_stranded.py, enrich_permutation.py and report_loci_per_cell.py.
Those run as `python src/<script>.py`, so `src/` is sys.path[0] and a plain
`from tecommon import ...` resolves without any packaging.
"""
import numpy as np


def read_bed(path):
    """BED6 or BED9 -> (chrom, start, end, strand) arrays, 0-based half-open.

    strand is +1/-1 (int8) so it can be compared directly against the TE index's
    strand column. The split cap of 6 stops us paying to parse the columns past
    strand, which nothing here uses.

    rstrip() on f[5] is load-bearing, not tidiness: in a BED6 file (the eCLIP
    peaks) strand is the LAST field, so f[5] is "+\\n"/"-\\n" and a bare == "+"
    silently labels every peak minus. That produces no error and no empty output
    -- it just inverts the sense/antisense split this whole analysis rests on.
    A BED9 input would have columns after strand and so would hide the bug
    entirely; do not drop the rstrip because "it works on my file".
    """
    ch, st, en, sd = [], [], [], []
    with open(path) as fh:
        for line in fh:
            f = line.split("\t", 6)
            ch.append(f[0]); st.append(int(f[1])); en.append(int(f[2]))
            sd.append(1 if f[5].rstrip() == "+" else -1)
    return (np.asarray(ch), np.asarray(st, np.int64), np.asarray(en, np.int64),
            np.asarray(sd, np.int8))


def by_chrom(ch, st, en, sd):
    """Group read_bed output into {chrom: (starts, ends, strands)} for the overlap query."""
    return {c: (st[ch == c], en[ch == c], sd[ch == c]) for c in np.unique(ch)}


def bh(p):
    """Benjamini-Hochberg q-values for a flat array of p-values."""
    n = p.size
    o = np.argsort(p, kind="stable")
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)
