#!/usr/bin/env python
"""Sanity checks on the m6A enrichment call, run before believing it.

Three things could produce a fake "young L1 is methylated" result:

  1. Ratio artefact. log2(IP/input) is a ratio with a noisy denominator, so
     categories with low input coverage drift upward on the median even with no
     real signal. Young L1s are poorly mappable and therefore low-coverage, which
     is exactly the wrong correlation to have. Checked by (a) correlating
     delta_vs_bg against input density across categories and (b) re-running the
     contrast *within* input-density deciles, where that confounder is held flat.

  2. Replicate-driven. Both replicates were pooled before testing, so one bad
     library could carry the whole effect. Checked by re-deriving each category's
     median log2fc per replicate and correlating.

  3. Length. L1s are ~6 kb and Alus ~300 bp. Density normalisation should make
     length cancel in the ratio; checked by correlating delta against copy length.

Usage: diagnostics_m6a.py [results_dir]
"""
import gzip
import os
import sys

import numpy as np
from scipy import stats

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJ, "results", "mintseq")

YOUNG_L1 = {"L1HS", "L1PA2", "L1PA3", "L1PA4", "L1PA5", "L1PA6", "L1PA7",
            "L1P1", "L1P2", "L1P3", "L1P4"}

totals = {}
with open(os.path.join(RES, "library_totals.tsv")) as fh:
    next(fh)
    for line in fh:
        s, strand, v = line.rstrip("\n").split("\t")
        totals[(s, strand)] = float(v)

subfam, strand_plus, length = [], [], []
cols = ("mint_r1_sense", "mint_r1_anti", "mint_r2_sense", "mint_r2_anti",
        "tt_r1_sense", "tt_r1_anti", "tt_r2_sense", "tt_r2_anti")
raw = {c: [] for c in cols}
with gzip.open(os.path.join(RES, "te_copy_signal.tsv.gz"), "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        subfam.append(f[ci["subfamily"]])
        strand_plus.append(f[ci["strand"]] == "+")
        length.append(int(f[ci["length"]]))
        for c in cols:
            raw[c].append(float(f[ci[c]]))

subfam = np.array(subfam)
strand_plus = np.array(strand_plus, bool)
length = np.array(length, np.float64)
for c in cols:
    raw[c] = np.array(raw[c], np.float64)
print("loaded %d copies\n" % len(subfam))

for orient in ("sense", "antisense"):
    key = "sense" if orient == "sense" else "anti"
    from_pos = strand_plus if orient == "sense" else ~strand_plus

    def dens(sample):
        v = raw["%s_%s" % (sample, key)]
        tot = np.where(from_pos, totals[(sample, "pos")], totals[(sample, "neg")])
        return v / (length * tot) * 1e9

    ip1, ip2 = dens("mint_r1"), dens("mint_r2")
    in1, in2 = dens("tt_r1"), dens("tt_r2")
    ip, inp = ip1 + ip2, in1 + in2

    covered = inp > 0
    floor = 0.25 * np.median(inp[covered])
    ok = inp >= floor
    eps = 0.5 * floor
    l2 = np.log2((ip + eps) / (inp + eps))

    print("=" * 72)
    print("%s  (%d copies tested)" % (orient.upper(), ok.sum()))
    print("=" * 72)

    names, delta, dens_med, len_med, n_cop = [], [], [], [], []
    bg = np.median(l2[ok])
    for nm in np.unique(subfam):
        sel = ok & (subfam == nm)
        if sel.sum() < 20:
            continue
        names.append(nm)
        delta.append(np.median(l2[sel]) - bg)
        dens_med.append(np.median(inp[sel]))
        len_med.append(np.median(length[sel]))
        n_cop.append(int(sel.sum()))
    names = np.array(names)
    delta = np.array(delta)
    dens_med = np.array(dens_med)
    len_med = np.array(len_med)

    # (1a) ratio artefact: is delta just a function of input coverage?
    r, p = stats.spearmanr(np.log2(dens_med), delta)
    print("1a. delta vs log2(input density) across %d categories:"
          "  Spearman rho = %+.3f  (p = %.3g)" % (len(names), r, p))

    # (3) length
    r3, p3 = stats.spearmanr(np.log2(len_med), delta)
    print("3.  delta vs log2(copy length):"
          "                        Spearman rho = %+.3f  (p = %.3g)" % (r3, p3))

    # (2) replicate concordance, per category
    d1, d2 = [], []
    for nm in names:
        sel = ok & (subfam == nm)
        a = np.log2((2 * ip1[sel] + eps) / (2 * in1[sel] + eps))
        b = np.log2((2 * ip2[sel] + eps) / (2 * in2[sel] + eps))
        d1.append(np.median(a))
        d2.append(np.median(b))
    r2, p2 = stats.pearsonr(d1, d2)
    print("2.  per-category median log2fc, rep1 vs rep2:"
          "           Pearson r = %+.3f  (p = %.3g)" % (r2, p2))

    # (1b) the decisive one: young L1 vs everything else, held within input-density
    # deciles. If the effect is a low-coverage artefact it vanishes here.
    young = ok & np.isin(subfam, list(YOUNG_L1))
    rest = ok & ~np.isin(subfam, list(YOUNG_L1))
    if young.sum() >= 20:
        edges = np.quantile(inp[ok], np.linspace(0, 1, 11))
        edges[-1] *= 1.001
        print("1b. young L1 vs rest, within input-density deciles:")
        print("      %-22s %8s %8s %9s %11s" %
              ("input-density decile", "n_young", "n_rest", "delta", "MWU p"))
        deltas_in_bin = []
        for i in range(10):
            lo, hi = edges[i], edges[i + 1]
            band = (inp >= lo) & (inp < hi)
            yb, rb = young & band, rest & band
            if yb.sum() < 20 or rb.sum() < 20:
                continue
            d = np.median(l2[yb]) - np.median(l2[rb])
            pv = stats.mannwhitneyu(l2[yb], l2[rb], alternative="two-sided").pvalue
            deltas_in_bin.append(d)
            print("      [%8.3f, %8.3f) %8d %8d %9.3f %11.3g"
                  % (lo, hi, yb.sum(), rb.sum(), d, pv))
        if deltas_in_bin:
            print("      -> median within-decile delta: %+.3f  (%d/%d deciles positive)"
                  % (np.median(deltas_in_bin),
                     sum(1 for d in deltas_in_bin if d > 0), len(deltas_in_bin)))
    print()
