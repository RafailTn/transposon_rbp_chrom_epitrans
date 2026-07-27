#!/usr/bin/env python
"""Separate "young L1" from "long L1".

diagnostics_m6a.py showed the sense-strand m6A signal survives stratification by
input density, but also that delta correlates with copy length at rho=+0.78. Young
L1s (L1HS, L1PA*) are the full-length ones and old L1s (L1M*, L1MA*) are mostly 5'
truncated fragments, so age and length are confounded and the enrichment could be a
length effect wearing an age label.

Two contrasts settle it:

  A. Young vs old L1, inside a length band AND an input-density decile. If age
     carries nothing, the difference goes to zero here.
  B. Within young L1 only, full-length vs truncated copies. L1 has a known m6A
     site cluster in its 5'UTR, which only full-length copies retain -- so a
     positive result here is the expected biology, not a confounder.

Usage: diagnostics_length.py [results_dir]
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
cols = ("mint_r1_sense", "mint_r2_sense", "tt_r1_sense", "tt_r2_sense")
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
from_pos = np.array(strand_plus, bool)          # sense only
length = np.array(length, np.float64)
for c in cols:
    raw[c] = np.array(raw[c], np.float64)


def dens(sample):
    v = raw["%s_sense" % sample]
    tot = np.where(from_pos, totals[(sample, "pos")], totals[(sample, "neg")])
    return v / (length * tot) * 1e9


ip = dens("mint_r1") + dens("mint_r2")
inp = dens("tt_r1") + dens("tt_r2")
floor = 0.25 * np.median(inp[inp > 0])
ok = inp >= floor
eps = 0.5 * floor
l2 = np.log2((ip + eps) / (inp + eps))

is_young = np.isin(subfam, list(YOUNG_L1))
is_l1 = np.char.startswith(subfam, "L1")
is_old_l1 = is_l1 & ~is_young

dec = np.quantile(inp[ok], np.linspace(0, 1, 6))
dec[-1] *= 1.001
BANDS = [(500, 1500), (1500, 3000), (3000, 5000), (5000, 6500), (6500, 1e9)]

print("A. YOUNG vs OLD L1, within a length band and an input-density quintile")
print("   (sense strand; both groups are L1, so mappability and composition match)\n")
print("   %-16s %-18s %7s %7s %8s %11s"
      % ("length band", "input-dens quintile", "n_young", "n_old", "delta", "MWU p"))
rows = []
for lo_l, hi_l in BANDS:
    lb = (length >= lo_l) & (length < hi_l)
    for i in range(5):
        band = lb & (inp >= dec[i]) & (inp < dec[i + 1]) & ok
        y, o = band & is_young, band & is_old_l1
        if y.sum() < 30 or o.sum() < 30:
            continue
        d = np.median(l2[y]) - np.median(l2[o])
        p = stats.mannwhitneyu(l2[y], l2[o], alternative="two-sided").pvalue
        rows.append(d)
        print("   %-16s %-18s %7d %7d %8.3f %11.3g"
              % ("%d-%d" % (lo_l, hi_l), "[%.2f, %.2f)" % (dec[i], dec[i + 1]),
                 y.sum(), o.sum(), d, p))
if rows:
    print("\n   -> median delta across %d matched cells: %+.3f  (%d/%d positive)"
          % (len(rows), np.median(rows), sum(1 for d in rows if d > 0), len(rows)))

print("\n\nB. WITHIN young L1: full-length (>=5.5 kb, retains 5'UTR) vs truncated")
print("   (a positive delta here is the known L1 5'UTR m6A cluster)\n")
print("   %-18s %7s %7s %8s %11s"
      % ("input-dens quintile", "n_full", "n_trunc", "delta", "MWU p"))
for i in range(5):
    band = (inp >= dec[i]) & (inp < dec[i + 1]) & ok & is_young
    full, trunc = band & (length >= 5500), band & (length < 5500)
    if full.sum() < 30 or trunc.sum() < 30:
        continue
    d = np.median(l2[full]) - np.median(l2[trunc])
    p = stats.mannwhitneyu(l2[full], l2[trunc], alternative="two-sided").pvalue
    print("   %-18s %7d %7d %8.3f %11.3g"
          % ("[%.2f, %.2f)" % (dec[i], dec[i + 1]), full.sum(), trunc.sum(), d, p))

print("\n\nC. Reference: length profile by group (median bp, copies passing floor)")
for nm, m in (("young L1", is_young), ("old L1", is_old_l1),
              ("Alu", np.char.startswith(subfam, "Alu")),
              ("SVA", np.char.startswith(subfam, "SVA"))):
    s = m & ok
    if s.sum():
        print("   %-10s n=%7d  median length %6.0f bp  median input dens %6.3f"
              % (nm, s.sum(), np.median(length[s]), np.median(inp[s])))
