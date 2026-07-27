#!/usr/bin/env python
"""Test which L1/Alu/SVA subfamilies carry more m6A than the nascent
transcriptome average, per strand orientation.

Input is te_copy_signal.tsv.gz from quantify_m6a.py. For each TE copy and each
orientation we form a per-copy log2 IP/input ratio:

    dens = signal_sum / (length * library_total)
    log2fc = log2((ip_dens + eps) / (in_dens + eps))

The library_total is taken per *file*, because a copy's sense signal comes from
.pos or .neg depending on its own strand and the two files have different
depths. Dividing by length makes short Alus and 6 kb full-length L1s comparable.

Only copies with enough input coverage are tested -- a copy that is not
transcribed has no nascent RNA to methylate, and its ratio is pseudocount noise
rather than evidence of no methylation. This is the analogue of the
loci_in_cat >= 10 floor in the eCLIP tree.

The test per (subfamily, orientation) is a two-sided Mann-Whitney U of that
category's per-copy log2fc against all other qualifying copies *in the same
orientation*, BH-corrected across categories. Rank-based rather than a
count-based model because bigWig values are CPM-scaled coverage, not raw counts,
so Poisson/binomial dispersion assumptions do not hold.

Usage: enrich_m6a.py [results_dir] [--min-input-frac F] [--min-copies N]
"""
import argparse
import gzip
import os
import sys

import numpy as np
from scipy import stats

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ap = argparse.ArgumentParser()
ap.add_argument("results_dir", nargs="?",
                default=os.path.join(PROJ, "results", "mintseq"))
ap.add_argument("--min-input-frac", type=float, default=0.25,
                help="a copy is tested only if its input density is at least "
                     "this fraction of the median input density over covered "
                     "copies (default 0.25)")
ap.add_argument("--min-copies", type=int, default=20,
                help="minimum qualifying copies for a category to be tested")
args = ap.parse_args()

RES = args.results_dir
SIGNAL = os.path.join(RES, "te_copy_signal.tsv.gz")
TOTALS = os.path.join(RES, "library_totals.tsv")


def bh(p):
    """Benjamini-Hochberg. Returns q in the original order."""
    p = np.asarray(p, float)
    n = p.size
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate(
        (p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main():
    for path in (SIGNAL, TOTALS):
        if not os.path.exists(path):
            sys.exit("missing %s -- run quantify_m6a.py first" % path)

    totals = {}
    with open(TOTALS) as fh:
        next(fh)
        for line in fh:
            sample, strand, v = line.rstrip("\n").split("\t")
            totals[(sample, strand)] = float(v)

    subfam, strand_plus, length = [], [], []
    raw = {c: [] for c in ("mint_r1_sense", "mint_r1_anti", "mint_r2_sense",
                           "mint_r2_anti", "tt_r1_sense", "tt_r1_anti",
                           "tt_r2_sense", "tt_r2_anti")}
    with gzip.open(SIGNAL, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            subfam.append(f[ci["subfamily"]])
            strand_plus.append(f[ci["strand"]] == "+")
            length.append(int(f[ci["length"]]))
            for c in raw:
                raw[c].append(float(f[ci[c]]))

    subfam = np.array(subfam)
    strand_plus = np.array(strand_plus, bool)
    length = np.array(length, np.float64)
    for c in raw:
        raw[c] = np.array(raw[c], np.float64)
    print("loaded %d TE copies" % len(subfam))

    rows = []
    for orient in ("sense", "antisense"):
        key = "sense" if orient == "sense" else "anti"
        # Which bigWig each copy's signal for this orientation came from.
        from_pos = strand_plus if orient == "sense" else ~strand_plus

        def density(sample):
            v = raw["%s_%s" % (sample, key)]
            tot = np.where(from_pos, totals[(sample, "pos")], totals[(sample, "neg")])
            return v / (length * tot) * 1e9

        ip = density("mint_r1") + density("mint_r2")
        inp = density("tt_r1") + density("tt_r2")

        covered = inp > 0
        if not covered.any():
            continue
        floor = args.min_input_frac * np.median(inp[covered])
        ok = inp >= floor
        eps = 0.5 * floor
        log2fc = np.log2((ip + eps) / (inp + eps))

        n_ok = int(ok.sum())
        print("%s: %d/%d copies pass input floor %.4g (%.1f%%)"
              % (orient, n_ok, len(ok), floor, 100.0 * n_ok / len(ok)))

        bg_median = float(np.median(log2fc[ok]))
        for name in np.unique(subfam):
            sel = ok & (subfam == name)
            n = int(sel.sum())
            if n < args.min_copies:
                continue
            other = ok & (subfam != name)
            x, y = log2fc[sel], log2fc[other]
            u = stats.mannwhitneyu(x, y, alternative="two-sided")
            fam = "Alu" if name.startswith("Alu") else ("L1" if name.startswith("L1") else "SVA")
            rows.append({
                "subfamily": name, "family": fam, "orientation": orient,
                "n_copies": n, "n_copies_total": int((subfam == name).sum()),
                "median_log2fc": float(np.median(x)),
                "bg_median_log2fc": bg_median,
                "delta_vs_bg": float(np.median(x)) - bg_median,
                "frac_positive": float((x > 0).mean()),
                "median_input_dens": float(np.median(inp[sel])),
                "p": float(u.pvalue),
            })

    if not rows:
        sys.exit("no category met --min-copies; nothing to test")

    q = bh([r["p"] for r in rows])
    for r, qq in zip(rows, q):
        r["q"] = float(qq)

    cols = ["subfamily", "family", "orientation", "n_copies", "n_copies_total",
            "median_log2fc", "bg_median_log2fc", "delta_vs_bg", "frac_positive",
            "median_input_dens", "p", "q"]
    rows.sort(key=lambda r: -r["delta_vs_bg"])

    os.makedirs(RES, exist_ok=True)
    out = os.path.join(RES, "te_m6a_enrichment.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                ("%.6g" % r[c]) if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")

    sig = [r for r in rows if r["q"] < 0.05 and r["delta_vs_bg"] > 0]
    out_sig = os.path.join(RES, "te_m6a_enrichment_significant.tsv")
    with open(out_sig, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sig:
            fh.write("\t".join(
                ("%.6g" % r[c]) if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")

    print("\n%d categories tested, %d enriched at q<0.05" % (len(rows), len(sig)))
    print("wrote %s\nwrote %s" % (out, out_sig))
    print("\ntop 15 by delta_vs_bg:")
    print("%-12s %-10s %7s %9s %9s %10s" %
          ("subfamily", "orient", "n", "med_l2fc", "delta", "q"))
    for r in rows[:15]:
        print("%-12s %-10s %7d %9.3f %9.3f %10.3g"
              % (r["subfamily"], r["orientation"], r["n_copies"],
                 r["median_log2fc"], r["delta_vs_bg"], r["q"]))


if __name__ == "__main__":
    main()
