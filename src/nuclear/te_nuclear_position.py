#!/usr/bin/env python
"""Where do L1 / Alu / SVA DNA copies sit relative to nuclear compartments?

H1-hESC TSA-seq, three axes, each the mean of two replicates:

    speckle    SON        4DNFI625PP2A + 4DNFIFKMOD1L
    lamina     Lamin B1   4DNFINXMR3OD + 4DNFILYDJU8T
    nucleolus  MKI67IP    4DNFIFNIK4HD + 4DNFIGU18ZPJ
               POL1RE     4DNFICLC27JC + 4DNFICYLNXF2   (kept separate: two
                          antibodies against the same compartment are a
                          reproducibility check, r~0.87, and are not averaged
                          into one column so a disagreement stays visible)

THIS IS A DNA-POSITION QUESTION. It asks where the copies are, not what their
RNA does. Nothing here bears on the sense/antisense hypothesis the RNA stages
test -- TSA-seq has no strand.

RESOLUTION, AND WHAT IT FORBIDS. The tracks are binned at 20 kb. An Alu is
~300 bp, SVA ~1-3 kb, a full L1 ~6 kb, so every copy is far sub-bin and inherits
the value of its neighbourhood. There is therefore NO per-copy nuclear position
here: two copies in one bin are indistinguishable by construction. The unit that
means anything is the subfamily aggregate, and the claim it supports is "the
neighbourhoods where these copies sit are compartment X-proximal", never "this
copy is at compartment X".

THE RESULT IS PREDETERMINED UNLESS GC IS CONTROLLED. Alu is concentrated in
GC-rich sequence, L1 in AT-rich, and the compartments themselves lie along that
same axis -- speckle-proximal chromatin is GC-rich, lamina-proximal is AT-rich.
So the raw ranking recovers 1980s compositional genomics with a 2021 assay. Both
columns are reported:

    <axis>_mean   raw mean over the subfamily's copies
    <axis>_resid  mean of (copy value - mean of all usable bins in the same GC
                  decile), i.e. the part NOT explained by base composition

Read the residual for anything load-bearing. The raw column is kept because it
is the descriptive answer to "where are they", and because seeing the two
disagree is the point.

lam_minus_nuc is reported for the same reason in a sharper form. The lamina and
the nucleolar periphery are BOTH heterochromatic destinations, so GC and gene
density predict that L1-rich DNA goes to "heterochromatin" but do not predict
WHICH. Their difference cancels the shared compartment-independent component and
is the one contrast here that is not a restatement of composition.

MASKING. cache/bin_covariates_20000.npz supplies gap (n_frac) and satellite
(sat_frac) fractions; bins over 10% of either are dropped, plus any bin where an
axis is exactly 0 (no reads -> the value is absent, not low). Unmasked, measured
nucleolar signal runs ~1.4 log2 units high across centromeres, which would land
on exactly the satellite-adjacent L1 copies and read as a finding.

COPY-WEIGHTED MEAN, BIN-CLUSTERED ERROR BAR. Alu averages 8.1 copies per
occupied 20 kb bin and L1 6.9, and every copy in a bin carries the same value,
so the two obvious aggregations answer different questions and one of them is
useless here:

  - Deduplicating to bins asks "what is a bin that contains >=1 copy like?".
    For Alu that is 97% of all usable bins and for L1 94%, so the answer is
    "like the genome" by construction. Tried it: Alu's speckle mean moved from
    +0.180 to -0.120 against a -0.126 background, i.e. the signal vanished --
    not because it was fake but because presence/absence throws away copy
    NUMBER, which is the whole compositional signal for a ubiquitous family.
  - Copy-weighting asks "what is the mean environment of a randomly chosen
    copy?". That is the question, and it is the point estimate reported.

The clustering problem is therefore about the ERROR BAR, not the estimate:
n_copies overstates precision by ~sqrt(copies per bin). So every SE here is
cluster-robust with the bin as the cluster, and n_eff is Kish's effective sample
size (sum w)^2 / sum w^2 over per-bin copy counts. Quote n_eff, not n_copies.
For SVA the distinction barely matters (1.2 copies/bin); for Alu it is a factor
of ~3.

GC IS NOT A SUFFICIENT CONTROL, AND top_chrom_frac IS WHY. L1M3a came out with
the highest speckle residual AND the lowest lamina residual of any subfamily
while sitting at GC 0.384, i.e. below average. The explanation is that 26% of its
copies are on chr19, the most gene-dense and speckle-proximal chromosome in the
genome, whose compartment behaviour is not captured by its bin-level GC. So
chromosome concentration is reported per subfamily, and a subfamily with a high
top_chrom_frac is describing a chromosome, not a repeat family.

Usage: te_nuclear_position.py [--min-bins N] [--bin N]
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pyBigWig

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TSA = os.path.join(PROJ, "data", "tsa")
KEY = os.path.join(PROJ, "data", "chrrna", "te_copy_key.tsv")
CLS = os.path.join(PROJ, "data", "chrrna", "te_copy_class.tsv")
CACHE = os.path.join(PROJ, "cache")
OUT = os.path.join(PROJ, "results", "nuclear")

AXES = [("speckle", ["4DNFI625PP2A", "4DNFIFKMOD1L"]),
        ("lamina", ["4DNFINXMR3OD", "4DNFILYDJU8T"]),
        ("nucleolus_MKI67IP", ["4DNFIFNIK4HD", "4DNFIGU18ZPJ"]),
        ("nucleolus_POL1RE", ["4DNFICLC27JC", "4DNFICYLNXF2"])]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin", type=int, default=20000)
    p.add_argument("--min-bins", type=int, default=50,
                   help="subfamilies occupying fewer usable BINS are dropped from "
                        "the per-subfamily table. Bins are the sample size here, "
                        "not copies -- see the docstring on pseudo-replication.")
    p.add_argument("--deciles", type=int, default=10,
                   help="number of GC strata for the residual")
    return p.parse_args()


def load_bins(binsize):
    p = os.path.join(CACHE, "bin_covariates_%d.npz" % binsize)
    if not os.path.exists(p):
        sys.exit("missing %s -- run src/nuclear/bin_covariates.py first" % p)
    z = np.load(p)
    if int(z["binsize"]) != binsize:
        sys.exit("binsize mismatch")
    chroms = [str(x) for x in z["chroms"]]
    return (chroms, dict(zip(chroms, z["offsets"])), dict(zip(chroms, z["nbins"])),
            z["gc_frac"], z["n_frac"], z["sat_frac"])


def load_axis(accs, chroms, off, nb, binsize, total):
    """Mean of the replicate bigWigs, on the common 20 kb grid.

    bw.stats(..., nBins=n) is used rather than values(): the track is already
    piecewise-constant at 20 kb, so a binned mean reproduces it exactly while
    reading a fraction of the data.
    """
    acc = np.zeros(total, np.float64)
    seen = np.zeros(total, np.int32)
    for a in accs:
        path = os.path.join(TSA, "%s.bw" % a)
        if not os.path.exists(path):
            sys.exit("missing %s" % path)
        bw = pyBigWig.open(path)
        for c in chroms:
            n = int(nb[c])
            L = bw.chroms().get(c)
            if L is None:
                continue
            m = min(n, L // binsize)
            if m <= 0:
                continue
            v = np.array(bw.stats(c, 0, m * binsize, type="mean", nBins=m),
                         dtype=np.float64)
            s = int(off[c])
            good = np.isfinite(v)
            acc[s:s + m][good] += v[good]
            seen[s:s + m][good] += 1
        bw.close()
    out = np.full(total, np.nan)
    ok = seen == len(accs)          # require BOTH replicates to cover the bin
    out[ok] = acc[ok] / len(accs)
    return out


def main():
    a = parse_args()
    os.makedirs(OUT, exist_ok=True)
    chroms, off, nb, gc, nfrac, sfrac = load_bins(a.bin)
    total = int(sum(nb[c] for c in chroms))

    print("loading %d TSA-seq tracks over %d bins ..." % (
        sum(len(v) for _, v in AXES), total), flush=True)
    sig = {}
    for name, accs in AXES:
        sig[name] = load_axis(name and accs, chroms, off, nb, a.bin, total)
        v = sig[name]
        print("  %-18s %d bins with signal  mean %+.3f" % (
            name, np.isfinite(v).sum(), np.nanmean(v)), flush=True)

    # ---- mask ----------------------------------------------------------------
    mask = (nfrac >= 0.10) | (sfrac >= 0.10) | ~np.isfinite(gc)
    for name in sig:
        mask |= ~np.isfinite(sig[name])
        mask |= sig[name] == 0.0
    usable = ~mask
    print("\nusable bins: %d of %d (%.1f%%)" % (usable.sum(), total,
                                                100.0 * usable.mean()))

    # ---- GC strata, on usable bins only -------------------------------------
    edges = np.nanpercentile(gc[usable], np.linspace(0, 100, a.deciles + 1))
    edges[0] -= 1e-9
    edges[-1] += 1e-9
    strat = np.full(total, -1, np.int8)
    strat[usable] = np.clip(np.searchsorted(edges, gc[usable], side="right") - 1,
                            0, a.deciles - 1)
    smean = {}
    for name in sig:
        m = np.full(a.deciles, np.nan)
        for d in range(a.deciles):
            k = usable & (strat == d)
            if k.any():
                m[d] = sig[name][k].mean()
        smean[name] = m
    print("GC deciles (usable bins): %s" % np.round(edges, 3))
    print("  speckle  by decile: %s" % np.round(smean["speckle"], 2))
    print("  lamina   by decile: %s" % np.round(smean["lamina"], 2))

    # ---- TE copies ----------------------------------------------------------
    print("\nreading TE copies ...", flush=True)
    kcls = {}
    for ln in open(CLS).readlines()[1:]:
        g, c = ln.rstrip("\n").split("\t")
        kcls[g] = c
    fam, sub, cls, binidx = [], [], [], []
    with open(KEY) as fh:
        fh.readline()
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            c = f[3]
            if c not in off:
                continue
            # midpoint decides the bin. te_copy_key is GTF 1-based inclusive, so
            # convert to 0-based before dividing; at 20 kb the off-by-one cannot
            # change a bin except exactly on a boundary, but be exact anyway.
            mid = ((int(f[4]) - 1) + int(f[5])) // 2
            b = mid // a.bin
            if b >= int(nb[c]):
                continue
            binidx.append(int(off[c]) + b)
            sub.append(f[1])
            fam.append(f[2])
            cls.append(kcls.get(f[0], "?"))
    binidx = np.array(binidx, np.int64)
    fam = np.array(fam)
    sub = np.array(sub)
    cls = np.array(cls)
    keep = usable[binidx]
    print("  %d copies, %d in usable bins (%.1f%%)"
          % (binidx.size, keep.sum(), 100.0 * keep.mean()))
    for f in ("Alu", "L1", "SVA"):
        k = fam == f
        nb_f = np.unique(binidx[k & keep]).size
        print("    %-4s %8d copies -> %8d usable (%.1f%% of this family) "
              "in %7d distinct bins (%.1f copies/bin)"
              % (f, k.sum(), (k & keep).sum(),
                 100.0 * (k & keep).sum() / max(k.sum(), 1), nb_f,
                 (k & keep).sum() / max(nb_f, 1)))

    binidx, fam, sub, cls = binidx[keep], fam[keep], sub[keep], cls[keep]

    # Bin -> chromosome, for the concentration diagnostic. offsets are ascending
    # by construction, so one searchsorted recovers the chromosome of any bin.
    ordered = sorted(chroms, key=lambda c: int(off[c]))
    obounds = np.array([int(off[c]) for c in ordered])

    names = [n for n, _ in AXES]

    def wstat(v, w):
        """Copy-weighted mean with a bin-clustered (Huber-White) SE.

        v and w are per-BIN: v the bin's value, w how many of this subfamily's
        copies fall in it. Var(m) = sum w^2 (v-m)^2 / (sum w)^2 treats each bin
        as one independent cluster, which is the assumption that fails if you
        divide by sqrt(n_copies) instead.
        """
        W = w.sum()
        m = float((w * v).sum() / W)
        if v.size < 2:
            return m, float("nan")
        se = float(np.sqrt((w ** 2 * (v - m) ** 2).sum()) / W)
        return m, se

    def rows(labels, tag):
        out = []
        for lab in sorted(set(labels)):
            k = labels == lab
            ub, w = np.unique(binidx[k], return_counts=True)
            w = w.astype(np.float64)
            n_bins = ub.size
            n_cop = int(k.sum())
            n_eff = float(w.sum() ** 2 / (w ** 2).sum())
            ci = np.searchsorted(obounds, ub, side="right") - 1
            _, cc = np.unique(ci, return_counts=True)
            r = {"level": tag, "name": lab, "n_copies": n_cop, "n_bins": n_bins,
                 "n_eff": n_eff, "copies_per_bin": n_cop / float(n_bins),
                 "top_chrom_frac": float(cc.max()) / n_bins,
                 "gc_mean": float((w * gc[ub]).sum() / w.sum())}
            if tag == "subfamily":
                fs = set(fam[k])
                r["family"] = fs.pop() if len(fs) == 1 else "|".join(sorted(fs))
            else:
                r["family"] = lab
            st = strat[ub]
            for nm in names:
                v = sig[nm][ub]
                m, se = wstat(v, w)
                mr, ser = wstat(v - smean[nm][st], w)
                r["%s_mean" % nm] = m
                r["%s_se" % nm] = se
                r["%s_resid" % nm] = mr
                r["%s_resid_se" % nm] = ser
            lm = sig["lamina"][ub] - sig["nucleolus_MKI67IP"][ub]
            lmr = ((sig["lamina"][ub] - smean["lamina"][st])
                   - (sig["nucleolus_MKI67IP"][ub] - smean["nucleolus_MKI67IP"][st]))
            r["lam_minus_nuc"], r["lam_minus_nuc_se"] = wstat(lm, w)
            r["lam_minus_nuc_resid"], r["lam_minus_nuc_resid_se"] = wstat(lmr, w)
            out.append(r)
        return out

    allrows = rows(fam, "family")
    subrows = [r for r in rows(sub, "subfamily") if r["n_bins"] >= a.min_bins]
    # genome background, for reference
    nb_u = int(usable.sum())
    bg = {"level": "background", "name": "all_usable_bins", "family": "-",
          "n_copies": nb_u, "n_bins": nb_u, "n_eff": float(nb_u),
          "copies_per_bin": 1.0, "top_chrom_frac": float("nan"),
          "gc_mean": float(gc[usable].mean()),
          "lam_minus_nuc": float((sig["lamina"][usable]
                                  - sig["nucleolus_MKI67IP"][usable]).mean()),
          "lam_minus_nuc_resid": 0.0, "lam_minus_nuc_se": float("nan"),
          "lam_minus_nuc_resid_se": float("nan")}
    for nm in names:
        bg["%s_mean" % nm] = float(sig[nm][usable].mean())
        bg["%s_resid" % nm] = 0.0
        bg["%s_se" % nm] = float("nan")
        bg["%s_resid_se" % nm] = float("nan")

    cols = (["level", "name", "family", "n_copies", "n_bins", "n_eff",
             "copies_per_bin", "top_chrom_frac", "gc_mean"]
            + ["%s_%s" % (n, s) for n in names
               for s in ("mean", "se", "resid", "resid_se")]
            + ["lam_minus_nuc", "lam_minus_nuc_se",
               "lam_minus_nuc_resid", "lam_minus_nuc_resid_se"])
    path = os.path.join(OUT, "te_nuclear_position.tsv")
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in [bg] + allrows + sorted(subrows, key=lambda x: (x["family"], x["name"])):
            fh.write("\t".join(
                ("%.5g" % r[c]) if isinstance(r[c], float) else str(r.get(c, ""))
                for c in cols) + "\n")

    # ---- report -------------------------------------------------------------
    print("\n=== family level: raw_mean / GC-residual (log2 TSA-seq) ===")
    print("%-6s %9s %8s %8s %6s  %s" % ("fam", "copies", "bins", "n_eff", "GC",
                                        "  ".join("%-15s" % n[:15] for n in names)))
    for r in [bg] + allrows:
        print("%-6s %9d %8d %8.0f %6.3f  %s" % (
            r["name"][:6], r["n_copies"], r["n_bins"], r["n_eff"], r["gc_mean"],
            "  ".join("%+6.3f/%+6.3f" % (r["%s_mean" % n], r["%s_resid" % n])
                      for n in names)))
    print("\n=== GC-residual with cluster-robust SE (bin = cluster) ===")
    for r in allrows:
        print("  %-4s %s" % (r["name"][:4], "  ".join(
            "%-9s %+6.3f+-%.3f" % (n[:9], r["%s_resid" % n], r["%s_resid_se" % n])
            for n in ("speckle", "lamina", "nucleolus_MKI67IP"))))

    print("\n=== lamina - nucleolus, the contrast GC does not predict ===")
    for r in [bg] + allrows:
        print("  %-6s raw %+.3f   GC-residual %+.3f  (SE %.4f)"
              % (r["name"][:6], r["lam_minus_nuc"], r["lam_minus_nuc_resid"],
                 r["lam_minus_nuc_se"]))

    print("\n=== nucleolar antibody agreement (should track each other) ===")
    for r in allrows:
        print("  %-4s MKI67IP %+.3f   POL1RE %+.3f   diff %+.3f"
              % (r["name"][:4], r["nucleolus_MKI67IP_mean"],
                 r["nucleolus_POL1RE_mean"],
                 r["nucleolus_MKI67IP_mean"] - r["nucleolus_POL1RE_mean"]))

    def show(title, keyf, sekey):
        top = sorted(subrows, key=lambda r: -r[keyf])[:8]
        bot = sorted(subrows, key=lambda r: r[keyf])[:8]
        print("\n--- %s: highest / lowest subfamilies (>=%d bins) ---"
              % (title, a.min_bins))
        print("     %8s %8s %5s %6s %6s  %s" % ("value", "SE", "bins", "GC",
                                                "topchr", "subfamily"))
        for grp in (top, ["..."], bot):
            for r in grp:
                if r == "...":
                    print("     ...")
                    continue
                flag = " <-- chrom-driven" if r["top_chrom_frac"] > 0.15 else ""
                print("     %+8.3f %8.3f %5d %6.3f %6.2f  %-12s %-4s%s"
                      % (r[keyf], r[sekey], r["n_bins"], r["gc_mean"],
                         r["top_chrom_frac"], r["name"][:12], r["family"], flag))

    for axis in ("speckle", "lamina"):
        show("%s GC-residual" % axis, "%s_resid" % axis, "%s_resid_se" % axis)
    show("lamina-minus-nucleolus GC-residual", "lam_minus_nuc_resid",
         "lam_minus_nuc_resid_se")

    conc = [r for r in subrows if r["top_chrom_frac"] > 0.15]
    if conc:
        print("\n%d of %d subfamilies have >15%% of their bins on one chromosome; "
              "those rows describe a chromosome, not a family." % (len(conc), len(subrows)))

    print("\nwrote %s  (%d rows: 1 background, %d family, %d subfamily)"
          % (path, 1 + len(allrows) + len(subrows), len(allrows), len(subrows)))


if __name__ == "__main__":
    main()
