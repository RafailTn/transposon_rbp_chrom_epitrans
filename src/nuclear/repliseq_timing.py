#!/usr/bin/env python
"""Build a per-20kb-bin replication-timing profile from 16-fraction Repli-seq.

WHY THIS REPLACES GC AS THE COVARIATE. The nuclear-position question is confounded
by large-scale genome compartmentalisation, and GC fraction is only a proxy for
it -- a single scalar that co-varies with gene density, replication timing, CpG
density and chromosome identity without holding any of them fixed. L1M3a proved
the gap concretely: below-average GC (0.384) but the most extreme compartment
residual of any subfamily, entirely because 26% of its copies are on chr19.

Replication timing is a direct measurement of the same compartmentalisation and
is the best-behaved summary of it available. It is measured, not inferred from
sequence, so it absorbs gene density and chromosome identity together.

THE DATA. 4DN multi-stage Repli-seq on H1-hESC (Tier 1): G1 plus 16 S-phase
fractions P2_S..P17_S, two replicates each, 5 kb bedGraph, GRCh38. The 17
"experiment sets" are one experiment split BY FRACTION, not 17 replicates.

THE STATISTIC. Per bin, the 16 S-fraction signals are normalised to sum to 1 and
RT is their weighted mean fraction index,

    RT = sum_i i * p_i        i = 1..16 over P2_S..P17_S

so RT is a continuous "average replication time" in fraction units. Direction is
NOT assumed: it is checked against GC downstream (early replication is GC-rich),
and asserted here against the sign of that correlation.

G1 is kept but not used in RT. It is the non-replicating control and its role is
copy-number / mappability correction; on a karyotypically normal line that
correction is small, and dividing by a noisy denominator per bin would add
variance. It is stored so the decision can be revisited.

COVERAGE MASK. s_total is the summed S-phase signal per bin. Bins with little
total signal have an RT that is a ratio of noise, so downstream code should
require s_total above a floor rather than trusting every bin.

Usage: repliseq_timing.py [--bin N] [--force]
"""
import argparse
import os
import subprocess
import sys

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MANIFEST = os.path.join(PROJ, "data", "nuclear_manifest.tsv")
CACHE = os.path.join(PROJ, "cache")

# S-phase fractions in temporal order as deposited. P2 -> P17 is 16 fractions.
S_FRACTIONS = ["P%d_S" % i for i in range(2, 18)]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin", type=int, default=20000)
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def load_grid(binsize):
    p = os.path.join(CACHE, "bin_covariates_%d.npz" % binsize)
    if not os.path.exists(p):
        sys.exit("missing %s -- run bin_covariates.py first" % p)
    z = np.load(p)
    chroms = [str(x) for x in z["chroms"]]
    return (chroms, dict(zip(chroms, z["offsets"])), dict(zip(chroms, z["nbins"])),
            z["gc_frac"])


def read_manifest():
    by = {}
    with open(MANIFEST) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ip, ia, il = hdr.index("path"), hdr.index("assay"), hdr.index("label")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[ia] != "repliseq":
                continue
            by.setdefault(f[il], []).append(f[ip])
    return by


def sum_into_bins(paths, chroms, off, nb, binsize, total):
    """Sum the replicate bedGraphs of one fraction onto the 20 kb grid.

    Aggregated with awk before it reaches python: 34 files x ~618k lines is 21 M
    lines, and awk does the bin reduction in a fraction of the time. Intervals
    are 5 kb and bin-aligned so each lands in one 20 kb bin, but the general
    overlap-weighted form is used anyway -- chromosome ends carry short bins.

    LC_ALL=C IS LOAD-BEARING, exactly as the project notes warn. Under this
    machine's comma-decimal locale awk parses "0.356" as 0, so every bedGraph
    value below 1 silently became zero: half the bins came out with no S-phase
    signal at all and the G1 total read 97.5 instead of ~2e4. The values here are
    normalised coverage and nearly all of them are <1, so this failure mode
    destroys the profile rather than perturbing it.
    """
    env = dict(os.environ, LC_ALL="C", LC_NUMERIC="C")
    out = np.zeros(total, np.float64)
    for p in paths:
        full = os.path.join(PROJ, p) if not os.path.isabs(p) else p
        if not os.path.exists(full):
            sys.exit("missing %s" % full)
        cmd = ("zcat '%s' | awk -F'\\t' -v B=%d '"
               "{ s=$2; e=$3; v=$4;"
               "  b0=int(s/B); b1=int((e-1)/B);"
               "  if (b0==b1) { a[$1\"\\t\"b0]+=v*(e-s) }"
               "  else { for (b=b0;b<=b1;b++) {"
               "      lo=(b*B>s)?b*B:s; hi=((b+1)*B<e)?(b+1)*B:e;"
               "      a[$1\"\\t\"b]+=v*(hi-lo) } } }"
               "END{ for (k in a) print k\"\\t\"a[k] }'" % (full, binsize))
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                           env=env)
        if r.returncode != 0:
            sys.exit("awk failed on %s: %s" % (full, r.stderr[:300]))
        for ln in r.stdout.splitlines():
            c, b, v = ln.split("\t")
            if c not in off:
                continue
            b = int(b)
            if b >= int(nb[c]):
                continue
            out[int(off[c]) + b] += float(v)
    return out / float(binsize)      # back to a mean-signal-per-bp scale


def main():
    a = parse_args()
    outp = os.path.join(CACHE, "repliseq_rt_%d.npz" % a.bin)
    if os.path.exists(outp) and not a.force:
        print("%s exists; --force to rebuild" % outp)
        return
    chroms, off, nb, gc = load_grid(a.bin)
    total = int(sum(int(nb[c]) for c in chroms))
    by = read_manifest()
    missing = [f for f in S_FRACTIONS + ["G1"] if f not in by]
    if missing:
        sys.exit("manifest lacks fractions %s" % missing)

    print("summing %d S-phase fractions (2 replicates each) onto %d bins ..."
          % (len(S_FRACTIONS), total), flush=True)
    F = np.zeros((len(S_FRACTIONS), total), np.float64)
    for i, fr in enumerate(S_FRACTIONS):
        F[i] = sum_into_bins(by[fr], chroms, off, nb, a.bin, total)
        nz = (F[i] > 0).mean()
        print("  %-6s %d files  total signal %.3g  nonzero bins %.1f%%"
              % (fr, len(by[fr]), F[i].sum(), 100 * nz), flush=True)
        # The bedGraphs are ~92% nonzero at 5 kb, so a fraction covering well
        # under half the bins means the values were mis-parsed, not that the
        # fraction is sparse. Guards the comma-decimal locale bug specifically.
        if nz < 0.50:
            sys.exit("%s covers only %.1f%% of bins -- expected ~90%%. Check that "
                     "awk is running under LC_ALL=C." % (fr, 100 * nz))
    g1 = sum_into_bins(by["G1"], chroms, off, nb, a.bin, total)
    print("  %-6s %d files  total signal %.3g" % ("G1", len(by["G1"]), g1.sum()))

    # Per-fraction depth normalisation before mixing them: the fractions were
    # sequenced to different depths, and RT is a ratio ACROSS fractions, so a
    # depth difference would masquerade as a timing shift.
    scale = F.sum(axis=1, keepdims=True)
    if (scale <= 0).any():
        sys.exit("a fraction has no signal at all")
    Fn = F / scale * scale.mean()

    s_total = Fn.sum(axis=0)
    idx = np.arange(1, len(S_FRACTIONS) + 1, dtype=np.float64)[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        rt = np.where(s_total > 0, (Fn * idx).sum(axis=0) / np.maximum(s_total, 1e-12),
                      np.nan)

    ok = np.isfinite(rt) & np.isfinite(gc) & (s_total > 0)
    r = np.corrcoef(rt[ok], gc[ok])[0, 1]
    print("\nRT range %.2f-%.2f (fraction units), median %.2f"
          % (np.nanmin(rt), np.nanmax(rt), np.nanmedian(rt)))
    print("corr(RT, GC) = %+.3f over %d bins" % (r, ok.sum()))
    # Early replication is GC-rich. A weighted-mean index that rises with GC
    # would mean P2 is LATE, i.e. the deposited order runs late->early.
    if r < -0.2:
        print("  => low RT = EARLY replicating (P2_S is earliest). "
              "As deposited, index runs early -> late.")
        direction = "low_is_early"
    elif r > 0.2:
        print("  => low RT = LATE replicating (P2_S is latest). "
              "As deposited, index runs late -> early.")
        direction = "low_is_late"
    else:
        sys.exit("corr(RT,GC) = %.3f is too weak to fix the direction; "
                 "the profile is probably wrong" % r)

    q = np.nanpercentile(s_total[np.isfinite(s_total)], [1, 5, 25, 50])
    print("s_total percentiles 1/5/25/50: %s" % np.round(q, 4))

    np.savez(outp, rt=rt.astype(np.float32), s_total=s_total.astype(np.float32),
             g1=g1.astype(np.float32), binsize=a.bin,
             fractions=np.array(S_FRACTIONS), direction=np.array(direction),
             corr_rt_gc=np.array(r))
    print("wrote %s" % outp)


if __name__ == "__main__":
    main()
