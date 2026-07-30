#!/usr/bin/env python
"""Per-20kb-bin covariates and mask for the TSA-seq nuclear-position stage.

Cached because it depends only on the assembly and the rmsk table, never on the
TSA-seq tracks or the TE subset under test. Rerun only if BIN changes.

WHY THESE THREE.

gc_frac -- THE confound. Alu is concentrated in GC-rich sequence and L1 in
    AT-rich sequence, which has been known since the 1980s, and nuclear
    compartments are themselves organised along that axis: speckle-proximal
    chromatin is GC-rich and gene-dense, lamina-proximal chromatin is not. So
    "Alu is near speckles" is predicted a priori by base composition alone and
    is not a finding about nuclear organisation. Any claim from this stage has
    to survive conditioning on GC, which is what this column exists for.

n_frac -- assembly gaps. GRCh38 leaves the acrocentric short arms and other
    heterochromatin unassembled. A bin that is mostly N has no sequence for
    reads to map to, so its TSA-seq value is meaningless rather than low.

sat_frac -- satellite density from the FULL rmsk table (repClass == Satellite),
    not the L1/Alu/SVA-filtered GTF. GRCh38 *models* centromeric alpha-satellite
    rather than N-masking it, so n_frac does not catch centromeres, and measured
    nucleolar TSA-seq runs ~1.4 log2 units high there -- plausibly real
    (pericentromeric heterochromatin does sit at the nucleolar periphery) but
    computed from few uniquely-mapping reads either way. Masked rather than
    interpreted.

Usage: bin_covariates.py [--bin N] [--force]
"""
import argparse
import gzip
import os
import subprocess
import sys

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FA = os.path.expanduser("~/Downloads/hg38/GRCh38.primary_assembly.genome.fa")
FAI = FA + ".fai"
RMSK = os.path.join(PROJ, "data", "hg38_rmsk.gtf.gz")
CACHE = os.path.join(PROJ, "cache")

CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bin", type=int, default=20000,
                   help="must match the TSA-seq track binning (20 kb)")
    p.add_argument("--force", action="store_true")
    return p.parse_args()


def chrom_sizes():
    d = {}
    for ln in open(FAI):
        f = ln.split("\t")
        if f[0] in CHROMS:
            d[f[0]] = int(f[1])
    missing = [c for c in CHROMS if c not in d]
    if missing:
        sys.exit("fai lacks %s" % missing)
    return d


def gc_and_n(sizes, binsize):
    """Stream the fasta once, accumulating per-bin G+C, A+T and N counts.

    Counted with str.count on the raw line rather than a per-base loop; the
    fasta is ~3 GB and anything per-base is unusably slow.
    """
    gc = {c: np.zeros(sizes[c] // binsize + 1, np.int64) for c in sizes}
    at = {c: np.zeros(sizes[c] // binsize + 1, np.int64) for c in sizes}
    nn = {c: np.zeros(sizes[c] // binsize + 1, np.int64) for c in sizes}
    cur, pos = None, 0
    with open(FA) as fh:
        for line in fh:
            if line[0] == ">":
                cur = line[1:].split()[0]
                pos = 0
                print("  %s%s" % (cur, "" if cur in sizes else " (skipped)"),
                      flush=True)
                continue
            if cur not in gc:
                continue
            s = line.strip()
            if not s:
                continue
            # A line (60 bp) can straddle a bin boundary. Split it there rather
            # than assigning the whole line to one bin: at 20 kb bins the error
            # would be small but it is free to be exact.
            start = pos
            pos += len(s)
            b0, b1 = start // binsize, (pos - 1) // binsize
            if b0 == b1:
                pieces = [(b0, s)]
            else:
                pieces, off = [], start
                while off < pos:
                    b = off // binsize
                    end = min((b + 1) * binsize, pos)
                    pieces.append((b, s[off - start:end - start]))
                    off = end
            for b, sub in pieces:
                u = sub.upper()
                g = u.count("G") + u.count("C")
                a = u.count("A") + u.count("T")
                gc[cur][b] += g
                at[cur][b] += a
                nn[cur][b] += len(u) - g - a
    return gc, at, nn


def satellite(sizes, binsize):
    """Satellite bp per bin from the raw rmsk table (0-based half-open)."""
    sat = {c: np.zeros(sizes[c] // binsize + 1, np.int64) for c in sizes}
    # awk is ~4x faster than parsing 152 MB of gzip in python for a 3-column
    # projection, and this file has no quoting to worry about.
    cmd = ("zcat %s | awk -F'\\t' 'NR>1 && $12==\"Satellite\" "
           "{print $6\"\\t\"$7\"\\t\"$8}'" % RMSK)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 22)
    n = 0
    for ln in p.stdout:
        c, s, e = ln.split("\t")
        if c not in sat:
            continue
        s, e = int(s), int(e)
        b0, b1 = s // binsize, (e - 1) // binsize
        if b0 == b1:
            sat[c][b0] += e - s
        else:
            for b in range(b0, min(b1 + 1, sat[c].size)):
                lo = max(s, b * binsize)
                hi = min(e, (b + 1) * binsize)
                sat[c][b] += hi - lo
        n += 1
    p.stdout.close()
    p.wait()
    print("  %d satellite intervals" % n, flush=True)
    return sat


def main():
    a = parse_args()
    os.makedirs(CACHE, exist_ok=True)
    out = os.path.join(CACHE, "bin_covariates_%d.npz" % a.bin)
    if os.path.exists(out) and not a.force:
        print("%s exists; --force to rebuild" % out)
        return
    for f in (FA, FAI, RMSK):
        if not os.path.exists(f):
            sys.exit("missing %s" % f)

    sizes = chrom_sizes()
    print("scanning fasta for GC/N ...", flush=True)
    gc, at, nn = gc_and_n(sizes, a.bin)
    print("scanning rmsk for satellite ...", flush=True)
    sat = satellite(sizes, a.bin)

    # One flat array per covariate, with an offset table, so downstream code can
    # index a (chrom, pos) pair with a single add -- the same layout the
    # permutation code uses for its flat genome coordinate.
    order = [c for c in CHROMS]
    nb = {c: sizes[c] // a.bin + 1 for c in order}
    off, t = {}, 0
    for c in order:
        off[c] = t
        t += nb[c]
    gcf = np.zeros(t, np.float32)
    nf = np.zeros(t, np.float32)
    sf = np.zeros(t, np.float32)
    for c in order:
        s = slice(off[c], off[c] + nb[c])
        tot = gc[c] + at[c] + nn[c]
        with np.errstate(invalid="ignore", divide="ignore"):
            gcf[s] = np.where(gc[c] + at[c] > 0, gc[c] / np.maximum(gc[c] + at[c], 1), np.nan)
            nf[s] = np.where(tot > 0, nn[c] / np.maximum(tot, 1), 1.0)
        sf[s] = sat[c] / float(a.bin)

    np.savez(out, chroms=np.array(order), nbins=np.array([nb[c] for c in order]),
             offsets=np.array([off[c] for c in order]), binsize=a.bin,
             gc_frac=gcf, n_frac=nf, sat_frac=sf)
    ok = np.isfinite(gcf) & (nf < 0.10) & (sf < 0.10)
    print("\n%d bins total" % t)
    print("  usable (n_frac<0.10, sat_frac<0.10, GC defined): %d (%.1f%%)"
          % (ok.sum(), 100.0 * ok.mean()))
    print("  gap-dominated  n_frac>=0.10 : %d" % (nf >= 0.10).sum())
    print("  satellite-rich sat_frac>=0.10: %d" % (sf >= 0.10).sum())
    print("  GC over usable bins: mean %.3f  range %.3f-%.3f"
          % (np.nanmean(gcf[ok]), np.nanmin(gcf[ok]), np.nanmax(gcf[ok])))
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
