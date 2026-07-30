#!/usr/bin/env python
"""Is intronic TE RNA more distal-acting than the rest of its own host gene?

WHY THIS EXISTS. trans_test.py runs on intergenic copies only, because a TE
inside a gene inherits the host's contact profile wholesale. That throws away
the 1.22 M intronic copies -- 2x the intergenic set, and in the SVA census the
intronic channel carried 116 k contacts against intergenic's 8 k. This recovers
them WITHOUT the confound, by changing the control rather than the filter.

THE CONTROL IS THE SAME GENE. For every protein-coding gene, its intron space
splits into TE-covered bases and TE-free bases. Both sit in the same transcript,
the same compartment, at the same expression level, in the same replication
timing domain. Comparing them is paired by construction, so host expression --
which drives distal fraction from 0.195 to 0.817 across deciles and is by far
the largest confound here -- cancels instead of needing to be matched.

That is why the toroidal shift null of trans_test.py CANNOT simply be pointed at
intronic copies. It moves anchors to random genomic positions, mostly outside the
host gene, so gene-body membership is not preserved. Gene bodies at high
expression sit at 0.817 against a genome-wide null of ~0.66: pooling classes
against that null would manufacture large positive z from membership alone.

THE CHANNEL. For an intronic copy the read channel running WITH the host strand
is host pre-mRNA passing through, co-directionally, swamping whatever the copy
itself produces. Only the HOST-OPPOSITE strand is interpretable, so only it is
counted, in both the TE and control regions. (src/chrrna quantified that
asymmetry at 49x, but in K562 and on a different assay -- the argument here is
mechanical and does not rest on that figure.) Since strand1 is inverted on this
library (check_pairs.py), host-opposite transcript means strand1 == the gene's
strand; that is asserted below on housekeeping genes, not assumed.

Which TE orientation that corresponds to follows from the copy's class:

    intronic_sense copy (TE strand == gene strand)  -> host-opposite is ANTISENSE
    intronic_anti  copy (TE strand != gene strand)  -> host-opposite is SENSE

so every copy contributes to exactly one (subfamily, orientation) cell and the
label is derived, never chosen.

THE CONTROL EXCLUDES ALL TEs, not just the subfamily under test. Introns are
TE-dense; a control that only removed the tested subfamily would be largely other
L1/Alu/SVA copies and the comparison would be TE-vs-TE. Copies of every indexed
subfamily are subtracted.

TWO REFERENCES, both within-gene:

1. THE GENE'S OWN TE-FREE RATE. Per gene i, expected TE-region numerator is
   n_i * p_i with p_i the control fraction in that same gene. Summed over genes
   this is a stratified (Mantel-Haenszel-style) statistic; z uses the binomial
   variance sum. Controls the host completely, but not WHERE in the gene the
   copies sit.

2. A SHIFT NULL INSIDE THE GENE'S OWN INTRON SPACE. The gene's introns are
   concatenated into a flat coordinate, the gene's TE footprint is shifted
   rigidly by one random offset modulo that length, and the count recomputed.
   Same construction as enrich_permutation.py, with the eligible set narrowed
   from the chromosome to one gene's introns. This adds control for position
   within the gene -- distance to TSS, intron rank -- which reference 1 lacks.

Contacts are placed in flat intron coordinates ONCE per gene, so a shift is two
searchsorted calls on that gene's own array; recomputing genomic overlaps per
shift would be hopeless.

Usage: within_intron_test.py [--shifts N] [--statistic S] [--min-mapq N]
"""
import argparse
import os
import sys
from collections import defaultdict
from math import erfc, sqrt

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trans_test import (ARR, CLS, GENCODE, HOUSEKEEPING, KEY, OUT,  # noqa: E402
                        STRAND1_INVERTED, merge, read_gencode_introns,
                        reduce_stat)

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Classes whose host context is a single protein-coding gene. exonic_* is mature
# mRNA containing the element, ambiguous has no single host to control for, and
# no_annotation was deliberately excluded upstream rather than folded in.
KEEP_CLASSES = ("intronic_sense", "intronic_anti")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--shifts", type=int, default=200)
    p.add_argument("--statistic", choices=("distal", "cisdistal", "trans"),
                   default="cisdistal",
                   help="default cisdistal: this is a spreading question, and "
                        "pooling trans in makes it a trans question -- see the "
                        "STATISTIC note in trans_test.py")
    p.add_argument("--min-mapq", type=int, default=0)
    p.add_argument("--min-contacts", type=int, default=200,
                   help="minimum TE-region contacts (statistic denominator) "
                        "for a subfamily x orientation to be reported")
    p.add_argument("--min-ctl", type=int, default=20,
                   help="minimum control contacts for a gene to contribute. "
                        "p_i from a handful of reads is noise, and it enters "
                        "the expectation directly.")
    p.add_argument("--arr", default=ARR)
    p.add_argument("--tag", default="")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def read_gene_introns_with_span():
    """gid -> (chrom, gene_start, gene_end, plus, [merged intron intervals]).

    read_gencode_introns() drops gene identity (it returns per-chromosome
    intervals for the pre-mRNA baseline), and identity is the whole point here,
    so the GTF is walked again keeping it.
    """
    import gzip
    span, ex = {}, defaultdict(list)
    with gzip.open(GENCODE, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            if 'gene_type "protein_coding"' not in f[8]:
                continue
            i = f[8].find('gene_id "')
            if i < 0:
                continue
            gid = f[8][i + 9:f[8].find('"', i + 9)]
            if f[2] == "gene":
                span[gid] = (f[0], int(f[3]), int(f[4]), f[6] == "+")
            else:
                ex[gid].append((int(f[3]), int(f[4])))
    out = {}
    for gid, (c, s, e, plus) in span.items():
        introns, prev = [], s
        for a, b in merge(ex.get(gid, [])):
            if a > prev:
                introns.append((prev, a - 1))
            prev = max(prev, b + 1)
        if prev <= e:
            introns.append((prev, e))
        if introns:
            out[gid] = (c, s, e, plus, introns)
    return out


def read_copies():
    """Intronic TE copies and, separately, ALL copies (for control masking)."""
    cls = {}
    with open(CLS) as fh:
        next(fh)
        for ln in fh:
            g, c = ln.rstrip("\n").split("\t")
            cls[g] = c
    per_chrom_all = defaultdict(list)          # chrom -> (start, end)
    per_chrom_te = defaultdict(list)           # chrom -> (start, end, sub, plus, klass)
    with open(KEY) as fh:
        next(fh)
        for ln in fh:
            g, sub, fam, c, s, e, st = ln.rstrip("\n").split("\t")
            s, e = int(s), int(e)
            per_chrom_all[c].append((s, e))
            k = cls.get(g)
            if k in KEEP_CLASSES:
                per_chrom_te[c].append((s, e, "%s|%s" % (fam, sub), st == "+", k))
    return per_chrom_all, per_chrom_te


def flat_map(introns):
    """Merged 1-based inclusive introns -> (starts, ends, cumulative offsets, L).

    Flat coordinate of a genomic position p inside intron j is
    off[j] + (p - starts[j]). Half-open in flat space, [0, L).
    """
    st = np.array([a for a, _ in introns], np.int64)
    en = np.array([b for _, b in introns], np.int64)
    ln = en - st + 1
    off = np.concatenate([[0], np.cumsum(ln)])
    return st, en, off, int(off[-1])


def to_flat(pos, st, en, off):
    """Genomic -> flat, for positions known to lie inside the intron set."""
    j = np.searchsorted(st, pos, side="right") - 1
    ok = (j >= 0) & (pos <= en[np.clip(j, 0, en.size - 1)])
    return off[np.clip(j, 0, st.size - 1)] + (pos - st[np.clip(j, 0, st.size - 1)]), ok


def segs_to_flat(s, e, st, en, off):
    """Clip genomic intervals to the intron set -> flat half-open segments.

    Returns (flat_start, flat_end, src) where src indexes back into s/e: an
    interval straddling an exon yields several flat segments, so the mapping is
    one-to-many and the caller needs to know which input each came from.
    """
    if s.size == 0:
        z = np.empty(0, np.int64)
        return z, z, z
    a = np.maximum(s[:, None], st[None, :])
    b = np.minimum(e[:, None], en[None, :])
    ii, jj = np.nonzero(a <= b)
    if ii.size == 0:
        z = np.empty(0, np.int64)
        return z, z, z
    fa = off[jj] + (a[ii, jj] - st[jj])
    fb = off[jj] + (b[ii, jj] - st[jj]) + 1
    return fa, fb, ii


def complement(segs, L):
    """[0,L) minus a set of half-open flat segments, merged first."""
    if segs[0].size == 0:
        return np.array([0], np.int64), np.array([L], np.int64)
    o = np.argsort(segs[0])
    s, e = segs[0][o], segs[1][o]
    ms, me = [s[0]], [e[0]]
    for a, b in zip(s[1:], e[1:]):
        if a <= me[-1]:
            me[-1] = max(me[-1], b)
        else:
            ms.append(a)
            me.append(b)
    ms, me = np.array(ms, np.int64), np.array(me, np.int64)
    cs = np.concatenate([[0], me])
    ce = np.concatenate([ms, [L]])
    m = cs < ce
    return cs[m], ce[m]


class FlatGene:
    """Prefix sums over dclass for one gene's host-opposite contacts, in flat
    intron coordinates. Counting any flat interval is then two searchsorted."""

    def __init__(self, fpos, dclass, L):
        o = np.argsort(fpos, kind="stable")
        self.pos = fpos[o]
        d = dclass[o]
        self.L = L
        self.cum = np.zeros((3, self.pos.size + 1), np.int32)
        for c in range(3):
            np.cumsum(d == c, dtype=np.int32, out=self.cum[c, 1:])

    def count_each(self, a, b):
        """(3, n) per-interval totals over half-open flat [a,b)."""
        lo = np.searchsorted(self.pos, a, side="left")
        hi = np.searchsorted(self.pos, b, side="left")
        return self.cum[:, hi] - self.cum[:, lo]

    def count(self, a, b):
        """(3,) totals summed over all intervals."""
        return self.count_each(a, b).sum(1)

    def count_wrap_each(self, a, ln):
        """(3, n) per-interval, for intervals that may run past L and wrap.

        Per-interval rather than summed so every shift can be counted in ONE
        searchsorted pass: the caller stacks (n_shifts x n_segments) starts,
        counts them all at once and reshapes. Looping shifts in Python instead
        made 200 shifts cost hours.
        """
        b = a + ln
        w = b > self.L
        c = self.count_each(a, np.where(w, self.L, b))
        if w.any():
            tail = self.count_each(np.zeros(int(w.sum()), np.int64),
                                   b[w] - self.L)
            c = c.copy()
            c[:, w] += tail
        return c


def main():
    a = parse_args()
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    genes = read_gene_introns_with_span()
    all_copies, te_copies = read_copies()
    print("%d protein-coding genes with introns" % len(genes), flush=True)

    by_chrom = defaultdict(list)
    for gid, (c, s, e, plus, introns) in genes.items():
        by_chrom[c].append((gid, s, e, plus, introns))

    subs = sorted({t[2] for v in te_copies.values() for t in v})
    sidx = {s: i for i, s in enumerate(subs)}
    ns = len(subs)
    print("%d subfamilies among intronic copies, %d shifts"
          % (ns, a.shifts), flush=True)

    # obs/null/ctl indexed [subfamily, orientation, dclass];
    # orientation 0 = sense, 1 = antisense (TE-relative, derived from class).
    obs = np.zeros((ns, 2, 3), np.int64)
    null = np.zeros((a.shifts, ns, 2, 3), np.int64)
    # Per (subfamily, orientation): accumulate the stratified expectation.
    exp_n = np.zeros((ns, 2))
    exp_v = np.zeros((ns, 2))
    ngene = defaultdict(set)
    ncopy = defaultdict(int)
    hk = defaultdict(lambda: [0, 0])

    chroms = sorted(f[:-4] for f in os.listdir(a.arr) if f.endswith(".npz"))
    for c in chroms:
        if c not in by_chrom:
            continue
        z = np.load(os.path.join(a.arr, "%s.npz" % c))
        pos, strand, dcl = z["pos"], z["strand"], z["dclass"]
        if a.min_mapq:
            if "mapq" not in z:
                sys.exit("--min-mapq needs mapq; re-run extract_rna_ends.py")
            k = z["mapq"] >= a.min_mapq
            pos, strand, dcl = pos[k], strand[k], dcl[k]
        pos = pos.astype(np.int64)

        for g, gc, gs, ge, gstr in HOUSEKEEPING:
            if gc != c:
                continue
            lo = np.searchsorted(pos, gs, "left")
            hi = np.searchsorted(pos, ge, "right")
            hk[g][0] = int((strand[lo:hi] == 1).sum())
            hk[g][1] = int((strand[lo:hi] == 0).sum())

        tes = te_copies.get(c, [])
        tstart = np.array([t[0] for t in tes], np.int64)
        torder = np.argsort(tstart)
        tstart = tstart[torder]
        tend = np.array([t[1] for t in tes], np.int64)[torder]
        tsub = [tes[i][2] for i in torder]
        tplus = np.array([t[3] for t in tes], bool)[torder]
        tcls = [tes[i][4] for i in torder]

        acp = all_copies.get(c, [])
        astart = np.array(sorted(x[0] for x in acp), np.int64)
        aend = np.array([x[1] for x in sorted(acp)], np.int64)

        for gid, gs2, ge2, plus, introns in by_chrom[c]:
            ist, ien, ioff, L = flat_map(introns)
            if L < 1000:
                continue
            # Host-opposite channel: strand1 == gene strand (library inverted).
            want = 1 if plus else 0
            lo = np.searchsorted(pos, ist[0], "left")
            hi = np.searchsorted(pos, ien[-1], "right")
            if hi - lo < a.min_ctl:
                continue
            p_, s_, d_ = pos[lo:hi], strand[lo:hi], dcl[lo:hi]
            sel = s_ == want
            if sel.sum() < a.min_ctl:
                continue
            fp, ok = to_flat(p_[sel], ist, ien, ioff)
            if ok.sum() < a.min_ctl:
                continue
            fg = FlatGene(fp[ok], d_[sel][ok], L)

            # TE copies of this gene, and ALL copies for the control mask.
            i0 = np.searchsorted(tstart, ist[0], "left")
            i1 = np.searchsorted(tstart, ien[-1], "right")
            j0 = np.searchsorted(astart, ist[0], "left")
            j1 = np.searchsorted(astart, ien[-1], "right")
            if i1 <= i0:
                continue
            mfa, mfb, _ = segs_to_flat(astart[j0:j1], aend[j0:j1],
                                       ist, ien, ioff[:-1])
            cs, ce = complement((mfa, mfb), L)
            ctl = fg.count(cs, ce)
            cn, cd = reduce_stat(ctl.astype(np.int64), a.statistic)
            if cd < a.min_ctl:
                continue
            p_i = cn / cd

            # All this gene's TE copies mapped to flat coords in one shot, then
            # grouped by (subfamily, derived orientation) via the src index.
            gfa, gfb, gsrc = segs_to_flat(tstart[i0:i1], tend[i0:i1],
                                          ist, ien, ioff[:-1])
            if gfa.size == 0:
                continue
            grp = defaultdict(list)
            for n_, k in enumerate(gsrc):
                oi = 1 if tcls[i0 + k] == "intronic_sense" else 0
                grp[(tsub[i0 + k], oi)].append(n_)
            for k in np.unique(gsrc):
                oi = 1 if tcls[i0 + k] == "intronic_sense" else 0
                ncopy[(tsub[i0 + k], oi)] += 1

            for (sub, oi), rowsn in grp.items():
                si = sidx[sub]
                sel_ = np.array(rowsn, np.int64)
                fs, fe = gfa[sel_], gfb[sel_]
                t = fg.count(fs, fe).astype(np.int64)
                obs[si, oi] += t
                tn, td = reduce_stat(t, a.statistic)
                if td:
                    ngene[(sub, oi)].add(gid)
                    exp_n[si, oi] += td * p_i
                    # Two variance terms, and dropping the second inflates z
                    # badly. The first is the binomial spread of the TE-region
                    # count given p_i. The second is the uncertainty in p_i
                    # ITSELF: it is estimated from that gene's control contacts,
                    # and at cd_i ~ 20 with p ~ 0.04 its standard error is about
                    # as large as p. Treating an estimated rate as known turns
                    # control noise into apparent signal.
                    exp_v[si, oi] += (td * p_i * (1 - p_i)
                                      + td * td * p_i * (1 - p_i) / cd)
                # One offset per gene per shift: preserves the spacing of this
                # gene's copies, same rigid construction as the genome-wide null.
                ln = fe - fs
                offs = rng.integers(0, L, size=a.shifts)
                nsx = ((fs[None, :] + offs[:, None]) % L).ravel()
                cw = fg.count_wrap_each(nsx, np.tile(ln, a.shifts))
                null[:, si, oi] += cw.reshape(3, a.shifts, ln.size).sum(2).T
        print("  %s done" % c, flush=True)

    GSTR = dict(GAPDH="+", RPL13A="+", ACTB="-", EEF1A1="-")
    covered = {g: v for g, v in hk.items() if sum(v) > 0}
    bad = [g for g, (gp, gm) in covered.items() if (gp > gm) == (GSTR[g] == "+")]
    if len(covered) < 2:
        sys.exit("strand assertion impossible: %d genes covered" % len(covered))
    if bad:
        sys.exit("strand assertion FAILED for %s -- host-opposite channel would "
                 "be the host channel; every number here would invert" % bad)
    print("strand assertion ok (strand1 inverted; %d/%d housekeeping genes)"
          % (len(covered), len(GSTR)))

    rows = []
    for sub, si in sidx.items():
        fam, name = sub.split("|")
        for oi in (0, 1):
            o = obs[si, oi]
            on, od = reduce_stat(o, a.statistic)
            if od < a.min_contacts:
                continue
            f = on / od
            e, v = exp_n[si, oi], exp_v[si, oi]
            z_gene = (on - e) / sqrt(v) if v > 0 else float("nan")
            nl = null[:, si, oi, :]
            nn, nd = reduce_stat(nl, a.statistic)
            ok = nd > 0
            nf = nn[ok] / nd[ok]
            mu = float(nf.mean()) if ok.any() else float("nan")
            sd = float(nf.std(ddof=1)) if ok.sum() > 1 else float("nan")
            z_shift = (f - mu) / sd if sd and sd == sd and sd > 0 else float("nan")
            rows.append(dict(
                subfamily=name, family=fam,
                orientation=("sense", "antisense")[oi],
                n_genes=len(ngene[(sub, oi)]), n_copies=ncopy[(sub, oi)],
                n_contacts=int(od), te_frac=f,
                ctl_exp_frac=(e / od) if od else float("nan"), z_gene=z_gene,
                shift_mean=mu, shift_sd=sd, z_shift=z_shift))

    for key in ("z_gene", "z_shift"):
        pk = "p_" + key.split("_")[1]
        for r in rows:
            zz = r[key]
            r[pk] = 0.5 * erfc(zz / sqrt(2)) if zz == zz else float("nan")
        ps = sorted((r[pk], i) for i, r in enumerate(rows) if r[pk] == r[pk])
        m, prev = len(ps), 1.0
        for rank in range(m - 1, -1, -1):
            p, i = ps[rank]
            prev = min(prev, p * m / (rank + 1))
            rows[i]["q_" + key.split("_")[1]] = prev

    rows.sort(key=lambda r: -r["z_gene"] if r["z_gene"] == r["z_gene"] else 1e9)
    cols = ["subfamily", "family", "orientation", "n_genes", "n_copies",
            "n_contacts", "te_frac", "ctl_exp_frac", "z_gene", "p_gene",
            "q_gene", "shift_mean", "shift_sd", "z_shift", "p_shift", "q_shift"]
    suffix = a.statistic + ("_mapq%d" % a.min_mapq if a.min_mapq else "") + a.tag
    path = os.path.join(OUT, "within_intron_%s.tsv" % suffix)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(("%.6g" % r[c]) if isinstance(r[c], float)
                               else str(r.get(c, "")) for c in cols) + "\n")
    print("\nwrote %s  (%d rows)" % (path, len(rows)))

    print("\ntop by z_gene (TE intron bases vs TE-free bases of the SAME gene):")
    print("  %-10s %-10s %-4s %7s %8s %9s %8s %8s %8s %9s %8s"
          % ("subfam", "orient", "fam", "genes", "copies", "contacts",
             "te_frac", "ctl_exp", "z_gene", "q_gene", "z_shift"))
    for r in rows[:20]:
        print("  %-10s %-10s %-4s %7d %8d %9d %8.3f %8.3f %8.2f %9.2g %8.2f"
              % (r["subfamily"], r["orientation"], r["family"], r["n_genes"],
                 r["n_copies"], r["n_contacts"], r["te_frac"], r["ctl_exp_frac"],
                 r["z_gene"], r["q_gene"], r["z_shift"]))


if __name__ == "__main__":
    main()
