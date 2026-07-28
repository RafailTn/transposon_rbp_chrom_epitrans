#!/usr/bin/env python
"""Do TE-derived RNAs contact DNA at distance more than their position predicts?

THE QUESTION. An architectural RNA acts somewhere other than where it was
transcribed. Chromatin retention (src/chrrna) and m6A competence (src/mintseq)
are necessary but very weak evidence, because pervasive nascent transcription has
both properties. This asks the discriminating question directly: of the contacts
made by RNA emanating from a TE subfamily, what fraction reach distal cis or
trans, and is that more than expected?

"More than expected" needs two separate references, because they control
different things and neither is sufficient alone:

1. THE SHIFT NULL controls for POSITION. TE copies are not randomly placed --
   Alu is isochore-biased into gene-dense open sequence, L1 the reverse -- and
   contact geometry varies enormously across the genome. A rigid toroidal shift
   of each subfamily's copies, one offset per chromosome, preserves copy spacing
   and asks: if these copies sat elsewhere, would RNA from those positions look
   equally distal? Same construction as enrich_permutation.py and
   te_tss_contacts.py.

2. THE PRE-mRNA BASELINE controls for the SCALE of the answer. Nascent
   pre-mRNA is the canonical locally-tethered RNA: transcribed, chromatin-bound,
   and doing nothing architectural. It defines what "acts locally" looks like in
   this library. Introns of protein-coding genes, sense strand.

THE CONFOUND THAT MATTERS MOST. A high trans fraction is also what NOISE looks
like -- ambient and random ligations are trans-dominated, so any low-abundance
source drifts toward trans for purely technical reasons. Comparing SVA antisense
(8 k contacts) against pre-mRNA (millions) unstratified would be measuring
abundance, not biology. The baseline is therefore reported per contact-count
decile, and each TE category is compared against the decile it belongs to.

INTERGENIC COPIES ONLY, by default. A TE inside a gene inherits the host's
contact profile wholesale; see 08_partition_copies.py.

STRAND. check_pairs.py measured that strand1 is INVERTED relative to transcript
strand on this library. That is re-derived here from the stored raw strand and
asserted, never assumed -- getting it backwards swaps sense and antisense, which
is the entire question.

Usage: trans_test.py [--shifts N] [--class CLASS] [--min-contacts N]
"""
import argparse
import gzip
import json
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ARR = os.path.join(PROJ, "data", "imargi")
PAIRS = os.path.join(PROJ, "data", "hic", "4DNFIVIHUHOE.pairs.gz")
KEY = os.path.join(PROJ, "data", "chrrna", "te_copy_key.tsv")
CLS = os.path.join(PROJ, "data", "chrrna", "te_copy_class.tsv")
GENCODE = os.path.expanduser(
    "~/Downloads/hg38/gencode.v47.primary_assembly.annotation.gtf.gz")
OUT = os.path.join(PROJ, "results", "imargi")

# Measured by check_pairs.py on this file: strand1 is inverted, so an RNA end
# whose alignment strand MATCHES a TE's annotated strand came from the antisense
# transcript. Asserted below against housekeeping genes.
STRAND1_INVERTED = True
HOUSEKEEPING = [("GAPDH", "chr12", 6534517, 6538371, "+"),
                ("RPL13A", "chr19", 49487608, 49492308, "+"),
                ("ACTB", "chr7", 5527151, 5530601, "-"),
                ("EEF1A1", "chr6", 73515750, 73523797, "-")]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--shifts", type=int, default=1000)
    p.add_argument("--class", dest="klass", default="intergenic")
    p.add_argument("--min-contacts", type=int, default=200)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def chromsizes():
    """From the pairs header -- same file the arrays came from, so they cannot
    disagree about chromosome length, which the toroidal shift depends on."""
    out = subprocess.run("zcat %s | head -400 | grep '^#chromsize'" % PAIRS,
                         shell=True, capture_output=True, text=True).stdout
    d = {}
    for ln in out.splitlines():
        f = ln.split()
        if len(f) >= 3:
            d[f[1]] = int(f[2])
    if not d:
        sys.exit("no #chromsize lines in the pairs header")
    return d


def merge(iv):
    if not iv:
        return []
    iv.sort()
    o = [list(iv[0])]
    for a, b in iv[1:]:
        if a <= o[-1][1] + 1:
            o[-1][1] = max(o[-1][1], b)
        else:
            o.append([a, b])
    return o


def read_gencode_introns():
    """Per-gene intron intervals for protein-coding genes: gene span minus that
    gene's own merged exons. Per gene, not globally merged, so the baseline keeps
    per-source contact counts and can be stratified by abundance."""
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
            a = f[8]
            i = a.find('gene_id "')
            if i < 0:
                continue
            gid = a[i + 9:a.find('"', i + 9)]
            if f[2] == "gene":
                span[gid] = (f[0], int(f[3]), int(f[4]), f[6])
            else:
                ex[gid].append((int(f[3]), int(f[4])))
    out = defaultdict(list)          # chrom -> (start, end, plus)
    for gid, (c, s, e, st) in span.items():
        introns, prev = [], s
        for a, b in merge(ex.get(gid, [])):
            if a > prev:
                introns.append((prev, a - 1))
            prev = max(prev, b + 1)
        if prev <= e:
            introns.append((prev, e))
        if introns:
            out[c].append((introns, st == "+"))
    return out


def read_anchors(klass):
    keep = {}
    with open(CLS) as fh:
        next(fh)
        for line in fh:
            g, c = line.rstrip("\n").split("\t")
            keep[g] = c
    per = defaultdict(lambda: ([], [], [], []))
    with open(KEY) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or keep.get(f[0]) != klass:
                continue
            a = per[f[3]]
            a[0].append(int(f[4])); a[1].append(int(f[5]))
            a[2].append(f[6] == "+"); a[3].append("%s|%s" % (f[2], f[1]))
    return per


class Chrom:
    """Prefix sums over (read_strand, dclass), so counting contacts inside any
    interval is two array lookups instead of a scan. This is what makes 1000
    shifts affordable."""

    def __init__(self, path):
        z = np.load(path)
        self.pos = z["pos"].astype(np.int64)
        code = z["strand"].astype(np.int64) * 3 + z["dclass"].astype(np.int64)
        # int32: counts per chromosome max out around 2e7, and int64 here would
        # cost ~1 GB on chr1 for no benefit.
        self.cum = np.zeros((6, self.pos.size + 1), np.int32)
        for c in range(6):
            np.cumsum(code == c, dtype=np.int32, out=self.cum[c, 1:])

    def counts(self, start, end):
        """-> (6, n_intervals) contact counts, one column per interval."""
        lo = np.searchsorted(self.pos, start, side="left")
        hi = np.searchsorted(self.pos, end, side="right")
        return self.cum[:, hi] - self.cum[:, lo]


def orient_rows(plus):
    """Row indices into the 6 codes giving (same-strand, opp-strand) triples.

    code = read_strand*3 + dclass, read_strand 1 == '+'. For a '+' anchor the
    same-strand reads are codes 3..5; for a '-' anchor they are 0..2.
    """
    same = np.where(plus, 3, 0)
    opp = np.where(plus, 0, 3)
    return same, opp


def tally(cnt, plus, groups, ngroup):
    """cnt (6,n) -> per-group (ngroup, 2, 3): [same|opp] x dclass."""
    same, opp = orient_rows(plus)
    out = np.zeros((ngroup, 2, 3), np.int64)
    idx = np.arange(cnt.shape[1])
    for d in range(3):
        s = cnt[same + d, idx]
        o = cnt[opp + d, idx]
        out[:, 0, d] = np.bincount(groups, weights=s, minlength=ngroup).astype(np.int64)
        out[:, 1, d] = np.bincount(groups, weights=o, minlength=ngroup).astype(np.int64)
    return out


def main():
    a = parse_args()
    os.makedirs(OUT, exist_ok=True)
    sizes = chromsizes()
    anchors = read_anchors(a.klass)
    if not anchors:
        sys.exit("no copies in class %s" % a.klass)
    introns = read_gencode_introns()

    subs = sorted({s for v in anchors.values() for s in v[3]})
    sidx = {s: i for i, s in enumerate(subs)}
    ns = len(subs)
    print("class=%s  %d subfamilies  %d shifts" % (a.klass, ns, a.shifts), flush=True)

    obs = np.zeros((ns, 2, 3), np.int64)
    null = np.zeros((a.shifts, ns, 2, 3), np.int64)
    ncopy = defaultdict(set)
    hk = defaultdict(lambda: [0, 0])
    base = []                       # (n_contacts, n_distal) per protein-coding gene
    rng = np.random.default_rng(a.seed)

    chroms = sorted(f[:-4] for f in os.listdir(ARR) if f.endswith(".npz"))
    for c in chroms:
        if c not in sizes:
            continue
        ch = Chrom(os.path.join(ARR, "%s.npz" % c))
        L = sizes[c]

        # --- strand assertion, on housekeeping genes ---
        for g, gc, gs, ge, gstr in HOUSEKEEPING:
            if gc != c:
                continue
            cnt = ch.counts(np.array([gs]), np.array([ge]))[:, 0]
            hk[g][0] = int(cnt[3:].sum())      # read '+'
            hk[g][1] = int(cnt[:3].sum())      # read '-'

        # --- pre-mRNA baseline: introns of protein-coding genes, sense ---
        for iv, plus in introns.get(c, []):
            s = np.array([x[0] for x in iv], np.int64)
            e = np.array([x[1] for x in iv], np.int64)
            cnt = ch.counts(s, e)
            same = 3 if plus else 0
            # Sense transcript = strand1 opposite (library is inverted).
            row = (0 if plus else 3) if STRAND1_INVERTED else same
            tot = int(cnt[row:row + 3].sum())
            dis = int(cnt[row + 1:row + 3].sum())
            if tot:
                base.append((tot, dis))

        st, en, plus, sub = anchors[c]
        st = np.array(st, np.int64); en = np.array(en, np.int64)
        plus = np.array(plus, bool)
        grp = np.array([sidx[s] for s in sub], np.int64)
        for s_, i_ in zip(sub, range(len(sub))):
            ncopy[s_].add((c, int(st[i_])))

        obs += tally(ch.counts(st, en), plus, grp, ns)

        ln = en - st
        for k in range(a.shifts):
            off = rng.integers(0, L)
            ns_ = (st - 1 + off) % L + 1
            ne_ = ns_ + ln
            wrap = ne_ > L
            ne_ = np.where(wrap, L, ne_)          # clip the few that wrap
            null[k] += tally(ch.counts(ns_, ne_), plus, grp, ns)
        del ch
        print("  %s done" % c, flush=True)

    # --- strand assertion ---
    GSTR = dict(GAPDH="+", RPL13A="+", ACTB="-", EEF1A1="-")
    # Genes with no coverage cannot vote; counting them would fail every
    # '-'-strand gene by default (0 > 0 is False).
    covered = {g: v for g, v in hk.items() if sum(v) > 0}
    bad = [g for g, (gp, gm) in covered.items()
           if (gp > gm) == (GSTR[g] == "+")]
    if len(covered) < 2:
        sys.exit("strand assertion impossible: %d housekeeping genes covered"
                 % len(covered))
    if bad:
        sys.exit("strand assertion FAILED for %s -- the convention differs from "
                 "what check_pairs.py measured; do not trust sense/antisense" % bad)
    print("strand assertion ok (strand1 inverted; %d/%d housekeeping genes agree)"
          % (len(covered), len(GSTR)))

    LBL = ("antisense", "sense") if STRAND1_INVERTED else ("sense", "antisense")

    # --- pre-mRNA baseline by abundance decile ---
    base = np.array(sorted(base), np.int64)
    q = np.quantile(base[:, 0], np.linspace(0, 1, 11))
    dec = []
    for i in range(10):
        m = (base[:, 0] >= q[i]) & (base[:, 0] <= q[i + 1])
        if m.sum():
            dec.append((q[i], q[i + 1], base[m, 1].sum() / base[m, 0].sum(), int(m.sum())))

    def ref_for(n):
        for lo, hi, f, _ in dec:
            if lo <= n <= hi:
                return f
        return dec[-1][2]

    rows = []
    for s, i in sidx.items():
        fam, name = s.split("|")
        for oi in (0, 1):
            o = obs[i, oi]
            n = int(o.sum())
            if n < a.min_contacts:
                continue
            f = (o[1] + o[2]) / n
            nl = null[:, i, oi, :]
            tot = nl.sum(1)
            ok = tot > 0
            nf = (nl[ok, 1] + nl[ok, 2]) / tot[ok]
            mu, sd = float(nf.mean()), float(nf.std(ddof=1))
            z = (f - mu) / sd if sd > 0 else float("nan")
            pe = (1 + int((nf >= f).sum())) / (len(nf) + 1)
            rows.append(dict(subfamily=name, family=fam, orientation=LBL[oi],
                             n_copies=len(ncopy[s]), n_contacts=n,
                             distal_frac=f, null_mean=mu, null_sd=sd, z=z,
                             p_emp=pe, premrna_ref=ref_for(n)))
    from math import erfc, sqrt
    for r in rows:
        r["p_z"] = 0.5 * erfc(r["z"] / sqrt(2)) if r["z"] == r["z"] else float("nan")
    ps = sorted((r["p_z"], i) for i, r in enumerate(rows) if r["p_z"] == r["p_z"])
    m = len(ps)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        p, i = ps[rank]
        prev = min(prev, p * m / (rank + 1))
        rows[i]["q"] = prev

    rows.sort(key=lambda r: -r["z"] if r["z"] == r["z"] else 1e9)
    cols = ["subfamily", "family", "orientation", "n_copies", "n_contacts",
            "distal_frac", "null_mean", "null_sd", "z", "p_emp", "p_z", "q",
            "premrna_ref"]
    path = os.path.join(OUT, "trans_test_%s.tsv" % a.klass)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                ("%.6g" % r[c]) if isinstance(r[c], float) else str(r.get(c, ""))
                for c in cols) + "\n")

    print("\npre-mRNA baseline, distal fraction by contact-count decile:")
    for lo, hi, f, k in dec:
        print("  n in [%8d,%8d]  distal=%.3f  (%d genes)" % (lo, hi, f, k))
    print("\nwrote %s  (%d rows)" % (path, len(rows)))
    print("\ntop by z (distal more than position predicts):")
    print("  %-10s %-10s %-5s %8s %9s %8s %8s %9s %9s"
          % ("subfam", "orient", "fam", "copies", "contacts", "distal",
             "null_mu", "z", "q"))
    for r in rows[:15]:
        print("  %-10s %-10s %-5s %8d %9d %8.3f %8.3f %9.2f %9.2g"
              % (r["subfamily"], r["orientation"], r["family"], r["n_copies"],
                 r["n_contacts"], r["distal_frac"], r["null_mean"], r["z"],
                 r.get("q", float("nan"))))


if __name__ == "__main__":
    main()
