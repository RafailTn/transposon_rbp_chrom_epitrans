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

STATISTIC -- AND WHY THE DEFAULT IS NOT THE TEST FOR SPREADING. --statistic
distal pools cis-distal AND trans into one fraction. That was the first thing
run here and it answers "is this RNA non-local at all", but it is the wrong
question for the architectural hypothesis people actually mean. Xist, Airn and
Kcnq1ot1 all act in CIS: they spread along the chromosome they were transcribed
from and do essentially nothing in trans. Trans-acting architectural RNA is the
rare exception, not the model.

That matters numerically, not just conceptually. Over TE-anchored contacts this
library is 56.9% trans, 9.6% cis-distal, 33.5% cis-proximal, so the pooled
statistic is roughly six parts trans to one part cis-distal: a cis-spreading
signal can move a long way
without shifting it, and the channel doing the shifting is the one where ambient
ligation lives. --statistic cisdistal conditions on the RNA having stayed on its
own chromosome -- cis-distal / (cis-proximal + cis-distal) -- which targets the
spreading hypothesis and drops the noise-dominated compartment from the
denominator entirely.

The shifts do not depend on the statistic, so the raw dclass tallies are dumped
to tallies_<class>.npz and any reduction can be re-derived without rerunning.

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
PAIRS = os.path.join(PROJ, "data", "hic", "4DNFIGDJIRV3.pairs.gz")
KEY = os.path.join(PROJ, "data", "chrrna", "te_copy_key.tsv")
CLS = os.path.join(PROJ, "data", "chrrna", "te_copy_class.tsv")
GENCODE = os.path.expanduser(
    "~/Downloads/hg38/gencode.v47.primary_assembly.annotation.gtf.gz")
OUT = os.path.join(PROJ, "results", "imargi")

# Measured by check_pairs.py on this file: strand1 is inverted, so an RNA end
# whose alignment strand MATCHES a TE's annotated strand came from the antisense
# transcript. Asserted below against housekeeping genes.
STRAND1_INVERTED = True
STAT = "distal"                 # set from --statistic in main()
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
    p.add_argument("--arr", default=ARR,
                   help="directory of per-chromosome RNA-end .npz. Point at a "
                        "per-replicate extraction to split the data.")
    p.add_argument("--pairs", default=PAIRS,
                   help="pairs file the --arr extraction came from. Read for its "
                        "#chromsize header ONLY, which the toroidal shift needs; "
                        "must be the same file, or the null shifts against "
                        "lengths the arrays do not share.")
    p.add_argument("--tag", default="",
                   help="suffix for output filenames, so runs on different "
                        "--arr do not overwrite each other.")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="drop RNA ends whose own MAPQ is below this. The 4DN "
                        "pipeline kept MAPQ>=1, so mismapped repeat reads are "
                        "present by design; raising this is the test for "
                        "whether a distal signal is really mismapping. Costs "
                        "depth on exactly the young subfamilies it protects.")
    p.add_argument("--statistic", choices=("distal", "cisdistal", "trans"),
                   default="distal",
                   help="which reduction of the dclass triple to test. "
                        "distal = (cis-distal+trans)/all, the pooled default. "
                        "cisdistal = cis-distal/(cis-prox+cis-distal) -- the "
                        "CIS-RESTRICTED statistic, which is the right one for "
                        "a Xist-like spreading hypothesis: it asks how far the "
                        "RNA reached GIVEN that it stayed on its own "
                        "chromosome, and drops the trans channel entirely. "
                        "trans = trans/all. See STATISTIC in the module "
                        "docstring for why the default is not the right test "
                        "for cis spreading.")
    return p.parse_args()


def reduce_stat(t, stat=None):
    """(..., 3) tally over dclass -> (numerator, denominator) for `stat`.

    `stat` defaults to the module-level STAT set from --statistic. It is an
    explicit argument so within_intron_test.py can import this rather than
    reimplement it -- the reduction must agree across every script here, for
    the same reason the searchsorted overlap idiom must.

    dclass is 0 cis-proximal (<=1 Mb), 1 cis-distal, 2 trans. Every statistic
    here is one ratio of sums over that axis, and NOTHING upstream of this
    function depends on which -- obs, null and the pre-mRNA baseline all carry
    the full triple. That is deliberate: the pooled default turned out to be
    the wrong question (see docstring), and re-deriving a different one must
    not cost another hour-long run.
    """
    stat = STAT if stat is None else stat
    if stat == "distal":
        return t[..., 1] + t[..., 2], t.sum(-1)
    if stat == "trans":
        return t[..., 2], t.sum(-1)
    # cisdistal: the trans channel is excluded from the denominator, not just
    # the numerator. Trans is 61% of this library and is where ambient
    # ligation lives, so conditioning on cis both targets the hypothesis and
    # removes the dominant noise compartment.
    return t[..., 1], t[..., 0] + t[..., 1]


def chromsizes(pairs=PAIRS):
    """From the pairs header -- same file the arrays came from, so they cannot
    disagree about chromosome length, which the toroidal shift depends on."""
    out = subprocess.run("zcat %s | head -400 | grep '^#chromsize'" % pairs,
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

    def __init__(self, path, min_mapq=0):
        z = np.load(path)
        pos = z["pos"].astype(np.int64)
        strand, dclass = z["strand"], z["dclass"]
        if min_mapq:
            if "mapq" not in z:
                sys.exit("--min-mapq needs mapq in the arrays; re-run "
                         "extract_rna_ends.py (it now stores mapq/mapq_dna)")
            # Filter on the RNA end's MAPQ. A mismapped RNA end puts the anchor
            # in the wrong place while its DNA end stays correct, which reads as
            # a distal contact -- see extract_rna_ends.py. Dropping rows before
            # the cumsum keeps `pos` sorted, which counts() depends on.
            keep = z["mapq"] >= min_mapq
            pos, strand, dclass = pos[keep], strand[keep], dclass[keep]
        self.pos = pos
        code = strand.astype(np.int64) * 3 + dclass.astype(np.int64)
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
    global STAT
    a = parse_args()
    STAT = a.statistic
    os.makedirs(OUT, exist_ok=True)
    sizes = chromsizes(a.pairs)
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

    chroms = sorted(f[:-4] for f in os.listdir(a.arr) if f.endswith(".npz"))
    for c in chroms:
        if c not in sizes:
            continue
        ch = Chrom(os.path.join(a.arr, "%s.npz" % c), a.min_mapq)
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
            # Stored as the full dclass triple, reduced at report time like
            # everything else, so --statistic re-selects the baseline too. A
            # baseline reduced on a different statistic than the observations
            # would be a silent unit mismatch.
            # cnt is (6, n_introns): sum over this gene's introns, keep the
            # dclass axis. One row per gene, so the decile stratification below
            # still has per-source abundance to cut on.
            tri3 = cnt[row:row + 3].sum(1).astype(np.int64)
            if tri3.sum():
                base.append(tri3)

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

    # The shifts are the whole cost of this script (~1 h at --shifts 1000) and
    # they do not depend on the statistic. Dump the raw dclass tallies so any
    # other reduction is a seconds-long re-read rather than a rerun. Written
    # per class, not per statistic -- one file serves all three.
    tally_path = os.path.join(OUT, "tallies_%s%s%s.npz" % (
        a.klass, "_mapq%d" % a.min_mapq if a.min_mapq else "", a.tag))
    np.savez_compressed(
        tally_path, obs=obs, null=null, base=np.array(base, np.int64),
        subs=np.array(subs), labels=np.array(LBL),
        ncopy=np.array([len(ncopy[s]) for s in subs]),
        shifts=a.shifts, seed=a.seed)
    print("wrote %s (obs/null dclass tallies, statistic-independent)"
          % tally_path, flush=True)

    # --- pre-mRNA baseline by abundance decile ---
    base = np.array(base, np.int64)
    bnum, bden = reduce_stat(base)
    # Deciles are cut on the statistic's OWN denominator, not on total
    # contacts: for --statistic cisdistal the trans contacts are not part of
    # the measurement, so a gene with a million trans and ten cis contacts
    # belongs in a low decile, not the top one.
    base, bnum, bden = base[bden > 0], bnum[bden > 0], bden[bden > 0]
    q = np.quantile(bden, np.linspace(0, 1, 11))
    dec = []
    for i in range(10):
        m = (bden >= q[i]) & (bden <= q[i + 1])
        if m.sum():
            dec.append((q[i], q[i + 1], bnum[m].sum() / bden[m].sum(), int(m.sum())))

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
            onum, oden = reduce_stat(o)
            n = int(oden)
            if n < a.min_contacts:
                continue
            f = onum / n
            nl = null[:, i, oi, :]
            nnum, nden = reduce_stat(nl)
            ok = nden > 0
            nf = nnum[ok] / nden[ok]
            mu, sd = float(nf.mean()), float(nf.std(ddof=1))
            z = (f - mu) / sd if sd > 0 else float("nan")
            pe = (1 + int((nf >= f).sum())) / (len(nf) + 1)
            rows.append(dict(subfamily=name, family=fam, orientation=LBL[oi],
                             n_copies=len(ncopy[s]), n_contacts=n,
                             n_all=int(o.sum()), statistic=STAT,
                             stat_frac=f, null_mean=mu, null_sd=sd, z=z,
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
            "n_all", "statistic", "stat_frac", "null_mean", "null_sd", "z",
            "p_emp", "p_z", "q", "premrna_ref"]
    # n_contacts is the STATISTIC's denominator, n_all every contact on the
    # copies; they differ for --statistic cisdistal, where trans is excluded.
    suffix = a.klass if STAT == "distal" else "%s_%s" % (a.klass, STAT)
    if a.min_mapq:
        suffix += "_mapq%d" % a.min_mapq
    suffix += a.tag
    path = os.path.join(OUT, "trans_test_%s.tsv" % suffix)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                ("%.6g" % r[c]) if isinstance(r[c], float) else str(r.get(c, ""))
                for c in cols) + "\n")

    print("\npre-mRNA baseline, %s fraction by decile of the statistic's own\n  denominator:" % STAT)
    for lo, hi, f, k in dec:
        print("  n in [%8d,%8d]  %s=%.3f  (%d genes)" % (lo, hi, STAT, f, k))
    print("\nwrote %s  (%d rows)" % (path, len(rows)))
    print("\ntop by z (%s more than position predicts):" % STAT)
    print("  %-10s %-10s %-5s %8s %9s %8s %8s %9s %9s"
          % ("subfam", "orient", "fam", "copies", "contacts", STAT[:8],
             "null_mu", "z", "q"))
    for r in rows[:15]:
        print("  %-10s %-10s %-5s %8d %9d %8.3f %8.3f %9.2f %9.2g"
              % (r["subfamily"], r["orientation"], r["family"], r["n_copies"],
                 r["n_contacts"], r["stat_frac"], r["null_mean"], r["z"],
                 r.get("q", float("nan"))))


if __name__ == "__main__":
    main()
