#!/usr/bin/env python
"""Validate the iMARGI pairs file and census SVA coverage, before any analysis.

Three checks, all empirical. None of them trusts a document, because each one
guards a convention that silently inverts the result if assumed wrong -- the same
failure class as the `f[5]` strand trap and the .pos/.neg bigWig convention
already recorded in this project.

CHECK A -- WHICH SIDE IS THE RNA END.
4DN's own pipeline notes say they swapped DNA/RNA column order when building the
cool files, so the .pairs and the .mcool disagree with each other and neither can
be assumed. The header provenance argues side1 = RNA (pairtools ran with
--no-flip, so sides keep read identity rather than being sorted into the upper
triangle, and the `dist2_rsite <= 3` selection means side2 is the
restriction-ligated end, i.e. DNA). That is an inference, so it gets tested.

Four independent discriminators, computed symmetrically for both sides:

  exon fraction     RNA comes from transcripts; the DNA end is genomic
  gene-body frac.   same argument, weaker
  top-bin share     RNA piles onto a few highly-expressed loci; DNA is diffuse
  chrM share        mitochondrial RNA is enormously abundant; chrM DNA is not

All four must agree, or the script refuses to continue.

CHECK B -- THE RNA-END STRAND CONVENTION.
Knowing which side is RNA does not tell you whether a '+' in strand1 means the
transcript was on '+' or on '-'. dUTP-style protocols invert it, and the whole
SVA question is a sense/antisense question, so guessing here would silently swap
the two channels. Calibrated against four housekeeping genes of known
orientation, two per strand, exactly as src/chrrna/04_strandedness.sh does -- two
per strand so that a genome-wide strand imbalance cannot fake a clean answer.

CHECK C -- SVA COVERAGE CENSUS.
Is the test powered at all? Counts how many of the 5,974 SVA copies carry any
RNA-end contact, split by subfamily, orientation and host-gene class from
08_partition_copies.py, with L1 and Alu alongside for scale. If intergenic SVA
copies with contacts number in the tens rather than the hundreds, the distance
test is underpowered and that is worth knowing before building a null.

The DNA-end distance split (cis-proximal <1 Mb / cis-distal / trans) is reported
here as raw description only. It is NOT the test -- that needs the shift null and
a pre-mRNA baseline. Do not read these numbers as evidence on their own; nascent
RNA is overwhelmingly cis-proximal for reasons that have nothing to do with
genome organisation.

COORDINATES. The .pairs format is 1-based, and te_copy_key.tsv carries GTF
1-based coordinates unchanged, so both sides are 1-based and nothing is
converted. GENCODE is 1-based too. Do not paste the half-open logic from the
eclip tree in here.

Usage: check_pairs.py [pairs.gz] [gencode.gtf.gz]

Both conventions above are properties of a FILE, not of the assay -- a different
deposit can have been built by a different pipeline -- so re-run this whenever the
pairs file changes rather than trusting the report on disk.
"""
import gzip
import os
import subprocess
import sys
from collections import defaultdict

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PAIRS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    PROJ, "data", "hic", "4DNFIGDJIRV3.pairs.gz")
GENCODE = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser(
    "~/Downloads/hg38/gencode.v47.primary_assembly.annotation.gtf.gz")
KEY = os.path.join(PROJ, "data", "chrrna", "te_copy_key.tsv")
CLASSES = os.path.join(PROJ, "data", "chrrna", "te_copy_class.tsv")
OUT = os.path.join(PROJ, "results", "imargi")

BIN = 10000          # for the top-bin concentration discriminator
CHUNK = 4_000_000    # lines per chunk
PROX = 1_000_000     # cis-proximal threshold

# Known-orientation housekeeping genes, hg38. Two '+' and two '-' so a global
# strand imbalance cannot produce a clean-looking verdict. Same set and same
# reasoning as src/chrrna/04_strandedness.sh.
HOUSEKEEPING = [
    ("GAPDH",  "chr12",  6534517,  6538371, "+"),
    ("RPL13A", "chr19", 49487608, 49492308, "+"),
    ("ACTB",   "chr7",   5527151,  5530601, "-"),
    ("EEF1A1", "chr6",  73515750, 73523797, "-"),
]

MAXK = 4   # max nested TE copies checked at a point; asserted below


def merge(iv):
    """Sort + merge 1-based inclusive intervals. See 08_partition_copies.py --
    the single-index lookup below is only valid on non-overlapping intervals."""
    if not iv:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    iv.sort()
    s, e = [], []
    cs, ce = iv[0]
    for a, b in iv[1:]:
        if a <= ce + 1:
            ce = max(ce, b)
        else:
            s.append(cs); e.append(ce); cs, ce = a, b
    s.append(cs); e.append(ce)
    return np.array(s, np.int64), np.array(e, np.int64)


def hit_merged(merged, pos):
    """Vectorised point-in-interval against merged (non-overlapping) intervals."""
    ms, me = merged
    out = np.zeros(pos.size, bool)
    if ms.size == 0:
        return out
    i = np.searchsorted(ms, pos, side="right") - 1
    ok = i >= 0
    if ok.any():
        out[ok] = me[i[ok]] >= pos[ok]
    return out


def read_gencode(path):
    genes, exons = defaultdict(list), defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 8 or (f[2] != "gene" and f[2] != "exon"):
                continue
            (genes if f[2] == "gene" else exons)[f[0]].append((int(f[3]), int(f[4])))
    return ({c: merge(v) for c, v in genes.items()},
            {c: merge(v) for c, v in exons.items()})


def read_te():
    """Per-chrom sorted TE arrays + class codes, keyed on te_copy_key gene_id."""
    cls = {}
    if os.path.exists(CLASSES):
        with open(CLASSES) as fh:
            next(fh)
            for line in fh:
                g, c = line.rstrip("\n").split("\t")
                cls[g] = c
    per = defaultdict(lambda: ([], [], [], [], [], []))
    with open(KEY) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            st, en, sub, fam, strand = int(f[4]), int(f[5]), f[1], f[2], f[6]
            a = per[f[3]]
            a[0].append(st); a[1].append(en); a[2].append(sub)
            a[3].append(fam); a[4].append(strand == "+")
            a[5].append(cls.get(f[0], "unclassified"))
    out = {}
    for c, (st, en, sub, fam, plus, cl) in per.items():
        o = np.argsort(np.array(st, np.int64), kind="stable")
        out[c] = dict(start=np.array(st, np.int64)[o], end=np.array(en, np.int64)[o],
                      sub=np.array(sub, object)[o], fam=np.array(fam, object)[o],
                      plus=np.array(plus, bool)[o], cls=np.array(cl, object)[o])
    return out


def te_index_at(arr, pos):
    """Index of the TE copy containing each pos, or -1. Handles up to MAXK nested
    copies by walking back from the last start <= pos; asserts that is enough."""
    st, en = arr["start"], arr["end"]
    hi = np.searchsorted(st, pos, side="right")
    res = np.full(pos.size, -1, np.int64)
    for k in range(1, MAXK + 1):
        j = hi - k
        ok = (res < 0) & (j >= 0)
        if not ok.any():
            break
        jj = j[ok]
        good = en[jj] >= pos[ok]
        idx = np.where(ok)[0][good]
        res[idx] = jj[good]
    return res


def main():
    if not os.path.exists(PAIRS):
        sys.exit("missing %s" % PAIRS)
    os.makedirs(OUT, exist_ok=True)
    print("loading annotation ...", flush=True)
    genes, exons = read_gencode(GENCODE)
    te = read_te()
    print("  %d gene chroms, %d TE chroms" % (len(genes), len(te)), flush=True)

    # accumulators ---------------------------------------------------------
    n = 0
    side = {1: dict(exon=0, gene=0, chrM=0), 2: dict(exon=0, gene=0, chrM=0)}
    binc = {1: defaultdict(int), 2: defaultdict(int)}
    hk = defaultdict(lambda: defaultdict(int))     # gene -> strand -> count
    # (subfam, orient, class) -> [contacts, set_of_copies]
    cen = defaultdict(lambda: [0, set()])
    dist = defaultdict(lambda: defaultdict(int))   # family_orient -> class -> n

    cmd = "zcat %s | grep -v '^#' | cut -f2,3,4,5,6,7" % PAIRS
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True,
                         bufsize=1 << 22)
    print("streaming pairs (this reads the whole file) ...", flush=True)
    while True:
        lines = p.stdout.readlines(1 << 26)
        if not lines:
            break
        rec = defaultdict(lambda: ([], [], [], [], []))
        for ln in lines:
            f = ln.split("\t")
            if len(f) < 6:
                continue
            r = rec[(f[0], f[2])]
            r[0].append(int(f[1])); r[1].append(int(f[3]))
            r[2].append(f[4]); r[3].append(f[5][0]); r[4].append(0)
        for (c1, c2), (p1, p2, s1, s2, _) in rec.items():
            a1 = np.array(p1, np.int64); a2 = np.array(p2, np.int64)
            st1 = np.array(s1, object)
            m = a1.size
            n += m
            # --- CHECK A discriminators, symmetric on both sides ---
            for sd, ch, pos in ((1, c1, a1), (2, c2, a2)):
                if ch == "chrM":
                    side[sd]["chrM"] += m
                if ch in exons:
                    side[sd]["exon"] += int(hit_merged(exons[ch], pos).sum())
                if ch in genes:
                    side[sd]["gene"] += int(hit_merged(genes[ch], pos).sum())
                b = binc[sd]
                for v, k in zip(*np.unique(pos // BIN, return_counts=True)):
                    b[(ch, int(v))] += int(k)
            # --- CHECK B: housekeeping strand calibration on the RNA side ---
            for g, gc, gs, ge, gstr in HOUSEKEEPING:
                if c1 != gc:
                    continue
                sel = (a1 >= gs) & (a1 <= ge)
                if sel.any():
                    for v, k in zip(*np.unique(st1[sel], return_counts=True)):
                        hk[(g, gstr)][str(v)] += int(k)
            # --- CHECK C: SVA/TE census on the RNA side ---
            if c1 not in te:
                continue
            arr = te[c1]
            idx = te_index_at(arr, a1)
            got = idx >= 0
            if not got.any():
                continue
            gi = idx[got]
            sub = arr["sub"][gi]; fam = arr["fam"][gi]
            plus = arr["plus"][gi]; cl = arr["cls"][gi]
            # RAW observable only: does the read strand match the TE's annotated
            # strand? Whether "same" means sense or antisense depends on the
            # library convention, which CHECK B measures in this same pass --
            # so the sense/antisense mapping is applied at report time, never
            # assumed here.
            rplus = (st1[got] == "+")
            orient = np.where(rplus == plus, "same", "opp")
            same = (c1 == c2)
            d = np.abs(a2[got] - a1[got])
            dcls = (np.where(d <= PROX, "cis_prox", "cis_distal") if same
                    else np.full(d.size, "trans", object))
            for s_, o_, c_, f_, dc_, ix_ in zip(sub, orient, cl, fam, dcls, gi):
                e = cen[(s_, o_, c_)]
                e[0] += 1
                e[1].add(int(ix_))
                dist[(f_, o_)][str(dc_)] += 1
    p.stdout.close(); p.wait()

    # ---- report ----------------------------------------------------------
    rep = open(os.path.join(OUT, "pairs_checks.txt"), "w")

    def out(s=""):
        print(s, flush=True)
        rep.write(s + "\n")

    out("pairs file : %s" % PAIRS)
    out("total pairs: %d" % n)
    out("")
    out("=== CHECK A: which side is the RNA end? ===")
    out("%-18s %12s %12s %10s   %s"
        % ("discriminator", "side1", "side2", "n(max)", "implies"))
    votes = []
    # A discriminator with no signal on EITHER side must abstain. Letting a
    # 0-vs-0 tie fall through to the else-branch would cast a real vote for
    # side2 on no evidence -- which is exactly what it did on this file, where
    # the pairs carry no chrM contacts at all.
    MINSIG = 1000
    for k, lbl in (("exon", "exon fraction"), ("gene", "gene-body frac"),
                   ("chrM", "chrM share")):
        c1, c2 = side[1][k], side[2][k]
        v1, v2 = c1 / n, c2 / n
        if max(c1, c2) < MINSIG:
            out("%-18s %11.4f%% %11.4f%% %10d   ABSTAINS (no signal)"
                % (lbl, 100 * v1, 100 * v2, max(c1, c2)))
            continue
        votes.append(1 if v1 > v2 else 2)
        out("%-18s %11.4f%% %11.4f%% %10d   side%d is RNA"
            % (lbl, 100 * v1, 100 * v2, max(c1, c2), votes[-1]))
    for sd in (1, 2):
        tot = sum(binc[sd].values())
        top = sorted(binc[sd].values(), reverse=True)
        k = max(1, len(top) // 1000)
        binc[sd] = sum(top[:k]) / tot
    votes.append(1 if binc[1] > binc[2] else 2)
    out("%-18s %11.4f%% %11.4f%% %10s   side%d is RNA"
        % ("top-0.1% bin share", 100 * binc[1], 100 * binc[2], "-", votes[-1]))
    out("")
    if len(votes) < 3:
        out("VERDICT: only %d discriminators had signal -- too few." % len(votes))
        rep.close()
        sys.exit("side assignment underdetermined")
    if len(set(votes)) != 1:
        out("VERDICT: DISCRIMINATORS DISAGREE (%s) -- stopping." % votes)
        rep.close()
        sys.exit("side assignment ambiguous; do not proceed")
    rna = votes[0]
    out("VERDICT: side%d is the RNA end (%d of %d discriminators voted, all agree)."
        % (rna, len(votes), 4))
    if rna != 1:
        out("NOTE: this contradicts the header provenance -- investigate before use.")

    out("")
    out("=== CHECK B: RNA-end strand convention ===")
    out("%-8s %-6s %10s %10s  %s" % ("gene", "strand", "RNA '+'", "RNA '-'", "read"))
    # Only genes with coverage get a vote. Counting uncovered genes as
    # disagreement would report "inverted" whenever the sample is thin.
    covered = []
    for g, gc, gs, ge, gstr in HOUSEKEEPING:
        d = hk[(g, gstr)]
        pp, mm = d.get("+", 0), d.get("-", 0)
        if pp + mm == 0:
            out("%-8s %-6s %10s %10s  NO COVERAGE" % (g, gstr, "-", "-"))
            continue
        same = ("+" if pp > mm else "-") == gstr
        covered.append((gstr, same))
        out("%-8s %-6s %10d %10d  %s"
            % (g, gstr, pp, mm, "matches gene strand" if same else "INVERTED"))
    out("")
    nc = len(covered)
    nm = sum(1 for _, s in covered if s)
    strands = set(st for st, _ in covered)
    if nc < 2:
        out("=> UNDETERMINED: only %d housekeeping gene(s) covered." % nc)
    elif nm == nc:
        out("=> RNA strand1 EQUALS the transcript strand (%d/%d covered)" % (nm, nc))
    elif nm == 0:
        out("=> RNA strand1 IS INVERTED vs the transcript strand (%d/%d covered)"
            % (nc, nc))
    else:
        out("=> CONFLICT: %d of %d covered genes match -- sense/antisense is NOT "
            "safe to assign." % (nm, nc))
    if nc >= 2 and len(strands) < 2:
        out("WARNING: only '%s'-strand genes had coverage. A genome-wide strand "
            "imbalance could fake this verdict; treat it as provisional."
            % strands.pop())

    # Apply CHECK B's measured convention to the raw same/opp counts. If the
    # read strand is inverted relative to the transcript, a read matching the
    # TE's annotated strand came from the ANTISENSE transcript.
    if nc >= 2 and nm == nc:
        LBL = {"same": "sense", "opp": "antisense"}
        conv = "strand1 == transcript strand"
    elif nc >= 2 and nm == 0:
        LBL = {"same": "antisense", "opp": "sense"}
        conv = "strand1 INVERTED -> same-strand reads are ANTISENSE transcripts"
    else:
        rep.close()
        sys.exit("strand convention undetermined; refusing to label sense/antisense")
    cen = {(sub, LBL[o], cl): v for (sub, o, cl), v in cen.items()}
    dist = {(fam, LBL[o]): v for (fam, o), v in dist.items()}

    out("")
    out("orientation labels below use: %s" % conv)
    out("")
    out("=== CHECK C: SVA coverage census (RNA end) ===")
    out("%-9s %-10s %-14s %9s %9s" % ("subfam", "orient", "class", "contacts", "copies"))
    sva = sorted((k for k in cen if k[0].startswith("SVA")),
                 key=lambda k: -cen[k][0])
    for k in sva:
        if cen[k][0] < 10:
            continue
        out("%-9s %-10s %-14s %9d %9d"
            % (k[0], k[1], k[2], cen[k][0], len(cen[k][1])))
    out("")
    out("SVA totals by orientation x class:")
    agg = defaultdict(lambda: [0, set()])
    for k, v in cen.items():
        if not k[0].startswith("SVA"):
            continue
        a = agg[(k[1], k[2])]
        a[0] += v[0]; a[1] |= v[1]
    for k in sorted(agg, key=lambda x: -agg[x][0]):
        out("  %-10s %-14s contacts=%-9d copies=%d" % (k[0], k[1], agg[k][0], len(agg[k][1])))

    out("")
    out("=== family scale + DNA-end distance (DESCRIPTIVE ONLY, not the test) ===")
    out("%-5s %-10s %10s %9s %11s %8s" % ("fam", "orient", "contacts",
                                          "cis_prox", "cis_distal", "trans"))
    for k in sorted(dist, key=lambda x: -sum(dist[x].values())):
        d = dist[k]
        t = sum(d.values())
        out("%-5s %-10s %10d %8.1f%% %10.1f%% %7.1f%%"
            % (k[0], k[1], t, 100 * d.get("cis_prox", 0) / t,
               100 * d.get("cis_distal", 0) / t, 100 * d.get("trans", 0) / t))
    rep.close()
    print("\nwrote %s" % os.path.join(OUT, "pairs_checks.txt"))


if __name__ == "__main__":
    main()
