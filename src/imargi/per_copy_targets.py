#!/usr/bin/env python
"""Does any single TE copy's RNA hit a RECURRENT set of target loci?

THE QUESTION THE AGGREGATE TESTS CANNOT ASK. trans_test.py and
within_intron_test.py both reduce a subfamily to one fraction, pooling hundreds
of copies and every target locus. An RNA acting on a few specific loci -- the
Xist-like case -- barely moves such a fraction. But the Xist analogy applied to a
repeat family only makes sense per COPY: Xist is one locus, one transcript, one
chromosome, whereas "SVA RNA" is thousands of loci each producing very little.
So the testable form is: does one particular copy behave like an architectural
lncRNA gene? That is this script.

THE STATISTIC IS CONCENTRATION, NOT REACH. Per copy and orientation, DNA ends are
binned at --bin (default 1 Mb) and summarised by the Herfindahl index
H = sum(share_i^2), whose reciprocal 1/H is the effective number of target bins.
H is near 1/n when contacts scatter and rises toward 1 when they pile onto few
loci. top1_share is reported alongside because H alone is hard to read.

Cis contacts within --prox (default 1 Mb) of the source are DROPPED. They are
dominated by the copy's own locus and by transcription in place, which is
tethering, not targeting -- and they are also where the local decay is steepest
and hardest to model.

THE NULL PRESERVES REACH AND RANDOMISES DESTINATION. For each copy the null
redraws the same number of contacts, keeping:

    - the cis/trans split
    - the full distance profile, stratified log2 (1-2, 2-4, 4-8 ... Mb)

and resampling only WHICH bin within each stratum, weighted by the marginal
DNA-end density of the whole TE-anchored library. So distance decay -- by far the
largest structure in the data, and entirely a property of the source's position
rather than of any targeting -- is held fixed by construction rather than
modelled, and accessibility enters through the weights. What is left to vary is
exactly the alternative hypothesis: recurrence on specific loci.

BOTH MAPQs ARE FILTERED. mapq (RNA end) fixes which copy the RNA is attributed
to; mapq_dna (DNA end) fixes where the target is. This screen depends on the
target, so the DNA end matters here in a way it did not for the distance tests.
Repeat-derived DNA ends mismapping onto a few loci would be indistinguishable
from targeting -- see the README for what low MAPQ did to AluYa5.

MASKING IS SYMMETRIC, AND THE SOURCE SIDE IS THE ONE THAT MATTERS. A bin excluded
as a target is also excluded as a source. Amplified regions are the reason: copies
inside one are copy-number inflated AND physically cis on the derivative chromosome
while the reference scores their contacts as trans, so they manufacture exactly the
recurrence this script looks for. Masking targets alone does not help -- the
artifact is in where the RNA came FROM.

H1 is karyotypically normal and this was verified on the data, not assumed: 1 Mb
DNA-end density tops out at 2.91x the median non-empty bin (99.9th percentile
2.55x), so --exclude defaults to empty and --mask-mult 3.5 masks nothing. Both
knobs are live for a rearranged line; see their help text.

DENSITY THRESHOLDING IS NOT A SUBSTITUTE FOR NAMING THE KARYOTYPE. The two failure
modes squeeze from both sides: a cutoff loose enough to spare ordinary
highly-transcribed bins can sit *above* an amplicon's actual density and mask
nothing at all, while a cutoff tight enough to catch it also eats gene-dense
regions. Karyotype is a property of the cell line, not of the data, so name the
coordinates in --exclude and leave --mask-mult as a backstop.

More importantly, on a line whose CNV is not fully characterised the uncertainty is
NOT confined to the bins you masked. Copy-number inflation and derivative-chromosome
topology act directly on this assay's readout -- a genomic distance between two
ligated ends -- so a coordinate mask handles only the part you can see. Prefer a
karyotypically normal line for per-locus work; this project uses H1 for exactly
that reason.

Usage: per_copy_targets.py [--min-contacts N] [--draws N] [--min-mapq N]
"""
import argparse
import os
import sys
from collections import defaultdict
from math import erfc, sqrt

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CON = os.path.join(PROJ, "data", "imargi_contacts")
KEY = os.path.join(PROJ, "data", "chrrna", "te_copy_key.tsv")
CLS = os.path.join(PROJ, "data", "chrrna", "te_copy_class.tsv")
OUT = os.path.join(PROJ, "results", "imargi")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--families", default="Alu,L1")
    p.add_argument("--bin", type=int, default=1_000_000)
    p.add_argument("--prox", type=int, default=1_000_000,
                   help="drop cis contacts within this of the source")
    p.add_argument("--min-contacts", type=int, default=50,
                   help="minimum contacts AFTER the proximal cut; below ~50 the "
                        "Herfindahl index is dominated by its own 1/n floor")
    p.add_argument("--min-mapq", type=int, default=30)
    p.add_argument("--draws", type=int, default=200)
    p.add_argument("--mask-mult", type=float, default=3.5,
                   help="mask bins whose DNA-end density exceeds this multiple of "
                        "the median non-empty bin. A backstop only -- on H1 it "
                        "fires on nothing, since the maximum bin is 2.91x. Do NOT "
                        "tighten it to H1's own 99th percentile (2.22x): the top "
                        "bins here are ordinary gene-dense regions (chr19, "
                        "chr3:47-49 Mb) and masking them removes real biology. "
                        "0 disables.")
    p.add_argument("--exclude", default="",
                   help="explicit bins to drop, as chrom:start-end, comma "
                        "separated, applied to sources and targets alike. Empty "
                        "by default because H1 is karyotypically normal -- "
                        "verified on the DNA-end density profile, where the "
                        "amplicon signature of a rearranged line (sharp multi-Mb "
                        "plateaus at several times median, abrupt drop-offs) is "
                        "absent. Karyotype is a property of the cell line, not of "
                        "the data, so a rearranged line needs its coordinates "
                        "NAMED here rather than caught by --mask-mult.")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--tag", default="")
    p.add_argument("--con", default=CON,
                   help="directory of extract_contacts.py output. Point at a "
                        "per-biorep extraction to split the data.")
    return p.parse_args()


def main():
    global CON
    a = parse_args()
    CON = a.con
    os.makedirs(OUT, exist_ok=True)
    fams = set(a.families.split(","))
    rng = np.random.default_rng(a.seed)
    chroms = [str(x) for x in np.load(os.path.join(CON, "chroms.npy"))]
    ci = {c: i for i, c in enumerate(chroms)}

    key = [l.rstrip("\n").split("\t") for l in open(KEY).readlines()[1:]]
    k_sub = np.array([r[1] for r in key])
    k_fam = np.array([r[2] for r in key])
    k_chr = np.array([r[3] for r in key])
    k_st = np.array([int(r[4]) for r in key], np.int64)
    k_plus = np.array([r[6] == "+" for r in key])
    kcls = {}
    for ln in open(CLS).readlines()[1:]:
        g, c = ln.rstrip("\n").split("\t")
        kcls[g] = c
    k_cls = np.array([kcls.get(r[0], "?") for r in key])

    # --- bin layout, and the marginal DNA-end weight map ---
    sizes = {}
    for c in chroms:
        p = os.path.join(CON, "%s.npz" % c)
        if os.path.exists(p):
            z = np.load(p)
            for j, cc in enumerate(chroms):
                m = z["dna_chrom"] == j
                if m.any():
                    sizes[cc] = max(sizes.get(cc, 0), int(z["dna_pos"][m].max()))
    nb = {c: sizes.get(c, 0) // a.bin + 1 for c in chroms}
    off = {}
    t = 0
    for c in chroms:
        off[c] = t
        t += nb[c]
    NB = t
    print("%d bins of %d bp over %d chromosomes" % (NB, a.bin, len(chroms)),
          flush=True)

    W = np.zeros(NB, np.float64)
    for c in chroms:
        p = os.path.join(CON, "%s.npz" % c)
        if not os.path.exists(p):
            continue
        z = np.load(p)
        m = (z["mapq"] >= a.min_mapq) & (z["mapq_dna"] >= a.min_mapq)
        gb = np.array([off[chroms[j]] for j in range(len(chroms))])[z["dna_chrom"][m]]
        W += np.bincount(gb + z["dna_pos"][m] // a.bin, minlength=NB)
    W /= W.sum()
    nz = W[W > 0]
    if a.mask_mult > 0:
        thr = a.mask_mult * float(np.median(nz))
        MASK = W > thr
        print("masking %d/%d bins above %.1fx median density (%.2f%% of DNA ends)"
              % (MASK.sum(), NB, a.mask_mult, 100 * W[MASK].sum()))
    else:
        MASK = np.zeros(NB, bool)
    if a.exclude.strip():
        nex = 0
        for spec in a.exclude.split(","):
            c, rng_ = spec.strip().split(":")
            lo_, hi_ = (int(x) for x in rng_.split("-"))
            if c not in off:
                print("  !! --exclude names %s, not in the contact map" % c,
                      file=sys.stderr)
                continue
            b0 = off[c] + lo_ // a.bin
            b1 = off[c] + min(hi_ // a.bin, nb[c] - 1)
            nex += int((~MASK[b0:b1 + 1]).sum())
            MASK[b0:b1 + 1] = True
        print("excluding %d further bins by --exclude (%s)" % (nex, a.exclude))
    print("masked %d/%d bins in total, holding %.2f%% of DNA ends"
          % (MASK.sum(), NB, 100 * W[MASK].sum()))
    # Renormalise over surviving bins: masked bins are removed from the null's
    # candidate set, so they must not keep probability mass.
    W = np.where(MASK, 0.0, W)
    W /= W.sum()

    # chromosome id per global bin, and bin centre, for stratified sampling
    binchr = np.zeros(NB, np.int32)
    bincen = np.zeros(NB, np.int64)
    for c in chroms:
        s = slice(off[c], off[c] + nb[c])
        binchr[s] = ci[c]
        bincen[s] = (np.arange(nb[c], dtype=np.int64) + 0.5) * a.bin
    EDGES = np.array([a.prox] + [a.prox * 2 ** k for k in range(1, 10)], np.int64)

    rows = []
    n_src_masked = 0
    for c in chroms:
        p = os.path.join(CON, "%s.npz" % c)
        if not os.path.exists(p):
            continue
        z = np.load(p)
        m = (z["mapq"] >= a.min_mapq) & (z["mapq_dna"] >= a.min_mapq)
        cp = z["copy"][m]
        keep = np.isin(k_fam[cp], list(fams))
        cp = cp[keep]
        if cp.size == 0:
            continue
        rs = z["rna_strand"][m][keep]
        dch = z["dna_chrom"][m][keep].astype(np.int64)
        dps = z["dna_pos"][m][keep].astype(np.int64)
        # strand1 == the copy's annotated strand -> ANTISENSE transcript
        anti = rs == k_plus[cp]
        gbin = np.array([off[chroms[j]] for j in range(len(chroms))])[dch] \
            + dps // a.bin
        this = ci[c]

        order = np.lexsort((anti, cp))
        cp, dch, dps, anti, gbin = (cp[order], dch[order], dps[order],
                                    anti[order], gbin[order])
        gid = cp.astype(np.int64) * 2 + anti
        bnd = np.flatnonzero(np.diff(gid)) + 1
        starts = np.concatenate([[0], bnd])
        ends = np.concatenate([bnd, [gid.size]])

        cbins = np.flatnonzero((binchr == this) & ~MASK)
        tbins = np.flatnonzero((binchr != this) & ~MASK)
        wt = W[tbins]
        wt_c = np.cumsum(wt)
        if wt_c[-1] <= 0:
            continue

        for s0, e0 in zip(starts, ends):
            n0 = e0 - s0
            if n0 < a.min_contacts:
                continue
            i = cp[s0]
            src = k_st[i]
            # SOURCE-SIDE MASK. An earlier version applied MASK to target bins
            # only, so a copy sitting INSIDE an excluded region was still
            # tested -- and that is the worst case, being both
            # copy-number-inflated and physically cis on the derivative
            # chromosome while the reference calls its contacts trans. It
            # dominated the top of the ranking outright. A masked bin cannot be
            # a target; it must not be a source either.
            if MASK[off[c] + min(src // a.bin, nb[c] - 1)]:
                n_src_masked += 1
                continue
            iscis = dch[s0:e0] == this
            d = np.abs(dps[s0:e0] - src)
            far = (~iscis) | (d > a.prox)
            if far.sum() < a.min_contacts:
                continue
            tb = gbin[s0:e0][far]
            unmasked = ~MASK[tb]
            far_idx = np.flatnonzero(far)[unmasked]
            tb = tb[unmasked]
            n = tb.size
            if n < a.min_contacts:
                continue
            u, cnt = np.unique(tb, return_counts=True)
            H = float((cnt.astype(np.float64) ** 2).sum()) / (n * n)
            top1 = int(cnt.max())

            # --- null: same reach, random destination ---
            cis_f = iscis[far_idx]
            dd = d[far_idx][cis_f]
            strat = np.searchsorted(EDGES, dd, side="right") - 1
            draws = np.zeros((a.draws, n), np.int64)
            col = 0
            for s_ in np.unique(strat):
                k = int((strat == s_).sum())
                lo = EDGES[s_]
                hi = EDGES[s_ + 1] if s_ + 1 < EDGES.size else np.inf
                dist = np.abs(bincen[cbins] - src)
                cand = cbins[(dist > lo) & (dist <= hi)]
                if cand.size == 0:
                    cand = cbins
                w = W[cand]
                if w.sum() <= 0:
                    w = np.ones(cand.size)
                cw = np.cumsum(w / w.sum())
                r = rng.random((a.draws, k))
                draws[:, col:col + k] = cand[np.searchsorted(cw, r)]
                col += k
            ktr = n - col
            if ktr:
                r = rng.random((a.draws, ktr))
                draws[:, col:] = tbins[np.searchsorted(wt_c / wt_c[-1], r)]

            did = np.repeat(np.arange(a.draws), n)
            flat = did * NB + draws.ravel()
            uu, cc2 = np.unique(flat, return_counts=True)
            dsq = np.bincount(uu // NB, weights=cc2.astype(np.float64) ** 2,
                              minlength=a.draws)
            Hn = dsq / (n * n)
            t1 = np.zeros(a.draws)
            np.maximum.at(t1, uu // NB, cc2)
            mu, sd = float(Hn.mean()), float(Hn.std(ddof=1))
            zsc = (H - mu) / sd if sd > 0 else float("nan")
            rows.append(dict(
                subfamily=k_sub[i], family=k_fam[i], chrom=c, start=int(src),
                klass=k_cls[i], orientation="antisense" if anti[s0] else "sense",
                n_contacts=int(n), n_bins=int(u.size), H=H, H_null=mu,
                H_null_sd=sd, z=zsc, top1=top1, top1_null=float(t1.mean()),
                eff_bins=1.0 / H if H > 0 else float("nan")))
        print("  %s done (%d cells so far)" % (c, len(rows)), flush=True)
    print("dropped %d (copy, orientation) cells whose SOURCE bin is masked"
          % n_src_masked)

    # SECOND REFERENCE: the empirical spread of H among copies with a SIMILAR
    # CONTACT COUNT. H carries a hard 1/n floor, so n must be matched. This
    # absorbs whatever nuclear-neighbourhood structure is common to all TE RNA
    # -- structure the global accessibility map demonstrably does not capture,
    # since against it essentially every copy looked concentrated. It asks the
    # narrower and more answerable question: is this copy unusual AMONG COPIES?
    if rows:
        ln = np.array([np.log2(r["n_contacts"]) for r in rows])
        hh = np.array([r["H"] for r in rows])
        lo = np.floor(ln * 2) / 2.0          # half-log2 strata
        for b in np.unique(lo):
            k = lo == b
            if k.sum() < 20:
                for i in np.flatnonzero(k):
                    rows[i]["z_emp"] = float("nan")
                continue
            mu, sd = hh[k].mean(), hh[k].std(ddof=1)
            for i in np.flatnonzero(k):
                rows[i]["z_emp"] = ((hh[i] - mu) / sd) if sd > 0 else float("nan")
                rows[i]["H_peer"] = float(mu)

    for r in rows:
        r.setdefault("z_emp", float("nan"))
        r.setdefault("H_peer", float("nan"))
        ze = r["z_emp"]
        r["p_emp"] = 0.5 * erfc(ze / sqrt(2)) if ze == ze else float("nan")
        r["p"] = 0.5 * erfc(r["z"] / sqrt(2)) if r["z"] == r["z"] else float("nan")
    ps = sorted((r["p"], i) for i, r in enumerate(rows) if r["p"] == r["p"])
    m2, prev = len(ps), 1.0
    for rank in range(m2 - 1, -1, -1):
        pv, i = ps[rank]
        prev = min(prev, pv * m2 / (rank + 1))
        rows[i]["q"] = prev

    rows.sort(key=lambda r: -r["z_emp"] if r["z_emp"] == r["z_emp"] else 1e9)
    ps2 = sorted((r["p_emp"], i) for i, r in enumerate(rows)
                 if r["p_emp"] == r["p_emp"])
    m3, prev2 = len(ps2), 1.0
    for rank in range(m3 - 1, -1, -1):
        pv, i = ps2[rank]
        prev2 = min(prev2, pv * m3 / (rank + 1))
        rows[i]["q_emp"] = prev2
    for r in rows:
        r.setdefault("q_emp", float("nan"))

    cols = ["subfamily", "family", "orientation", "klass", "chrom", "start",
            "n_contacts", "n_bins", "eff_bins", "H", "H_null", "H_null_sd",
            "z", "p", "q", "H_peer", "z_emp", "p_emp", "q_emp",
            "top1", "top1_null"]
    path = os.path.join(OUT, "per_copy_targets%s.tsv" % a.tag)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(("%.6g" % r[c]) if isinstance(r[c], float)
                               else str(r.get(c, "")) for c in cols) + "\n")
    print("\nwrote %s  (%d copy x orientation cells)" % (path, len(rows)))
    print("cells q < 0.05 vs reach-matched null : %d" %
          sum(1 for r in rows if r.get("q", 1) < 0.05))
    print("cells q_emp < 0.05 vs matched-depth peers : %d" %
          sum(1 for r in rows if r.get("q_emp", 1) == r.get("q_emp", 1)
              and r.get("q_emp", 1) < 0.05))
    print("\ntop by z_emp (more concentrated than equally-deep peer copies):")
    print("  %-9s %-10s %-14s %-7s %8s %7s %8s %7s %6s %8s"
          % ("subfam", "orient", "locus", "class", "contacts", "effbin",
             "H", "H_peer", "z_emp", "q_emp"))
    for r in rows[:20]:
        print("  %-9s %-10s %-14s %-7s %8d %7.1f %8.4f %7.4f %6.2f %8.2g"
              % (r["subfamily"], r["orientation"],
                 "%s:%d" % (r["chrom"], r["start"]), r["klass"][:7],
                 r["n_contacts"], r["eff_bins"], r["H"], r["H_peer"],
                 r["z_emp"], r.get("q_emp", float("nan"))))


if __name__ == "__main__":
    main()
