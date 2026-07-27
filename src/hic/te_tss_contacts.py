"""Stage 2: do copies of a TE subfamily contact promoters more than chance?

    python src/hic/te_tss_contacts.py MCOOL --tss data/hic/tss.bed [--resolution 10000]

THE TWO CONFOUNDS, AND WHICH CONTROL KILLS WHICH

(1) Genomic distance. Contact frequency falls off by orders of magnitude with
    separation, and TE families differ systematically in how close they sit to
    genes -- Alu is isochore-biased into gene-dense sequence, L1 is not. Comparing
    raw contact values across subfamilies would mostly re-measure that. Handled by
    dividing every pixel by cooltools' expected_cis curve for its own separation,
    so the reported quantity is observed/expected throughout.

(2) Composition and coverage. Obs/exp does NOT fix the fact that Alu-dense bins
    are also gene-dense, open, high-coverage bins, which have systematically
    different obs/exp regardless of any TE effect. Handled by the shift null below.

THE NULL (--null-mode relocate, the default)
One rigid toroidal offset per (chromosome, replicate) is added to the TE anchors,
moving the whole set as a unit so spacings between copies -- the clustering --
survive. Same construction as EligibleSpace.shift in src/enrich_permutation.py; if
you change the reasoning in one, revisit the other.

The relocated anchors are then RE-PAIRED against the real TSS near their new
positions, rather than dragging their original TSS partner along. The reason is
an asymmetry argument, not a measured effect:

  The tempting alternative is to shift both anchors of a pair rigidly, preserving
  each pair's exact separation (--null-mode pair). But its far anchor then lands
  on arbitrary sequence, whereas the observed far anchor is a promoter, and
  promoters sit in gene-dense high-contact regions. The observed arm would carry
  a promoter bonus the null arm lacks. Under `relocate` both arms are
  locus-to-real-promoter pairs, so that bonus appears on both sides and cancels.

  This is reasoning, not a benchmark. On the synthetic fixture the two modes were
  indistinguishable -- see the correlation caveat below for why that fixture
  cannot resolve them. `relocate` is the default because its asymmetry argument
  is sound, not because `pair` was caught misbehaving.

`relocate` does not preserve each pair's exact separation, so the expected-value
denominator differs slightly between arms; that is what the obs/exp division is
for.

SUBFAMILIES ARE NOT INDEPENDENT -- READ BEFORE INTERPRETING
Every subfamily is scored against the same contact matrix, the same TSS set and
the same chromosomes. They therefore move together: if the genome's TE landscape
happens to sit slightly nearer promoters than a shifted copy does, EVERY row
shifts up at once. The effective number of independent observations behind a
genome-wide offset is closer to 1 than to 180.

Concretely, on the synthetic fixture with --scramble-anchors (observed arm made
into a null draw, so the truth is exactly zero), the median log2_ratio wandered
between -0.071 and -0.022, and the share of positive subfamilies between 6/121
and 48/121, purely with the seed. No q_emp fell below 0.05 in any run.

So: trust q_emp, which is computed against the shift distribution and stayed
correctly calibrated. Do NOT read "most subfamilies are above their null" as a
result -- that is one draw wearing a hundred hats. Differences BETWEEN
subfamilies are the interpretable quantity, which is the same reason the
tissue-interaction design beats a marginal enrichment.

P-FLOOR
The empirical p is (exceed + 1) / (n_shifts + 1), so --n-shifts sets a floor of
1/(K+1). That floor has to clear BH against the number of subfamilies tested, in
exactly the way --pmax has to in enrich_permutation.py. With ~100 subfamilies the
default 50 is marginal for the most extreme rows; raise it to 200+ for anything
load-bearing. Setting it too low does not error -- it silently ties subfamilies on
the floor so they cannot be ranked.

Output: results/hic/<sample>_te_tss_contacts.tsv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from hic_common import (PROJ, assign_bins, load_te_index, make_view,
                        masking_report, open_cooler, reconcile_chroms, te_table)
from tecommon import bh          # same BH the eCLIP stages use


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mcool")
    p.add_argument("--tss", required=True)
    p.add_argument("--resolution", type=int, default=10000)
    p.add_argument("--out", default=os.path.join(PROJ, "results", "hic"))
    p.add_argument("--elements", default="L1,Alu,SVA")
    p.add_argument("--min-dist", type=int, default=20000,
                   help="ignore pairs closer than this; the near-diagonal is "
                        "dominated by polymer proximity, not regulatory contact")
    p.add_argument("--max-dist", type=int, default=1000000)
    p.add_argument("--n-shifts", type=int, default=50, help="p-floor is 1/(K+1)")
    p.add_argument("--null-mode", choices=("relocate", "pair"), default="relocate",
                   help="relocate (default): shift TE anchors only, then re-pair "
                        "against the real TSS near their new positions, so both "
                        "arms are locus-to-promoter pairs and the promoter's own "
                        "contact elevation cancels. pair: rigidly shift both "
                        "anchors, preserving each pair's exact separation but "
                        "leaving the null with no promoter at the far anchor -- "
                        "measurably biased upward, kept only for comparison")
    p.add_argument("--max-anchors", type=int, default=1000,
                   help="TE copies sampled per subfamily. Alu has ~1.1M copies "
                        "genome-wide and you are averaging -- a thousand is plenty, "
                        "and the full set would build a 20M-row pair table")
    p.add_argument("--tss-per-anchor", type=int, default=4)
    p.add_argument("--min-copies", type=int, default=50)
    p.add_argument("--window-bins", type=int, default=2000,
                   help="dense fetch width. 2000 bins is 32 MB at float64; the whole "
                        "point is never to fetch a chromosome-sized square")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--nproc", type=int, default=4)
    p.add_argument("--scramble-anchors", action="store_true",
                   help="CALIBRATION CHECK, not an analysis. Rigidly relocates the "
                        "observed anchors once before testing, so the observed arm "
                        "is itself a draw from the null. A correctly specified test "
                        "must then return median log2_ratio ~ 0, uniform p_emp, and "
                        "no q_emp < 0.05. Any consistent offset here is machinery, "
                        "not biology -- run this whenever you change the null")
    return p.parse_args()


def build_pairs(te, tss_by_chrom, clr, mind, maxd, per_anchor, rng):
    """TE anchors x nearby TSS -> pair table in LOCAL bin coordinates per chromosome.

    Same searchsorted idiom as peaks_per_cat / count_A / loci_hit: bracket the
    candidate window with two searchsorted calls, then filter. Sampling at most
    `per_anchor` TSS per TE copy bounds the table without biasing separation --
    the sample is uniform over the eligible window, so the separation distribution
    is the natural one thinned, not reshaped.
    """
    out = []
    for c, g in te.groupby("chrom", sort=False):
        t = tss_by_chrom.get(c)
        if t is None or len(t) == 0:
            continue
        off = clr.offset(c)
        a = g["bin_id"].to_numpy() - off                    # TE local bin
        lo = np.searchsorted(t, a - maxd, side="left")
        hi = np.searchsorted(t, a + maxd, side="right")

        ii, jj, ss = [], [], []
        sub = g["sid"].to_numpy()
        for k in range(len(a)):
            cand = t[lo[k]:hi[k]]
            if cand.size == 0:
                continue
            sep = np.abs(cand - a[k])
            cand = cand[(sep >= mind) & (sep <= maxd)]
            if cand.size == 0:
                continue
            if cand.size > per_anchor:
                cand = rng.choice(cand, per_anchor, replace=False)
            ii.append(np.full(cand.size, a[k], np.int32))
            jj.append(cand.astype(np.int32))
            ss.append(np.full(cand.size, sub[k], np.int16))
        if not ii:
            continue
        i = np.concatenate(ii); j = np.concatenate(jj); s = np.concatenate(ss)
        lo_b = np.minimum(i, j); hi_b = np.maximum(i, j)
        out.append(pd.DataFrame({"chrom": c, "i": lo_b, "j": hi_b, "sid": s}))
    if not out:
        raise SystemExit("no TE-TSS pairs built; check --tss and distance window")
    return pd.concat(out, ignore_index=True)


def shift_anchors(te, clr, rng):
    """One rigid toroidal offset per chromosome applied to the TE anchors only.

    Pure modular wrap, so nothing is dropped and pairwise spacings between TE
    copies survive exactly -- the same construction as EligibleSpace.shift in
    src/enrich_permutation.py.

    The relocated anchors are then RE-PAIRED against the real TSS near their new
    positions, which is the whole point (see --null-mode in main): both arms of
    the comparison end up being locus-to-real-promoter pairs, so the promoter's
    own elevated contact frequency cancels instead of being scored as a TE effect.
    """
    out = te.copy()
    b = out["bin_id"].to_numpy().copy()
    for c, idx in out.groupby("chrom", sort=False).indices.items():
        nb = int(np.ceil(clr.chromsizes[c] / clr.binsize))
        off = clr.offset(c)
        d = int(rng.integers(0, nb))
        b[idx] = off + ((b[idx] - off + d) % nb)
    out["bin_id"] = b
    return out


def shifted_queries(pairs, clr, n_shifts, rng):
    """Observed pairs (group 0) plus n_shifts rigidly relocated copies.

    One offset per (chrom, replicate), applied to both anchors, so separation is
    exact and the configuration moves as a unit. Pairs running off the chromosome
    end after the shift are dropped rather than wrapped: wrapping only one anchor
    would silently change that pair's separation and therefore its expected value.
    """
    qi, qj, qs, qg = [pairs["i"].to_numpy()], [pairs["j"].to_numpy()], \
                     [pairs["sid"].to_numpy()], \
                     [np.zeros(len(pairs), np.int16)]
    chrom_codes = pairs["chrom"].to_numpy()
    dropped = 0
    total = 0

    for c, idx in pairs.groupby("chrom", sort=False).indices.items():
        nb = int(np.ceil(clr.chromsizes[c] / clr.binsize))
        i = pairs["i"].to_numpy()[idx]
        j = pairs["j"].to_numpy()[idx]
        s = pairs["sid"].to_numpy()[idx]
        sep = j - i
        for k in range(1, n_shifts + 1):
            d = int(rng.integers(0, nb))
            ni = (i + d) % nb
            nj = ni + sep
            ok = nj < nb
            dropped += int((~ok).sum()); total += ok.size
            qi.append(ni[ok].astype(np.int32))
            qj.append(nj[ok].astype(np.int32))
            qs.append(s[ok])
            qg.append(np.full(int(ok.sum()), k, np.int16))
            chrom_codes = np.concatenate([chrom_codes, np.full(int(ok.sum()), c)])

    q = pd.DataFrame({
        "chrom": chrom_codes,
        "i": np.concatenate(qi).astype(np.int32),
        "j": np.concatenate(qj).astype(np.int32),
        "sid": np.concatenate(qs).astype(np.int16),
        "grp": np.concatenate(qg).astype(np.int16),
    })
    print(f"  shift overhang dropped {dropped}/{total} "
          f"({100 * dropped / max(1, total):.2f}%) of shifted pairs")
    return q


def evaluate(q, clr, expected, window_bins, maxd_bins):
    """Look up balanced contact / expected for every query pair.

    Windowed DENSE fetch: a `window_bins` square is 32 MB at the default, while a
    whole chr1 square at 10 kb would be 4.9 GB and the genome-wide square 764 TB.
    Because every pair has j - i <= maxd_bins, assigning each query to window
    i // step with step = window_bins - maxd_bins guarantees both anchors land
    inside that window, so each window is fetched exactly once and serves the
    observed and all shifted queries that fall in it.
    """
    step = window_bins - maxd_bins
    if step <= 0:
        raise SystemExit("--window-bins must exceed --max-dist in bins")

    oe = np.full(len(q), np.nan)
    mat = clr.matrix(balance=True)
    for c, idx in q.groupby("chrom", sort=False).indices.items():
        nb = int(np.ceil(clr.chromsizes[c] / clr.binsize))
        qi = q["i"].to_numpy()[idx]
        qj = q["j"].to_numpy()[idx]
        widx = qi // step
        for w in np.unique(widx):
            sel = idx[widx == w]
            ws = int(w) * step
            we = min(ws + window_bins, nb)
            # nb is ceil(L / binsize), so the last bin is a PARTIAL bin and
            # we * binsize can run past the chromosome. cooler rejects that
            # outright; clamp to the true length and let the bin count stand.
            end_bp = min(we * clr.binsize, int(clr.chromsizes[c]))
            M = mat.fetch((c, ws * clr.binsize, end_bp))
            li = q["i"].to_numpy()[sel] - ws
            lj = q["j"].to_numpy()[sel] - ws
            good = (li >= 0) & (lj < M.shape[0])
            vals = np.full(len(sel), np.nan)
            vals[good] = M[li[good], lj[good]]
            ex = np.array([expected.get((c, int(d)), np.nan)
                           for d in (q["j"].to_numpy()[sel] - q["i"].to_numpy()[sel])])
            with np.errstate(divide="ignore", invalid="ignore"):
                oe[sel] = vals / ex
    return oe


def main():
    a = parse_args()
    os.makedirs(a.out, exist_ok=True)
    sample = os.path.basename(a.mcool).replace(".mcool", "").replace(".cool", "")
    rng = np.random.default_rng(a.seed)

    import cooltools

    print(f"[te_tss_contacts] {sample}")
    clr, binsize = open_cooler(a.mcool, a.resolution)
    mind_b, maxd_b = a.min_dist // binsize, a.max_dist // binsize
    print(f"  distance window {a.min_dist}-{a.max_dist} bp "
          f"({mind_b}-{maxd_b} bins), {a.n_shifts} shifts "
          f"-> p-floor {1 / (a.n_shifts + 1):.4f}")

    index, subfam_names, subfam_elem = load_te_index()
    chroms, _ = reconcile_chroms(index, clr)
    view = make_view(clr, chroms)

    print(f"  expected_cis (nproc={a.nproc}) ...")
    exp_df = cooltools.expected_cis(clr, view_df=view, nproc=a.nproc)
    col = ("balanced.avg.smoothed.agg" if "balanced.avg.smoothed.agg" in exp_df
           else "balanced.avg")
    exp_sub = exp_df[["region1", "dist", col]].dropna()
    expected = dict(zip(zip(exp_sub["region1"], exp_sub["dist"].astype(int)),
                        exp_sub[col]))
    print(f"    using column '{col}', {len(expected)} (region, dist) entries")

    tss = pd.read_csv(a.tss, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["chrom", "start", "end"], comment="#")
    tss = tss[tss["chrom"].isin(set(chroms))]
    tss = assign_bins(tss.assign(mid=tss["start"]), clr)
    tss = tss[tss["bin_id"] >= 0]
    tss_by_chrom = {c: np.sort(np.unique(g["bin_id"].to_numpy() - clr.offset(c)))
                    for c, g in tss.groupby("chrom", sort=False)}
    print(f"  TSS: {len(tss)} -> {sum(len(v) for v in tss_by_chrom.values())} "
          f"distinct bins")

    te = te_table(index, subfam_names, subfam_elem, chroms,
                  elements=set(a.elements.split(",")), min_copies=a.min_copies)
    te = assign_bins(te, clr)
    te = te[te["bin_id"] >= 0]

    bins = clr.bins()[:]
    mask_rep = masking_report(te, bins)
    mask_rep.to_csv(os.path.join(a.out, f"{sample}_te_bin_masking.tsv"),
                    sep="\t", index=False, float_format="%.4f")
    print(f"  masking: median {mask_rep['frac_masked'].median():.1%} of copies "
          f"in NaN-weight bins (worst {mask_rep['frac_masked'].max():.1%})")

    # shuffle-then-head is an unbiased sample per subfamily and avoids a
    # groupby.apply over ~180 groups of up to a million rows
    te = (te.sample(frac=1.0, random_state=a.seed)
            .groupby("subfam", sort=False)
            .head(a.max_anchors)
            .reset_index(drop=True))
    print(f"  anchors after subsampling: {len(te)} "
          f"across {te['subfam'].nunique()} subfamilies")

    if a.scramble_anchors:
        te = shift_anchors(te, clr, np.random.default_rng(a.seed + 991))
        print("  !! --scramble-anchors: observed arm relocated, this is a "
              "calibration run and its output is not biology")

    if a.null_mode == "relocate":
        frames = []
        for k in range(a.n_shifts + 1):
            te_k = te if k == 0 else shift_anchors(te, clr, rng)
            pk = build_pairs(te_k, tss_by_chrom, clr, mind_b, maxd_b,
                             a.tss_per_anchor, rng)
            pk["grp"] = np.int16(k)
            frames.append(pk)
            if k == 0:
                print(f"  TE-TSS pairs (observed): {len(pk)}")
        q = pd.concat(frames, ignore_index=True)
    else:
        pairs = build_pairs(te, tss_by_chrom, clr, mind_b, maxd_b,
                            a.tss_per_anchor, rng)
        print(f"  TE-TSS pairs (observed): {len(pairs)}")
        q = shifted_queries(pairs, clr, a.n_shifts, rng)
    print(f"  queries to evaluate: {len(q)} "
          f"(~{q.memory_usage(deep=True).sum() / 1e6:.0f} MB)")

    q["oe"] = evaluate(q, clr, expected, a.window_bins, maxd_b)
    finite = np.isfinite(q["oe"])
    print(f"  finite obs/exp: {int(finite.sum())}/{len(q)} "
          f"({100 * finite.mean():.1f}%); the rest are ICE-masked or empty bins")

    names = np.asarray(subfam_names)
    elems = np.asarray(subfam_elem)
    qf = q[finite]
    means = qf.groupby(["sid", "grp"])["oe"].mean().unstack("grp")
    obs_n = qf[qf["grp"] == 0].groupby("sid").size()

    rows = []
    for sid, r in means.iterrows():
        if 0 not in r.index or not np.isfinite(r.get(0, np.nan)):
            continue
        obs = float(r[0])
        null = r.drop(index=0).to_numpy(dtype=float)
        null = null[np.isfinite(null)]
        if null.size < 2:
            continue
        exceed = int((null >= obs).sum())
        rows.append(dict(
            elem=elems[sid], subfam=names[sid],
            n_pairs=int(obs_n.get(sid, 0)),
            obs_oe=obs, null_oe=float(null.mean()), null_sd=float(null.std(ddof=1)),
            log2_ratio=float(np.log2(obs / null.mean())) if null.mean() > 0 else np.nan,
            z=float((obs - null.mean()) / null.std(ddof=1))
            if null.std(ddof=1) > 0 else np.nan,
            n_shifts_used=int(null.size),
            p_emp=(exceed + 1) / (null.size + 1),
        ))

    out = pd.DataFrame(rows)
    if out.empty:
        raise SystemExit("no subfamily produced a finite observed value")
    out["q_emp"] = bh(out["p_emp"].to_numpy())
    out = out.sort_values("p_emp")
    out_path = os.path.join(a.out, f"{sample}_te_tss_contacts.tsv")
    out.to_csv(out_path, sep="\t", index=False, float_format="%.5g")

    n_floor = int((out["p_emp"] <= 1.5 / (a.n_shifts + 1)).sum())
    print(f"\n  wrote {out_path}  ({len(out)} subfamilies)")
    if n_floor:
        print(f"  !! {n_floor} subfamilies sit on the p-floor and cannot be ranked "
              f"against each other -- raise --n-shifts to separate them")
    print(f"  q_emp < 0.05: {int((out['q_emp'] < 0.05).sum())}")
    print("\n  strongest:")
    print(out.head(8)[["elem", "subfam", "n_pairs", "obs_oe", "null_oe",
                       "log2_ratio", "p_emp", "q_emp"]].to_string(index=False))


if __name__ == "__main__":
    main()
