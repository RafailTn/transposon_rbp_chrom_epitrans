"""Permutation null for the RBP x TE-subfamily hits (the claim, not the filter).

The hypergeometric in enrich_stranded.py is a fast screen: it assumes peaks land
independently at the pooled-peak rate, which they do not (peaks cluster on genes,
TEs cluster in the genome). This re-tests every hit it flagged against a null that
keeps both structures: circularly shift each RBP's peak set -- same count, same
widths, same strands, same relative spacing -- by one random toroidal offset in
the cell line's expressed/mappable space, re-run the identical locus-deduped
overlap, and read the p-value off the empirical distribution of loci-in-category.

Rigid circular shift (not independent per-peak placement) is the point: it moves
the whole configuration by a common offset, so the pairwise spacings between an
RBP's peaks -- its gene-level clustering -- survive, and only the absolute phase
is randomized. Independent placement would destroy exactly the clustering the
hypergeometric already ignores.

Edge effect -- known, measured, and deliberately left in. A shifted peak's start
lands in flat eligible space but its body extends in GENOMIC space, so it can
overhang its segment into a gap. An earlier version of this docstring claimed
segments greatly exceed peak widths and the effect was negligible and unbiased
w.r.t. TE identity. Both halves are false. Measured on the pooled footprints:

    footprint   segments   median seg   median peak   P(overhang)
    24.0 Mb     341,218       57 bp        48 bp         59%

The footprint is fragmented into segments barely wider than the peaks themselves,
which is what drives the overhang rate that high -- it is a property of this peak
set's geometry, so re-measure it rather than assuming 59% on another one.

It is also not TE-neutral: the 100 bp flanks just outside K562 segments are
13.8% L1/Alu/SVA against 9.8% inside the footprint, ~1.4x denser. So overhanging
peaks sample TE-richer sequence than the observed peaks occupy, which INFLATES
the null's TE hit rate. The bias is therefore conservative: exp_loci_perm too
high, fold_perm understated, p_perm overstated. Surviving hits are safe; hits
that just miss q_perm < 0.05 may be false negatives.

Do not "fix" this by requiring shifted peaks to sit wholly inside a segment. On
eCLIP K562 that makes 37% of segments unable to host a median-width peak at all
(87% at p95 width 119) and shrinks usable start-space to 41% of the footprint,
concentrating placement in the long segments -- which are the most heavily
co-bound, exon/UTR-biased, TE-POOR regions. That trades a bounded conservative
bias for an unbounded anti-conservative one. Removing the effect properly needs a
different null, not a patch.

Eligible/mappable space: no expression or mappability track ships with this
project, so the operational proxy is the pooled peak footprint -- the merged union
of all that cell line's RBP peaks (116 in K562, 97 in HepG2). That is exactly "where this assay sees
RNA-bound signal at all in this cell", and is the permutation analog of the
pooled-peak background. Override with an external BED via --eligible if you have a
real expressed+mappable region set. Peaks are only ever relocated inside this
footprint, so the null preserves the transcriptome/TE bias the hypergeometric
ignores; what is left is positional enrichment.

Cost control:
  * only hits that pass the hypergeometric screen are permuted (the fast filter);
  * permutations for one RBP are vectorised in chunks (all its peaks x B perms at
    once);
  * sequential stopping (Besag & Clifford 1991): once every candidate category of
    an RBP has accumulated H_STOP exceedances it is decided non-significant and we
    stop permuting it, so only genuinely-enriched categories pay the full budget.
Empirical p = (exceed + 1) / (n_perm + 1); BH within cell line, matching the
hypergeometric stage. The floor 1/(P_MAX+1) must clear BH, so P_MAX is set so that
1/(P_MAX+1) * (max candidates in a cell) < 0.05. That bound scales with the
CANDIDATE COUNT, so it must be re-derived for any new peak set rather than
inherited. Here the largest candidate set is HepG2's 1,037, needing P_MAX >
20,740; src/eclip/run.sh therefore passes --pmax 25000, which satisfies the rule
outright (1/25001 * 1037 = 0.042). This module's own default is deliberately far
lower and is NOT safe for a screen this size.

Beware the failure mode: an undersized budget does not error, and its q-values can
still come out significant, because when many hits tie at the floor BH gives them
q = floor * n / k with k = the tie count, which lands well below 0.05. That is the
tie structure rescuing the budget, not the rule being met.

Budget also sets p-value RESOLUTION, which matters for ranking: at --pmax 2000,
56% of candidates sat on the floor and were mutually unorderable; at 25000 that is
43%, across 623 distinct p_perm values. Raising it barely moves the verdict (1789
-> 1792 survivors, measured on the superseded hg38_fixed annotation) but it does
sharpen fold_perm, whose denominator exp_loci_perm is a Monte Carlo mean -- e.g.
HNRNPC/AluSg7/K562 moves 79.4 -> 74.0, a better estimate rather than a different
result.

Two folds, two baselines -- do not mix them up:
  fold_pooled  obs / expected-from-the-pooled-peak-background. The screen's effect
               size, computed in enrich_stranded.py: "this RBP concentrates on the
               subfamily more than the OTHER RBPs in this cell line do".
  fold_perm    obs / mean(loci under the shift null), from exp_loci_perm. The effect
               size matched to q_perm: "this RBP concentrates on the subfamily more
               than a random relocation of ITSELF inside the RBP-bound footprint".

They do NOT agree, and the gap is structural, not noise:

    median fold_perm / fold_pooled   (loci >= 10)
           0.64          (HepG2 0.75, K562 0.55)

Below 1 means exp_loci_perm sits ABOVE the hypergeometric's expected count. Two
properties of this peak set plausibly drive that, both raising the null's locus
yield: peaks are wide (median 48 bp), so a relocated peak touches several copies;
and the segment overhang documented above is large (59% of placements), leaking
the null into TE-denser gap sequence. That is the likely mechanism, not a verified
one -- it has not been isolated experimentally.

Do NOT memorise the sign. The competing effect runs the other way: relocate a
tightly clustered set of narrow peaks and they keep piling onto the same few
copies, covering FEWER distinct loci than the pooled loci-per-peak rate predicts,
which puts exp_loci_perm below the hypergeometric expectation and the ratio above
1. Which effect wins is a property of peak width and clustering, so a peak set
with different geometry can invert this ratio.

What does NOT change: the gap is not a constant rescaling, so ranking by one fold is
not ranking by the other, and a cross-cell "top N by fold" table compares hits whose
denominators differ -- read it as a within-cell statement, not a league table. Here
the top-15 by fold_perm are 12 K562 / 3 HepG2 in both orientations.

The hypergeometric OVERSTATES significance either way (it ignores the variance that
clustering adds), which is why 269 of the 2,063 screened hits still lose their
q -- 195 of them sense-strand L1.

Quote fold_perm next to q_perm; fold_pooled belongs next to q_hyper.

Usage: python enrich_permutation.py [data_dir] [results_dir] [--eligible BED] [--pmax N]
"""
import numpy as np, pickle, os, sys, csv, glob, time, zlib
from collections import defaultdict
from tecommon import read_bed, bh


def _fmt(v):
    return f"{v:.3e}" if isinstance(v, float) else str(v)


def _seed(*parts):
    return zlib.crc32(("|".join(map(str, parts))).encode()) & 0xFFFFFFFF

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARGV = sys.argv[1:]


def take_flag(name, default=None):
    """Pull `--name VALUE` out of ARGV so it is not mistaken for a positional arg."""
    if name in ARGV:
        i = ARGV.index(name)
        val = ARGV[i + 1]
        del ARGV[i:i + 2]
        return val
    return default


ELIGIBLE_BED = take_flag("--eligible")                # optional external region set
# Default P_MAX clears within-cell BH at the floor: 1/(P+1) < 0.05 / max-candidates.
# The largest candidate set here is HepG2's 1,037, so P must exceed ~20,740; 25000
# gives margin and resolves p to ~4e-5. Re-derive this if the peak set changes --
# see the budget section of the module docstring for why too-low fails silently.
P_MAX = int(take_flag("--pmax", 25000))               # perms before we give up resolving
POS = [a for a in ARGV if not a.startswith("--")]
DATA = POS[0] if len(POS) > 0 else f"{PROJ}/data/eclip"
RESULTS = POS[1] if len(POS) > 1 else f"{PROJ}/results/eclip"
H_STOP = 25                                            # exceedances -> decided non-significant
CHUNK_CELLS = 6_000_000                                # ~ samples per vectorised chunk
QCUT = 0.05
SEED = 0
SIG = f"{RESULTS}/rbp_te_enrichment_stranded_significant.tsv"

with open(f"{PROJ}/cache/te_index.pkl", "rb") as fh:
    blob = pickle.load(fh)
INDEX, SUBFAM, SELEM = blob["index"], blob["subfam_names"], blob["subfam_elem"]
SID = {s: i for i, s in enumerate(SUBFAM)}
NSUB = len(SUBFAM)
NCAT = NSUB * 2
ORIENT = ("sense", "antisense")


def merge_intervals(s, e):
    """Merge overlapping [s,e) intervals. Returns sorted, disjoint starts/ends."""
    if len(s) == 0:
        return s, e
    o = np.argsort(s, kind="stable")
    s, e = s[o], e[o]
    cummax = np.maximum.accumulate(e)
    new = np.ones(len(s), bool)
    new[1:] = s[1:] > cummax[:-1]                     # a gap opens a new interval
    idx = np.where(new)[0]
    return s[idx], np.maximum.reduceat(e, idx)


class Eligible:
    """Concatenated eligible-region coordinate space for uniform position sampling."""

    def __init__(self, per_chrom):
        self.chroms = sorted(per_chrom)
        self.cid = {c: i for i, c in enumerate(self.chroms)}
        seg_c, seg_s, seg_l = [], [], []
        for c in self.chroms:
            es, ee = per_chrom[c]
            seg_c.append(np.full(len(es), self.cid[c], np.int32))
            seg_s.append(es); seg_l.append(ee - es)
        # seg_c is non-decreasing (chroms visited in sorted cid order) and seg_s is
        # ascending within each chrom block -- both relied on by to_flat below.
        self.seg_c = np.concatenate(seg_c)
        self.seg_s = np.concatenate(seg_s)
        self.seg_l = np.concatenate(seg_l)
        self.cum = np.cumsum(self.seg_l)              # cum[k] = end offset of seg k
        self.seg_off = self.cum - self.seg_l          # cum[k-1] = flat start offset of seg k
        self.total = int(self.cum[-1])

    def _flat_to_genomic(self, flat):
        """Flat concatenated-eligible coordinate(s) -> (chrom_id, genomic_start)."""
        k = np.searchsorted(self.cum, flat, side="right")   # segment containing flat
        return self.seg_c[k], self.seg_s[k] + (flat - self.seg_off[k])

    def _locate(self, cids, starts):
        """Per position: segment index, and whether the position is inside it.

        Per chrom the candidate segment is the rightmost one whose start <= the
        position. That candidate only actually CONTAINS the position if the
        position also falls before the segment's end -- a position in a gap
        between segments, or before the first segment on its chrom, has no
        containing segment at all. Callers must not use gk where inside is False.
        """
        gk = np.zeros(len(starts), np.int64)
        inside = np.zeros(len(starts), bool)
        for ci in np.unique(cids):
            m = cids == ci
            lo = np.searchsorted(self.seg_c, ci, side="left")
            hi = np.searchsorted(self.seg_c, ci, side="right")
            if lo == hi:                       # chrom carries no eligible segment
                continue
            segs = self.seg_s[lo:hi]
            j = np.searchsorted(segs, starts[m], side="right") - 1
            ok = j >= 0                        # not before the first segment
            k = lo + np.clip(j, 0, len(segs) - 1)
            ok &= starts[m] < self.seg_s[k] + self.seg_l[k]   # not in a gap
            gk[m], inside[m] = k, ok
        return gk, inside

    def contains(self, cids, starts):
        """Bool mask: which (chrom_id, start) positions lie inside an eligible segment.

        With the default eligible space -- the merged union of all the cell's peaks
        -- every peak start is inside by construction, so this is all-True and
        filtering on it is a no-op. It can be False only under --eligible with an
        external region set that does not cover the peaks, which is exactly the
        case to_flat cannot represent.
        """
        return self._locate(cids, starts)[1]

    def to_flat(self, cids, starts):
        """Genomic (chrom_id, start) -> flat concatenated-eligible coordinate.

        Inverse of _flat_to_genomic, and defined ONLY for positions inside
        eligible space. Callers must filter with contains() first; passing an
        outside position raises rather than returning a silently wrong
        coordinate. An earlier version clipped instead, which mapped a position
        in a gap onto the preceding segment -- no error, no warning, just a peak
        permuted from the wrong place. That was unreachable on the default
        eligible space but live under --eligible.
        """
        gk, inside = self._locate(cids, starts)
        if not inside.all():
            raise ValueError(
                f"{int((~inside).sum())} of {len(starts)} positions lie outside "
                "eligible space; filter with contains() before calling to_flat")
        return self.seg_off[gk] + (starts - self.seg_s[gk])

    def shift(self, rng, flat_starts, b):
        """b circular shifts of a fixed peak configuration.

        One random toroidal offset per permutation is added (mod total) to every
        peak's flat coordinate, so the whole set moves rigidly: pairwise spacings --
        i.e. the RBP's clustering -- are preserved, only the global phase is random.
        Returns (chrom_id, start) for b*n samples in perm-major order (perm = i // n).
        """
        deltas = rng.integers(0, self.total, size=b)                 # one offset per perm
        flat = ((flat_starts[None, :] + deltas[:, None]) % self.total).reshape(-1)
        return self._flat_to_genomic(flat)


# TE index re-keyed by eligible chrom-id, only for chroms that carry TEs
def build_te_by_cid(elig):
    te = {}
    for c, ci in elig.cid.items():
        ix = INDEX.get(c)
        if ix is None:
            continue
        te[ci] = (ix["starts"], ix["ends"], ix["sid"].astype(np.int64),
                  ix["strand"], int(ix["maxlen"]), len(ix["starts"]))
    return te


def count_A(TE, gchrom, gstart, width, strand, permid, B):
    """Distinct-locus count per (perm, category). Returns (B, NCAT) int64.

    Same locus-deduped overlap as enrich_stranded, but batched over B permutations:
    a copy counts once per (perm, category). permid[i] in [0,B) tags each sample's
    permutation; samples are laid out perm-major so permid = index // n_per_perm.
    """
    A = np.zeros(B * NCAT, np.int64)
    spe = gstart + width
    for ci in np.unique(gchrom):
        t = TE.get(int(ci))
        if t is None:
            continue
        ts, te, sid, tstr, L, nflat = t
        m = gchrom == ci
        sps, spe_, sstr, pm = gstart[m], spe[m], strand[m], permid[m]
        hi = np.searchsorted(ts, spe_, side="left")
        lo = np.searchsorted(ts, sps - L, side="left")
        cnt = hi - lo
        tot = int(cnt.sum())
        if tot == 0:
            continue
        loc = np.repeat(np.arange(len(sps), dtype=np.int64), cnt)   # local sample idx
        flat = np.repeat(lo, cnt) + (np.arange(tot, dtype=np.int64)
                                     - np.repeat(np.cumsum(cnt) - cnt, cnt))
        keep = (te[flat] > sps[loc]) & (ts[flat] < spe_[loc])
        loc, flat = loc[keep], flat[keep]
        if loc.size == 0:
            continue
        anti = (sstr[loc] != tstr[flat]).astype(np.int64)
        perm = pm[loc].astype(np.int64)
        # dedupe (perm, locus, orientation): composite key, unique, then bincount
        comp = (perm * nflat + flat) * 2 + anti
        u = np.unique(comp)
        anti_u = u & 1
        rest = u >> 1
        flat_u = rest % nflat
        perm_u = rest // nflat
        binidx = perm_u * NCAT + (sid[flat_u] * 2 + anti_u)
        A += np.bincount(binidx, minlength=B * NCAT)
    return A.reshape(B, NCAT)


# ---- candidates: every hit the hypergeometric screen flagged
cand = defaultdict(list)     # (cell, rbp) -> [(cat, row_dict), ...]
with open(SIG) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cat = SID[r["subfamily"]] * 2 + (0 if r["orientation"] == "sense" else 1)
        cand[(r["cell_line"], r["rbp"])].append((cat, r))
cells = sorted({c for c, _ in cand})
n_cand_total = sum(len(v) for v in cand.values())
print(f"{n_cand_total} hypergeometric hits to defend, across {len(cand)} (cell, RBP) pairs, "
      f"{len(cells)} cell lines; P_MAX={P_MAX}", flush=True)

out_rows = []
checked_shift = False
for cell in cells:
    t0 = time.time()
    beds = sorted(glob.glob(f"{DATA}/*_{cell}.bed"))
    rbp_of = {os.path.basename(b)[:-4].rsplit("_", 1)[0]: b for b in beds}
    cand_rbps = {rbp for (c, rbp) in cand if c == cell}

    # eligible space = merged union of ALL peaks in this cell (or external BED)
    acc = defaultdict(lambda: ([], []))
    keep_peaks = {}
    if ELIGIBLE_BED:
        ech, est, een, _ = read_bed(ELIGIBLE_BED)
        for c in np.unique(ech):
            acc[c][0].append(est[ech == c]); acc[c][1].append(een[ech == c])
    for rbp, path in rbp_of.items():
        ch, st, en, sd = read_bed(path)
        if not ELIGIBLE_BED:
            for c in np.unique(ch):
                acc[c][0].append(st[ch == c]); acc[c][1].append(en[ch == c])
        if rbp in cand_rbps:
            keep_peaks[rbp] = (ch, st, en, sd)
    elig = Eligible({c: merge_intervals(np.concatenate(a[0]), np.concatenate(a[1]))
                     for c, a in acc.items()})
    TE = build_te_by_cid(elig)

    for rbp in sorted(cand_rbps):
        ch, st, en, sd = keep_peaks[rbp]
        width = (en - st).astype(np.int64)
        cats = [c for c, _ in cand[(cell, rbp)]]

        # Restrict to peaks on eligible chroms, then to peaks whose start lies
        # inside an eligible SEGMENT, and use this SAME subset for the observed
        # count and every permutation, so they are strictly comparable. A peak the
        # null cannot place must not count toward the observed total either.
        # Both filters are no-ops on the default eligible space (union of all the
        # cell's peaks, so every peak start is inside by construction); they bite
        # only under --eligible with an external region set.
        gc = np.array([elig.cid.get(c, -1) for c in ch], np.int32)
        km = gc >= 0
        gc_e, st_e, w_e, sd_e = gc[km], st[km], width[km], sd[km]
        inside = elig.contains(gc_e, st_e)
        if not inside.all():
            print(f"  {rbp}_{cell}: {int((~inside).sum()):,} of {len(inside):,} peaks start "
                  f"outside eligible space -- dropped from observed and null alike", flush=True)
            gc_e, st_e, w_e, sd_e = gc_e[inside], st_e[inside], w_e[inside], sd_e[inside]
        n = len(st_e)
        if n == 0:
            for c, row in cand[(cell, rbp)]:
                out_rows.append({**row, "obs_loci": 0, "n_perm": 0, "exceed": 0, "p_perm": 1.0,
                                 "exp_loci_perm": "nan", "fold_perm": "nan"})
            continue

        # Map the observed peak starts into concatenated-eligible coordinates once;
        # each permutation is a single rigid toroidal offset of this configuration.
        flat_starts = elig.to_flat(gc_e, st_e)
        if not checked_shift:
            rc, rs = elig._flat_to_genomic(flat_starts)
            ok = bool(np.array_equal(rc, gc_e) and np.array_equal(rs, st_e))
            print(f"  flat round-trip on {rbp}_{cell}: {'PASS' if ok else 'FAIL'}", flush=True)
            checked_shift = True

        # observed loci per candidate category, via the SAME overlap routine
        A_obs = count_A(TE, gc_e, st_e, w_e, sd_e, np.zeros(n, np.int64), 1)[0]
        for c, row in cand[(cell, rbp)]:
            if A_obs[c] != int(row["loci_in_cat"]):
                print(f"  !! observed mismatch {rbp} {cell} cat {c}: "
                      f"{A_obs[c]} vs table {row['loci_in_cat']}", file=sys.stderr)

        exceed = {c: 0 for c in cats}
        null_sum = {c: 0 for c in cats}     # running sum of null loci, for the perm fold
        done = 0
        B = max(1, CHUNK_CELLS // max(n, 1))
        rng = np.random.default_rng(_seed(SEED, cell, rbp))
        while done < P_MAX and any(exceed[c] < H_STOP for c in cats):
            b = min(B, P_MAX - done)
            gchrom, gstart = elig.shift(rng, flat_starts, b)   # rigid circular shift
            permid = np.repeat(np.arange(b, dtype=np.int64), n)
            widths = np.tile(w_e, b)
            strands = np.tile(sd_e, b)
            A = count_A(TE, gchrom, gstart, widths, strands, permid, b)
            for c in cats:
                exceed[c] += int((A[:, c] >= A_obs[c]).sum())
                null_sum[c] += int(A[:, c].sum())
            done += b

        for c, row in cand[(cell, rbp)]:
            p = (exceed[c] + 1) / (done + 1)
            # Effect size against the SAME null the q-value comes from. Sequential
            # stopping shares one `done` across an RBP's categories, so this mean is
            # over the same permutations that produced `exceed`.
            exp = null_sum[c] / done
            fperm = (A_obs[c] / exp) if exp > 0 else float("inf")
            out_rows.append({**row, "obs_loci": int(A_obs[c]), "n_perm": done,
                             "exceed": exceed[c], "p_perm": p,
                             "exp_loci_perm": f"{exp:.4f}", "fold_perm": f"{fperm:.4f}"})
    print(f"{cell}: {len(cand_rbps)} RBPs, {sum(1 for k in cand if k[0]==cell)} hits "
          f"defended  [{time.time()-t0:.0f}s]", flush=True)

# ---- BH within cell line, matching the hypergeometric stage
by_cell = defaultdict(list)
for r in out_rows:
    by_cell[r["cell_line"]].append(r)
for cell, rows in by_cell.items():
    q = bh(np.array([r["p_perm"] for r in rows]))
    for r, qq in zip(rows, q):
        r["q_perm"] = qq

out_rows.sort(key=lambda r: (r["cell_line"], r["p_perm"], -float(r["fold_enrichment"])))
COLS = ["cell_line", "rbp", "element", "subfamily", "orientation", "obs_loci",
        "pooled_loci_in_cat", "fold_enrichment", "exp_loci_perm", "fold_perm",
        "q_BH", "n_perm", "exceed", "p_perm", "q_perm"]
HDRMAP = {"q_BH": "q_hyper",                           # rename for clarity in output
          "fold_enrichment": "fold_pooled"}            # the screen's fold, vs other RBPs
with open(f"{RESULTS}/rbp_te_enrichment_permutation.tsv", "w") as fh:
    fh.write("\t".join(HDRMAP.get(c, c) for c in COLS) + "\n")
    for r in out_rows:
        fh.write("\t".join(_fmt(r[c]) for c in COLS) + "\n")

sig = [r for r in out_rows if r["q_perm"] < QCUT]
with open(f"{RESULTS}/rbp_te_enrichment_permutation_significant.tsv", "w") as fh:
    fh.write("\t".join(HDRMAP.get(c, c) for c in COLS) + "\n")
    for r in sig:
        fh.write("\t".join(_fmt(r[c]) for c in COLS) + "\n")

ns = sum(1 for r in sig if r["orientation"] == "sense")
print(f"\n{len(sig)} / {len(out_rows)} hypergeometric hits survive the permutation null "
      f"at q < {QCUT}  ({ns} sense / {len(sig)-ns} antisense)")
print(f"wrote {RESULTS}/rbp_te_enrichment_permutation.tsv (+ _significant.tsv)")
