"""Collapse the two eCLIP replicates per RBP into one reproducible-peak BED.

Input: data/eCLIP_Peaks_K562_HepG2_{K562,HepG2}/<RBP>_<CELL>_{1,2}_chr.bed, BED6
       (chrom, start, end, name, signal, strand), 0-based half-open, hg38.
Output: data/eclip/<RBP>_<CELL>.bed, BED6, in the `<RBP>_<CELL>.bed` naming the
       rest of the pipeline parses with basename[:-4].rsplit("_", 1).

Why intersection and not union: the downstream background is pooled peaks over
all RBPs in the cell line, so a peak present in only one replicate does not just
add noise to its own row, it also inflates the shared denominator every other
RBP is tested against. Keeping the intersected REGION (not "rep1 peaks that have
rep2 support") makes the operation symmetric -- there is no primary replicate --
and it is the strictest honest reading of "reproducible peak" available without
IDR, which does not ship with these files.

Consequences, all reported at the end of the run:
  * an RBP whose replicate is empty intersects to nothing and is dropped, rather
    than silently inheriting the surviving replicate's peaks;
  * SFPQ_HepG2 has only one replicate at all and is dropped for the same reason.

Intervals are merged within each replicate before intersecting, so the output is
sorted and disjoint per (chrom, strand). That matters because `n_peaks_rbp` is a
peak count and is the denominator of every fold: overlapping duplicates there
would deflate the fold of the RBP that has them.

Column 5 is written as "." -- after merging, no single input signal value maps to
an output interval. Nothing downstream reads column 5.

Usage: python src/eclip/prep_replicates.py [data_dir] [out_dir]
"""
import numpy as np, os, sys, glob, pickle, re
from collections import defaultdict

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = sys.argv[1] if len(sys.argv) > 1 else f"{PROJ}/data"
OUT = sys.argv[2] if len(sys.argv) > 2 else f"{PROJ}/data/eclip"
SRC_DIRS = [f"{DATA}/eCLIP_Peaks_K562_HepG2_K562", f"{DATA}/eCLIP_Peaks_K562_HepG2_HepG2"]
CACHE = f"{PROJ}/cache/te_index.pkl"

# <RBP>_<CELL>_<rep>_chr.bed -- the rep digit is what we collapse over.
NAME_RE = re.compile(r"^(?P<rbp>.+)_(?P<cell>[^_]+)_(?P<rep>[12])_chr\.bed$")


def read_bed6(path):
    """BED6 -> {(chrom, strand): (starts, ends)}, strand kept as the raw +/- char."""
    by = defaultdict(lambda: ([], []))
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            s, e = by[(f[0], f[5])]
            s.append(int(f[1])); e.append(int(f[2]))
    return {k: (np.asarray(v[0], np.int64), np.asarray(v[1], np.int64))
            for k, v in by.items()}


def merge(starts, ends):
    """Sort and merge overlapping/bookended intervals -> sorted, disjoint arrays."""
    if len(starts) == 0:
        return starts, ends
    o = np.argsort(starts, kind="stable")
    s, e = starts[o], ends[o]
    # a new block begins wherever this start is beyond every end seen so far
    keep = np.empty(len(s), bool)
    keep[0] = True
    running = np.maximum.accumulate(e)
    keep[1:] = s[1:] > running[:-1]
    idx = np.flatnonzero(keep)
    # each block ends at the running max just before the next block starts
    ends_at = np.append(idx[1:] - 1, len(s) - 1)
    return s[idx], running[ends_at]


def intersect(a, b):
    """Pairwise intersection of two sorted DISJOINT interval sets.

    Same searchsorted candidate-window idiom as the overlap query in
    enrich_stranded.peaks_per_cat. Because b is disjoint and sorted, b_end is
    sorted too, so the window [lo, hi) contains exactly the genuine overlaps
    (b_start < a_end and b_end > a_start) with no filtering pass needed.
    """
    (as_, ae), (bs, be) = a, b
    if len(as_) == 0 or len(bs) == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    hi = np.searchsorted(bs, ae, side="left")
    lo = np.searchsorted(be, as_, side="right")
    cnt = np.maximum(hi - lo, 0)
    tot = int(cnt.sum())
    if tot == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    pid = np.repeat(np.arange(len(as_), dtype=np.int64), cnt)
    flat = np.repeat(lo, cnt) + (np.arange(tot, dtype=np.int64)
                                 - np.repeat(np.cumsum(cnt) - cnt, cnt))
    return np.maximum(as_[pid], bs[flat]), np.minimum(ae[pid], be[flat])


def brute_check(a, b, seed=0):
    """O(n*m) intersection on a small random slice, against the vectorised path."""
    rng = np.random.default_rng(seed)
    n = min(200, len(a[0]))
    if n == 0 or len(b[0]) == 0:
        return "skip"
    pick = np.sort(rng.choice(len(a[0]), n, replace=False))
    sub = (a[0][pick], a[1][pick])
    ref = []
    for s, e in zip(*sub):
        for bs, be in zip(*b):
            if bs < e and be > s:
                ref.append((max(s, bs), min(e, be)))
    got = list(zip(*intersect(sub, b)))
    return "PASS" if sorted(ref) == sorted((int(x), int(y)) for x, y in got) else "FAIL"


with open(CACHE, "rb") as fh:
    INDEX_CHROMS = set(pickle.load(fh)["index"].keys())

os.makedirs(OUT, exist_ok=True)

reps = defaultdict(dict)                 # (rbp, cell) -> {rep: path}
for d in SRC_DIRS:
    for p in sorted(glob.glob(f"{d}/*.bed")):
        m = NAME_RE.match(os.path.basename(p))
        if not m:
            print(f"  ! unparsed filename, skipped: {os.path.basename(p)}")
            continue
        reps[(m["rbp"], m["cell"])][m["rep"]] = p

written, dropped, off_index = [], [], defaultdict(int)
checked = False
for (rbp, cell), paths in sorted(reps.items()):
    if set(paths) != {"1", "2"}:
        dropped.append((rbp, cell, f"only replicate {'/'.join(sorted(paths))}"))
        continue
    r1, r2 = read_bed6(paths["1"]), read_bed6(paths["2"])
    if not r1 or not r2:
        dropped.append((rbp, cell, "a replicate is empty"))
        continue

    rows, npk = [], 0
    for key in sorted(set(r1) & set(r2)):
        a, b = merge(*r1[key]), merge(*r2[key])
        if not checked and key[0] == "chr21":
            print(f"  intersection brute-force check on {rbp}_{cell} chr21{key[1]}: "
                  f"{brute_check(a, b)}", flush=True)
            checked = True
        s, e = intersect(a, b)
        if len(s) == 0:
            continue
        chrom, strand = key
        if chrom not in INDEX_CHROMS:
            off_index[chrom] += len(s)
        npk += len(s)
        rows.extend((chrom, int(x), int(y), strand) for x, y in zip(s, e))

    if npk == 0:
        dropped.append((rbp, cell, "replicates share no peak"))
        continue
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    name = f"{rbp}_{cell}"
    with open(f"{OUT}/{name}.bed", "w") as fh:
        for c, s, e, sd in rows:
            fh.write(f"{c}\t{s}\t{e}\t{name}\t.\t{sd}\n")
    written.append((name, cell, npk))

for cell in sorted({c for _, c, _ in written}):
    n = [k for k in written if k[1] == cell]
    print(f"{cell}: {len(n)} RBPs, {sum(k[2] for k in n):,} reproducible peaks")
print(f"\nwrote {len(written)} files to {OUT}")

if off_index:
    tot = sum(off_index.values())
    print(f"\n{tot:,} peaks sit on {len(off_index)} contigs absent from the TE index "
          f"and will be silently ignored downstream:")
    for c, n in sorted(off_index.items(), key=lambda kv: -kv[1]):
        print(f"    {c:28s} {n:,}")

if dropped:
    print(f"\ndropped {len(dropped)} (RBP, cell) pairs:")
    for rbp, cell, why in dropped:
        print(f"    {rbp}_{cell:6s}  {why}")
