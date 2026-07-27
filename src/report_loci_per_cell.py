"""Report the distinct TE loci that drive the locus-deduped enrichment, per cell.

`obs_loci` in the enrichment table is distinct copies per (RBP, category); it
double-counts a copy bound by two different RBPs. This unions the *actual* locus
identities (chrom, element-index) across every significant hit in a cell, so the
per-cell total is genuinely distinct copies -- how many independent TE insertions
carry the enriched signal, not how many binding events.

Reads results/rbp_te_enrichment_permutation_significant.tsv -- the hits that survive
the circular-shift null, not the raw hypergeometric screen -- re-opens only the BED
files of RBPs that have a significant hit, and recollects loci for exactly the
categories they were significant in. Writes results/distinct_loci_per_cell.tsv.

Usage: python report_loci_per_cell.py [data_dir] [results_dir]
"""
import numpy as np, pickle, os, sys, csv, glob
from collections import defaultdict
from tecommon import read_bed, by_chrom

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = sys.argv[1] if len(sys.argv) > 1 else f"{PROJ}/data/eclip"
RESULTS = sys.argv[2] if len(sys.argv) > 2 else f"{PROJ}/results/eclip"
SIG = f"{RESULTS}/rbp_te_enrichment_permutation_significant.tsv"

with open(f"{PROJ}/cache/te_index.pkl", "rb") as fh:
    blob = pickle.load(fh)
INDEX, SUBFAM, SELEM = blob["index"], blob["subfam_names"], blob["subfam_elem"]
SID = {s: i for i, s in enumerate(SUBFAM)}


def loci_hit(bed, want):
    """Distinct (chrom, element-index) touched by this RBP in categories `want`.

    want: set of cat = sid*2 + anti. Returns set of (chrom, flat) global locus ids.
    Mirrors enrich_stranded.peaks_per_cat's overlap, but keeps locus identities.
    """
    out = set()
    for c, (ps, pe, pstr) in bed.items():
        ix = INDEX.get(c)
        if ix is None:
            continue
        ts, te, sid, tstr, L = ix["starts"], ix["ends"], ix["sid"], ix["strand"], ix["maxlen"]
        hi = np.searchsorted(ts, pe, side="left")
        lo = np.searchsorted(ts, ps - L, side="left")
        cnt = hi - lo
        tot = int(cnt.sum())
        if tot == 0:
            continue
        pid = np.repeat(np.arange(len(ps), dtype=np.int64), cnt)
        flat = np.repeat(lo, cnt) + (np.arange(tot, dtype=np.int64)
                                     - np.repeat(np.cumsum(cnt) - cnt, cnt))
        keep = (te[flat] > ps[pid]) & (ts[flat] < pe[pid])
        pid, flat = pid[keep], flat[keep]
        cat = sid[flat].astype(np.int64) * 2 + (pstr[pid] != tstr[flat]).astype(np.int64)
        for f, k in zip(flat.tolist(), cat.tolist()):
            if k in want:
                out.add((c, int(f)))
    return out


# ---- significant categories per (cell, rbp)
sig = defaultdict(lambda: defaultdict(set))          # cell -> rbp -> {cat}
n_hits = defaultdict(int)
with open(SIG) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        cat = SID[r["subfamily"]] * 2 + (0 if r["orientation"] == "sense" else 1)
        sig[r["cell_line"]][r["rbp"]].add(cat)
        n_hits[r["cell_line"]] += 1

rows = []
for cell in sorted(sig):
    per_elem = defaultdict(set)                       # element -> {global locus}
    all_loci = set()
    for rbp, want in sig[cell].items():
        path = glob.glob(f"{DATA}/{rbp}_{cell}.bed")
        if not path:
            print(f"  !! missing BED for {rbp}_{cell}", file=sys.stderr); continue
        loci = loci_hit(by_chrom(*read_bed(path[0])), want)
        all_loci |= loci
        for c, f in loci:
            # element class from the subfamily of that locus copy
            per_elem[SELEM[int(INDEX[c]["sid"][f])]] |= {(c, f)}
    rows.append((cell, len(sig[cell]), n_hits[cell], len(all_loci),
                 len(per_elem.get("Alu", ())), len(per_elem.get("L1", ())),
                 len(per_elem.get("SVA", ()))))
    print(f"{cell:9s}  {len(sig[cell]):3d} RBPs  {n_hits[cell]:3d} hits  "
          f"{len(all_loci):6d} distinct loci  "
          f"(Alu {len(per_elem.get('Alu', ())):5d}  L1 {len(per_elem.get('L1', ())):4d}  "
          f"SVA {len(per_elem.get('SVA', ())):4d})")

out = f"{RESULTS}/distinct_loci_per_cell.tsv"
with open(out, "w") as fh:
    fh.write("cell_line\tn_sig_rbps\tn_sig_hits\tdistinct_loci\tloci_Alu\tloci_L1\tloci_SVA\n")
    for r in rows:
        fh.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % r)
print(f"\nwrote {out}")
