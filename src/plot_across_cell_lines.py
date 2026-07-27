"""Plot the most reproducible RBP x subfamily hits across all 7 cell lines.

Form: the job is "compare magnitude over a grid" (hit x cell line) -> heatmap with
a sequential one-hue ramp. Sequential, not ordinal: per the palette reference the
full 100->700 blue range is allowed for continuous magnitude, and the lightest
step is *meant* to recede toward the surface ("near zero"). The ordinal 2:1
light-end floor applies to discrete ordered marks, not to a heatmap.

Two variables per cell, so neither rests on color alone:
  fill  = fold enrichment (sequential blue, binned)
  hatch = NOT a robust hit in that cell line (did not survive the permutation null,
          or fold/copies below cut)
  text  = the fold value, printed in every cell

Robustness is judged on the permutation null, not the hypergeometric screen: a cell
counts as a hit only if (rbp, subfamily, orientation, cell) appears in
rbp_te_enrichment_permutation_significant.tsv. The fold/copy numbers still come from
the per-cell stranded tables, because the grid must show a value for every cell --
including the ones where the hit does not replicate, which were never permuted at
all (only screened hits are). Those cells are correctly hatched.

Rows are ranked by how many cell lines replicate the hit, then by median fold,
so the plot answers "what holds up everywhere?" rather than "what is biggest once?".

Rows are faceted into Alu / L1 / SVA blocks and ranked within each. Under the
per-copy dedupe of enrich_stranded.py the surviving hits are all Alu (L1 and SVA
rest on too few distinct copies), so in practice only the Alu block renders; the
faceting is kept so the figure stays correct if the inputs ever change.

Usage: python plot_across_cell_lines.py [results_dir]
"""
import csv, glob, math, os, sys
from statistics import median
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = sys.argv[1] if len(sys.argv) > 1 else f"{PROJ}/results/eclip"

PER_ELEMENT, MIN_LOCI, MIN_FOLD, QCUT = 7, 10, 2.0, 0.05
ELEMENTS = ("Alu", "L1", "SVA")
GAP = 0.7                    # blank rows between element blocks

SURFACE = "#fcfcfb"
INK, INK_2, INK_MUTED = "#0b0b0b", "#52514e", "#8a8880"
GRID = "#e6e5e1"

# sequential blue, light -> dark (palette steps 100/200/300/400/500/650)
BINS = [(0, 1), (1, 2), (2, 4), (4, 6), (6, 10), (10, 1e9)]
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#104281"]
BIN_LABEL = ["<1", "1–2", "2–4", "4–6", "6–10", "≥10"]
DARK_TEXT_FROM = 3          # ramp index at which cell text flips to white
NODATA = "#e9e8e4"          # off-ramp neutral: the category has no pooled peaks at all


def bin_of(f):
    """Ramp index, or None when fold is not a finite number.

    fold is NaN when the (subfamily, orientation) category has zero pooled peaks
    across all RBPs (0/0), and +inf when only this RBP binds it (background rate 0).
    Both must be handled explicitly: every `lo <= nan < hi` is False, so a naive
    loop falls through and paints "no data" as the darkest, most-enriched step.
    """
    if math.isnan(f):
        return None
    if math.isinf(f):
        return len(BINS) - 1
    for i, (lo, hi) in enumerate(BINS):
        if lo <= f < hi:
            return i
    return len(BINS) - 1


# ---- hits that survived the permutation null: (rbp, subfam, orient, cell)
PSIG = set()
with open(f"{RESULTS}/rbp_te_enrichment_permutation_significant.tsv") as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if float(r["q_perm"]) < QCUT:
            PSIG.add((r["rbp"], r["subfamily"], r["orientation"], r["cell_line"]))


def robust(r):
    """Survived the permutation null in this cell, on enough copies, at enough fold.

    Rows absent from PSIG include both permuted-and-rejected hits and the far larger
    set that never passed the hypergeometric screen, so were never permuted; neither
    is a robust hit, and both hatch.
    """
    return ((r["rbp"], r["subfamily"], r["orientation"], r["cell_line"]) in PSIG
            and int(r["loci_in_cat"]) >= MIN_LOCI
            and float(r["fold_enrichment"]) >= MIN_FOLD)


# ---- load every stranded per-cell table into grid[(rbp, subfam, orient)][cell]
grid, cells = {}, []
for path in sorted(glob.glob(f"{RESULTS}/rbp_te_enrichment_stranded_*.tsv")):
    if path.endswith("_significant.tsv"):
        continue
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cell = r["cell_line"]
            key = (r["rbp"], r["subfamily"], r["orientation"], r["element"])
            grid.setdefault(key, {})[cell] = (
                float(r["fold_enrichment"]), float(r["q_BH"]), int(r["loci_in_cat"]), robust(r))
    cells.append(os.path.basename(path)[len("rbp_te_enrichment_stranded_"):-4])
cells = sorted(set(cells))

# ---- rank by replication (n cell lines robust), then median fold, WITHIN each element
ranked = []
for key, per in grid.items():
    hits = [v for v in per.values() if v[3]]
    if not hits:
        continue
    ranked.append((len(hits), median(h[0] for h in hits), key))
ranked.sort(key=lambda t: (t[0], t[1]), reverse=True)
blocks = [(e, [r for r in ranked if r[2][3] == e][:PER_ELEMENT]) for e in ELEMENTS]
blocks = [(e, b) for e, b in blocks if b]

total = sum(len(b) for _, b in blocks) + GAP * (len(blocks) - 1)

# ---- draw
fig, ax = plt.subplots(figsize=(13.2, 10.2))
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)

yticks, ylabels, y = [], [], total
for bi, (elem, block) in enumerate(blocks):
    ytop = y
    for nhit, medf, key in block:
        y -= 1
        for ci, cell in enumerate(cells):
            v = grid[key].get(cell)
            fold, q, n, ok = v if v else (float("nan"), 1.0, 0, False)
            b = bin_of(fold)
            fc = NODATA if b is None else RAMP[b]
            ax.add_patch(Rectangle((ci, y), 1, 1, fc=fc, ec=SURFACE, lw=2, zorder=2))
            if not ok:   # texture: replication failure is not encoded by color alone
                # hatch must contrast with the cell it sits on, not always be white
                hc = "#ffffff" if (b is not None and b >= DARK_TEXT_FROM - 1) else INK_MUTED
                ax.add_patch(Rectangle((ci, y), 1, 1, fc="none", ec=hc, lw=0,
                                       hatch="////", alpha=0.5, zorder=3))
            if b is None:
                txt, tc = "n/a", INK_MUTED
            elif math.isinf(fold):
                txt, tc = "∞", "#ffffff"
            else:
                txt, tc = f"{fold:.1f}", ("#ffffff" if b >= DARK_TEXT_FROM else INK)
            ax.text(ci + 0.5, y + 0.5, txt, ha="center", va="center", fontsize=8.2,
                    zorder=4, color=tc, fontweight="bold" if ok else "normal")
        yticks.append(y + 0.5)
        ylabels.append(f"{key[0]} · {key[1]}")
        ax.text(len(cells) + 0.15, y + 0.5, f"{key[2]:9s}  {nhit}/{len(cells)}",
                fontsize=7.8, color=INK_2, va="center", family="monospace")
    # element label + block separator
    ax.text(-0.34, (ytop + y) / 2, elem, transform=ax.get_yaxis_transform(),
            fontsize=11, color=INK, fontweight="bold", rotation=90,
            va="center", ha="center")
    if bi < len(blocks) - 1:
        ax.axhline(y - GAP / 2, color=GRID, lw=1, zorder=1)
        y -= GAP

ax.set_xlim(0, len(cells)); ax.set_ylim(0, total)
ax.set_xticks([c + 0.5 for c in range(len(cells))])
ax.set_xticklabels(cells, fontsize=9, color=INK)
ax.xaxis.set_ticks_position("top")
ax.tick_params(axis="both", length=0)
ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=8.8, color=INK)
for s in ax.spines.values():
    s.set_visible(False)
ax.text(len(cells) + 0.15, total + 0.20, "orientation  replicates",
        fontsize=7.4, color=INK_MUTED, va="center", family="monospace")

legend = [Patch(fc=RAMP[i], ec=SURFACE, label=BIN_LABEL[i]) for i in range(len(RAMP))]
legend.append(Patch(fc=NODATA, ec=SURFACE, hatch="////", label="not a robust hit"))
legend.append(Patch(fc=NODATA, ec=GRID, label="n/a — no pooled peaks"))
fig.legend(handles=legend, loc="lower center", ncol=len(legend), frameon=False,
           fontsize=8.6, bbox_to_anchor=(0.5, 0.005), labelcolor=INK,
           title="fold enrichment", title_fontsize=8.6)

fig.suptitle("Reproducibility of the top RBP–subfamily hits across cell lines",
             fontsize=13.5, color=INK, x=0.008, ha="left", y=0.985)
fig.text(0.008, 0.945,
         f"Top {PER_ELEMENT} per element, ranked by how many of the {len(cells)} cell lines "
         f"replicate the hit (permutation q < {QCUT}, fold ≥ {MIN_FOLD:g}, ≥ {MIN_LOCI} copies), then by "
         "median fold. Bold = robust in that cell line; hatched = not.",
         fontsize=8.6, color=INK_2, ha="left")

fig.tight_layout(rect=[0, 0.055, 0.995, 0.905])
for ext in ("png", "pdf"):
    fig.savefig(f"{RESULTS}/top_hits_across_cell_lines.{ext}", dpi=200, facecolor=SURFACE)
print(f"wrote {RESULTS}/top_hits_across_cell_lines.png and .pdf")
print(f"  {len(grid):,} (rbp, subfamily, orientation) combos; {len(ranked):,} robust in >=1 cell line")
for elem, block in blocks:
    print(f"  {elem}: {len(block)} rows, best = {block[0][2][0]} {block[0][2][1]} "
          f"({block[0][2][2]}) replicates {block[0][0]}/{len(cells)}, median fold {block[0][1]:.1f}")
