"""Plot the top RBP x TE-subfamily enrichments that survive the permutation null.

Form: the job is "compare magnitude across discrete items", with element class as
identity (3 categories) and orientation as a facet -> horizontal bars, sorted by
fold, small-multiple panels for sense / antisense.

Source is the permutation-significant table, not the hypergeometric one: the
hypergeometric is only a screen (it assumes independent peak placement), so every
bar here has already been defended against the circular-shift null at q_perm < 0.05.
On top of that we require fold >= MIN_FOLD and at least MIN_LOCI distinct supporting
copies. Counting is per TE copy (see enrich_stranded.py), so after this filter only
high-copy Alu subfamilies survive.

The bar length is fold_perm (obs / mean loci under the shift null), NOT fold_pooled
(obs vs the other 170 RBPs). fold_perm is the effect size matched to the q_perm this
plot filters on -- pairing a q from one null with a fold from another was the bug this
fixes. The two differ by a median 1.62x (see enrich_permutation.py), so bars here run
longer than the pooled folds quoted in the stranded tables, and the ranking differs.

Palette: categorical slots blue / orange / violet, validated with the data-viz
validator (all checks pass on the light surface; worst adjacent CVD dE 96.7).
This is a static figure for print, so it commits to the light surface only.

Usage: python plot_enrichment.py [results_dir]
"""
import csv, math, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = sys.argv[1] if len(sys.argv) > 1 else f"{PROJ}/results/eclip"
SRC = f"{RESULTS}/rbp_te_enrichment_permutation_significant.tsv"

TOP_N, MIN_LOCI, MIN_FOLD = 15, 10, 2.0

SURFACE = "#fcfcfb"
INK, INK_2, INK_MUTED = "#0b0b0b", "#52514e", "#8a8880"
GRID = "#e6e5e1"
COLOR = {"Alu": "#2a78d6", "L1": "#eb6834", "SVA": "#4a3aa7"}   # fixed slot order


def load():
    with open(SRC) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    # fold_perm is inf when the shift null never touches the category (exp = 0). Such
    # a bar has no length to draw and would take xmax to inf, so drop it here; with
    # MIN_LOCI >= 10 copies backing every hit, none currently occur.
    keep = [r for r in rows
            if int(r["obs_loci"]) >= MIN_LOCI and math.isfinite(float(r["fold_perm"]))
            and float(r["fold_perm"]) >= MIN_FOLD]
    out, totals = {}, {}
    for o in ("sense", "antisense"):
        sel = [r for r in keep if r["orientation"] == o]
        sel.sort(key=lambda r: float(r["fold_perm"]), reverse=True)
        out[o] = sel[:TOP_N]
        # The top-N by fold is not representative of the hit *counts*; the panel
        # title reports the per-element totals. (After locus dedup these are all Alu.)
        totals[o] = {e: sum(1 for r in sel if r["element"] == e) for e in ("Alu", "L1", "SVA")}
    return out, totals


def rounded_barh(ax, y, w, h, color, r):
    """Bar with 4px-rounded data-end, square against the x=0 baseline."""
    r = min(r, h / 2, w / 2)
    v = [(0, y - h / 2), (w - r, y - h / 2), (w, y - h / 2), (w, y - h / 2 + r),
         (w, y + h / 2 - r), (w, y + h / 2), (w - r, y + h / 2), (0, y + h / 2), (0, y - h / 2)]
    c = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3,
         Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(v, c), fc=color, ec="none", zorder=3))


data, totals = load()
xmax = max(float(r["fold_perm"]) for o in data for r in data[o]) * 1.18

fig, axes = plt.subplots(1, 2, figsize=(14.5, 7.2), sharex=True)
fig.patch.set_facecolor(SURFACE)

for ax, orient in zip(axes, ("sense", "antisense")):
    rows = data[orient]
    ax.set_facecolor(SURFACE)
    labels = []
    for i, r in enumerate(rows):
        y = len(rows) - 1 - i
        fold = float(r["fold_perm"])
        rounded_barh(ax, y, fold, 0.68, COLOR[r["element"]], r=xmax * 0.006)
        ax.text(fold + xmax * 0.012, y, f"{fold:.1f}×", va="center", ha="left",
                fontsize=8.5, color=INK_2)                      # direct label, text token
        labels.append(f"{r['rbp']} · {r['subfamily']}")
        ax.text(-xmax * 0.012, y - 0.30, f"{r['cell_line']}  n={r['obs_loci']} copies",
                va="center", ha="right", fontsize=6.6, color=INK_MUTED)

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels[::-1], fontsize=8.6, color=INK)
    ax.set_ylim(-0.75, len(rows) - 0.25)
    ax.set_xlim(0, xmax)
    ax.set_xlabel("Fold enrichment vs the circular-shift null", fontsize=9, color=INK_2)
    tot = totals[orient]
    ax.set_title(f"{orient}   —   top {len(rows)} of "
                 f"{sum(tot.values())} robust hits "
                 f"(Alu {tot['Alu']}, L1 {tot['L1']}, SVA {tot['SVA']})",
                 fontsize=10, color=INK, loc="left", pad=10)
    ax.xaxis.grid(True, color=GRID, lw=0.8, zorder=0)           # recessive grid
    ax.set_axisbelow(True)
    ax.yaxis.grid(False)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(axis="x", colors=INK_2, labelsize=8.5, length=0)
    ax.tick_params(axis="y", length=0)
    # above the bars (zorder 4) so the annotated null line is actually visible
    ax.axvline(1, color=SURFACE, lw=2.2, zorder=4)
    ax.axvline(1, color=INK_MUTED, lw=0.9, ls=(0, (3, 3)), zorder=5)
    ax.text(1, len(rows) - 0.45, " fold = 1 (no enrichment)", fontsize=7,
            color=INK_MUTED, va="center")

present = [e for e in ("Alu", "L1", "SVA") if any(r["element"] == e for o in data for r in data[o])]
fig.legend(handles=[Patch(fc=COLOR[e], ec="none", label=e) for e in present],
           loc="lower center", ncol=len(present), frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, -0.005), labelcolor=INK)

fig.suptitle("Top RBP enrichments on L1 / Alu / SVA subfamilies", fontsize=13.5,
             color=INK, x=0.011, ha="left", y=0.985)
fig.text(0.011, 0.938,
         f"Hits that survive the circular-shift permutation null: q_perm < 0.05, fold ≥ {MIN_FOLD:g}, "
         f"≥ {MIN_LOCI} distinct copies. Fold is observed copies ÷ mean copies hit when the RBP's whole "
         "peak set is rigidly shifted to a random phase inside the cell's peak footprint.",
         fontsize=8.6, color=INK_2, ha="left")

fig.tight_layout(rect=[0, 0.045, 1, 0.925])
fig.subplots_adjust(wspace=0.42)
for ext in ("png", "pdf"):
    fig.savefig(f"{RESULTS}/top_enrichments.{ext}", dpi=200, facecolor=SURFACE)
print(f"wrote {RESULTS}/top_enrichments.png and .pdf")
for o in ("sense", "antisense"):
    n = len(data[o])
    top = data[o][0]
    print(f"  {o:9s} {n} bars, top = {top['rbp']} {top['subfamily']} "
          f"{float(top['fold_perm']):.2f}x ({top['cell_line']})")
