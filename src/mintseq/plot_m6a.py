"""Plot the L1 / Alu / SVA subfamilies whose nascent transcripts carry more m6A
than the transcriptome average.

Form: "compare magnitude across discrete items", element class as identity
(3 categories), orientation as a facet -> horizontal bars, small-multiple panels
for sense / antisense. Matches plot_enrichment.py in the eCLIP tree.

The measure is signed (a subfamily can be under-methylated), so bars run both
ways off a zero baseline. Polarity is carried by *direction*, not hue, which
leaves color free to carry element identity -- mixing a diverging ramp with a
categorical one would put two meanings on the same channel.

Palette: the same validated categorical slots as plot_enrichment.py
(blue / orange / violet). Static figure for print, light surface only.

Usage: python plot_m6a.py [results_dir]
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Patch

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJ, "results", "mintseq")
SRC = os.path.join(RESULTS, "te_m6a_enrichment.tsv")

TOP_N, MIN_COPIES, QMAX = 15, 20, 0.05

SURFACE = "#fcfcfb"
INK, INK_2, INK_MUTED = "#0b0b0b", "#52514e", "#8a8880"
GRID = "#e6e5e1"
COLOR = {"Alu": "#2a78d6", "L1": "#eb6834", "SVA": "#4a3aa7"}   # fixed slot order


def load():
    with open(SRC) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    keep = [r for r in rows
            if int(r["n_copies"]) >= MIN_COPIES and float(r["q"]) < QMAX]
    out, totals = {}, {}
    for o in ("sense", "antisense"):
        sel = [r for r in keep if r["orientation"] == o]
        sel.sort(key=lambda r: float(r["delta_vs_bg"]), reverse=True)
        out[o] = sel[:TOP_N]
        totals[o] = {
            "up": sum(1 for r in sel if float(r["delta_vs_bg"]) > 0),
            "down": sum(1 for r in sel if float(r["delta_vs_bg"]) < 0),
        }
    return out, totals


def rounded_barh(ax, y, w, h, color, r):
    """Bar with a 4px-rounded data-end, square against the x=0 baseline.

    Handles w < 0: the rounded end follows the data, so a depleted subfamily is
    rounded on its left edge and square where it meets zero.
    """
    s = 1.0 if w >= 0 else -1.0
    aw = abs(w)
    r = min(r, h / 2, aw / 2) * s
    v = [(0, y - h / 2), (w - r, y - h / 2), (w, y - h / 2), (w, y - h / 2 + abs(r)),
         (w, y + h / 2 - abs(r)), (w, y + h / 2), (w - r, y + h / 2),
         (0, y + h / 2), (0, y - h / 2)]
    c = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3,
         Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(v, c), fc=color, ec="none", zorder=3))


data, totals = load()
if not any(data.values()):
    sys.exit("no categories passed q < %g and n_copies >= %d" % (QMAX, MIN_COPIES))

span = max(abs(float(r["delta_vs_bg"])) for o in data for r in data[o])
# Right margin holds the per-bar n / %-positive column, so it has to be wide
# enough that the text never runs back over the bars.
xlo, xhi = -span * 0.12, span * 1.62
XCOL = span * 1.60          # right-aligned annotation column

fig, axes = plt.subplots(1, 2, figsize=(14.5, 7.2), sharex=True)
fig.patch.set_facecolor(SURFACE)

for ax, orient in zip(axes, ("sense", "antisense")):
    rows = data[orient]
    ax.set_facecolor(SURFACE)
    labels = []
    for i, r in enumerate(rows):
        y = len(rows) - 1 - i
        d = float(r["delta_vs_bg"])
        rounded_barh(ax, y, d, 0.68, COLOR[r["family"]], r=span * 0.008)
        off = span * 0.014
        ax.text(d + (off if d >= 0 else -off), y, "%+.2f" % d, va="center",
                ha="left" if d >= 0 else "right", fontsize=8.5, color=INK_2)
        labels.append(r["subfamily"])
        ax.text(XCOL, y, "n=%s  ·  %.0f%% pos"
                % (r["n_copies"], 100 * float(r["frac_positive"])),
                va="center", ha="right", fontsize=6.8, color=INK_MUTED)

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels[::-1], fontsize=8.6, color=INK)
    ax.set_ylim(-0.75, len(rows) - 0.25)
    ax.set_xlim(xlo, xhi)
    ax.set_xlabel("Δ median log2(m6A IP / nascent input) vs all TE copies",
                  fontsize=9, color=INK_2)
    t = totals[orient]
    ax.set_title("%s   —   top %d of %d enriched subfamilies (%d depleted)"
                 % (orient, len(rows), t["up"], t["down"]),
                 fontsize=10, color=INK, loc="left", pad=10)
    ax.xaxis.grid(True, color=GRID, lw=0.8, zorder=0)           # recessive grid
    ax.set_axisbelow(True)
    ax.yaxis.grid(False)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(axis="x", colors=INK_2, labelsize=8.5, length=0)
    ax.tick_params(axis="y", length=0)
    ax.axvline(0, color=SURFACE, lw=2.2, zorder=4)
    ax.axvline(0, color=INK_MUTED, lw=0.9, ls=(0, (3, 3)), zorder=5)
    ax.text(0, len(rows) - 0.45, " transcriptome average", fontsize=7,
            color=INK_MUTED, va="center")

present = [e for e in ("Alu", "L1", "SVA")
           if any(r["family"] == e for o in data for r in data[o])]
fig.legend(handles=[Patch(fc=COLOR[e], ec="none", label=e) for e in present],
           loc="lower center", ncol=len(present), frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, -0.005), labelcolor=INK)

fig.suptitle("m6A on nascent transcripts of L1 / Alu / SVA subfamilies (K562)",
             fontsize=13.5, color=INK, x=0.011, ha="left", y=0.985)
fig.text(0.011, 0.945,
         "MINT-Seq (m6A IP on nascent RNA) over matched TT-Seq input, 2 replicates pooled, "
         "per TE copy and length-normalised.",
         fontsize=8.6, color=INK_2, ha="left")
fig.text(0.011, 0.917,
         "Bars are the shift in median per-copy log2 ratio relative to all tested copies "
         "in the same orientation; Mann-Whitney q < %g, ≥ %d copies over the input-coverage floor."
         % (QMAX, MIN_COPIES),
         fontsize=8.6, color=INK_2, ha="left")

fig.tight_layout(rect=[0, 0.045, 1, 0.905])
fig.subplots_adjust(wspace=0.42)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(RESULTS, "te_m6a_enrichment.%s" % ext),
                dpi=200, facecolor=SURFACE)
print("wrote %s/te_m6a_enrichment.png and .pdf" % RESULTS)
for o in ("sense", "antisense"):
    if data[o]:
        top = data[o][0]
        print("  %-9s %d bars, top = %s %+.3f (n=%s, q=%.3g)"
              % (o, len(data[o]), top["subfamily"], float(top["delta_vs_bg"]),
                 top["n_copies"], float(top["q"])))
