"""Stage 1 of the 3D analysis: which compartment does each TE subfamily live in?

    python src/hic/te_compartments.py MCOOL --tss data/hic/tss.bed [--resolution 250000]

This is the cheap, mappability-tolerant question. Bins are 100-250 kb, so the
unique-alignment suppression that hides young subfamilies from peak-level assays
(CLAUDE.md: L1H at 0.36% coverage, zero hits) barely bites -- a 250 kb bin gets
its coverage from the 99.9% of sequence that is not the TE.

WHY THE PHASING TRACK IS REQUIRED
The cis eigendecomposition returns eigenvectors whose SIGN is arbitrary and
independently arbitrary per chromosome. Without phasing, "E1 > 0" means the A
compartment on some chromosomes and the B compartment on others, and a
genome-wide A-fraction per subfamily is then meaningless noise. cooltools phases
against a track known to track activity; we use TSS density rather than GC.

That choice is not cosmetic here. GC phasing would be actively dangerous for this
particular analysis: Alu is GC-rich and L1 is AT-rich BY COMPOSITION, so phasing
the compartment call on GC and then asking whether Alu prefers the A compartment
partially tests the phasing track against itself. TSS density is not free of that
worry either (Alu is gene-dense by isochore), but it is a downstream biological
quantity rather than the same nucleotide statistic, so the circularity is weaker.
Interpret a positive Alu result with that caveat live; the honest test of
subfamily effects is the tissue-interaction design, not this marginal one.

Output: results/hic/<sample>_te_compartments.tsv plus a masking report.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hic_common import (PROJ, assign_bins, load_te_index, make_view,
                        masking_report, open_cooler, reconcile_chroms, te_table)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("mcool")
    p.add_argument("--tss", required=True,
                   help="TSS BED from fetch_gencode_tss.sh; phases the eigenvector")
    p.add_argument("--resolution", type=int, default=250000)
    p.add_argument("--out", default=os.path.join(PROJ, "results", "hic"))
    p.add_argument("--elements", default="L1,Alu,SVA")
    p.add_argument("--min-copies", type=int, default=50,
                   help="drop subfamilies with fewer copies; extreme fractions on a "
                        "handful of copies are noise, same rationale as the "
                        "loci_in_cat >= 10 floor in the eCLIP plots")
    p.add_argument("--nproc", type=int, default=4,
                   help="cooltools workers. Do NOT set this to your core count: "
                        "each worker holds a chromosome-scale working set")
    return p.parse_args()


def load_tss(path, clr, chroms):
    """TSS BED -> per-bin count, aligned to clr.bins() exactly."""
    tss = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["chrom", "start", "end"], comment="#")
    tss = tss[tss["chrom"].isin(chroms)]
    tss = assign_bins(tss.assign(mid=tss["start"]), clr)
    tss = tss[tss["bin_id"] >= 0]
    counts = np.bincount(tss["bin_id"].to_numpy(),
                         minlength=clr.info["nbins"]).astype(float)
    print(f"  TSS: {len(tss)} in {len(chroms)} chroms -> phasing track")
    return counts


def main():
    a = parse_args()
    os.makedirs(a.out, exist_ok=True)
    sample = os.path.basename(a.mcool).replace(".mcool", "").replace(".cool", "")

    import cooltools

    print(f"[te_compartments] {sample}")
    clr, binsize = open_cooler(a.mcool, a.resolution)

    index, subfam_names, subfam_elem = load_te_index()
    chroms, _ = reconcile_chroms(index, clr)
    if not chroms:
        raise SystemExit("no chromosomes shared between TE index and cooler")
    view = make_view(clr, chroms)

    bins = clr.bins()[:]
    phasing = bins[["chrom", "start", "end"]].copy()
    phasing["tss_density"] = load_tss(a.tss, clr, set(chroms))

    print(f"  eigendecomposition at {binsize} bp (nproc={a.nproc}) ...")
    _, eigvecs = cooltools.eigs_cis(
        clr, phasing_track=phasing, view_df=view, n_eigs=3)

    # eigs_cis returns rows only for bins inside the view, so a positional
    # assignment onto the full bin table would be off by every excluded scaffold.
    # Merge on coordinates; bins outside the view get E1 = NaN and drop out below.
    bins = clr.bins()[:].merge(
        eigvecs[["chrom", "start", "end", "E1"]], on=["chrom", "start", "end"],
        how="left")

    te = te_table(index, subfam_names, subfam_elem, chroms,
                  elements=set(a.elements.split(",")), min_copies=a.min_copies)
    te = assign_bins(te, clr)
    te = te[te["bin_id"] >= 0]
    print(f"  TE copies placed in bins: {len(te)} "
          f"across {te['subfam'].nunique()} subfamilies")

    mask_rep = masking_report(te, bins)
    mask_path = os.path.join(a.out, f"{sample}_te_bin_masking.tsv")
    mask_rep.to_csv(mask_path, sep="\t", index=False, float_format="%.4f")
    worst = mask_rep.iloc[0]
    print(f"  masking: worst subfamily {worst['subfam']} "
          f"{worst['frac_masked']:.1%} of copies in NaN-weight bins "
          f"-> {os.path.basename(mask_path)}")

    e1_of_te = bins["E1"].to_numpy()[te["bin_id"].to_numpy()]
    te = te.assign(E1=e1_of_te)
    scored = te[np.isfinite(te["E1"])]

    valid_bins = bins[np.isfinite(bins["E1"])]
    bg_a = float((valid_bins["E1"] > 0).mean())
    print(f"  background: {bg_a:.1%} of valid bins are A")

    rows = []
    for (elem, subfam), g in scored.groupby(["elem", "subfam"]):
        n = len(g)
        rows.append(dict(
            elem=elem, subfam=subfam,
            n_copies=int((te["subfam"] == subfam).sum()),
            n_scored=n,
            frac_masked=float(mask_rep.loc[mask_rep["subfam"] == subfam,
                                           "frac_masked"].iloc[0]),
            frac_A=float((g["E1"] > 0).mean()),
            frac_A_minus_bg=float((g["E1"] > 0).mean()) - bg_a,
            mean_E1=float(g["E1"].mean()),
            median_E1=float(g["E1"].median()),
        ))
    out = pd.DataFrame(rows).sort_values("frac_A_minus_bg", ascending=False)
    out_path = os.path.join(a.out, f"{sample}_te_compartments.tsv")
    out.to_csv(out_path, sep="\t", index=False, float_format="%.4f")

    print(f"\n  wrote {out_path}  ({len(out)} subfamilies)")
    print("\n  most A-shifted:")
    print(out.head(5)[["elem", "subfam", "n_scored", "frac_A",
                       "frac_A_minus_bg"]].to_string(index=False))
    print("\n  most B-shifted:")
    print(out.tail(5)[["elem", "subfam", "n_scored", "frac_A",
                       "frac_A_minus_bg"]].to_string(index=False))


if __name__ == "__main__":
    main()
