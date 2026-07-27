"""Primitives shared by the 3D-genome stage: cooler/TE-index reconciliation.

Imported by te_compartments.py and te_tss_contacts.py. Those run as
`python src/hic/<script>.py`, so `src/hic/` is sys.path[0] and a plain
`from hic_common import ...` resolves without packaging -- same convention as
src/tecommon.py.

RUNS UNDER deps-3d, NOT deps. cooltools caps at numpy <2 in every packaged
version and ../deps pins numpy ==2.4.6 for the eCLIP permutation stream. See
deps-3d/pixi.toml.
"""
import os
import pickle

import numpy as np
import pandas as pd

PROJ = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TE_INDEX = os.path.join(PROJ, "cache", "te_index.pkl")


def load_te_index(path=TE_INDEX):
    """cache/te_index.pkl -> (index, subfam_names, subfam_elem).

    index is {chrom: {starts, ends, sid, strand, maxlen}} with BED 0-based
    half-open coordinates -- build_te_index.py already did the GTF conversion, so
    nothing here re-converts. subfam_elem[i] is 'L1' | 'Alu' | 'SVA' for
    subfam_names[i].
    """
    with open(path, "rb") as fh:
        d = pickle.load(fh)
    return d["index"], d["subfam_names"], d["subfam_elem"]


def reconcile_chroms(index, clr, verbose=True):
    """Chromosomes usable by BOTH the TE index and the cooler, plus a loss report.

    This is the mirror image of the INDEX.get(c) trap recorded in CLAUDE.md. There,
    peaks on chromosomes absent from the TE index inflated the denominator but could
    never hit a TE. Here the loss runs the other way: the TE annotation carries 631
    sequences (deliberately -- "canonical" in the filename means canonical
    SUBFAMILIES, not canonical chromosomes), while a 4DN cooler is built against a
    fixed chromsizes file carrying ~24 plus a handful of scaffolds. Every TE copy on
    a sequence the cooler does not know silently vanishes from the numerator.

    That loss is real but it is NOT uniform across subfamilies -- young elements are
    over-represented on unplaced scaffolds -- so it has to be a printed number, not
    an assumption. Returns (chroms, report_df) and leaves the decision to the caller.
    """
    cool_chroms = set(clr.chromnames)
    te_chroms = set(index.keys())
    keep = sorted(te_chroms & cool_chroms, key=lambda c: clr.chromnames.index(c))
    dropped = te_chroms - cool_chroms

    n_kept = sum(len(index[c]["starts"]) for c in keep)
    n_drop = sum(len(index[c]["starts"]) for c in dropped)
    if verbose:
        print(f"  chroms: {len(keep)} shared, {len(dropped)} TE-only "
              f"(cooler has {len(cool_chroms)}, TE index has {len(te_chroms)})")
        print(f"  elements: {n_kept} usable, {n_drop} dropped "
              f"({100 * n_drop / max(1, n_kept + n_drop):.2f}%)")
    return keep, dropped


def te_table(index, subfam_names, subfam_elem, chroms, elements=None,
             subfamilies=None, min_copies=0):
    """Flatten the per-chromosome index into one tidy frame of TE loci.

    Columns: chrom, start, end, mid, strand, sid, subfam, elem. Coordinates stay
    BED 0-based half-open. `elements` filters on L1/Alu/SVA, `subfamilies` on exact
    subfamily name, `min_copies` drops subfamilies too small to average over --
    the 3D analogue of the loci_in_cat >= 10 floor the eCLIP plots apply.
    """
    names = np.asarray(subfam_names)
    elems = np.asarray(subfam_elem)

    frames = []
    for c in chroms:
        ix = index[c]
        n = len(ix["starts"])
        if n == 0:
            continue
        frames.append(pd.DataFrame({
            "chrom": np.repeat(c, n),
            "start": ix["starts"],
            "end": ix["ends"],
            "strand": ix["strand"],
            "sid": ix["sid"],
        }))
    df = pd.concat(frames, ignore_index=True)
    df["subfam"] = names[df["sid"].to_numpy()]
    df["elem"] = elems[df["sid"].to_numpy()]
    df["mid"] = (df["start"] + df["end"]) // 2

    if elements:
        df = df[df["elem"].isin(elements)]
    if subfamilies:
        df = df[df["subfam"].isin(subfamilies)]
    if min_copies > 0:
        keep = df["subfam"].value_counts()
        keep = keep[keep >= min_copies].index
        df = df[df["subfam"].isin(keep)]

    return df.reset_index(drop=True)[
        ["chrom", "start", "end", "mid", "strand", "sid", "subfam", "elem"]]


def assign_bins(df, clr, pos_col="mid"):
    """Add a genome-wide bin_id for each row, using the cooler's own bin layout.

    Uses clr.offset(chrom) rather than recomputing cumulative sizes so the ids line
    up with clr.bins()[:] exactly. Rows whose position falls past the end of the
    chromosome (possible when a TE annotation and a cooler disagree on a scaffold
    length) get bin_id -1 and should be dropped by the caller.
    """
    binsize = clr.binsize
    if binsize is None:
        raise ValueError("variable-binsize cooler; this pipeline assumes fixed bins")

    sizes = clr.chromsizes
    df = df.reset_index(drop=True)
    pos = df[pos_col].to_numpy()
    bin_id = np.full(len(df), -1, dtype=np.int64)
    for c, idx in df.groupby("chrom", sort=False).indices.items():
        p = pos[idx]
        ok = (p >= 0) & (p < sizes[c])
        bin_id[idx[ok]] = clr.offset(c) + p[ok] // binsize
    out = df.copy()
    out["bin_id"] = bin_id
    return out


def masking_report(df, bins, group_col="subfam"):
    """Per-group count of copies landing in NaN-weight (ICE-masked) bins.

    ICE assigns NaN weights to low-coverage bins, and in unique-alignment data
    repeat-dense bins are disproportionately low-coverage. If subfamily A loses 5%
    of copies to masking and subfamily B loses 20%, then comparing their contact
    averages compares masking rates as much as biology. Report this before
    comparing anything; it is the 3D analogue of the mappability caveat in
    CLAUDE.md (L1H at 0.36% coverage vs L1M at 1.85%).
    """
    w = bins["weight"].to_numpy()
    valid = np.isfinite(w[df["bin_id"].to_numpy()])
    rep = (df.assign(_valid=valid)
             .groupby(group_col)
             .agg(n_copies=("_valid", "size"), n_valid=("_valid", "sum")))
    rep["n_masked"] = rep["n_copies"] - rep["n_valid"]
    rep["frac_masked"] = rep["n_masked"] / rep["n_copies"]
    return rep.reset_index().sort_values("frac_masked", ascending=False)


def make_view(clr, chroms):
    """Viewframe over `chroms` only, for cooltools calls that take view_df."""
    sizes = clr.chromsizes
    return pd.DataFrame({
        "chrom": chroms,
        "start": 0,
        "end": [int(sizes[c]) for c in chroms],
        "name": chroms,
    })


def open_cooler(path, resolution=None):
    """Open a .cool or .mcool. For .mcool a resolution is required.

    Returns (clr, resolution). Kept here so every script reports the same thing
    about what it opened -- resolution silently differing between the compartment
    stage and the contact stage would be very hard to notice downstream.
    """
    import cooler

    if path.endswith(".mcool") or "::" in path:
        if "::" in path:
            clr = cooler.Cooler(path)
        else:
            avail = cooler.fileops.list_coolers(path)
            res_avail = sorted(int(u.rsplit("/", 1)[-1]) for u in avail
                               if u.rsplit("/", 1)[-1].isdigit())
            if resolution is None:
                raise SystemExit(
                    f"{path} is multi-resolution; pass --resolution. "
                    f"available: {res_avail}")
            if resolution not in res_avail:
                raise SystemExit(
                    f"resolution {resolution} not in {path}. available: {res_avail}")
            clr = cooler.Cooler(f"{path}::resolutions/{resolution}")
    else:
        clr = cooler.Cooler(path)

    print(f"  cooler: {os.path.basename(path)} binsize={clr.binsize} "
          f"chroms={len(clr.chromnames)} bins={clr.info['nbins']}")
    return clr, clr.binsize
