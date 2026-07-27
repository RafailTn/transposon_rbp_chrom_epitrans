# RBP x transposable-element enrichment (eCLIP)

Tests whether each RBP's **observed eCLIP binding sites** are over-represented on
each L1 / Alu / SVA subfamily, per cell line, strand-aware and counting each TE
copy once (not once per overlapping peak).

Significance is two-stage: a fast hypergeometric **screen** (`enrich_stranded.py`)
proposes hits, and a circular-shift permutation null (`enrich_permutation.py`)
defends them. The permutation q-value is the claim; the hypergeometric one is the
filter that keeps the permutation affordable.

## Pipeline

Run from the project root. Steps 1 and 2 are cached artefacts; rerun them only if
the TE annotation changes.

```bash
bash   src/rmsk_to_canonical_gtf.sh  # data/hg38_rmsk.gtf.gz -> hg38_L1_Alu_SVA_canonical.gtf
python src/build_te_index.py         # -> cache/te_index.pkl

bash src/eclip/run.sh                # -> results/eclip/   (~1 h; --no-prep reuses data/eclip/)
```

Operational detail — install, input provenance, checksums, how to verify a rerun — is in the
top-level `README.md`. This file is the design document: it covers *why*, not *how to run*.

`run.sh` is the five analysis stages in order, preceded by the replicate prep.
Every script defaults to `data/eclip` and `results/eclip`, so each also runs bare:

```bash
python src/eclip/prep_replicates.py            # -> data/eclip/<RBP>_<CELL>.bed
python src/enrich_stranded.py                  # -> rbp_te_enrichment_stranded_*.tsv
python src/enrich_permutation.py               # -> rbp_te_enrichment_permutation*.tsv
python src/report_loci_per_cell.py             # -> distinct_loci_per_cell.tsv
python src/plot_enrichment.py                  # -> top_enrichments.{png,pdf}
python src/plot_across_cell_lines.py           # -> top_hits_across_cell_lines.{png,pdf}
```

All of them still take optional positional dirs (`data_dir`, `results_dir`), so a
scratch run is `python src/enrich_stranded.py data/eclip /tmp/out`.
`enrich_permutation.py` also takes `--eligible BED` (an external expressed/mappable
region set; without it the eligible space is the merged union of all peaks in the
cell line) and `--pmax N` (permutation budget, **default 25000** — see below).

Requires `numpy` and `scipy` (plus `matplotlib` for the plots). Prep and the screen
take about a minute each; the permutation stage dominates at ~54 min for the two
cell lines. Use the project's own env, `deps/.pixi/envs/default/bin/python` — the
system `python`/`python3` and the `rna-seq` pixi env carry neither library. The
versions there are pinned exactly (numpy 2.4.6, scipy 1.17.1, matplotlib 3.10.9)
because the tables in `results/eclip/` are byte-reproducible only against them.

## Inputs

- `data/hg38_rmsk.gtf.gz` — UCSC Table Browser dump of the hg38 `rmsk` table. Not a GTF despite
  the name: the raw 17-column RepeatMasker table, BED 0-based half-open.
  `rmsk_to_canonical_gtf.sh` converts to GTF 1-based inclusive. It replaced `hg38_fixed.gtf`, a
  2022 download whose provenance was lost; see the top-level `README.md` §7 for the diff between
  the two annotations.
- `data/eCLIP_Peaks_K562_HepG2_{K562,HepG2}/<RBP>_<CELL>_{1,2}_chr.bed` — raw
  eCLIP peaks, BED6 `(chrom, start, end, name, signal, strand)`, 0-based
  half-open, hg38, two replicates per RBP. 240 K562 + 205 HepG2 files, 8.36 M
  peaks. Despite the folder names, each folder holds exactly one cell line. Not
  read directly — `prep_replicates.py` turns them into `data/eclip/`.
- `data/eclip/<RBP>_<CELL>.bed` — the prepped form and the pipeline's actual
  input. 116 K562 + 97 HepG2 = 213 files, 1.82 M reproducible peaks. Derived and
  regenerable; safe to delete.

### The eCLIP input

`src/eclip/prep_replicates.py` does two things and nothing else.

**Rename.** Every script parses `(RBP, cell)` with `basename[:-4].rsplit("_", 1)`.
On `AARS_K562_1_chr.bed` that yields RBP `AARS_K562_1` and cell line `chr`, so the
files must become `<RBP>_<CELL>.bed` before anything will run.

**Intersect the replicates.** Reproducible peaks only: the intersected *region* of
rep1 and rep2, per `(chrom, strand)`, with each replicate merged first. Union was
the alternative and was rejected — the background here is pooled peaks over all
RBPs in the cell line, so a peak seen in only one replicate does not merely add
noise to its own row, it inflates the shared denominator every other RBP is tested
against. Taking the region (rather than "rep1 peaks with rep2 support") keeps the
operation symmetric: there is no primary replicate. Merging first matters because
`n_peaks_rbp` is a peak count and is the denominator of every fold — overlapping
duplicates would deflate the fold of whichever RBP has them.

This is strict: 8.36 M raw peaks become 1.82 M reproducible ones (~22%). Ten
(RBP, cell) pairs drop out entirely — nine because one replicate is empty in the
download, plus `SFPQ_HepG2`, which was only ever deposited with one replicate.
Under an intersection rule those drop honestly instead of silently inheriting the
surviving replicate.

Column 5 is written as `.`: after merging, no single input signal value maps to an
output interval. Nothing downstream reads column 5. 4,603 peaks (0.25%) sit on
contigs the TE index does not carry — 4,073 of them on chrM, which has no TEs —
and the prep script reports them rather than letting the index lookup drop them
silently.

## Outputs

In `results/eclip/`: one table per cell line, plus a combined file of significant
hits (q < 0.05, fold > 1). The `*_in_cat` columns are distinct TE copies, not
peaks; `n_peaks_rbp`/`n_peaks_bg` stay peak counts (see the counting note below):

- `rbp_te_enrichment_stranded_<CELL>.tsv` — one row per (RBP, subfamily, orientation)
- `rbp_te_enrichment_stranded_significant.tsv` — the hypergeometric screen's hits
- `rbp_te_enrichment_permutation.tsv` — every screened hit re-tested against the shift null, with `p_perm` / `q_perm` alongside the original `q_hyper`
- `rbp_te_enrichment_permutation_significant.tsv` — those surviving at `q_perm < 0.05`
- `distinct_loci_per_cell.tsv` — distinct TE copies driving the hits, per cell
- `top_enrichments.{png,pdf}` — top 15 robust hits per orientation
- `top_hits_across_cell_lines.{png,pdf}` — replication heatmap

**2,063 screened → 1,794 defended** (461 sense / 1,333 antisense), from 1,037
HepG2 and 1,026 K562 candidates. `report_loci_per_cell.py` and both plot scripts
key off the *permutation-significant* table, so everything downstream shows the
defended hits. The screened table is only an intermediate: the input
`enrich_permutation.py` consumes.

Headline hits, all `loci >= 10`:

- `HNRNPC` antisense Alu — AluSg7 74.0x in K562 (q_perm 9.9e-05), median fold 109.4
  across both cell lines, replicates 2/2. HNRNPC owns 11 of the top 12 antisense rows.
- `EXOSC5` antisense L1PB1 — median fold 50.3, replicates 2/2.
- `KHDRBS1` sense L1MEi — 27.9x in K562, the top sense hit.
- `PRPF4` antisense SVA_D — median fold 16.8, but replicates only 1/2.

## Decisions worth knowing

**Canonical elements only.** Family-level filtering sweeps in the ancestral Alu
monomers (FLAM_A/FLAM_C/FRAM/FAM), the half-L1s (HAL1*), and a mis-filed
`X9_LINE`. `rmsk_to_canonical_gtf.sh` requires the subfamily name to start with
`Alu` or `L1`, dropping 96,213 of 2,320,319 family-level elements and leaving
2,224,106.

**Parsing the attribute field.** `subfamily_id` contains `family_id` as a
substring, so an unanchored regex on `family_id` silently matches the subfamily
slot. Always split on `|` and key off exact field names.

**Pooled-peak background, not a uniform-genome null.** eCLIP peaks sit on bound
RNA and are globally *depleted* in TEs relative to genomic coverage:

| | Alu | L1 | SVA |
|---|---|---|---|
| % of K562 peaks | 1.77% | 4.70% | 0.01% |
| % of HepG2 peaks | 2.79% | 4.04% | 0.02% |
| % of the genome | 10.17% | 17.00% | 0.14% |

A uniform null reports sweeping depletion and buries the per-RBP signal. Fold here
means "relative to where RBPs bind in general", not "relative to random genomic
placement". An RBP that is TE-biased exactly like the average RBP scores fold ~ 1.

**The hypergeometric is a screen, not the claim.** It assumes peaks land
independently at the pooled-peak rate, and they do not: peaks cluster on genes and
TEs cluster in the genome, so its p-values are anticonservative.
`enrich_permutation.py` re-tests every hit it flags against a null that keeps both
structures — circularly shift the RBP's whole peak set by one random toroidal
offset inside the cell's eligible space (same count, widths, strands, and
*relative spacing*), re-run the identical locus-deduped overlap, read p off the
empirical distribution. The shift is rigid on purpose: independent per-peak
placement would destroy exactly the clustering the hypergeometric already ignores.
Eligible space defaults to the merged union of all that cell's RBP peaks — the
permutation analog of the pooled-peak background, i.e. "where this assay sees
RNA-bound signal at all".

**The screen's parameterization mixes units, deliberately.** Written out, the call
in `enrich_stranded.py` is `hypergeom.sf(A-1, Ntot, C, NX)` — i.e. P(X >= A) with

| scipy arg | value | meaning | column |
|---|---|---|---|
| `M` = `Ntot` | total peaks over all RBPs in the cell | population | — |
| `n` = `C[j]` | pooled loci in category j, summed over RBPs | successes in population | `pooled_loci_in_cat` |
| `N` = `NX` | this RBP's peak count | draws | `n_peaks_rbp` |
| `k` = `A[i,j]` | this RBP's distinct loci in category j | observed successes | `loci_in_cat` |

The population and the draw size are counted in **peaks**; the successes are
counted in **distinct TE copies**. Those are not the same unit, so this is not a
genuine 2x2 contingency table over one set of items — a copy is not a peak, and
the "successes" are not a subset of the "population". It is a stylized
approximation, and it is the reason the hypergeometric is labelled a screen
everywhere in this repo rather than a result.

Two consequences worth stating explicitly rather than leaving to be rediscovered:

- **`C` is a sum over RBPs of per-RBP deduped counts, not the number of distinct
  copies in the pooled peak set.** A copy bound by three RBPs contributes 3. That
  is intentional and consistent with `fold` reading as loci-per-peak of this RBP
  against loci-per-peak of everyone else, but it means `pooled_loci_in_cat` is
  *not* a locus census — `report_loci_per_cell.py` exists precisely because
  summing these columns double-counts (see "Distinct copies per cell" below).
- **The support could in principle be violated** — one wide peak can touch two
  copies of the same subfamily and orientation, so nothing in the arithmetic
  forbids `A > NX`, which would put `k` outside the hypergeometric's support and
  return p = 0. Measured on the current tables it never happens: over 77,958 rows
  (13,037 with a nonzero count) not one has `A > NX`; the largest
  `loci_in_cat / n_peaks_rbp` ratio is **0.091** in K562 and **0.102** in HepG2,
  an order of magnitude clear. The failure mode is theoretical here, but it is a
  property of this peak set's width, not a guarantee — re-check it with
  `awk '$7 > $6'` on a peak set with wider intervals.

Why this is tolerable: every hit the screen flags is re-tested by
`enrich_permutation.py` against a null with no such assumption — the shift null
compares an observed locus count to the distribution of the *same statistic*
computed the *same way* on relocated peaks, so units cancel and the mixing above
never enters the claim. The screen only has to be a sensible ordering that
over-selects; it is not required to be calibrated. What it must not do is *miss*
real hits, and its anticonservatism is the safe direction for that. Do not quote
`p_fisher_greater` / `q_BH` as evidence — quote `q_perm`, as the Outputs section
says.

(The column is named `p_fisher_greater` because the one-sided Fisher exact test
and the hypergeometric upper-tail probability are the same number; the name
records the interpretation, `scipy.stats.hypergeom.sf` is the implementation.)

**The 269 casualties are concentrated, and it is informative where.** Of the 269
screened hits that fail the shift null, **195 are *sense* L1** and 65 more are
sense Alu (144 HepG2, 125 K562). Antisense loses almost nothing: 9 of 1,342. A hit that
dissolves under a rigid shift was never positional — it was the RBP's clustering
coinciding with where the element sits, and sense-strand L1 is where that
coincidence is most common (L1 is long, so an intronic peak cluster lands on one
by chance far more easily than on a 300 bp Alu). The antisense signal, which is
the biologically interesting one here, is essentially untouched by the null.

**`fold_perm` is SMALLER than `fold_pooled`.** Median ratio 0.64 overall (HepG2
0.75, K562 0.55) — that is, `exp_loci_perm` sits *above* the hypergeometric's
expectation, so the fold measured against it comes out lower. Two properties of
eCLIP plausibly drive that, both pushing the null's locus yield up: peaks are wide
(median 48 bp), so a relocated peak touches several copies; and the segment
overhang is large (59% of placements), leaking the null into TE-denser gap
sequence. Treat that as the likely mechanism, not a verified one — and do not
assume the direction is a property of the method. It is a property of this peak
set's width and clustering, and a set with different geometry can invert it.

The consequence is unchanged and is the point: **quote `fold_perm` next to
`q_perm`, and `fold_pooled` next to `q_hyper`.** The two do not agree, the gap is
not a constant rescaling, and ranking by one is not ranking by the other.

**Budget scales with candidate count.** `P_MAX` must be large enough that the
p-value floor `1/(P_MAX+1)` still clears BH against the largest candidate set in a
cell. eCLIP has 1,037 candidates in HepG2, needing `P_MAX > 20,740` — hence the
**default 25000**. Getting this
wrong does not error: at 2000, 56% of hits tied on the p-floor and were mutually
unorderable, and BH cleared only because the ties rescued it. At 25000 that is
43%, across 623 distinct `p_perm` values, and the rule is satisfied outright
(`1/25001 * 1037 = 0.042`). The budget barely moves the verdict (1,789 → 1,792
survivors) but it sharpens `fold_perm`, whose denominator is a Monte Carlo mean —
HNRNPC/AluSg7/K562 moved 79.4 → 74.0, a better estimate of the same result.
(That 2000-vs-25000 comparison was measured on the retired `hg38_fixed.gtf`
annotation and has not been repeated since the switch to the UCSC `rmsk` dump;
the 25000 run on the current annotation gives 1,794. The conclusion — the budget
sharpens `fold_perm` far more than it moves the verdict — is unaffected.)

Two knobs keep it affordable: only screened hits are permuted, and sequential
stopping (Besag & Clifford 1991) abandons a category once it has accumulated 25
exceedances. The run is seeded per (cell, RBP) and reproduces exactly; a rerun that
disagrees means the script changed, not that the null is noisy.

**The shift null has a known, deliberate, conservative bias.** Shifted peaks
overhang their eligible segment 59% of the time, into gap sequence ~1.4x TE-denser
than the footprint, which inflates the null's TE hit rate. So `fold_perm` is
understated and `p_perm` overstated: survivors are safe, near-misses may be false
negatives. Do **not** patch it by requiring full containment — that would confine
placement to long, heavily co-bound, TE-poor segments and trade a bounded
conservative bias for an unbounded anti-conservative one. The full measurement and
the reasoning are in `enrich_permutation.py`'s docstring.

**Strand matters.** Each overlap is split into sense and antisense (peak strand vs
element strand); BH runs across both orientations within a cell line. The eCLIP
result is overwhelmingly antisense — 1,333 of 1,794 defended hits — and it is not
one protein: **40 distinct RBPs** carry antisense hits, led by MATR3 (173),
HNRNPM (170), SUGP2 and EXOSC5 (117), PTBP1 (115), KHSRP (100), QKI (94) and
HNRNPC (71). Those are polypyrimidine-tract and intronic binders, and antisense Alu
is poly-U rich; HNRNPC's antisense Alu signal matches Zarnack et al. 2013 and is the
biological positive control. Collapsing strands would average this against sense
noise.

**Counting is per TE copy, not per peak.** A category count could be "overlapping
peaks" or "distinct TE copies touched". Counting peaks lets *one* heavily-tiled
copy carry a whole "subfamily enrichment", so the pipeline dedupes each peak to the
copy it hits (`pid` -> locus before the `bincount`): a copy counts once per RBP per
orientation, the pooled background is deduped the same way, and
`n_peaks_rbp`/`n_peaks_bg` stay peak totals so `fold` reads as loci-per-peak of the
RBP vs background.

**On eCLIP the dedup is a safeguard, not the load-bearing filter**: reproducible
eCLIP peaks are already close to non-redundant per copy (HNRNPC averages 1.04
peaks per copy), so it changes little here. Keep it anyway — it costs nothing, the
three overlap routines must agree, and on a peak set that tiles a bound region
more heavily it would be decisive rather than cosmetic.

**L1 is not a minor contributor here.** L1 supplies 13,383 of K562's 19,889
distinct hit loci — outnumbering Alu's 6,492 — and 12,118 of HepG2's 27,560. Do
not assume Alu dominates just because it is the higher-copy element. Fraction of
copies covered by any pooled K562 peak:

| | covered |
|---|---|
| L1P | 1.35% |
| L1M | 1.85% |
| AluS | 0.70% |

**Mappability biases against young subfamilies.** eCLIP peak calling uses unique
alignments, so high-identity young copies are systematically harder to see. L1H
sits at 0.36% coverage against L1M's 1.85%, and there are zero L1H hits at all.
Because the shift null relocates peaks *inside the same footprint*,
observed and null are equally blind to unmappable copies, so the **test stays fair**
— the bias limits what is testable, not what is valid. But absence of AluY/L1HS
hits is **not** evidence of absence.

**Interpret extreme folds with care.** Filter on `loci_in_cat >= 10` (`obs_loci` in
the permutation table) for anything load-bearing; the plots apply exactly that
floor. Deduping to copies lowers counts, so power drops, and a consistent
sub-threshold fold across both cell lines is still signal.

**Ranking by fold is not ranking by count.** Fold ranks a handful of copies above
broad subfamilies, so `plot_enrichment.py` prints per-element totals in each panel
title and `plot_across_cell_lines.py` facets by element — with eCLIP all three
blocks render (Alu 7 rows, L1 7, SVA 1, out of 990 combos robust in >= 1 cell
line, from 52,704 tested).
Don't read the top-N bars as a census. Ranked by `fold_perm` the top 15 are
K562-heavy in both orientations (12 K562 / 3 HepG2), so read a cross-cell "top N by
fold" table as a within-cell statement.

**Distinct copies per cell.** `report_loci_per_cell.py` unions the real copy
identities across a cell's permutation-significant hits (not `sum(obs_loci)`, which
would double-count a copy bound by two RBPs):

| cell | RBPs | hits | distinct loci | Alu | L1 | SVA |
|---|---|---|---|---|---|---|
| HepG2 | 36 | 893 | 27,560 | 15,392 | 12,118 | 50 |
| K562 | 42 | 901 | 19,889 | 6,492 | 13,383 | 14 |

**`fold_enrichment` is NaN or inf in some rows.** Of 77,958 stranded rows, NaN in
7,221 (the (subfamily, orientation) category has zero pooled copies — nobody binds
it — so the fold is 0/0) and inf in 32 (only that one RBP binds it, so the
background rate is zero). Both must be special-cased: every `lo <= nan < hi` test is
False, so a naive binning loop falls through to the last bin and paints "no data" as
"maximally enriched". `plot_across_cell_lines.py` renders NaN as a neutral `n/a`
cell.

## Correctness checks

`enrich_stranded.py` recounts a random chr21 peak subset with a brute-force O(n*m)
overlap (deduped to distinct copies) and asserts it matches the vectorised
`searchsorted` path. The index builder asserts its element count equals the GTF
record count. `enrich_permutation.py` checks the flat-to-genomic coordinate
round-trip (peaks mapped into concatenated-eligible space and back must land where
they started) and recomputes each hit's observed locus count through the permutation
code path, warning to stderr if it disagrees with the screened table — that catch
fires if the two overlap routines ever drift apart. `prep_replicates.py`
brute-forces the replicate intersection on a random chr21 slice against its
vectorised path, the same way.

The same `searchsorted` overlap idiom is written three times (`peaks_per_cat`,
`count_A`, `loci_hit`), each keeping something different (counts / batched counts /
locus identities). They must agree; change one, change all three.
`prep_replicates.intersect` is a fourth use of the *window* trick, but it computes
interval overlaps between two peak sets rather than peak-to-TE hits, and is not
bound by that agreement.

Three traps that cost real time here, in case you edit these:

- `read_bed` splits with `line.split("\t", 6)` and tests `f[5]` for strand. In a
  BED6 file strand is the last field, so `f[5]` is `"+\n"` and a bare `== "+"`
  silently labels **every** peak minus — no error, no empty file, just an inverted
  sense/antisense split, which is the one thing this analysis is about. Hence the
  `.rstrip()`.
- `sub` is an awk builtin and cannot be used as a variable name. It fails with a
  syntax error and, in a pipeline, yields an empty output file.
- Sanity-checking these TSVs with `awk` under a comma-decimal locale parses
  `2.107e-03` as `2` and `1.651` as `1`, so numeric filters silently reject
  everything. Always `export LC_ALL=C` first.
