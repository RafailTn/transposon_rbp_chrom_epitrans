# RBP × transposable-element enrichment — procedures and reproduction

Tests whether each RBP's **observed eCLIP binding sites** are over-represented on each L1 / Alu / SVA
subfamily, per cell line, strand-aware, counting each TE copy once rather than once per overlapping
peak. Two cell lines (K562, HepG2), 213 (RBP, cell) pairs, 1.82 M reproducible peaks against 2.22 M
TE copies.

This file is the **operational** document: what to install, what inputs are needed, which commands to
run in which order, how long they take, and how to tell whether a rerun reproduced.

**The reasoning is not here.** `src/README.md` is the design document — why the background is
pooled-peak rather than uniform-genome, why counting is per copy, why the permutation budget scales
with candidate count, why the shift null's conservative bias is deliberate. Read it before changing
anything statistical. 

A separate m6A / MinT-seq analysis lives in `src/mintseq/` with its own README. It is not covered
here and shares only the RepeatMasker filtering idea.

---

## 1. Repository layout

```
hg38_L1_Alu_SVA_canonical.gtf      stage-1 product: filtered TE annotation (2,224,106 elements)
cache/te_index.pkl                 stage-2 product: interval index consumed by every analysis script
data/
  hg38_rmsk.gtf.gz                 INPUT  UCSC hg38 rmsk table dump (see §3)
  eCLIP_Peaks_K562_HepG2_K562/     INPUT  raw eCLIP peaks, 240 files
  eCLIP_Peaks_K562_HepG2_HepG2/    INPUT  raw eCLIP peaks, 205 files
  eclip/                           derived: 213 prepped <RBP>_<CELL>.bed files
deps/                              self-contained pixi environment
results/eclip/                     all output tables and figures
backup_hg38_fixed_provenance/      the superseded annotation, index, and results (see §7)
src/                               scripts; src/README.md is the design document
```

---

## 2. Environment

`deps/` is a self-contained [pixi](https://pixi.sh) environment. Nothing outside the project
directory is required.

```bash
cd deps && pixi install && cd ..
```

Then use the interpreter directly — every command below assumes it:

```
deps/.pixi/envs/default/bin/python     # python 3.12.13, numpy 2.4.6, scipy 1.17.1, matplotlib 3.10.9
deps/.pixi/envs/default/bin/Rscript    # R 4.5.3, Bioconductor SummarizedExperiment 1.40.0
```

or equivalently `cd deps && pixi run python ...`.

Two things that will waste your time if ignored:

- **The system `python`/`python3` and the `rna-seq` pixi env carry neither numpy nor scipy.** Use the
  path above.
- **The python versions are pinned exactly** (`==`, not `>=`) in `deps/pixi.toml`. scipy's
  `hypergeom.sf` and numpy's PCG64 stream both feed the committed tables; loosen the pins and
  bit-level reproducibility is gone.

Set `export LC_ALL=C` in any shell where you will `awk` over the output TSVs. Under a comma-decimal
locale awk parses `2.107e-03` as `2` and `1.651` as `1`, so numeric filters silently reject
everything. `run.sh` and `rmsk_to_canonical_gtf.sh` set it internally.

---

## 3. Inputs

### TE annotation — `data/hg38_rmsk.gtf.gz`

A UCSC Table Browser dump of the hg38 `rmsk` table: genome **hg38**, group **Repeats**, track
**RepeatMasker**, table **rmsk**, output format **all fields from selected table**.

**It is not a GTF**, despite the filename the Table Browser hands out. It is the raw 17-column
RepeatMasker table with a `#`-prefixed header:

```
bin  swScore  milliDiv  milliDel  milliIns  genoName  genoStart  genoEnd  genoLeft  strand
repName  repClass  repFamily  repStart  repEnd  repLeft  id
```

and its coordinates are **BED 0-based half-open**, not GTF 1-based inclusive. The converter handles
both facts; see §4.

Equivalently, and without the header line:

```bash
curl -o data/hg38_rmsk.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
```

`rmsk_to_canonical_gtf.sh` accepts either form — its header skip is
`NR > 1 || $1 !~ /^#/`, which passes a headerless file through untouched. Gzipped or plain both work.

### eCLIP peaks — `data/eCLIP_Peaks_K562_HepG2_{K562,HepG2}/`

BED6 `(chrom, start, end, name, signal, strand)`, hg38, 0-based half-open, two replicates per RBP,
8.36 M peaks total. Supplied as the two zips in `data/`; unzip in place.

**Downloaded from the PrismNet web server**, <http://prismnetweb.zhanglab.net/download/>, as the
entries *"RBP binding sites in K562 from eCLIP"* and *"RBP binding sites in HepG2 from eCLIP"* — one
zip per cell line, which is why the two folders in `data/` share a stem and differ only in the
trailing cell name. Cite PrismNet if you use these files; see §9.

The peaks are experimental eCLIP calls redistributed by PrismNet, **not** PrismNet model
predictions — the same server also offers predicted binding sites, and those are a different
download that this analysis does not use. The underlying experiments are ENCODE's eCLIP set, so
ENCODE is cited in §9 as well; note that the download page itself carries no provenance or assembly
statement, so that attribution comes from the PrismNet papers rather than from the files.

Two traps in the naming:

- **Despite the folder names, each folder holds exactly one cell line.**
- The filenames are `AARS_K562_1_chr.bed`. Every analysis script parses `(RBP, cell)` with
  `basename[:-4].rsplit("_", 1)`, which on that name yields RBP `AARS_K562_1` and cell line `chr`.
  **This is why the prep step must rename** — the raw folders can never be read directly.

### Checksums of the inputs this tree was last built from

```
c49fc30435fbc705acc26f2e9344528b  data/hg38_rmsk.gtf.gz
979797974ea0b9377b780a065908ef4e  data/eCLIP_Peaks_K562_HepG2_K562.zip
1749b54d3ba3270fe6511b0f20b17e70  data/eCLIP_Peaks_K562_HepG2_HepG2.zip
```

UCSC updates `rmsk` in place across releases, so a fresh download may not match that first hash. It
will still run; expect small count differences and rerun the §6 comparison rather than assuming a
bug.

---

## 4. Reproducing the analysis

Run from the project root, in order. Stages 1–2 are cached artefacts — rerun them **only if the TE
annotation changes**; stage 3 is the analysis.

```bash
bash src/rmsk_to_canonical_gtf.sh    # data/hg38_rmsk.gtf.gz -> hg38_L1_Alu_SVA_canonical.gtf   ~11 s
deps/.pixi/envs/default/bin/python src/build_te_index.py         # -> cache/te_index.pkl         ~4 s
bash src/eclip/run.sh                # -> results/eclip/                                         ~1 h
```

Expected console output, in order:

| stage | expects |
|---|---|
| `rmsk_to_canonical_gtf.sh` | `wrote …/hg38_L1_Alu_SVA_canonical.gtf (2224106 elements)` |
| `build_te_index.py` | `chroms=631 elements=2224106 subfamilies=183` then `element count matches GTF: PASS` |
| `run.sh` | prep summary, then `HepG2: 937,552 peaks … 1,037 enriched`, `K562: 885,395 peaks … 1,026 enriched`, `TOTAL significant: 2,063`, and finally `1794 / 2063 hypergeometric hits survive the permutation null at q < 0.05  (461 sense / 1333 antisense)` |

The full console record of the last run is kept at `results/eclip/_perm_run.log`.

### Stage 1 — `src/rmsk_to_canonical_gtf.sh [in] [out]`

Filters the rmsk table to canonical L1 / Alu / SVA and emits the 9-column GTF the rest of the
pipeline reads. Adds 1 to `genoStart` for GTF 1-based inclusive; `build_te_index.py` converts
straight back to BED. The round trip is deliberate — the GTF stays the one interchange format
downstream expects.

Selection rules, all three required:

```
SINE       + Alu + repName starting "Alu"     # excludes ancestral monomers FLAM_A/FLAM_C/FRAM/FAM
LINE       + L1  + repName starting "L1"      # excludes half-L1s HAL1* and a mis-filed X9_LINE
Retroposon + SVA
```

**"Canonical" means canonical *subfamilies*, not canonical chromosomes.** There is deliberately no
chromosome filter: eCLIP peaks land on 10 non-canonical scaffolds, and TE annotation there is real
signal. TE loci on scaffolds carrying no peaks are inert anyway — every count in the analysis is
driven off the peak files, and `peaks_per_cat` does `INDEX.get(c)`, skipping unindexed chromosomes.

### Stage 2 — `src/build_te_index.py [gtf] [pkl]`

Builds per-chromosome sorted `starts/ends/sid/strand/maxlen` arrays plus `subfam_names` /
`subfam_elem`. `maxlen` bounds the left edge of the candidate window in every overlap query.
Category encoding is `cat = subfamily_id * 2 + (0 sense | 1 antisense)`, used as an index into
length-`NCAT` arrays throughout.

### Stage 3 — `src/eclip/run.sh [--no-prep] [--pmax N]`

The replicate prep plus five analysis stages. Every script defaults to `data/eclip` and
`results/eclip`, so each also runs bare:

```bash
python src/eclip/prep_replicates.py    # -> data/eclip/<RBP>_<CELL>.bed          ~1 min
python src/enrich_stranded.py          # -> rbp_te_enrichment_stranded_*.tsv      ~1 min
python src/enrich_permutation.py       # -> rbp_te_enrichment_permutation*.tsv   ~54 min
python src/report_loci_per_cell.py     # -> distinct_loci_per_cell.tsv
python src/plot_enrichment.py          # -> top_enrichments.{png,pdf}
python src/plot_across_cell_lines.py   # -> top_hits_across_cell_lines.{png,pdf}
```

All take optional positional dirs (`data_dir`, `results_dir`), so a scratch run is
`python src/enrich_stranded.py data/eclip /tmp/out`. The plot scripts and `report_loci_per_cell.py`
take `results_dir` only (the report takes both).

Order constraints:

- `enrich_permutation.py` reads `rbp_te_enrichment_stranded_significant.tsv`, so
  **`enrich_stranded.py` must run first.**
- The report and both plot scripts read `rbp_te_enrichment_permutation_significant.tsv` —
  everything downstream is filtered on `q_perm`, and the screened table is only the permutation's
  input. Note the column is named `obs_loci` there, not `loci_in_cat`.

Useful flags:

- `--no-prep` skips `prep_replicates.py`. Safe whenever `data/eclip/` already holds its 213 files —
  the prep depends only on the raw peaks, never on the TE annotation, so **a TE-annotation rerun
  should always use it.**
- `--pmax N` sets the permutation budget (**default 25000**). `enrich_permutation.py` also takes
  `--eligible BED`, an external expressed/mappable region set; without it the eligible space is the
  merged union of all peaks in that cell line.

**Do not lower `--pmax` casually.** The p-floor `1/(P_MAX+1)` must still clear BH against the largest
candidate set in a cell — HepG2 has 1,037, so `P_MAX > 20,740`. Setting it too low does not error: it
silently ties hits on the floor so they cannot be ranked.

---

## 5. Outputs

All in `results/eclip/`. The `*_in_cat` columns are distinct TE copies, not peaks; `n_peaks_rbp` /
`n_peaks_bg` stay peak counts.

| file | contents |
|---|---|
| `rbp_te_enrichment_stranded_<CELL>.tsv` | one row per (RBP, subfamily, orientation) |
| `rbp_te_enrichment_stranded_significant.tsv` | the hypergeometric screen's hits |
| `rbp_te_enrichment_permutation.tsv` | every screened hit re-tested against the shift null |
| `rbp_te_enrichment_permutation_significant.tsv` | those surviving at `q_perm < 0.05` — **the claim** |
| `distinct_loci_per_cell.tsv` | distinct TE copies driving the hits, per cell |
| `top_enrichments.{png,pdf}` | top 15 robust hits per orientation |
| `top_hits_across_cell_lines.{png,pdf}` | replication heatmap |
| `_perm_run.log` | permutation stage console record |

**Quote `q_perm`, not `q_hyper`, for anything load-bearing** — and quote `fold_perm` next to `q_perm`,
`fold_pooled` next to `q_hyper`. The two folds do not agree, the gap is not a constant rescaling, and
ranking by one is not ranking by the other.

Filter on `loci_in_cat >= 10` (`obs_loci` in the permutation table) for anything load-bearing;
extreme folds on a handful of copies are noise. The plots apply exactly that floor.

---

## 6. Verifying a reproduction

### Built-in assertions

These are the only tests in the project. They run automatically and must all pass; keep them.

| script | check |
|---|---|
| `build_te_index.py` | indexed element count == GTF record count |
| `enrich_stranded.py` | brute-forces a random chr21 peak subset, asserts it matches the vectorised path |
| `enrich_permutation.py` | flat↔genomic coordinate round-trip; warns if a recomputed observed count disagrees with the screened table |
| `prep_replicates.py` | brute-forces the replicate intersection on a random chr21 slice |

### Determinism

The permutation is seeded per (cell, RBP) and reproduces exactly. **Two runs of the same code on the
same annotation give identical tables.** A rerun that disagrees means the script or the inputs
changed — not that the null is noisy.

### Diffing two result trees

Compare on the natural key rather than by line, since row order can shift:

```bash
export LC_ALL=C
# key = cell_line, rbp, element, subfamily, orientation
cut -f1-5 results/eclip/rbp_te_enrichment_permutation_significant.tsv | sort > /tmp/new.key
cut -f1-5 OTHER_TREE/rbp_te_enrichment_permutation_significant.tsv    | sort > /tmp/old.key
comm -23 /tmp/old.key /tmp/new.key | wc -l   # lost
comm -13 /tmp/old.key /tmp/new.key | wc -l   # gained
```

---

## 7. Annotation provenance — read before comparing to older results

Stage 1 used to be `src/filter_te_gtf.sh`, reading `hg38_fixed.gtf` — a 2022 download whose
provenance was lost, so the annotation could not be regenerated from a named source. It was replaced
by `src/rmsk_to_canonical_gtf.sh` reading a re-downloadable UCSC dump. `filter_te_gtf.sh` and
`hg38_fixed.gtf` are kept as a record; **nothing runs them.** The superseded annotation, index, and
result tables are archived in `backup_hg38_fixed_provenance/`.

The two annotations were diffed on `(chrom, start, end, strand, subfamily)` before the switch:

- Restricted to chr1–22/X/Y/M they agree on **2,059,498** elements. The new dump is a **strict subset**
  there — 94 rows only in the old file, **none** only in the new.
- Those 94 are 11–35 bp `+`-strand Alu micro-fragments across 92 `gene_id`s: gap-filling
  sub-fragments the old GTF conversion inserted between an inserted element and the resumption of its
  host, which the raw table leaves unannotated. Every one lies inside the span of siblings present in
  both files.
- Genome-wide the new annotation is **larger**: 631 sequences vs 375, 2,224,106 elements vs
  2,153,399, 183 subfamilies vs 180 (gaining AluYb10, AluYb11, AluYb8a1; losing none).

**The rebuild is not purely cosmetic.** eCLIP peaks land on 10 non-canonical scaffolds and the old
GTF carried TE annotation for only 4 of them, so peaks on the other 6 counted toward the denominator
but could never hit a TE. That is a real, small correctness fix. Post-switch tables are therefore
**not** expected to match `backup_hg38_fixed_provenance/results_eclip/`, and a difference there is
not evidence of a bug.

### Measured effect of the switch

The rerun is complete. Both trees, compared on `(cell_line, rbp, element, subfamily, orientation)`:

| | pre-switch | post-switch |
|---|---|---|
| screened (hypergeometric) | 2,063 | 2,063 |
| defended (`q_perm < 0.05`) | 1,792 | **1,794** (461 sense / 1,333 antisense) |
| HepG2 distinct loci | 27,538 | 27,560 |
| K562 distinct loci | 19,866 | 19,889 |

**1,791 of the 1,792 pre-switch hits are retained.** One is lost (`K562 EIF4G2 Alu AluY antisense`)
and three are gained (`HepG2/K562 HNRNPL L1 L1MCa antisense`, `K562 PTBP1 L1 L1MB3 antisense`) — all
four sit at the q < 0.05 boundary. Of the 1,791 shared hits, `obs_loci` is **identical in 1,761
(98.3%)**; the 30 that move all move *up*, by 1 to 5 loci, as expected from a strictly larger
annotation. The largest is `HepG2 SUGP2 L1MB3 antisense`, 72 → 77 (`fold_perm` 4.08 → 4.30).

Every headline result is unchanged: `HNRNPC` antisense AluSg7 74.00× in K562, `EXOSC5` antisense
L1PB1 median fold 50.3, `KHDRBS1` sense L1MEi 27.92×, `PRPF4` antisense SVA_D 16.8.

---

## 8. Traps

- **`tecommon.read_bed` does `line.split("\t", 6)` and tests `f[5]` for strand.** On BED6 strand is
  the last field, so `f[5]` is `"+\n"` and a bare `== "+"` silently marks **every** peak minus,
  inverting the sense/antisense split the whole analysis rests on. The `.rstrip()` is load-bearing.
- **The same `searchsorted` overlap idiom is written three times** — `enrich_stranded.peaks_per_cat`,
  `enrich_permutation.count_A`, `report_loci_per_cell.loci_hit`. They must agree; change one, change
  all three. The dedupe of `(locus, category)` before the `bincount` is what makes counts per-copy.
- **`fold_enrichment` is NaN or inf in many rows** — NaN when the category has no pooled peaks, inf
  when only that RBP binds it. Every `lo <= nan < hi` test is False, so a naive binning loop paints
  "no data" as "maximally enriched". Special-case both, as `plot_across_cell_lines.bin_of` does.
- **`export LC_ALL=C` before any `awk` over the TSVs.** See §2.
- **`sub` is an awk builtin and cannot be a variable name.** It fails with a syntax error and, in a
  pipeline, yields an empty output file.
- **The GTF attribute field packs everything into one pipe-delimited string, and `subfamily_id`
  contains `family_id` as a substring.** Split on `|` and key off exact field names; never use an
  unanchored regex.
- **The shift null has a known conservative bias** — shifted peaks overhang their eligible segment
  59% of the time into gap sequence ~1.4× TE-denser than the footprint. Survivors are safe;
  near-misses may be false negatives. Do **not** "fix" it by requiring full containment; see
  `enrich_permutation.py`.
- **Mappability suppresses young subfamilies.** eCLIP peak calling uses unique alignments, so
  high-identity young copies are harder to see: L1H sits at 0.36% coverage against L1M's 1.85%, with
  zero L1H hits. The test stays fair because the null lives in the same footprint, but absence of
  AluY/L1HS hits is **not** evidence of absence.
- **L1 is not a minor contributor.** It supplies 13,383 of K562's hit loci against Alu's 6,492. Don't
  assume Alu dominates just because it is the higher-copy element.
- **`fold_perm` runs smaller than `fold_pooled` here** (median 0.64). That sign is a property of this
  peak set's width and clustering, not of the method — don't memorise it as a rule.

---

## 9. Citations

### eCLIP peaks — PrismNet

The eCLIP peak files in `data/eCLIP_Peaks_K562_HepG2_{K562,HepG2}/` were obtained from the PrismNet
web server (§3). Cite the web server paper for the download, and the original method paper for
PrismNet itself:

> Xu, Y., Zhu, J., Huang, W., Xu, K., Yang, R., Zhang, Q.C., & Sun, L. (2023). PrismNet: predicting
> protein–RNA interaction using in vivo RNA structural information. *Nucleic Acids Research*,
> **51**(W1), W468–W477. <https://doi.org/10.1093/nar/gkad353>

> Sun, L., Xu, K., Huang, W., Yang, Y.T., Li, P., Tang, L., Xiong, T., & Zhang, Q.C. (2021).
> Predicting dynamic cellular protein–RNA interactions by deep learning using in vivo RNA structures.
> *Cell Research*, **31**(5), 495–516. <https://doi.org/10.1038/s41422-021-00476-y>

The eCLIP experiments themselves were generated by the ENCODE consortium:

> Van Nostrand, E.L., Freese, P., Pratt, G.A., et al. (2020). A large-scale binding and functional map
> of human RNA-binding proteins. *Nature*, **583**, 711–719.
> <https://doi.org/10.1038/s41586-020-2077-3>

### TE annotation — RepeatMasker via UCSC

`data/hg38_rmsk.gtf.gz` is the hg38 `rmsk` table (§3), i.e. RepeatMasker output distributed by the
UCSC Genome Browser:

> Smit, A.F.A., Hubley, R., & Green, P. RepeatMasker Open-4.0. <https://www.repeatmasker.org>

> Navarro Gonzalez, J., Zweig, A.S., Speir, M.L., et al. (2021). The UCSC Genome Browser database:
> 2021 update. *Nucleic Acids Research*, **49**(D1), D1046–D1057. <https://doi.org/10.1093/nar/gkaa1070>
