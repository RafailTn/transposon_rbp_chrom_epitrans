# m6A on TE nascent transcripts (MINT-Seq, K562, hg19)

Which L1/Alu/SVA subfamilies carry more m6A on their **nascent** transcripts than the
transcriptome average, per strand orientation. Independent of the eCLIP tree: different
assembly (hg19), different assay, different statistics. Nothing here reads `cache/te_index.pkl`.

## Why hg19 and not liftOver

The MINT-Seq bigWigs (GSE137752) are hg19. Lifting the hg38 TE annotation over was rejected:
liftOver chains are built from *unique* alignments, which is exactly what repeats lack, so the
dropout would be systematic and biased against the young subfamilies that turn out to carry the
signal. Instead the hg19 `rmsk` table is filtered natively by `filter_rmsk_hg19.sh` — same
selection rules as `rmsk_to_canonical_gtf.sh`, no projection step.

Result: 2,144,124 canonical copies (Alu 1,189,771 / L1 950,326 / SVA 4,027), within 0.4% of the
hg38 set, so the two annotations are comparable in aggregate.

## Data

Only the 4 MINT-Seq bigWigs were in `GSE137752_RAW.tar`. The matched **inputs** were downloaded
separately, and are not optional:

| | IP | input |
|---|---|---|
| rep1 | GSM4086534 MINT-Seq | GSM4086533 TT-Seq |
| rep2 | GSM4086536 MINT-Seq | GSM4086535 TT-Seq rep2 |

MINT-Seq is an m6A IP on nascent RNA; TT-Seq is the same nascent RNA without the IP. Raw IP signal
over a TE is (transcription x methylation), and intronic Alus sit in highly transcribed genes, so
IP alone ranks **expression**. Every number here is a per-copy IP/input ratio.

All 8 libraries are pre-scaled to sumData ~5e8, so library normalisation is nearly a no-op; it is
still applied per file because the .pos/.neg pair within a sample differ by ~4%.

## Pipeline

```bash
bash   src/mintseq/filter_rmsk_hg19.sh   # data/hg19_rmsk.gtf -> hg19_L1_Alu_SVA_canonical.bed
python src/mintseq/quantify_m6a.py       # -> results/mintseq/te_copy_signal.tsv.gz   (~8 min)
python src/mintseq/enrich_m6a.py         # -> results/mintseq/te_m6a_enrichment*.tsv
python src/mintseq/plot_m6a.py           # -> results/mintseq/te_m6a_enrichment.{png,pdf}
python src/mintseq/diagnostics_m6a.py    # confounder checks, prints only
python src/mintseq/diagnostics_length.py # age-vs-length deconfounding, prints only
```

## Decisions that were expensive to make

**Strand convention was verified, not assumed.** `.pos`/`.neg` denote the *transcript* strand, not
the read strand — GAPDH(+) is 10x higher in `.pos`, ACTB(-) is 32x higher in `.neg`, in both
replicates. Under the dUTP convention this would have been inverted and the entire sense/antisense
split with it. This is the same class of silent error as the `f[5]` strand bug in `tecommon.read_bed`.

**SVA is classed differently in hg19.** 3733 of 4027 copies are `Other/Other`, only 294 are
`Retroposon/SVA`. The hg38 rule `repClass == "Retroposon"` is exact on hg38 (all 5974 there) but
keeps 7% on hg19. See the header comment in `filter_rmsk_hg19.sh` — the two scripts differ because
the two tables differ, so do not harmonise them.

**Only transcribed copies are tested.** A copy with no nascent RNA has nothing to methylate, and
its ratio is pseudocount noise rather than evidence of no methylation. The floor is 25% of the
median input density over covered copies; ~25% of copies pass. This is the analogue of the
`loci_in_cat >= 10` floor in the eCLIP tree.

**Rank test, not a count model.** bigWig values are CPM-scaled coverage, not raw counts, so
Poisson/binomial dispersion assumptions do not hold. Per-copy log2 ratios are compared with
Mann-Whitney against all other qualifying copies **in the same orientation**, BH-corrected.

**Length is normalised away before the ratio.** Both IP and input are divided by copy length, so
length cancels; L1s at ~6 kb and Alus at ~300 bp are on the same scale.

## The result survived three confounders

`delta_vs_bg` correlates with input density at rho = -0.43 (sense) — the classic noisy-denominator
artefact, and young L1s are low-coverage, which is the wrong direction to be comfortable with. It is
not the explanation:

- **Coverage.** Young L1 vs everything else, *within* input-density deciles: enriched in 10/10
  deciles, median delta +1.69, still p = 7e-13 in the highest-coverage decile.
- **Replicates.** Per-category median log2fc, rep1 vs rep2: Pearson r = +0.88 (sense), +0.80
  (antisense). Not one bad library.
- **Length.** delta correlates with copy length at rho = +0.78 (sense), and young L1s are the
  full-length ones, so age and length are confounded. Splitting them (`diagnostics_length.py`):
  young vs old L1 matched on *both* length band and density quintile is positive in 18/18 cells,
  median +1.48. Independently, within young L1, full-length (>=5.5 kb) beats truncated by +1.6 to
  +2.2 in every quintile — that is the known L1 5'UTR m6A cluster, which only full-length copies
  retain.

## Traps

- **Antisense deltas are relative to a strongly depleted baseline.** The antisense background median
  log2fc is -1.76 (the IP does not retain the mostly-spurious antisense signal in the input), vs
  -0.33 sense. So SVA_B antisense at delta +1.76 has an *absolute* log2fc of ~0.005: it is
  "not depleted", not "actively enriched". Quote `median_log2fc` alongside `delta_vs_bg` for
  anything antisense, or the direction will be overstated.
- Mappability suppresses young subfamilies under unique alignment, same as the eCLIP tree. Here it
  works *against* the reported result rather than for it, so the L1 finding is conservative — but
  absence of a hit on a young subfamily is still not evidence of absence.
- 113,445 copies on alt/random contigs are dropped. hg19 `_hap` contigs are alternate copies of
  chr6/chr4/chr17 sequence, so counting them would double-weight those subfamilies.
- `export LC_ALL=C` before any `awk` over these TSVs, for the same comma-decimal reason as the
  eCLIP tree.
