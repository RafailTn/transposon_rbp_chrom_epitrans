# Chromatin retention of L1 / Alu / SVA RNA in K562 (GSE296212, hg38)

Re-quantification of the GSE296212 chromatin and cytoplasmic RNA-seq from raw reads, to
ask whether TE RNA stays on chromatin or is exported — for **Alu and SVA** as well as L1.

Full write-up of the results: [Chromatin retention of L1, Alu and SVA RNA in K562](https://claude.ai/artifact/HjaW12Gc5ymympWN3CbniU).

## Why this exists

GSE296212 (SubSeries of GSE296213, "LINE-1 locus transcription nucleates oncogenic
chromatin architecture") deposits one `.txt.gz` count table per sample and nothing else:
~28,300 gene symbols plus 10,171 L1 loci, **0 Alu and 0 SVA rows**. The tables are
htseq-count output (read level, multi-mappers dropped), and the L1 rows are full-length
copies (they match our L1 copies >= 5 kb per subfamily). No sample ships a bigWig, BED
or BAM, so raw reads are the only route.

The assay fits: ribo-depleted (non-polyadenylated TE transcripts survive), stranded,
paired-end, already hg38.

## Samples

| sample | fraction | kit | run | role |
|---|---|---|---|---|
| K562_chr_Takara_rep1/2 | chromatin | Takara (SMARTer Pico v2-type) | SRR33381451 / 450 | contrast numerator |
| K562_cyto_Takara_rep1/2 | cytoplasm | Takara | SRR33381449 / 448 | contrast denominator |
| K562_chrRNA_rep1/2 | chromatin | NEB Ultra II Directional | SRR33381486 / 484 | independent chromatin check |

Sample metadata lives in `config/runs.tsv`. 17.3 GB of SRA, ~100 GB of FASTQ; the STAR index
needs ~32 GB RAM. The cytoplasm libraries are ~35% rRNA.

```
retention = log2( chromatin_takara / cytoplasm_takara )
```

chrRNA alone has no denominator ("Alu is abundant in chromatin RNA" mostly restates Alu
copy number); the Takara cytoplasm pair supplies it. NEB has no cytoplasm partner, so it
cannot replicate a retention result.

## Environment

`deps_chrrna/` is a self-contained pixi project; nothing else needs to be installed.

```bash
cd deps_chrrna && pixi install && pixi install -e em
export PATH=/mnt/raid1/rafail/chrrna/deps_chrrna/.pixi/envs/default/bin:$PATH LC_ALL=C THREADS=32
```

- `default`: STAR **2.7.2b** (pinned exactly — the GSE296212 version; multi-mapper
  handling is the measurement for repeats), sra-tools, cutadapt, subread, samtools,
  python 3.12/numpy, R/DESeq2, fgsea, msigdbr.
- `em`: TEtranscripts (`TEcount`), which needs python 3.10 and cannot share the default solve.

It is deliberately separate from the sibling project's `deps/`
(`/mnt/raid1/rafail/aluaca/deps/`), whose `==` pins must not move.

## Layout

```
chrrna/
├── config/runs.tsv     sample sheet: every step iterates over it
├── src/                the pipeline, 00–19 + make_saf.py
├── deps_chrrna/        pixi environment (pixi.toml, pixi.lock)
├── data/               references: STAR index, TE and GENCODE-derived annotations, GEO tables
├── work/               intermediates: fastq/, fastq_trimmed/, bam/, bam_multi/, counts*/, strand_check/
├── results/            final tables; results/em/ is the EM arm
├── logs/               run logs and step stdout
└── archive/v1/         the first, buggy run, kept for comparison
```

`data/`, `work/`, `archive/` and `logs/` are git-ignored. Every script resolves its
default paths from the project root, so steps run from any directory with no arguments;
arguments override the defaults.

## Pipeline

Every step skips a sample whose output exists — delete the output to redo it.

```bash
# primary arm (authors' alignment settings)
python src/make_saf.py <hg38_L1_Alu_SVA_canonical.gtf>     # done: data/te_{sense,antisense}.saf, te_copy_key.tsv
bash   src/00_fetch.sh                                     # -> work/fastq
bash   src/01_trim.sh                                      # -> work/fastq_trimmed (TruSeq adapters + 3 nt off Takara R2)
bash   src/02_star_index.sh <hg38.fa>                      # -> data/ (STAR index)
bash   src/03_align.sh                                     # -> work/bam
bash   src/04_strandedness.sh                              # -> work/strand_check; read the verdicts (-s 2 for all six)
bash   src/05_count.sh                                     # -> work/counts
python src/06_enrich_chrrna.py                             # -> results/subfamily_counts.tsv, chromatin_retention.tsv
Rscript src/07_deseq2.R                                    # TE-only normalisation: a ranking, not retention calls

# genomic context + gene-based normalisation (the retention calls)
python src/08_context.py                                   # -> data/te_context.tsv, data/gencode_intron.saf
bash   src/09_gene_count.sh                                # -> work/gene_counts
python src/10_context.py
Rscript src/11_context_deseq2.R                            # -> results/chromatin_retention_context*.tsv

# EM arm (TEtranscripts over all alignments)
bash   src/12_align_multi.sh                               # -> work/bam_multi
bash   src/13_em_count.sh                                  # -> work/counts_em
python src/14_em_matrix.py                                 # -> results/em/
Rscript src/07_deseq2.R results/em

# cross-checks: authors' GEO tables, random vs EM arm, NEB vs Takara
python src/15_crosscheck.py

# Alu boundary test (own transcript vs read-through)
python src/16_boundary.py                                  # -> results/boundary_test.tsv

# host-gene function test for intronic Alu retention
python src/17_gene_alu.py                                  # -> results/gene_alu_counts.tsv
Rscript src/18_gene_alu_gsea.R                             # per-gene score, raw + regression-adjusted GSEA
Rscript src/19_gene_alu_matched.R                          # expression-matched GSEA + matched controls
FAMILY=L1 python src/17_gene_alu.py                        # same for intronic L1, + per-lineage counts
Rscript src/18_gene_alu_gsea.R results L1                  # -> results/gene_l1_*
Rscript src/19_gene_alu_matched.R results L1               # + per-lineage matched GSEA
```

References outside the project: GENCODE v47 and GRCh38 in `/mnt/raid1/rafail/aluaca/ref/`.

## Results in brief

| Question | Answer | Output |
|---|---|---|
| Does the pipeline reproduce the authors? | Genes r = 0.97; full-length L1 abundance ρ 0.95–1.00; L1 retention ρ 0.93, 37/37 same direction | `results/crosscheck_summary.tsv` |
| L1 retention | Chromatin-retained in every context (median +5.5); same-strand intronic L1 +3.8 above host introns; L1HS not special | `results/chromatin_retention_context*.tsv` |
| Alu retention | At or below host introns inside genes; ~0 in intergenic, opposite-strand and 3′UTR copies | same |
| Alu subfamilies | Sense-only age gradient: AluJ chromatin-bound, AluY cytoplasm-leaning, also away from genes | same |
| Is cytoplasmic young-Alu RNA its own transcript? | Ends at the first downstream T-run (Pol III-like) but starts upstream; Pol III vs Pol II undecided | `results/boundary_test.tsv` |
| Host-gene functions | Intronic Alu RNA ~1.6x more retained than its introns overall; after expression matching only OXPHOS (+) and TGF-β (−) hold | `results/gene_alu_matched_sets.tsv` |
| L1 host-gene functions | Intronic L1 RNA ~3.9x more retained than its introns; after matching only TGF-β (−0.91) and TNF-α/NF-κB (−0.88) hold, both lower; no lineage (L1PA…L1ME) has its own pathway theme | `results/gene_l1_matched_sets.tsv` |

## Decisions that will bite if forgotten

**Retention calls come from the gene normalisation.** Cytoplasm carries ~20x fewer TE
fragments per mapped fragment than chromatin. TE-only size factors (`06`, `07`,
`norm == te` in `11`) erase that and put host introns on the cytoplasm side (−0.8),
which is impossible; gene median-of-ratios factors put them at +2.5. No spike-ins, so
"cytoplasmic" means relative to mRNA.

**Context before family.** Sense Alus in 3′UTRs leave with mRNA; intronic TEs stay with
pre-mRNA. Pooled retention largely restates exon vs intron. `11` fits each context
separately with FIXED size factors (re-estimating per subset would make each subset
compositional) and reports `delta_vs_intron` for intronic TEs.

**The counting unit is the subfamily, not the copy.** STAR keeps one alignment per
multi-mapper, chosen at random, so per-copy counts are a random partition of the
subfamily total. Copies are used only as a support filter, a context label, the
full-length L1 subset, through their unique flanks (`16`), or with uniquely mapped
reads only (`17`).

**Count fragments.** featureCounts >= 2.0.2 needs `--countReadPairs` beside `-p`, or
each mate is counted. The first run had this bug.

**Takara R2 carries 3 template-switch bases** (`GGG` in ~50% of reads). Under the
authors' `--alignEndsType EndToEnd` they cannot be soft-clipped, so `01_trim.sh` clips
them (`CLIP_R2_<sample>` to override). Adapter trimming is mandatory for the same reason.

**Orientation comes from two SAF/GTF files** (TE strand as annotated vs flipped), counted
with the same `-s` / `--stranded`. Never flip `-s` to get antisense.

**Strandedness is checked per library** (`04`), from four housekeeping genes (two per
strand); the GEO `-s 2` note is SuperSeries boilerplate.

**No GTF in STAR `genomeGenerate`**: annotated junctions bias alignment away from the
intronic and intergenic reads being measured.

**The intron baseline masks rDNA and the exons of all gene types.** Unmasked, rRNA in
lncRNA "introns" was ~75% of cytoplasmic intronic fragments.

**Sense and antisense are separate DESeq2 models**; their baselines differ.

**`06` produces no p-values.** Quote DESeq2 `padj`.

**Two quantification arms.** Subfamilies whose log2FC moves by >1 between random
placement and EM (`results/crosscheck_arms.tsv`) are quantification-sensitive.

## Open issues

- NEB chromatin is 60% antisense vs 36% for Takara chromatin, unchanged by the fixes — a
  kit property, unexplained.
- Clipping Takara R2 raised unique mapping but not the "too short" rate.
- Cytoplasmic young-Alu RNA: Pol III with upstream starts, or Pol II ending at an Alu
  polyadenylation signal? Next test: read 3′-end positions vs TTTT vs AATAAA, and
  poly(A)-tailed reads among the unmapped.
- `export LC_ALL=C` before any `awk` over these TSVs.
