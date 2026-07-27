# Chromatin retention of L1 / Alu / SVA in K562 (GSE296212, hg38)

Re-quantification of the GSE296212 chromatin RNA-seq from raw reads, to get **Alu and
SVA** numbers that the deposited files do not contain.

## Why this exists

GSE296212 (SubSeries of GSE296213, "LINE-1 locus transcription nucleates oncogenic
chromatin architecture") deposits one `.txt.gz` count table per sample and nothing else.
That table was checked directly, not assumed: 38,476 rows = ~28,300 gene symbols plus
10,171 individual L1 loci (332 `L1HS_*`, 8,381 `L1PA*`, rest `L1M*`). Rows matching
`^Alu`: **0**. Rows matching `^SVA`: **0**. The authors ran featureCounts against an
L1-only SAF.

None of the 80 samples in the series ships a bigWig, BED or BAM, so there is no coverage
shortcut of the kind that made the `mintseq` analysis cheap. Raw reads are the only route.

The assay itself is a good fit: ribo-depleted (not polyA-selected, so non-polyadenylated
TE transcripts survive), directionally stranded, and already **hg38** — so
`hg38_L1_Alu_SVA_canonical.gtf` applies unchanged, no liftOver anywhere.

## Scope

K562 only, which is the cell line shared with the eCLIP tree and the MINT-Seq m6A work.

| sample | fraction | run | SRA |
|---|---|---|---|
| K562_chrRNA_rep1 | chromatin (NEB) | SRR33381486 | 4.0 GB |
| K562_chrRNA_rep2 | chromatin (NEB) | SRR33381484 | 3.0 GB |
| K562_chr_Takara_rep1 | chromatin (Takara) | SRR33381451 | 2.8 GB |
| K562_chr_Takara_rep2 | chromatin (Takara) | SRR33381450 | 3.4 GB |
| K562_cyto_Takara_rep1 | cytoplasm (Takara) | SRR33381449 | 2.7 GB |
| K562_cyto_Takara_rep2 | cytoplasm (Takara) | SRR33381448 | 1.5 GB |

17.3 GB of SRA, ~100 GB of FASTQ. STAR index needs ~32 GB RAM.

## The contrast, and why it is not raw chrRNA abundance

chrRNA on its own has no denominator. "Alu is highly represented in chromatin RNA" is
close to a restatement of Alu copy number, and the same trap that made a matched TT-Seq
input mandatory for the m6A analysis applies here. The Takara chromatin/cytoplasm pair
restores the denominator:

```
retention = log2( chromatin_takara / cytoplasm_takara )
```

The NEB `chrRNA` libraries ride along as an independent check that the Takara chromatin
fraction agrees with the primary chromatin libraries.

## Interpreter

`deps_chrrna/` is a self-contained pixi env holding the whole toolchain — STAR, sra-tools,
cutadapt, subread, samtools, python/numpy, R and DESeq2. Every step below runs inside it,
so the server needs nothing preinstalled.

```bash
cd deps_chrrna && pixi install     # then use pixi run, or the env's bin/ directly
```

**It is deliberately separate from `deps/`.** `deps/` pins python/numpy/scipy/matplotlib
with `==` because `results/eclip/` must reproduce byte-for-byte, and pulling STAR,
sra-tools and DESeq2 into that solve would force those pins to move. The two envs share
nothing.

`star` is pinned to **exactly 2.7.2b**, the version in the GSE296212 methods. That pin
carries scientific weight rather than being hygiene: with `--outFilterMultimapNmax 10000`
and `--outMultimapperOrder Random`, multi-mapper handling *is* the measurement for repeat
elements, and that behaviour is not guaranteed stable across STAR releases. Everything
else floats — nothing here feeds a committed table.

Resolved on lock: STAR 2.7.2b, subread 2.1.1, sra-tools 3.4.1, cutadapt 5.2,
samtools 1.24, python 3.12.13, numpy 2.5.1, r-base 4.5.3, DESeq2 1.50.2.

## Pipeline

Steps 00–05 run on the server; 06–07 are cheap and run anywhere. `make_saf.py` is
annotation-side and independent of the read pipeline — it is not numbered because it can
run at any point before 05.

```bash
python src/chrrna/make_saf.py                          # -> data/chrrna/te_{sense,antisense}.saf + key
bash   src/chrrna/00_fetch.sh          fastq/
bash   src/chrrna/01_trim.sh           fastq/ fastq_trimmed/
bash   src/chrrna/02_star_index.sh     /path/hg38.fa  star_index_hg38/
bash   src/chrrna/03_align.sh          star_index_hg38/ fastq_trimmed/ bam/
bash   src/chrrna/04_strandedness.sh   bam/            # READ THE OUTPUT
STRAND_DEFAULT=2 bash src/chrrna/05_count.sh bam/ data/chrrna counts/
python src/chrrna/06_enrich_chrrna.py  counts/         # -> results/chrrna/
Rscript src/chrrna/07_deseq2.R         results/chrrna
```

`make_saf.py` has already been run: 2,224,106 copies (Alu 1,228,449 / L1 989,683 /
SVA 5,974)

**Do not skip `01_trim.sh`.** `03_align.sh` uses the authors' `--alignEndsType EndToEnd`,
which forbids soft-clipping, so an untrimmed adapter cannot be clipped — the read just
fails to align, silently. The loss falls hardest on short inserts, which is a large
fraction of a fragmented chromatin RNA library. Skipping it does not error; it quietly
deletes data.

## Decisions that will bite if forgotten

**The counting unit is the subfamily here, not the copy.** Everywhere else in this
project counting is deliberately per TE copy. That is impossible with this data. STAR is
run with the authors' `--outSAMmultNmax 1 --outMultimapperOrder Random`, so a read
matching 500 AluY copies equally well is credited to one of them **at random**. Per-copy
counts are a random partition of a subfamily total, not locus-level expression.
Aggregating over copies undoes the randomisation; reading the per-copy column does not.
The same caveat applies to the authors' own published `L1HS_1..332` per-locus columns.

**Orientation comes from two SAF files, not from flipping `-s`.** `make_saf.py`
writes `te_sense.saf` and `te_antisense.saf` differing only in the Strand column, and
both are counted with the same `-s`. Reasoning about what `-s 1` means for a reversely
stranded paired library is easy to get backwards, and getting it backwards silently
inverts the sense/antisense split the whole project rests on.

**Verify strandedness per library; do not trust the one `-s 2` in the GEO record.** That
note is boilerplate for the whole SuperSeries, but the `chrRNA` samples used the NEB Ultra
II Directional kit and the `Takara` samples a different kit. `04_strandedness.sh` decides
it empirically from four housekeeping genes of known orientation (two `+`, two `-`, so a
genome-wide strand imbalance cannot fake a clean answer). If the two groups disagree, that
is a real kit difference — set `STRAND_<sample>` separately rather than forcing one value.

**No GTF is given to STAR `genomeGenerate`.** Splice junctions from a gene model bias
alignment toward annotated exons, and most TE copies are intronic or intergenic — the
reads being measured are exactly the ones a gene-model-guided aligner down-weights.

**Size factors are computed over TE-assigned reads only**, not genome-wide totals. A
genome-wide factor would import gene-level composition differences between the chromatin
and cytoplasmic fractions into a TE-only comparison.

**Sense and antisense are fitted as separate DESeq2 models.** In the MINT-Seq analysis the
antisense baseline sat ~1.4 log2 units below sense; pooling would let the sense baseline
set the dispersion prior for antisense.

**`06_enrich_chrrna.py` does not produce a p-value, on purpose.** 2v2 replicates over a few
hundred count features is a negative-binomial problem. Quote the DESeq2 `padj` from
`07_deseq2.R`, not the descriptive `log2_retention` — that column is an effect size with a
replicate-concordance flag beside it, nothing more.

## Traps carried over

- `export LC_ALL=C` before any `awk` over these TSVs.
- Mappability still suppresses young subfamilies under unique alignment — but far less
  here than in the eCLIP/MINT-Seq trees, because `--outFilterMultimapNmax 10000` keeps
  multi-mapping reads instead of discarding them. This is the one dataset in the project
  where AluY and L1HS are not structurally under-observed.
