# iMARGI stage: does SVA/Alu antisense RNA contact DNA in trans?

Tests the hypothesis that Alu/SVA RNA participates in 3D genome organisation, using
RNA–DNA contact maps rather than the RNA-only assays in the other stages.

## Why this assay and not the others

`src/chrrna` established that Alu and SVA antisense RNA is chromatin-retained, and
`src/mintseq` that SVA antisense is the only antisense repeat class where m6A is not
depleted. Both are **necessary but very weak** evidence: pervasive nascent transcription
is chromatin-retained by definition and looks identical in those assays. Neither can show
the defining property of an architectural RNA — acting in *trans*, at loci other than
where it was transcribed. iMARGI measures exactly that.

`src/hic` is a different question again: it asks where the TE **DNA copies** sit relative
to compartments and promoters. DNA position is not RNA function.

## Data

4DN `4DNESIKCVASO`, "iMARGI on K562 cells", 2 replicates, GRCh38, released.

| file | what | size |
|---|---|---|
| `4DNFIVIHUHOE` | **pairs, contact list-combined** | 7.35 GB |
| `4DNFIOHY9ZX7` | mcool | 2.04 GB |
| 4 × `4DNFI…` | pairs, per replicate | 2.1–2.4 GB |
| 2 × `4DNFI…` | bam | 68–73 GB |

Take the **combined pairs**. The mcool bins both ends, and at any resolution a cooler
offers an SVA is sub-bin — binning discards which copy the RNA came from, which is the
entire measurement. `fetch_4dn.py`'s docstring recommends the mcool; that advice is
correct for `src/hic` and wrong here.

Downloaded md5 verified against the portal: `3f53f1c4779159ad10b2439f6fe4ebd8`.

## Pipeline provenance, read off the pairs header

```
bwa mem -SP5M GRCh38_no_alt_analysis_set (+chrEBV decoy)
pairtools parse --no-flip --min-mapq 1 --report-alignment-end 5 --walks-policy 5any
pairtools dedup ; imargi_restrict.py
pairtools select regex_match(pair_type,"[UuR][UuR]") and abs(dist2_rsite)<=3 and <self-ligation filter>
```

**`--min-mapq 1`, not a unique filter.** Only MAPQ 0 — genuinely tied alignments — is
dropped; rows with MAPQ 4–9 pass. Repeat-derived reads with a distinguishable best hit
survive, which is why SVA coverage is far better than a strict-unique pipeline would give.
The residual bias still runs against near-identical young copies, so absence at SVA_E/F is
not evidence of absence.

**Do not realign from FASTQ to "recover" multimappers.** The readout is a distance between
where the RNA was transcribed and where the DNA end landed, so it needs a unique RNA
origin. Random multimapper assignment (as STAR does in `src/chrrna`, where it is correct
for counting) would randomise that origin and destroy the measurement. Unique alignment is
a requirement here, not a limitation.

## The two conventions, both measured not assumed

`check_pairs.py` settles both in one streaming pass over all 248,414,406 pairs (~5 min).
4DN's pipeline notes say DNA/RNA column order was swapped when building the cool files, so
the `.pairs` and `.mcool` disagree and neither can be taken on faith.

**Side 1 is the RNA end.** Three of four discriminators voted, unanimously:

| discriminator | side1 | side2 |
|---|---|---|
| exon fraction | 26.31% | 11.30% |
| gene-body fraction | 93.53% | 77.50% |
| top-0.1% bin share | 8.95% | 1.09% |
| chrM share | 0% | 0% — **abstains** |

The file contains **zero chrM contacts**, so that discriminator has no signal. It abstains
rather than voting; a 0-vs-0 tie falling through to an `else` branch cast a spurious vote
for side 2 in the first version and halted the run.

**`strand1` is INVERTED relative to the transcript strand** — 4/4 housekeeping genes, by
factors of 300–800×:

| gene | strand | RNA '+' | RNA '-' |
|---|---|---|---|
| GAPDH | + | 128 | 39,191 |
| RPL13A | + | 123 | 100,980 |
| ACTB | − | 54,101 | 132 |
| EEF1A1 | − | 25,544 | 542 |

**A read whose strand matches a TE's annotated strand therefore came from the ANTISENSE
transcript.** Assuming the natural convention inverts every sense/antisense claim — and
this whole stage is a sense/antisense question. `check_pairs.py` accumulates only the raw
`same`/`opp` observable and applies the mapping at report time, so the labels can never
drift from the measurement.

## Coverage: the test is powered

SVA copies carrying at least one RNA-end contact, from the census:

| orientation | class | contacts | copies |
|---|---|---|---|
| antisense | intronic_anti | 116,459 | 1,468 |
| sense | intronic_sense | 23,130 | 722 |
| sense | intergenic | 9,416 | 972 |
| **antisense** | **intergenic** | **8,055** | **959** |

959 of the 1,624 intergenic SVA copies (59%) carry antisense contacts. That is hundreds,
not the tens that a strict-unique pipeline would have left — the `--min-mapq 1` finding
explains the difference.

## Result: the hypothesis is not supported for SVA

`extract_rna_ends.py` then `trans_test.py --shifts 1000` (intergenic copies, 329
categories). Statistic = fraction of a subfamily's contacts reaching cis-distal (>1 Mb)
or trans, against a 1000-shift toroidal null and an abundance-matched pre-mRNA baseline.

**Every SVA category sits at or below its position-matched null.** z ranges from −0.08
(SVA_E antisense) to −5.65 (SVA_A sense); all q = 1. Nothing is distal-enriched.

The two references disagree in an instructive way, which is why both exist:

| family | orient | median (obs − pre-mRNA) | median (obs − shift null) |
|---|---|---|---|
| SVA | antisense | **+0.123** | −0.087 |
| SVA | sense | −0.160 | −0.330 |
| Alu | antisense | −0.179 | −0.193 |
| Alu | sense | −0.207 | −0.196 |
| L1 | antisense | −0.369 | −0.373 |
| L1 | sense | −0.427 | −0.465 |

SVA antisense is the least locally-tethered TE category and is the **only** one that beats
abundance-matched pre-mRNA. Against the shift null it does not. The difference between
those two statements is the whole point of the null: "SVA RNA acts at a distance" and "SVA
copies sit in places where RNA is generally distal" predict the same pre-mRNA comparison
and are separated only by position matching. The pre-mRNA comparison alone would have been
read as a positive.

The broader finding is a negative that points the other way: **TE-derived RNA is
systematically *more* locally tethered than its position predicts**, most extremely L1
sense (−0.465) — the same category `src/chrrna` identified as dominated by host pre-mRNA
passengers. Consistent story across two unrelated assays.

**One hit survives BH: AluYa5 sense** (z = 5.13, q = 4.7e-05; distal 0.896 vs null 0.714
vs pre-mRNA 0.642; 3,397 contacts over 1,201 copies). Treat it cautiously. AluYa5 is the
youngest and most self-similar Alu subfamily, so it is the biggest magnet for reads
mismapped from other Alu copies genome-wide — and such reads carry DNA ends scattered
across the genome, which is exactly a distal-looking signature. The shift null does not
control for this, because it moves anchors to positions that attract no such mismapping.
AluYe5 sense is second (z = 2.26, n.s. after BH), which is the same youth ordering.

### What this test cannot see

The statistic is a genome-wide **aggregate fraction** per subfamily. An architectural RNA
acting on a few hundred specific loci — the Xist-like case — would barely move it. This
result rules out *broad* distal action by SVA/Alu RNA; it does not rule out focused
targeting. The natural follow-up is to look for recurrent specific DNA targets of SVA
antisense RNA rather than an aggregate distal fraction.

The 1 Mb cis-distal threshold is a convention, not a measurement.

## The references, and why both

### An abundance effect large enough to swamp everything

The pre-mRNA baseline's distal fraction rises monotonically with contact count, from
**0.165 in the lowest decile to 0.845 in the highest** — a 5× range. Highly transcribed
genes are the distal-heavy ones, not the noise-dominated low-abundance ones. Any
unstratified comparison of distal fractions here measures expression level and nothing
else. Every comparison above is decile-matched.

## Not yet the test (superseded -- kept for the reasoning)

The DNA-end distance split in `pairs_checks.txt` is descriptive only. **A high trans
fraction is also what noise looks like** — ambient and random ligations are
trans-dominated, so the low-abundance antisense channel having more trans contacts may
only mean a worse signal-to-noise ratio. Nothing can be concluded without:

1. a **pre-mRNA baseline** — protein-coding nascent RNA, which is strongly cis-tethered
   and defines what "locally acting" looks like in this library, and
2. the **rigid toroidal shift null** used in `enrich_permutation.py` and
   `te_tss_contacts.py`, preserving copy spacing.

Restrict to `intergenic` copies from `08_partition_copies.py` before any of it; copies
inside genes inherit the host's contact profile.

One internal consistency check already passes: **L1 sense is the most cis-proximal
category (40.1%)** of any measured, which is what host pre-mRNA passengers predict — and
`src/chrrna` independently found 60% of L1-sense reads come from same-strand intronic
copies. Two different assays, same artefact, same direction.
