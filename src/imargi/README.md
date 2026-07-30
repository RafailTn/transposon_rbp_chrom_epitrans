# iMARGI stage: does SVA/Alu antisense RNA contact DNA in trans?

Tests the hypothesis that Alu/SVA RNA participates in 3D genome organisation, using
RNA–DNA contact maps rather than the RNA-only assays in the other stages.

## Provenance rule: this stage rests on H1 alone

**Every number, null and conclusion below is measured on H1-hESC and nothing else.** An
earlier K562 iMARGI dataset was analysed and then deleted, and none of it is carried
forward as evidence — not as corroboration, not as a cross-check, not as a
"replicates in a second line" argument.

The reason is specific and not merely tidiness: **K562's copy-number variation was never
resolved.** The BCR-ABL amplicon was the visible part and was handled by naming
coordinates, but K562 is broadly aneuploid, the rest of its CNV landscape was never
characterised here, and the effect of that on an assay whose readout is *genomic distance
between two ligated ends* is unknown. Copy-number inflation and derivative-chromosome
topology both act directly on that readout. An unquantified bias cannot be used to
support a conclusion, so it is not used at all.

Where another stage is cited below it is flagged, because **`src/chrrna`, `src/mintseq`
and the eCLIP stage are all K562/HepG2**. Those citations are design rationale or
context, never evidence for the H1 result.

## Why this assay and not the others

`src/chrrna` established that Alu and SVA antisense RNA is chromatin-retained, and
`src/mintseq` that SVA antisense is the only antisense repeat class where m6A is not
depleted. Both are **necessary but very weak** evidence: pervasive nascent transcription
is chromatin-retained by definition and looks identical in those assays. Neither can show
the defining property of an architectural RNA — acting in *trans*, at loci other than
where it was transcribed. iMARGI measures exactly that.

(Those stages are K562. They motivate the question; they cannot confirm the answer, both
because they measure something else and because they are a different cell line.)

`src/hic` is a different question again: it asks where the TE **DNA copies** sit relative
to compartments and promoters. DNA position is not RNA function.

## Data

4DN `4DNESNOJ7HY7`, "iMARGI on H1 cells" — H1-hESC (Tier 1), no-treatment control,
Sheng Zhong lab, 4 bioreps, GRCh38, released.

| file | what | size |
|---|---|---|
| `4DNFIGDJIRV3` | **pairs, contact list-combined** | 8.84 GB |
| `4DNFIKWK6NSD` | mcool | 2.34 GB |

Take the **combined pairs**. The mcool bins both ends, and at any resolution a cooler
offers an SVA is sub-bin — binning discards which copy the RNA came from, which is the
entire measurement. `fetch_4dn.py`'s docstring recommends the mcool; that advice is
correct for `src/hic` and wrong here. (`fetch_4dn.py` also needs the `deps-3d` env, not
`deps` — it imports `requests`, which `deps` does not carry.)

Downloaded md5 verified against the portal: `f5d1f95a763065592987fc123981bd01`.
**285,576,525 pairs.**

**H1, not a cancer line, is the deliberate choice.** The per-copy screen below is a
per-locus test, and on a rearranged karyotype copies inside an amplified region are
copy-number inflated *and* physically cis on the derivative chromosome while the
reference scores their contacts as trans — which manufactures precisely the recurrence
that screen looks for. H1 has no such region (measured, not assumed — see the density
profile below).

## Pipeline provenance, read off the pairs header

```
bwa mem -t 4 -SP5M GRCh38_no_alt_analysis_set_GCA_000001405.15 (+chrEBV decoy)
pairtools parse --add-columns mapq,cigar --no-flip --min-mapq 1
               --max-inter-align-gap 20 --report-alignment-end 5 --walks-policy 5any
pairtools sort ; pairtools dedup --mark-dups ; imargi_restrict.py
pairtools select regex_match(pair_type,"[UuR][UuR]") and abs(dist2_rsite)<=3
               and <self-ligation filter> --chrom-subset <hg38 chrom.sizes>
```

`--chrom-subset` is why **100%** of pairs land on chr1–22/X/Y: there are no scaffold or
chrM rows to discard, so `extract_rna_ends.py` reporting "100.0% stored" is correct and
not a filter that failed to fire.

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

`check_pairs.py` settles both in one streaming pass over all 285,576,525 pairs (~7 min).
4DN's pipeline notes say DNA/RNA column order was swapped when building the cool files, so
the `.pairs` and `.mcool` disagree and neither can be taken on faith.

**These are properties of a FILE, not of the assay** — a different deposit may have been
built by a different pipeline — so `check_pairs.py` must be re-run on any new pairs file
rather than trusting the report on disk. Both scripts re-derive the strand mapping and
assert it, so a silent drift is not possible.

**Side 1 is the RNA end.** Three of four discriminators voted, unanimously:

| discriminator | side1 | side2 |
|---|---|---|
| exon fraction | 26.16% | 10.43% |
| gene-body fraction | 91.02% | 76.58% |
| top-0.1% bin share | 6.73% | 0.42% |
| chrM share | 0% | 0% — **abstains** |

The file contains **zero chrM contacts**, so that discriminator has no signal. It abstains
rather than voting; a 0-vs-0 tie falling through to an `else` branch cast a spurious vote
for side 2 in the first version and halted the run.

**`strand1` is INVERTED relative to the transcript strand** — 4/4 housekeeping genes, by
factors of 170–820×:

| gene | strand | RNA '+' | RNA '-' |
|---|---|---|---|
| GAPDH | + | 363 | 62,082 |
| RPL13A | + | 227 | 185,625 |
| ACTB | − | 109,392 | 175 |
| EEF1A1 | − | 36,361 | 513 |

**A read whose strand matches a TE's annotated strand therefore came from the ANTISENSE
transcript.** Assuming the natural convention inverts every sense/antisense claim — and
this whole stage is a sense/antisense question. `check_pairs.py` accumulates only the raw
`same`/`opp` observable and applies the mapping at report time, so the labels can never
drift from the measurement.

## Coverage: the test is powered

50,697,140 of 285,576,525 contacts (17.75%) have their RNA end inside an indexed TE copy.
SVA copies carrying at least one RNA-end contact, from the census:

| orientation | class | contacts | copies |
|---|---|---|---|
| antisense | intronic_anti | 73,160 | 1,553 |
| sense | intronic_sense | 18,673 | 800 |
| sense | intergenic | 10,279 | 1,038 |
| **antisense** | **intergenic** | **8,051** | **1,047** |

1,047 of the 1,624 intergenic SVA copies (64%) carry antisense contacts. That is hundreds,
not the tens that a strict-unique pipeline would have left — the `--min-mapq 1` finding
explains the difference.

## Result: the hypothesis is not supported for SVA

`extract_rna_ends.py` then `trans_test.py --shifts 1000` (intergenic copies, 337
categories). Statistic = fraction of a subfamily's contacts reaching cis-distal (>1 Mb)
or trans, against a 1000-shift toroidal null and an abundance-matched pre-mRNA baseline.

**Every SVA category sits at or below its position-matched null.** z ranges from −6.39
(SVA_D sense) to +1.04 (SVA_E sense); all q = 1. Nothing is distal-enriched. **No category
of any family survives BH** — the maximum z over all 337 is +2.99.

The two references behave differently, which is why both exist:

| family | orient | median (obs − pre-mRNA) | median (obs − shift null) |
|---|---|---|---|
| SVA | antisense | **−0.023** | −0.155 |
| SVA | sense | −0.275 | −0.366 |
| Alu | antisense | −0.313 | −0.279 |
| Alu | sense | −0.288 | −0.268 |
| L1 | antisense | −0.401 | −0.343 |
| L1 | sense | −0.472 | −0.418 |

SVA antisense is the least locally-tethered TE category and the only one anywhere near
parity with abundance-matched pre-mRNA (−0.023). Against the shift null it is clearly
below. That gap is the whole point of the null: "SVA RNA acts at a distance" and "SVA
copies sit in places where RNA is generally distal" predict the same pre-mRNA comparison
and are separated only by position matching. **Read on the pre-mRNA baseline alone, SVA
antisense would look like a near-positive; it is not one.**

The broader finding points the other way: **TE-derived RNA is systematically *more*
locally tethered than its position predicts**, most extremely L1 sense (−0.418). That
ordering is internally consistent within H1 — L1 sense is also the most cis-proximal
category measured here (45.1%), which is what host pre-mRNA passengers predict. (`src/chrrna`
reached a compatible conclusion about L1 sense, but it is K562 and a different assay, so it
is noted rather than counted.)

**The top category is AluYa5 sense** (z = 2.99, q = 0.47 — nominally highest, not
significant; distal 0.784 vs null 0.666 vs pre-mRNA 0.601; 2,378 contacts over 1,201
copies). Do not chase it. AluYa5 is the youngest and most self-similar Alu subfamily, so
it is the biggest magnet for reads mismapped from other Alu copies genome-wide — and such
reads carry DNA ends scattered across the genome, which is exactly a distal-looking
signature. The shift null does not control for this, because it moves anchors to positions
that attract no such mismapping. AluYe5 sense follows the same youth ordering. The MAPQ
section below settles that this is mismapping and not signal.

### What this test cannot see

The statistic is a genome-wide **aggregate fraction** per subfamily. An architectural RNA
acting on a few hundred specific loci would barely move it. This rules out *broad* distal
action; it does not rule out focused targeting.

The 1 Mb cis-distal threshold is a convention, not a measurement.

## The cis-restricted rerun, which is the better-aimed test

`distal_frac` above pools cis-distal **and** trans. That answers "is this RNA non-local at
all", but it is not the test for the architectural hypothesis people actually mean. **Xist,
Airn and Kcnq1ot1 all act in cis** — they spread along the chromosome they were transcribed
from and do essentially nothing in trans. Trans-acting architectural RNA (Firre) is the rare
exception, not the model.

That is a numerical problem, not only a conceptual one. Over TE-anchored contacts this
library is **56.9% trans, 9.6% cis-distal, 33.5% cis-proximal**, so the pooled statistic
runs about six parts trans to one part cis-distal: a cis-spreading signal can move a long
way without shifting it, and the channel that dominates it is where ambient ligation lives.

`--statistic cisdistal` conditions on the RNA having stayed on its own chromosome —
`cis-distal / (cis-proximal + cis-distal)` — dropping trans from the denominator entirely.
The dclass triple was already carried through `obs`, `null` and the baseline; only the last
line collapsed it. Tallies now dump to `tallies_<class>.npz`, so a further reduction never
costs another hour-long run.

```bash
python src/imargi/trans_test.py --statistic cisdistal --shifts 1000
# -> results/imargi/trans_test_intergenic_cisdistal.tsv
```

### Result: not one category exceeds its null

**The maximum z over all 331 categories is +0.03.** Not a single category exceeds its
position-matched null by any meaningful margin, let alone reaches significance; every
q = 1. Family medians stay negative against both references, so the "TE RNA is more locally
tethered than its position predicts" finding holds on the cis channel too:

| family | orient | median (obs − pre-mRNA) | median (obs − shift null) |
|---|---|---|---|
| SVA | antisense | **−0.091** | −0.108 |
| SVA | sense | −0.159 | −0.167 |
| Alu | antisense | −0.256 | −0.150 |
| Alu | sense | −0.262 | −0.155 |
| L1 | antisense | −0.314 | −0.154 |
| L1 | sense | −0.344 | −0.175 |

Two things change materially against the pooled statistic, and both remove a near-positive:

**SVA antisense vs pre-mRNA goes from −0.023 to −0.091.** Its near-parity with
abundance-matched pre-mRNA lived in the trans channel. On the cis channel — the one the
spreading hypothesis actually predicts — it is clearly below nascent pre-mRNA, not at
parity with it.

**AluYa5 sense collapses from z = 2.99 to z = −0.76 (q = 1).** The pooled test's top
category does not survive cis-restriction. This is exactly what the mismapping explanation
predicts: reads mismapped from Alu copies genome-wide carry DNA ends scattered across
*chromosomes*, inflating trans and not cis-distal. The cis test is the one to quote.

Every SVA category sits at or below its null (z −4.32 to +0.03), with the shift null saying
RNA from comparable positions reaches ~21% cis-distal while SVA sense manages far less.

## The AluYa5 hit was mismapping — confirmed by MAPQ, not by replicates

The pairs carry **`mapq1` and `mapq2`** (columns 9–10; `pairtools parse --add-columns
mapq,cigar`). `extract_rna_ends.py` now stores both, and `trans_test.py --min-mapq N`
filters on the RNA end's own MAPQ before the cumsum.

**Replicate splits cannot test this and it was a mistake to suggest they could.**
Mismapping is systematic, not stochastic: both bioreps sequence the same repeat families,
align to the same reference with the same aligner, and mismap identically. An artifact
reproduces across replicates perfectly. Replicates test *sampling* stability, a different
question.

MAPQ does test it, because the mechanism is specific. A mismapped RNA end puts the
**anchor** in the wrong place while its DNA end stays correctly mapped near the RNA's true
origin — so measured against the wrong anchor it looks far away. Mismapping manufactures a
distal signature, and `mapq1` is what discriminates.

Rerunning the *same* pooled `distal` statistic that produced the top ranking, one variable
changed:

| category | statistic | n | frac | null | z | q |
|---|---|---|---|---|---|---|
| AluYa5 sense | distal, MAPQ≥1 | 2,378 | 0.784 | 0.666 | **+2.99** | 0.47 |
| AluYa5 sense | distal, MAPQ≥30 | 546 | **0.291** | 0.657 | **−9.53** | 1 |
| AluYe5 sense | distal, MAPQ≥1 | 1,728 | 0.628 | 0.660 | −0.55 | 1 |
| AluYe5 sense | distal, MAPQ≥30 | 708 | 0.360 | 0.653 | −5.30 | 1 |

The signal does not merely weaken — it **inverts**, from well above the null to far below it.
High-confidence AluYa5 RNA is strongly locally tethered like every other TE category; the
0.784 was carried entirely by ambiguous alignments.

Three checks make this conclusive rather than a power failure:

**The depth cost tracks subfamily age**, which is the mismapping prediction and nothing
else. Pooled over all categories 86.2% of RNA ends survive MAPQ 30. AluYa5 sense keeps 23%:

| subfamily | kept at MAPQ≥30 |
|---|---|
| AluYa5 sense (youngest) | 23.0% |
| SVA_D antisense | 26.4% |
| L1PA2 sense | 26.8% |
| L1PA2 antisense | 29.2% |
| L1HS antisense | 39.2% |
| L1M3d antisense (old) | 99.5% |
| L1M2a1 sense (old) | 99.7% |
| L1M2a1 antisense (old) | 100.0% |

**The null is unmoved**: median shift −0.0090 across shared categories (range −0.0113 to
−0.0063). A 49-point collapse in the observation against a 1-point shift in the reference is
not a filtering artifact.

**The surviving 546 contacts still estimate the fraction precisely** — SE ≈ 0.019, putting
the observed 0.291 about **19 SE below** the null. Depth collapsed, but the remaining
estimate is sharp and points the other way. That is a confirmed artifact, not an
inconclusive one.

Two further H1-internal results point the same way, and neither shares the depth-loss
weakness:

- **Cis-restriction alone kills it**, at full MAPQ≥1 depth: AluYa5 sense goes from z = 2.99
  pooled to **z = −0.76** on `cisdistal` with no MAPQ filter at all. Mismapped reads scatter
  DNA ends across *chromosomes*, so they inflate trans and not cis-distal. No contacts were
  discarded to get this.
- **The within-intron test reproduces the collapse where depth does not collapse** — Alu
  categories there keep 69–83% of their contacts and their `z_gene` still evaporates
  (AluSc5 sense +6.38 → +0.02 on 83% retention). See that section.

So the case rests on four independent H1 observations: the inversion, the age-ordered depth
cost, the unmoved null, and two full-depth reproductions that spend no reads. **This is why
no second cell line is invoked.** Mismapping is a property of the *reference* and of a
subfamily's self-similarity, so it is demonstrable within one library — and the alternative
corroboration is unavailable anyway, since the only other dataset here had unresolved CNV
(see the provenance rule). Do not re-litigate the AluYa5 hit on H1 evidence; it is settled.

**At MAPQ≥30 the maximum z over all 335 categories is 1.84 and nothing has q < 0.05.** The
pooled test agrees with the cis-restricted test: no TE subfamily in either orientation
shows distal enrichment over its position-matched null, on either statistic.

### What this costs

Cis-restriction roughly halves an already thin channel — SVA_D antisense keeps 1,199 cis
contacts of 2,652 total. Aggregate and feature-level tests are fine at that depth; a
per-source genome-wide bin scan is not (see the power table in the follow-up notes).

## Recovering the intronic copies: the within-intron control

`trans_test.py` runs on intergenic copies only, discarding 1.22 M intronic copies — 2× the
intergenic set, and in the SVA census the intronic channel carried ~107 k contacts against
intergenic's ~18 k. `within_intron_test.py` recovers them by changing the **control** rather
than the filter.

**The control is the same gene.** Each protein-coding gene's intron space splits into
TE-covered and TE-free bases. Both sit in the same transcript, compartment, expression
level and replication-timing domain, so host expression — which drives distal fraction
0.195 → 0.817 across deciles, by far the largest confound — cancels by construction instead
of needing to be matched.

This is why the genome-wide toroidal null **cannot** simply be pointed at intronic copies:
it moves anchors mostly outside the host gene, and gene bodies at high expression sit at
0.817 against a genome-wide null of ~0.66. Pooling classes against that null would
manufacture large positive z from gene-body membership alone.

**Only the host-opposite strand is counted.** For an intronic copy the channel running with
the host strand is pre-mRNA passing through, co-directionally, and swamping whatever the copy
itself produces. (`src/chrrna` quantified that asymmetry at 49× in K562; the argument here is
mechanical and does not depend on that figure, which is why no H1 equivalent is needed to
justify the design.) Since strand1 is inverted,
host-opposite means `strand1 == gene strand`, asserted on housekeeping genes. Which TE
orientation that is follows from the copy's class, so the label is derived, never chosen:

    intronic_sense copy (TE strand == gene strand)  ->  host-opposite is ANTISENSE
    intronic_anti  copy (TE strand != gene strand)  ->  host-opposite is SENSE

**The control subtracts every indexed TE copy**, not just the tested subfamily. Introns are
TE-dense; otherwise the "TE-free" control would largely be other L1/Alu/SVA and the
comparison would be TE-vs-TE.

Two within-gene references, and they do different jobs:

1. **`z_gene`** — the gene's own TE-free rate. Controls the host completely, but not *where
   in the gene* the copies sit.
2. **`z_shift`** — the gene's TE footprint shifted rigidly inside its own concatenated
   intron space, 200 offsets. Adds control for position within the gene (distance to TSS,
   intron rank). Same construction as `enrich_permutation.py`, eligible set narrowed from
   the chromosome to one gene's introns.

```bash
python src/imargi/within_intron_test.py --shifts 200 --min-ctl 100 --min-mapq 30
```

### Result: nothing survives either reference

| | MAPQ≥1 | MAPQ≥30 |
|---|---|---|
| rows | 314 | 306 |
| categories with q_gene < 0.05 | **10** | **0** |
| max z_gene | 16.37 | 3.53 (q = 0.063) |
| categories with q_shift < 0.05 | **0** | **0** |
| **passing BOTH** | **0** | **0** |

At MAPQ≥30 nothing is significant on *either* reference, not even nominally, and the
position control (`z_shift`) already rejected everything at MAPQ≥1 — the ten `q_gene`
survivors there are all explained by *where in the introns* the copies sit.

**This is the strongest of the mismapping demonstrations, because depth does not collapse.**
In the intergenic run AluYa5 lost 77% of its contacts, leaving room to argue power. Here
categories keep 69–83% and the effect vanishes anyway:

| subfamily | orient | n Q1→Q30 | kept | te_frac | z_gene |
|---|---|---|---|---|---|
| AluSc5 | sense | 4,019 → 3,344 | 83.2% | 0.042 → 0.025 | +6.38 → **+0.02** |
| AluSx1 | antisense | 72,068 → 59,694 | 82.8% | 0.034 → 0.026 | −3.53 → −12.14 |
| AluY | sense | 45,609 → 31,646 | 69.4% | 0.035 → 0.025 | +1.42 → **−8.74** |
| L1PA2 | sense | 2,556 → 985 | 38.5% | 0.074 → 0.031 | **+16.37 → +1.10** |
| L1PA6 (old) | sense | 16,168 → 11,507 | 71.2% | 0.017 → 0.016 | −13.32 → −10.50 |
| L1MA1 (old) | sense | 4,667 → 4,453 | 95.4% | 0.009 → 0.010 | −7.65 → −5.52 |

AluSc5 sense retains 83% of its contacts and lands exactly on its control rate. Old L1,
whose `te_frac` barely moves at all, stays put — the test is stable where mismapping is not
in play. That is a controlled comparison, not a power artifact. L1PA2 sense is the extreme
case: the single largest `z_gene` in the whole MAPQ≥1 table, and it evaporates.

**The dominant signal is negative**, in every family and both orientations — TE intron bases
are *less* distal than TE-free bases of the same host gene (MAPQ≥30):

| family | orient | n | median z_gene | median (te − ctl) |
|---|---|---|---|---|
| L1 | sense | 121 | −5.95 | −0.0130 |
| L1 | antisense | 113 | −4.12 | −0.0118 |
| SVA | antisense | 1 | −3.21 | −0.0283 |
| SVA | sense | 5 | −4.36 | −0.0215 |
| Alu | antisense | 33 | −3.98 | −0.0091 |
| Alu | sense | 33 | −4.04 | −0.0097 |

Same conclusion as the intergenic tests, now with the host controlled by construction rather
than by matching: TE-derived RNA is more locally tethered than the sequence around it.

## Sampling stability: the biorep split is NOT yet run here

The deposit has 4 bioreps and the machinery is in place — extract each biorep's pairs to
its own directory and pass `--arr`/`--pairs`/`--con` plus a distinct `--tag`. **It has not
been done on this dataset**, so nothing below is defended against sampling noise by
replication; the defence is depth (8 k–24 k contacts per tested category) and the shift
null. If you run it, compare `stat_frac` across shared categories.

Note what a replicate split would **not** show. Replicate agreement cannot test mismapping
— that is systematic, not stochastic, and reproduces perfectly across replicates aligned to
the same reference. MAPQ, cis-restriction and the within-gene control are what test it, all
three within this one dataset. Replicates test sampling stability, a different question.

## The karyotype is normal, and that is measured, not assumed

This matters because the per-copy screen below is a per-locus test, and an amplified region
manufactures exactly the recurrence it looks for. DNA-end density, 1 Mb bins, MAPQ≥30 on
both ends, as a multiple of the median non-empty bin (13,526 DNA ends over 2,950 bins):

| | value |
|---|---|
| max bin | **2.91×** (chr19:10–11 Mb) |
| 99th / 99.5th / 99.9th percentile | 2.22× / 2.36× / 2.55× |

No multi-Mb plateau, no abrupt drop-off — none of the amplicon signature. So
**`--exclude` is empty** and `--mask-mult 3.5` masks zero bins and drops zero source cells,
which is the correct outcome rather than a misconfiguration.

Do **not** retune `--mask-mult` down to this line's own 99th percentile (2.22×): the top
bins here are ordinary gene-dense regions (chr19, chr3:47–49 Mb), and masking them removes
real biology.

**If this stage is ever pointed at a rearranged line, name the coordinates — do not trust a
density threshold to find them.** The two failure modes squeeze from both sides: a cutoff
loose enough to spare gene-dense bins can sit *above* an amplicon's actual density and mask
nothing, while a cutoff tight enough to catch it starts eating the most transcribed regions
in the genome. Karyotype is a property of the cell line, not of the data, so it has to be
asserted from outside rather than inferred from a quantile. Masking must also be applied to
the **source** side, not just to targets: a copy inside an amplified region is the thing
generating the spurious recurrence, so excluding it as a target while still testing it as a
source does not help.

**And note what that costs in trust, not just in bins.** Copy-number inflation and
derivative-chromosome topology act directly on the readout of this assay — a genomic
distance between two ligated ends — so on a line whose CNV is not fully characterised, the
uncertainty is not confined to the masked regions. That is why this stage does not carry any
K562-derived result forward (see the provenance rule); a coordinate mask handles the part
you can see.

## Per-copy screen: clean of the karyotype, but not of the host

`per_copy_targets.py` asks the one question the aggregate tests cannot: does a *single* copy
hit a recurrent target set, as an architectural lncRNA gene would? Statistic is the
Herfindahl concentration of its DNA ends at 1 Mb, against a null preserving the copy's
cis/trans split and full distance profile and randomising only destination.

`per_copy_targets.tsv` — 123,547 cells, **162 at q_emp < 0.05 (0.13%)**, hits scattered
(most-loaded 10 Mb window holds 4 of the top 100; only 5 of the top 100 lie within 10 kb of
another). Against the *global* reach-matched null 123,062 of 123,547 cells "hit", which is
the script's own documented reason for preferring the matched-depth peer reference.

**Those 162 are not architectural-RNA candidates, because this screen has no host control.**
142 of them (88%) are intronic and only 2 intergenic — and both spatial clusters in the top
100 sit in the most highly expressed loci available to this cell type: **JARID2** introns
(chr6:15.51–15.52 Mb, a PRC2 subunit essential in hESC) and **PVT1** (chr8:128.02 Mb, the
MYC-adjacent super-enhancer lncRNA). `within_intron_test.py` exists precisely because host
expression drives distal fraction 0.195 → 0.817; the Herfindahl index is exposed to the same
confound, and the peer reference matches **depth, not host gene**. Any follow-up needs the
within-gene control ported across, not a random-copy null.

`per_copy_targets_sva.tsv` adds SVA (`--families Alu,L1,SVA`, not in the default `Alu,L1`):
**1 of 191 SVA cells** clears q_emp < 0.05 — SVA_F sense, chr6:76,743,673, 182 contacts,
z_emp 6.46 — and it too is intronic, hosted in an unnamed lncRNA (ENSG00000302959). Note
`z_emp` is computed over all rows pooled into half-log2 depth strata **regardless of
family**, so adding SVA shifts every Alu/L1 peer mean: the two tables are not row-comparable,
and the default-families one is the reference table.

## Running it

```bash
python src/hic/fetch_4dn.py get 4DNFIGDJIRV3 --out data/hic   # needs deps-3d (requests)
python src/imargi/check_pairs.py            # -> results/imargi/pairs_checks.txt
python src/imargi/extract_rna_ends.py       # -> data/imargi/
python src/imargi/extract_contacts.py       # -> data/imargi_contacts/
python src/imargi/trans_test.py --statistic cisdistal --shifts 1000
python src/imargi/trans_test.py --statistic distal --shifts 1000 --min-mapq 30
python src/imargi/within_intron_test.py --shifts 200 --min-ctl 100 --min-mapq 30
python src/imargi/per_copy_targets.py
```

Every script now defaults to this dataset, so each runs bare. `--arr`, `--pairs`, `--con`
and `--tag` exist for biorep splits; `--pairs` is read for its `#chromsize` header only,
which the toroidal shift depends on, so it must match whatever `--arr` was built from.

**There is no `--out`, so a run always overwrites the canonical table.** Do not
"smoke-test" with a reduced `--shifts` against the default paths — it silently replaces a
1000-shift null with a 25-shift one, leaving the observed columns right and the inference
wrong. Also note the tallies filename is keyed on class + MAPQ + tag but **not** on
statistic (they are statistic-independent by design, which is what makes re-reducing free),
so two *concurrent* runs differing only in `--statistic` race on the same `.npz` — give one
a distinct `--tag`.

## The references, and why both

### An abundance effect large enough to swamp everything

The pre-mRNA baseline's distal fraction rises monotonically with contact count, from
**0.195 in the lowest decile to 0.817 in the highest** — a 4× range (on the cis-restricted
statistic, 0.074 → 0.421, a 6× range). Highly transcribed genes are the distal-heavy ones,
not the noise-dominated low-abundance ones. Any unstratified comparison of distal fractions
here measures expression level and nothing else. Every comparison above is decile-matched.

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
category (45.1%)** of any measured, which is what host pre-mRNA passengers predict, and
**SVA antisense is the least (25.6%)** — the same ordering the tests above recover from a
different statistic. Both come from this dataset, so the check is genuinely internal.
(`src/chrrna` reported a compatible L1-sense result in K562; noted, not counted.)
