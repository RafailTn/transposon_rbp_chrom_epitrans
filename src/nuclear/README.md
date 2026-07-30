# nuclear stage: where do L1 / Alu / SVA DNA copies sit in the nucleus?

A **DNA-position** question, deliberately separate from the RNA stages. TSA-seq has no
strand, so nothing here bears on the sense/antisense hypothesis `src/imargi` and
`src/chrrna` test. This is the better-instrumented version of what `src/hic` asks: actual
distance to nuclear bodies rather than A/B compartment.

## Data

H1-hESC, 4DN, GRCh38 — **same cell line as the iMARGI stage**, which is why it can later be
joined to it. Nothing has been joined yet; this stage stands alone.

| axis | target | files | set |
|---|---|---|---|
| speckle | SON | `4DNFI625PP2A`, `4DNFIFKMOD1L` | `4DNESC3D6NGQ` |
| lamina | Lamin B1 | `4DNFINXMR3OD`, `4DNFILYDJU8T` | `4DNESGGXKI1H` |
| nucleolus | MKI67IP (NIFK) | `4DNFIFNIK4HD`, `4DNFIGU18ZPJ` | `4DNESO6HFSAD` |
| nucleolus | POL1RE (POLR1E) | `4DNFICLC27JC`, `4DNFICYLNXF2` | `4DNESIH5OOWT` |

All md5-verified; see `data/nuclear_manifest.tsv`. Also pulled but **not yet used**:
16-fraction Repli-seq (`data/repliseq/`, 17 fractions × 2 reps, bedGraph, 5 kb bins).

**The tracks are log2 enrichment scores, not counts**, despite the portal labelling them
`normalized counts` / "20kb bin (read normalized) count track": 70.7% of bins are negative,
range −1.99…+2.26. Whether the denominator is the matched input could not be determined
from portal metadata — the input sets (`4DNES6PANOF4`, `4DNESBRNMFK6`) carry **raw fastq
only**, and the linked workflow is a provenance tracker. Treat the values as relative, not
absolute distances.

**Two nucleolar antibodies were taken deliberately, and they agree.** Replicate r = 0.952
(MKI67IP) and 0.921 (POL1RE); cross-antibody r = 0.855–0.885. Replicates agree more than
antibodies do, which is the right ordering, and two different nucleolar proteins converging
at r ≈ 0.87 means the axis is real. They are **not averaged into one column**, so a
family-specific disagreement would stay visible. None appeared: the MKI67IP−POL1RE offset is
−0.13/−0.10/−0.11 for Alu/L1/SVA, i.e. a constant, not a family effect.

## Three things that had to be got right first

**20 kb bins forbid per-copy claims.** An Alu is ~300 bp against a 20 kb bin, so every copy
inherits its neighbourhood's value and two copies in one bin are indistinguishable by
construction. The subfamily aggregate is the only meaningful unit. Same objection as taking
an mcool over pairs in `src/imargi`.

**Masking, because the tracks mask nothing themselves** — they carry a value for 100% of
bins. Unmasked, nucleolar signal reads **+1.55 across the chr1 centromere against +0.16 on
a normal arm**, and chr21's p-arm is 86.5% exact zeros. `bin_covariates.py` drops bins over
10% assembly-N or 10% satellite (satellite from the *full* rmsk table — GRCh38 models
centromeric alpha-satellite rather than N-masking it, so N alone does not catch centromeres).
92.3% of bins survive; 99.6% of TE copies land in a surviving bin.

Note what the chr21 row costs: the most repeat-dense nucleolar material is rDNA-adjacent, on
acrocentric short arms that GRCh38 does not assemble. The nucleolar axis here is the
*assembled* fraction of the compartment, systematically missing its most repeat-rich part.

**Copy-weighted mean, bin-clustered SE.** Both obvious aggregations were tried:

- Deduplicating to bins asks "what is a bin containing ≥1 copy like?" — for Alu that is 97%
  of all usable bins, so the answer is "like the genome" by construction. Alu's speckle mean
  collapsed from +0.180 to −0.120 against a −0.126 background. The signal did not turn out
  to be fake; presence/absence throws away copy *number*, which is the entire compositional
  signal for a ubiquitous family.
- Copy-weighting asks "what is the mean environment of a randomly chosen copy?" — the actual
  question, and the point estimate reported.

So the clustering problem is about the **error bar**, not the estimate. Every SE is
cluster-robust with the bin as cluster; `n_eff` (Kish) is ~74 k for Alu against 1.13 M
copies. Quote `n_eff`.

## Result 1: the raw ranking is 1980s compositional genomics

GC is not a nuisance here, it is nearly the whole signal. Over usable bins, mean TSA-seq by
GC decile:

| GC decile | 1 (0.26–0.35) | … | 10 (0.49–0.68) |
|---|---|---|---|
| speckle | −0.67 | monotone | **+0.92** |
| lamina | +0.56 | monotone | **−1.32** |

Both axes traverse ~1.6–1.9 log2 units monotonically across GC. Against that, the raw family
means:

| family | GC | speckle | lamina | nucleolus |
|---|---|---|---|---|
| background | 0.409 | −0.126 | −0.256 | −0.196 |
| Alu | 0.430 | **+0.180** | **−0.703** | −0.081 |
| L1 | 0.398 | −0.168 | −0.267 | −0.202 |
| SVA | 0.434 | +0.164 | −0.738 | −0.058 |

Alu speckle-proximal and lamina-depleted, L1 the reverse. This is real and it is also exactly
what Alu's GC-richness and L1's AT-richness predict without any nuclear biology. **Quoting
this table as a nuclear-organisation finding would be a mistake.**

## Result 2: after GC control, all three families are lamina-*depleted*

The GC-decile residual, with cluster-robust SE:

| family | speckle | lamina | nucleolus (MKI67IP) | lamina − nucleolus |
|---|---|---|---|---|
| Alu | +0.116 ± 0.002 | **−0.182 ± 0.003** | +0.010 ± 0.001 | −0.192 ± 0.005 |
| L1 | +0.045 ± 0.001 | **−0.121 ± 0.003** | +0.038 ± 0.001 | −0.159 ± 0.005 |
| SVA | +0.071 ± 0.006 | **−0.152 ± 0.012** | +0.009 ± 0.008 | −0.161 ± 0.021 |

**The naive picture inverts for L1.** "L1 lives at the nuclear lamina" is entirely a
composition effect: at matched GC, L1-containing neighbourhoods are *less* lamina-proximal
than average (−0.121, ~40 SE). All three families sit slightly speckle-ward and
nucleolus-ward of GC-matched sequence, and they barely differ from each other on the
lamina−nucleolus contrast (−0.16 to −0.19) — the contrast chosen precisely because GC does
not predict which heterochromatic destination a locus takes. So the large apparent family
differences are composition; what remains is small and shared.

Effect sizes are 0.05–0.19 log2 units against a GC-driven range of ~1.9, i.e. **composition
explains roughly an order of magnitude more than family identity does.** The SEs are tiny
only because n is enormous; statistical significance here is not the same as importance.

SVA tracks Alu, not L1, on every axis — consistent with its Alu-derived, GC-rich composition
rather than with its being a young hominid-specific element.

## Result 3: the L1 age gradient does not survive testing

The top and bottom of the lamina-residual ranking look like a clean age story — young L1 at
the lamina (L1PBb +0.200, L1PA8A +0.098, L1HS +0.095, L1PA7 +0.070, L1PA2 +0.051) and old
L1M\* far from it (L1M4c −0.467, L1MB2 −0.366, L1MB4 −0.300, L1M4 −0.297). **Tested against
the known L1HS→L1PA17 age ordering, it is not significant:**

| axis | Spearman(age, residual) | p | n |
|---|---|---|---|
| lamina | −0.456 | 0.076 | 16 |
| speckle | +0.553 | 0.026 | 16 |
| nucleolus | −0.326 | 0.22 | 16 |

Only speckle reaches nominal significance, on 16 points, and would not survive correction
across three axes. Within the series the values are non-monotone (L1PA4 −0.004, L1PA6
−0.018, L1PA12 −0.037, L1PA17 −0.061) and span only ±0.05. SVA_A→SVA_F shows nothing
(ρ = +0.486, p = 0.33).

What does hold is coarser: **L1M (oldest) vs L1P/L1HS differ on the lamina axis**, with L1M
strongly depleted. That is a two-group contrast, not a graded gradient — do not describe it
as one.

## Trap: GC matching does not absorb chromosome identity

**L1M3a is the largest outlier on every axis** (speckle residual +0.237, lamina −0.930,
lamina−nucleolus −1.445, each 2–3× any other subfamily) at a below-average GC of 0.384 —
which makes no sense compositionally. The reason is that **26% of its copies (17.7% of its
bins) are on chr19**, the most gene-dense, speckle-proximal chromosome in the genome, whose
compartment behaviour is not captured by bin-level GC.

`top_chrom_frac` is reported for this. Median across 179 subfamilies is 0.090 (≈ what
chromosome sizes predict) and the max is 0.177, so the flag threshold is 0.15 — 8 subfamilies
trip it. **A tripped row describes a chromosome, not a repeat family.** A chromosome-matched
or replication-timing-matched control would be the fix; the Repli-seq is downloaded for it
but not yet used.

## Running it

```bash
python src/nuclear/bin_covariates.py        # -> cache/bin_covariates_20000.npz (~4 min)
python src/nuclear/te_nuclear_position.py   # -> results/nuclear/te_nuclear_position.tsv
```

`bin_covariates.py` is cached and depends only on the assembly and the rmsk table; rerun
with `--force` only if `--bin` changes. The output table has 183 rows: 1 background,
3 family, 179 subfamily (≥50 usable bins). Read `*_resid` for anything load-bearing and
check `top_chrom_frac` before believing a subfamily.

## What this stage cannot say

- Nothing per-copy. 20 kb bins.
- Nothing about RNA, and nothing about orientation. TSA-seq is DNA.
- Nothing about absolute nuclear distance — the values are relative enrichment scores whose
  input normalisation could not be confirmed.
- Nothing about the compartment's most repeat-dense regions: acrocentric p-arms and
  centromeric satellite are unassembled or masked out.
