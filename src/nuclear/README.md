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

All md5-verified; see `data/nuclear_manifest.tsv`. Plus **16-fraction Repli-seq**
(`data/repliseq/`, G1 + P2_S…P17_S × 2 reps, bedGraph, 5 kb bins), used as the
replication-timing covariate.

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

## How GC is controlled, and the check that it works

**Stratified mean-subtraction** — a piecewise-constant regression of TSA signal on GC:

1. GC fraction per 20 kb bin, `G+C / (G+C+A+T)`, N excluded from the denominator.
2. Equal-count GC strata cut on percentiles **of the usable bins only**, so the masked
   centromeric and gap bins cannot shift the strata.
3. Per stratum and per axis, the mean signal over all usable bins in it.
4. `resid` for a copy = (its bin's value) − (its stratum's mean).
5. Copy-weighted average of those residuals, bin-clustered SE.

So the residual is the part of a subfamily's signal not attributable to the GC band its
copies sit in.

**The failure mode this has, quantified.** Equal-*count* strata are wildly unequal in
*width*: at 10 deciles the middle ones span ~0.013 GC units but the top spans **0.194**
(0.486→0.680) and the bottom 0.089. The tails are exactly where Alu (GC-rich) and L1
(AT-rich) concentrate, and where the GC→signal slope is steepest — so coarse strata leave
residual within-stratum GC in the worst possible place.

**Tested by refining the stratification until the estimate stops moving:**

| strata | Alu speckle | L1 speckle | SVA speckle | Alu lamina | L1 lamina |
|---|---|---|---|---|---|
| 10 | +0.1162 | +0.0451 | +0.0705 | −0.1816 | −0.1214 |
| 50 | +0.1209 | +0.0495 | +0.0749 | −0.1816 | −0.1236 |
| 100 | +0.1214 | +0.0497 | +0.0755 | −0.1816 | −0.1236 |
| 500 | +0.1214 | +0.0496 | +0.0748 | −0.1812 | −0.1233 |

Converged by ~50 and flat to 500. Ten deciles understated the speckle residual by ~0.005
(4–10% of the estimate) and the lamina residual by ~0.002; no conclusion changes. The
default is therefore **`--gc-strata 50`**, and the numbers below are the converged ones.

Residual within-stratum GC imbalance at 10 deciles confirms the diagnosis: copy-weighted
minus bin-mean GC is ≤0.0008 in deciles 0–8 but −0.006 to −0.009 in decile 9, the wide one.

**What this does not control.** GC is a proxy. Gene density, replication timing, CpG-island
density and chromosome identity all co-vary with it and are *not* held fixed — the L1M3a /
chr19 case below is the proof, a subfamily at below-average GC whose extreme residual is
pure chromosome identity. The strata count was never the real limitation; the choice of a
single scalar covariate is.

## Replication timing as the covariate: weaker alone, better jointly

`repliseq_timing.py` builds a per-bin replication-timing profile from the 16-fraction
Repli-seq (G1 + P2_S…P17_S, 2 reps each; the 17 "sets" are one experiment split by
fraction). Per bin the 16 S-phase signals are depth-normalised, converted to proportions and
reduced to `RT = Σ i·p_i`, a continuous average replication time in fraction units.

**Direction is measured, not assumed:** corr(RT, GC) = **−0.542**, and early replication is
GC-rich, so low RT = early and the deposited order P2→P17 runs early→late. The script exits
rather than guess if that correlation is too weak to fix the sign.

RT was expected to beat GC, being a direct measurement of compartmentalisation rather than a
sequence proxy. **It does not** — variance explained in the TSA axes, 50 strata:

| axis | GC | RT | GC × RT jointly |
|---|---|---|---|
| speckle | **0.700** | 0.517 | **0.771** |
| lamina | **0.426** | 0.406 | **0.530** |

GC is the stronger single covariate on both axes. But RT is not redundant: adding it lifts R²
by 0.07 (speckle) and 0.10 (lamina), so it carries compartment information GC misses. Hence
`--covariate {gc,rt,both}`, **default `both`** (15×15 = 225 cells; cells under `--min-cell`
fall back to the marginal GC-stratum mean, which affected 1 cell of 225).

Requiring finite RT above `--min-s-total` costs 2,506 bins, leaving 139,980 (90.6%).

### Better control roughly halves every residual

| | GC only | RT only | GC × RT |
|---|---|---|---|
| Alu speckle | +0.121 | +0.072 | **+0.049** |
| L1 speckle | +0.050 | −0.014 | **+0.038** |
| SVA speckle | +0.075 | +0.084 | **+0.023** |
| Alu lamina | −0.182 | −0.122 | **−0.085** |
| L1 lamina | −0.124 | −0.042 | **−0.102** |
| SVA lamina | −0.150 | −0.175 | **−0.086** |

**More than half of what the GC-only analysis called "residual" was compartment structure GC
could not see.** All SEs are ≤0.011, so these shifts are far outside the error bars. The
GC-only numbers in the sections above are superseded by the `both` column; they are kept
because the size of the gap is the point.

Note L1 under RT-only: speckle **−0.014** and lamina **−0.042**, i.e. essentially at its
control. L1 is the family whose apparent position is most fully explained by replication
timing — unsurprising, since L1 density and late replication are close to the same variable.

## Result 2: after full control, all three families are mildly lamina-*depleted*

Conditioned on GC × RT jointly (the `both` default), with cluster-robust SE:

| family | speckle | lamina | nucleolus (MKI67IP) |
|---|---|---|---|
| Alu | +0.049 ± 0.001 | **−0.085 ± 0.002** | +0.034 ± 0.001 |
| L1 | +0.038 ± 0.001 | **−0.102 ± 0.002** | +0.042 ± 0.001 |
| SVA | +0.023 ± 0.005 | **−0.086 ± 0.011** | +0.030 ± 0.006 |

**The naive picture inverts for L1.** "L1 lives at the nuclear lamina" is not a nuclear
finding: at matched composition and replication timing, L1-containing neighbourhoods are
*less* lamina-proximal than average (−0.102, ~50 SE). All three families sit slightly
speckle-ward and nucleolus-ward of matched sequence.

**And the three families have become nearly indistinguishable.** Under GC only the lamina
residuals were −0.182 / −0.124 / −0.150 (Alu/L1/SVA); jointly they are −0.085 / −0.102 /
−0.086. The remaining spread is comparable to the SEs. So essentially all of the apparent
family difference in nuclear position is composition plus replication timing, and what is
left is a small, *shared* offset rather than a family-specific one.

Effect sizes are 0.02–0.10 log2 units against a covariate-driven range of ~1.9, i.e.
**compartment structure explains one to two orders of magnitude more than family identity
does.** The SEs are tiny only because n is enormous; significance here is not importance.

SVA tracks Alu rather than L1 — consistent with its Alu-derived, GC-rich composition rather
than with its being a young hominid-specific element.

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
(ρ = +0.486, p = 0.33). Those figures are GC-only.

**Under joint GC × RT control the pattern disappears altogether.** The positive tail of the
lamina ranking collapses: the maximum subfamily lamina residual falls from **+0.202** (GC) to
**+0.030** (joint), and the young L1s that headed the GC-only list land on their control —
L1PBb +0.030 ± 0.077, L1PA8A +0.028 ± 0.019, L1PA7 +0.013 ± 0.007, L1PA5 +0.001 ± 0.007.
**No subfamily is lamina-enriched once composition and replication timing are held fixed.**
The subfamily spread also narrows (sd 0.114 → 0.088). Do not report a young-L1-at-the-lamina
result.

What survives is one-sided: a group of **old L1M subfamilies remains distinctly
lamina-depleted** beyond both covariates — L1M4c −0.398 ± 0.029, L1M3de −0.325 ± 0.066,
L1MDb −0.301 ± 0.042, L1MB2 −0.264 ± 0.022, all with `top_chrom_frac` 0.08–0.14, so not
chromosome artifacts. That is the only subfamily-level signal in this stage that survives
full control.

## Trap: neither GC nor RT absorbs chromosome identity

**L1M3a is the largest outlier on every axis and every control** — lamina residual −0.933
(GC), −0.770 (RT), **−0.910 (joint)**; lamina−nucleolus −1.365 jointly, 2–3× any other
subfamily — at a below-average GC of 0.384, which makes no compositional sense. The reason is
that **26% of its copies (18% of its bins) are on chr19**, the most gene-dense,
speckle-proximal chromosome in the genome.

**Replication timing does not fix this**, which was the hope when the Repli-seq was pulled:
the residual barely moves between covariates. Chromosome identity is a distinct confound from
both composition and timing, and the fix is to stratify on chromosome explicitly (not yet
done).

`top_chrom_frac` is reported so the failure mode is visible. Median across 179 subfamilies is
0.090 (≈ what chromosome sizes predict) and the max is 0.18, so the flag threshold is 0.15 —
8 subfamilies trip it. **A tripped row describes a chromosome, not a repeat family.**

## Running it

```bash
python src/nuclear/bin_covariates.py        # -> cache/bin_covariates_20000.npz (~4 min)
python src/nuclear/repliseq_timing.py       # -> cache/repliseq_rt_20000.npz     (~6 min)
python src/nuclear/te_nuclear_position.py   # -> results/nuclear/te_nuclear_position_both.tsv
```

`--covariate {gc,rt,both}` selects the control and is part of the output filename, so the
three runs cannot overwrite each other. Default `both`.

`--gc-strata` defaults to 50 (see the convergence check above). `bin_covariates.py` is cached
and depends only on the assembly and the rmsk table; rerun
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
