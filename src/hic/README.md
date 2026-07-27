# 3D genome stage: TE subfamilies vs compartments and promoter contacts

Asks, per TE subfamily, two questions of a Hi-C/HiChIP/Micro-C contact matrix:

1. **Which compartment do its copies live in?** (`te_compartments.py`)
2. **Do its copies contact promoters more than the same loci relocated?**
   (`te_tss_contacts.py`)

This is a **separate stage** from the eCLIP analysis. It shares `cache/te_index.pkl`,
the `searchsorted` overlap idiom and `tecommon.bh`, and nothing else.

## Interpreter — NOT `deps/`

Runs under **`deps-3d/`**, which exists solely because the two cannot be merged:

```
deps-3d/.pixi/envs/default/bin/python   # python 3.12.13, numpy 1.26.4, cooler 0.10.4, cooltools 0.7.1
```

Every packaged `cooltools` caps at `numpy <2.0`. `deps/` pins `numpy ==2.4.6` because
`enrich_permutation.py` draws its circular shifts from `Generator.integers` and
`enrich_stranded.py` calls `hypergeom.sf`, and both feed `results/eclip/`. NumPy freezes
`RandomState` but explicitly does **not** guarantee `Generator` method streams across
feature releases, so moving `deps/` to numpy 1.x would silently respecify the permutation
null and permanently break the "rerun reproduces the tables" property.

`pixi add cooltools` into `deps/` fails to solve rather than downgrading anything — the
`==` pins protect themselves. That failure is the correct outcome. Do not "fix" it.

## Pipeline

```bash
bash src/hic/fetch_gencode_tss.sh                       # -> data/hic/tss_protein_coding.bed
python src/hic/fetch_4dn.py search --experiment-type HiChIP --biosource K562
python src/hic/fetch_4dn.py get 4DNFIxxxxxxx --out data/hic

python src/hic/te_compartments.py  data/hic/<f>.mcool --tss data/hic/tss_protein_coding.bed --resolution 250000
python src/hic/te_tss_contacts.py  data/hic/<f>.mcool --tss data/hic/tss_protein_coding.bed --resolution 10000
```

Start with `te_compartments.py`. It is minutes of compute, tolerates mappability loss
because its bins are 250 kb, and if its masking report shows subfamilies losing wildly
different fractions of copies then `te_tss_contacts.py` will not be interpretable either.

## What the design decisions were, and why

**Resolution is not a free parameter.** Alu is ~300 bp and there is one Alu per ~3 kb of
genome; at any resolution a real cooler offers, an Alu is sub-bin and a 10 kb bin holds
about three. Everything here is therefore an **aggregate over copies**, never a per-copy
call. L1 full-length (~6 kb) and SVA are the elements where a bin means something.

**Distance is the dominant confound**, and obs/exp against `cooltools.expected_cis`
removes it. Without it, comparing subfamilies mostly re-measures the fact that Alu is
isochore-biased into gene-dense sequence and L1 is not.

**Obs/exp does not remove composition.** Alu-dense bins are also gene-dense, open,
high-coverage bins. That is what the shift null is for.

**The null is a rigid toroidal shift** of the TE anchors, one offset per (chromosome,
replicate), preserving the spacing between copies. Deliberately the same construction as
`EligibleSpace.shift` in `src/enrich_permutation.py`. It inherits the same known
conservative bias discussed there. If you change the reasoning in one, revisit the other.

**ICE masking is a per-subfamily confound, not a nuisance.** ICE gives NaN weights to
low-coverage bins, and in unique-alignment data repeat-dense bins are disproportionately
low-coverage — the 3D form of the mappability caveat in `CLAUDE.md` (L1H at 0.36%
coverage against L1M's 1.85%, zero L1H hits). If subfamily A loses 5% of copies to masking
and B loses 20%, comparing their averages compares masking rates. Both scripts write
`<sample>_te_bin_masking.tsv` before computing anything; read it first.

**Eigenvector phasing is required, and GC would be the wrong track.** Cis eigenvector signs
are arbitrary and independently arbitrary per chromosome, so unphased "E1 > 0" means A on
some chromosomes and B on others. We phase on TSS density rather than GC because Alu is
GC-rich and L1 is AT-rich *by composition* — phasing on GC and then asking whether Alu
prefers A partly tests the track against itself. TSS density is not fully free of this
(Alu is gene-dense by isochore) so a positive Alu compartment result still needs that
caveat stated.

## Traps

- **The cooler knows ~24 chromosomes; the TE index knows 631.** "Canonical" in
  `hg38_L1_Alu_SVA_canonical.gtf` means canonical *subfamilies*, not chromosomes. Every TE
  copy on a sequence absent from the cooler silently leaves the numerator — the mirror of
  the `INDEX.get(c)` trap in `CLAUDE.md`, running the other way. `reconcile_chroms` prints
  the loss; it is not uniform across subfamilies, so read the number rather than assuming
  it is small.
- **Never fetch a chromosome-sized dense square.** chr1 at 10 kb is 4.9 GB; the
  genome-wide square is 764 TB. `evaluate()` fetches `--window-bins` squares (32 MB at the
  default) and relies on `j - i <= maxd` to guarantee each pair lands wholly inside one
  window. This is the only way this script OOMs, and it is a one-line mistake.
- **`nb = ceil(L / binsize)` means the last bin is partial**, so `nb * binsize` runs past
  the chromosome and cooler rejects the region outright. Clamp to `chromsizes`.
- **Do not pass `--nproc` your core count.** Each cooltools worker holds a
  chromosome-scale working set; you are I/O-bound on HDF5 reads well before core-bound.
  4–8 is right on a 20-core box.
- **`--n-shifts` sets a p-floor of `1/(K+1)`** that must clear BH against the number of
  subfamilies tested, exactly as `--pmax` must in `enrich_permutation.py`. Too low does not
  error; it silently ties subfamilies on the floor so they cannot be ranked. The script
  warns when rows land there.
- **Subfamilies are not independent.** One matrix, one TSS set, shared chromosomes: a
  genome-wide offset moves every row at once, so "most subfamilies are above their null"
  is one observation, not a hundred. Differences *between* subfamilies are the
  interpretable quantity.
- **`export LC_ALL=C` before any `awk` over these TSVs**, per `CLAUDE.md` — under a
  comma-decimal locale `2.107e-03` parses as `2`.

## Built-in correctness check

`te_tss_contacts.py --scramble-anchors` relocates the observed anchors once before
testing, so the observed arm is itself a draw from the null and the true answer is exactly
zero. A correctly specified test returns no `q_emp < 0.05`. **Run it after any change to
the null.**

Note what it does *not* promise: because subfamilies are correlated, the median
`log2_ratio` under scramble wanders substantially with the seed (−0.071 to −0.022 on the
synthetic fixture, with 6/121 to 48/121 subfamilies positive). Judge calibration on
`q_emp`, not on the sign of the median.

## What this cannot show

A contact is not regulation. These scripts can support "copies of subfamily X are in
closer-than-chance contact with promoters" — a real claim. They cannot show the element
*modulates* expression: that needs perturbation (CRISPRi at the element) or natural
variation (polymorphic MEI genotypes crossed with eQTLs). The MEI route is complementary
rather than redundant, since polymorphic insertions concentrate in exactly the young
subfamilies that unique-alignment coverage hides here.

The strong design across tissues is an **interaction**, not a marginal enrichment: is
subfamily X enriched at promoter contacts specifically in the tissues where the linked
genes are expressed? Alu's GC bias and its mappability profile do not change between
tissues, so an interaction self-controls for both. A marginal enrichment does not.
