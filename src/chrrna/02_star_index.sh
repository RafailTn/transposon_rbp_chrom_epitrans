#!/usr/bin/env bash
# Build the hg38 STAR index. Run once on the server; ~32 GB RAM, ~1 h.
#
# No GTF is passed to genomeGenerate on purpose. Splice junctions from a gene
# annotation would bias alignment toward annotated exons, and most TE copies are
# intronic or intergenic -- the thing being measured here is precisely the reads
# that a gene-model-guided aligner would down-weight. Quantification happens
# afterwards with featureCounts against the TE SAF.
set -euo pipefail

GENOME_FA="${1:?usage: 02_star_index.sh <hg38.fa> [outdir]}"
IDX="${2:-$PWD/star_index_hg38}"
THREADS="${THREADS:-16}"

command -v STAR >/dev/null || { echo "missing STAR" >&2; exit 1; }

mkdir -p "$IDX"
STAR --runMode genomeGenerate \
     --runThreadN "$THREADS" \
     --genomeDir "$IDX" \
     --genomeFastaFiles "$GENOME_FA"

echo "index written to $IDX"
