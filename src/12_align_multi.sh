#!/usr/bin/env bash
# Second alignment for the EM arm: identical to 03_align.sh except that STAR
# reports EVERY alignment of a multi-mapping read instead of one at random.
#
# 03_align.sh keeps the authors' --outSAMmultNmax 1 --outMultimapperOrder Random,
# which spreads a repeat read over equivalent loci at random -- not only within a
# subfamily but across subfamilies that share sequence (AluY vs AluSx). That is
# kept as the primary arm for comparability with the published L1 numbers. This
# arm hands all alignments to TEcount (13_em_count.sh), whose EM assigns them in
# proportion to the evidence instead. Subfamilies whose retention differs much
# between the two arms are quantification-sensitive and should be reported so.
#
# Multi-mapping limits follow the TEtranscripts recommendation (100 / 100)
# rather than the authors' 10000; reads with >100 equally good loci carry almost
# no subfamily information. Everything else -- EndToEnd, the index, the trimmed
# FASTQs -- is shared with 03, so the arms differ only in multi-mapper handling.
# Output is unsorted (mates adjacent), as TEcount expects.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IDX="${1:-$ROOT/data}"
FQ="${2:-$ROOT/work/fastq_trimmed}"
OUT="${3:-$ROOT/work/bam_multi}"
THREADS="${THREADS:-16}"

command -v STAR >/dev/null || { echo "missing STAR" >&2; exit 1; }
mkdir -p "$OUT"

tail -n +2 "$ROOT/config/runs.tsv" | while IFS=$'\t' read -r sample fraction rep run srx layout; do
    [ -z "${sample:-}" ] && continue
    bam="$OUT/${sample}.Aligned.out.bam"
    if [ -s "$bam" ]; then
        echo "== $sample already aligned, skipping"
        continue
    fi
    r1="$FQ/${sample}_1.fastq.gz"
    r2="$FQ/${sample}_2.fastq.gz"
    for f in "$r1" "$r2"; do
        [ -s "$f" ] || { echo "missing $f -- run 01_trim.sh first" >&2; exit 1; }
    done
    echo "== aligning $sample (all multi-mapper alignments)"
    STAR --runThreadN "$THREADS" \
         --genomeDir "$IDX" \
         --readFilesIn "$r1" "$r2" \
         --readFilesCommand zcat \
         --runMode alignReads \
         --outSAMtype BAM Unsorted \
         --outFilterMultimapNmax 100 \
         --winAnchorMultimapNmax 100 \
         --alignEndsType EndToEnd \
         --outSAMattributes NH HI AS nM \
         --outFileNamePrefix "$OUT/${sample}."
done
