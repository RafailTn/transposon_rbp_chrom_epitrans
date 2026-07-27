#!/usr/bin/env bash
# Align the 6 K562 runs with STAR, reproducing GSE296212's published parameters.
#
# The parameters are the authors' own, copied verbatim from the series
# !Sample_data_processing so our Alu/SVA numbers sit on the same footing as their
# published L1 numbers:
#
#   --outFilterMultimapNmax 10000  --winAnchorMultimapNmax 20000
#   --alignTranscriptsPerReadNmax 1000000  --alignEndsType EndToEnd
#   --outSAMmultNmax 1  --outMultimapperOrder Random
#
# READ THIS BEFORE INTERPRETING ANYTHING DOWNSTREAM. The last two are the
# consequential pair: STAR finds up to 10000 alignments for a repeat-derived read
# but reports exactly ONE, chosen at random among them. So:
#
#   * Subfamily-level totals are approximately unbiased -- random assignment
#     spreads a read evenly over the equivalent copies it could have come from,
#     and aggregating puts them back together.
#   * PER-COPY counts are close to meaningless. A read is credited to an
#     arbitrary member of its subfamily, not the one it came from. Do not read
#     the per-copy column as locus-level expression; it exists only so the
#     downstream filter can drop copies with no support at all.
#
# This is also why the authors' own L1HS_1..L1HS_332 per-locus columns should be
# treated with the same caution.
#
# --outSAMattributes NH HI AS nM keeps NH, so multi-mapping depth stays
# inspectable after the fact even though only one alignment is kept.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IDX="${1:?usage: 03_align.sh <star_index> [fastq_dir] [outdir]}"
FQ="${2:-$PWD/fastq}"
OUT="${3:-$PWD/bam}"
THREADS="${THREADS:-16}"

command -v STAR >/dev/null || { echo "missing STAR" >&2; exit 1; }
mkdir -p "$OUT"

tail -n +2 "$HERE/runs.tsv" | while IFS=$'\t' read -r sample fraction rep run srx layout; do
    [ -z "${sample:-}" ] && continue
    bam="$OUT/${sample}.Aligned.sortedByCoord.out.bam"
    if [ -s "$bam" ]; then
        echo "== $sample already aligned, skipping"
        continue
    fi
    r1="$FQ/${sample}_1.fastq.gz"
    r2="$FQ/${sample}_2.fastq.gz"
    for f in "$r1" "$r2"; do
        [ -s "$f" ] || { echo "missing $f -- run 00_fetch.sh first" >&2; exit 1; }
    done
    echo "== aligning $sample"
    STAR --runThreadN "$THREADS" \
         --genomeDir "$IDX" \
         --readFilesIn "$r1" "$r2" \
         --readFilesCommand zcat \
         --runMode alignReads \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 10000 \
         --alignTranscriptsPerReadNmax 1000000 \
         --outSAMmultNmax 1 \
         --outMultimapperOrder Random \
         --winAnchorMultimapNmax 20000 \
         --alignEndsType EndToEnd \
         --outSAMattributes NH HI AS nM \
         --outFileNamePrefix "$OUT/${sample}."
done

echo
echo "alignment summary (uniquely + multi-mapped %):"
for log in "$OUT"/*.Log.final.out; do
    s=$(basename "$log" .Log.final.out)
    u=$(awk -F'\t' '/Uniquely mapped reads %/{print $2}' "$log")
    m=$(awk -F'\t' '/% of reads mapped to multiple loci/{print $2}' "$log")
    printf "  %-26s unique %-8s multi %s\n" "$s" "$u" "$m"
done
