#!/usr/bin/env bash
# Count GENCODE genes (exons) and host introns, to give the TE analysis a
# gene-based normalisation and an intronic baseline.
#
# Two outputs per sample:
#   <sample>.genes.txt    featureCounts over GENCODE exons, by gene_id, unique
#                         alignments only (no -M: standard gene counting).
#                         Used for gene-based size factors (10_context.py), and
#                         by 15_compare_geo.py against the authors' gene rows.
#   <sample>.intron.txt   featureCounts over data/gencode_intron.saf (gene body
#                         minus all exons, per strand; from 08_context.py), WITH
#                         -M like the TE counts, so intronic TE retention can be
#                         read against the retention of the introns around it.
#
# Same strand handling as 05_count.sh: STRAND_<sample> or STRAND_DEFAULT, and
# --countReadPairs so the unit is the fragment.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAMDIR="${1:-$ROOT/work/bam}"
OUT="${2:-$ROOT/work/gene_counts}"
GTF="${GTF:-/mnt/raid1/rafail/aluaca/ref/gencode.v47.primary_assembly.annotation.gtf.gz}"
INTRON_SAF="${INTRON_SAF:-$ROOT/data/gencode_intron.saf}"
THREADS="${THREADS:-16}"
STRAND_DEFAULT="${STRAND_DEFAULT:-2}"

command -v featureCounts >/dev/null || { echo "missing featureCounts (subread)" >&2; exit 1; }
for f in "$GTF" "$INTRON_SAF"; do
    [ -s "$f" ] || { echo "missing $f (intron SAF: run 08_context.py)" >&2; exit 1; }
done
mkdir -p "$OUT"

tail -n +2 "$ROOT/config/runs.tsv" | while IFS=$'\t' read -r sample fraction rep run srx layout; do
    [ -z "${sample:-}" ] && continue
    bam="$BAMDIR/${sample}.Aligned.sortedByCoord.out.bam"
    [ -s "$bam" ] || { echo "missing $bam -- run 03_align.sh first" >&2; exit 1; }
    var="STRAND_${sample}"
    strand="${!var:-$STRAND_DEFAULT}"

    if [ ! -s "$OUT/${sample}.genes.txt" ]; then
        echo "== $sample genes (-s $strand)"
        featureCounts -T "$THREADS" -p --countReadPairs -s "$strand" \
            -t exon -g gene_id -a "$GTF" \
            -o "$OUT/${sample}.genes.txt" "$bam"
    fi
    if [ ! -s "$OUT/${sample}.intron.txt" ]; then
        echo "== $sample introns (-s $strand)"
        featureCounts -M -F SAF -T "$THREADS" -p --countReadPairs -s "$strand" \
            -a "$INTRON_SAF" \
            -o "$OUT/${sample}.intron.txt" "$bam"
    fi
done

echo
echo "gene assignment rates:"
for f in "$OUT"/*.genes.txt.summary; do
    s=$(basename "$f" .genes.txt.summary)
    awk -F'\t' -v s="$s" '
        $1=="Assigned"{a=$2} {if(NR>1) t+=$2}
        END{printf "  %-26s %12d / %12d  (%.1f%%)\n", s, a, t, 100*a/t}' "$f"
done
