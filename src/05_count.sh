#!/usr/bin/env bash
# Count reads over TE copies, sense and antisense, with featureCounts.
#
# Matches the authors' quantification flags (-M -F SAF -p, plus their -s) so the
# Alu/SVA numbers are produced the same way as the L1 numbers they published --
# the only difference is that the SAF has Alu and SVA in it.
#
# --countReadPairs is required on top of -p. Since subread 2.0.2, -p alone only
# declares the input paired and still counts READS, so every fragment is counted
# twice and DESeq2 treats two mates of one fragment as independent evidence. (The
# first run of this pipeline had exactly that bug.) The authors' deposited tables
# are htseq-count output whose totals equal read counts, not fragments -- so this
# is a deliberate departure from them, and our totals are ~half of theirs.
#
# -M counts multi-mapping reads. Note this is nearly a no-op given STAR was told
# --outSAMmultNmax 1: only one alignment per read reaches the BAM. It is kept for
# fidelity to the published pipeline, and it still matters for the handful of
# reads STAR reports more than once.
#
# Orientation comes from counting two SAF files that differ only in the Strand
# column, with the SAME -s value, rather than from flipping -s. See make_saf.py.
#
# Per-sample strand override (from 04_strandedness.sh), e.g.:
#   STRAND_DEFAULT=2 STRAND_K562_cyto_Takara_rep1=1 bash 05_count.sh
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAMDIR="${1:-$ROOT/work/bam}"
SAFDIR="${2:-$ROOT/data}"
OUT="${3:-$ROOT/work/counts}"
THREADS="${THREADS:-16}"
STRAND_DEFAULT="${STRAND_DEFAULT:-2}"

command -v featureCounts >/dev/null || { echo "missing featureCounts (subread)" >&2; exit 1; }
for f in "$SAFDIR/te_sense.saf" "$SAFDIR/te_antisense.saf"; do
    [ -s "$f" ] || { echo "missing $f -- run make_saf.py first" >&2; exit 1; }
done
mkdir -p "$OUT"

tail -n +2 "$ROOT/config/runs.tsv" | while IFS=$'\t' read -r sample fraction rep run srx layout; do
    [ -z "${sample:-}" ] && continue
    bam="$BAMDIR/${sample}.Aligned.sortedByCoord.out.bam"
    [ -s "$bam" ] || { echo "missing $bam -- run 03_align.sh first" >&2; exit 1; }

    var="STRAND_${sample}"
    strand="${!var:-$STRAND_DEFAULT}"

    for orient in sense antisense; do
        dest="$OUT/${sample}.${orient}.txt"
        if [ -s "$dest" ]; then
            echo "== $sample $orient already counted, skipping"
            continue
        fi
        echo "== $sample $orient (-s $strand)"
        featureCounts -M -F SAF -T "$THREADS" -p --countReadPairs -s "$strand" \
            -a "$SAFDIR/te_${orient}.saf" \
            -o "$dest" "$bam"
    done
done

echo
echo "assignment rates:"
for f in "$OUT"/*.summary; do
    s=$(basename "$f" .txt.summary)
    awk -F'\t' -v s="$s" '
        $1=="Assigned"{a=$2} {if(NR>1) t+=$2}
        END{printf "  %-38s %12d / %12d  (%.1f%%)\n", s, a, t, 100*a/t}' "$f"
done
echo
echo "Next:"
echo "  python src/06_enrich_chrrna.py"
