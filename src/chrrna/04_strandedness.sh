#!/usr/bin/env bash
# Determine each library's strand convention empirically. Do not skip this.
#
# The series !Sample_data_processing quotes one featureCounts setting, "-s 2",
# for the whole SuperSeries. But the chrRNA samples were built with the NEB
# Ultra II Directional kit and the *_Takara_* samples with a different Takara
# kit, and a single boilerplate processing note is not evidence that both
# libraries share a convention. Getting this backwards silently swaps sense and
# antisense -- the same failure mode as assuming .pos/.neg in the mintseq
# bigWigs, which turned out fine only because it was checked.
#
# Method: count four housekeeping genes of known orientation under -s 1 and
# -s 2. The correct setting is the one that puts the reads on the annotated
# strand. Two genes are '+' and two are '-' so a result that merely reflects a
# genome-wide strand imbalance cannot masquerade as a clean answer.
#
# Expect a lopsided ratio (>5x). Anything near 1:1 means the library is
# effectively unstranded, and the sense/antisense split cannot be made at all.
set -euo pipefail

BAMDIR="${1:-$PWD/bam}"
WORK="${2:-$PWD/strand_check}"
THREADS="${THREADS:-8}"

command -v featureCounts >/dev/null || { echo "missing featureCounts (subread)" >&2; exit 1; }
mkdir -p "$WORK"

# hg38, GENCODE spans. Approximate ends are fine: this is a ratio between two
# strand settings over the same intervals, so a few hundred bp either way cancels.
cat > "$WORK/housekeeping.saf" <<'EOF'
GeneID	Chr	Start	End	Strand
GAPDH_plus	chr12	6534517	6538371	+
RPL13A_plus	chr19	49487608	49492308	+
ACTB_minus	chr7	5527151	5530601	-
EEF1A1_minus	chr6	73515750	73523797	-
EOF

printf "%-26s %10s %10s   %s\n" "sample" "s=1" "s=2" "verdict"
for bam in "$BAMDIR"/*.Aligned.sortedByCoord.out.bam; do
    s=$(basename "$bam" .Aligned.sortedByCoord.out.bam)
    for strand in 1 2; do
        featureCounts -M -F SAF -T "$THREADS" -p -s "$strand" \
            -a "$WORK/housekeeping.saf" \
            -o "$WORK/${s}.s${strand}.txt" "$bam" >/dev/null 2>&1
    done
    c1=$(awk 'NR>2{n+=$NF} END{print n+0}' "$WORK/${s}.s1.txt")
    c2=$(awk 'NR>2{n+=$NF} END{print n+0}' "$WORK/${s}.s2.txt")
    verdict=$(awk -v a="$c1" -v b="$c2" 'BEGIN{
        if (a==0 && b==0) { print "NO READS - check the BAM"; exit }
        hi = (a>b)? a : b; lo = (a>b)? b : a;
        r = (lo>0)? hi/lo : hi;
        if (r < 5) { printf "UNSTRANDED? ratio %.1fx - do NOT split by orientation", r }
        else if (a>b) { printf "use -s 1  (%.1fx)", r }
        else          { printf "use -s 2  (%.1fx)", r }
    }')
    printf "%-26s %10s %10s   %s\n" "$s" "$c1" "$c2" "$verdict"
done

echo
echo "Feed the winning value per sample into 05_count.sh via STRAND_<sample>."
echo "If the chrRNA and Takara groups disagree, that is a real kit difference,"
echo "not an error -- set them separately rather than forcing one value."
