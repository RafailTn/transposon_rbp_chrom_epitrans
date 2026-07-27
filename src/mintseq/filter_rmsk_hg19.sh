#!/usr/bin/env bash
# Filter the UCSC hg19 RepeatMasker table dump down to canonical L1, Alu and SVA
# elements, emitting BED6.
#
# This is the hg19 analogue of filter_te_gtf.sh, but the input is a different
# format and the class labels do not agree with hg38:
#
#   * Input is the rmsk MySQL table dump (tab-separated, one "#bin ..." header
#     line), not a GTF. Fields used:
#       $6 genoName  $7 genoStart  $8 genoEnd  $10 strand
#       $11 repName  $12 repClass  $13 repFamily
#     genoStart/genoEnd are already 0-based half-open, i.e. BED convention, so
#     unlike the GTF path there is no coordinate conversion anywhere.
#
#   * SVA is class "Other"/family "Other" for 3733 copies and only
#     "Retroposon"/"SVA" for 294. Keying off repClass=="Retroposon" as
#     rmsk_to_canonical_gtf.sh does would drop 93% of SVAs without erroring, so
#     this script keys off repName instead.
#
#     That is an assembly difference, NOT a bug in the hg38 script: the hg38
#     dump classes all 5974 SVAs as Retroposon/SVA, so its rule is exact there.
#     Do not "harmonise" the two by copying this rule back -- on hg38, repName
#     ~ /^SVA/ is equivalent, but the check that documents the hg38 table is the
#     class one. The rules differ because the tables differ.
#
# "Canonical" excludes family members that are not the element proper, matching
# filter_te_gtf.sh: Alu ancestral monomers FLAM_A/FLAM_C/FRAM/FAM, and the
# half-L1s HAL1* plus the mis-filed X9_LINE. The ^Alu / ^L1 repName anchors do
# that on their own.
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IN="${1:-$PROJ/data/hg19_rmsk.gtf}"
OUT="${2:-$PROJ/hg19_L1_Alu_SVA_canonical.bed}"

# LC_ALL=C: awk under a comma-decimal locale mangles numeric fields.
export LC_ALL=C

awk -F'\t' -v OFS='\t' '
NR == 1 && /^#/ { next }
{
  keep = 0
  if      ($12 == "SINE" && $13 == "Alu" && $11 ~ /^Alu/) keep = 1
  else if ($12 == "LINE" && $13 == "L1"  && $11 ~ /^L1/)  keep = 1
  else if ($11 ~ /^SVA/)                                  keep = 1
  if (keep) print $6, $7, $8, $11, ".", $10
}' "$IN" | sort -k1,1 -k2,2n > "$OUT"

echo "wrote $OUT ($(wc -l < "$OUT") elements)"
awk -F'\t' '{ f = ($4 ~ /^Alu/) ? "Alu" : ($4 ~ /^L1/) ? "L1" : "SVA"; n[f]++ }
            END { for (k in n) printf "  %-4s %8d\n", k, n[k] }' "$OUT"
