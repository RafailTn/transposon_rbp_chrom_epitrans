#!/usr/bin/env bash
# Build hg38_L1_Alu_SVA_canonical.gtf from the UCSC RepeatMasker table dump.
#
# This REPLACES filter_te_gtf.sh as stage 1. That script read hg38_fixed.gtf, a
# 2022 download whose provenance was lost; this one starts from a re-downloadable
# UCSC Table Browser dump of the hg38 `rmsk` table, so the annotation is
# reproducible from a named source.
#
# Input is the raw table, NOT a GTF despite the .gtf.gz filename UCSC hands out:
#   1 bin  2 swScore  3 milliDiv  4 milliDel  5 milliIns  6 genoName  7 genoStart
#   8 genoEnd  9 genoLeft  10 strand  11 repName  12 repClass  13 repFamily
#   14 repStart  15 repEnd  16 repLeft  17 id
# with a leading '#'-prefixed header line. Coordinates are BED 0-based half-open;
# GTF is 1-based inclusive, so start becomes genoStart+1 and end stays genoEnd.
# (build_te_index.py converts straight back -- the round trip is deliberate, it
# keeps the GTF the single interchange format every downstream script expects.)
#
# Selection is identical to filter_te_gtf.sh, just keyed off table columns rather
# than the packed attribute string. "Canonical" means canonical *subfamilies*, not
# canonical chromosomes -- there is deliberately no chromosome filter here:
#   Alu family also contains ancestral monomers FLAM_A/FLAM_C/FRAM/FAM
#   L1  family also contains half-L1s HAL1* and a mis-filed X9_LINE
#
# Output reproduces the 9-column layout the old file had, so nothing downstream
# changes:
#   chrom  RepeatMasker  repeat_element  start  end  swScore  strand  0  attrs
# gene_id is the genome-wide positional rank over ALL repeats (not just the
# filtered ones), matching the old file's numbering scheme. Nothing reads it --
# build_te_index.py keys off the exact field names class_id/family_id/subfamily_id
# and never touches a[2] -- it is kept only so the two files stay diffable.
set -euo pipefail
export LC_ALL=C          # comma-decimal locales break numeric sort and awk here

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN="${1:-$PROJ/data/hg38_rmsk.gtf.gz}"
OUT="${2:-$PROJ/hg38_L1_Alu_SVA_canonical.gtf}"

case "$IN" in
    *.gz) CAT=zcat ;;
    *)    CAT=cat  ;;
esac

# Sort every repeat by position first so the gene_id counter is a positional rank,
# then filter. Numbering before filtering is what makes the ids comparable to the
# old file's.
"$CAT" "$IN" \
  | awk -F'\t' 'NR > 1 || $1 !~ /^#/' \
  | sort -k6,6 -k7,7n -k8,8n \
  | awk -F'\t' -v OFS='\t' '
{
  gid++
  if (($12 == "SINE"       && $13 == "Alu" && $11 ~ /^Alu/) ||
      ($12 == "LINE"       && $13 == "L1"  && $11 ~ /^L1/)  ||
      ($12 == "Retroposon" && $13 == "SVA"))
    print $6, "RepeatMasker", "repeat_element", $7 + 1, $8, $2, $10, 0,
          "gene_id |" gid "|subfamily_id|" $11 "|class_id|" $12 "|family_id|" $13 "|"
}' > "$OUT"

echo "wrote $OUT ($(wc -l < "$OUT") elements)"
