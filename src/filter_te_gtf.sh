#!/usr/bin/env bash
# SUPERSEDED by src/rmsk_to_canonical_gtf.sh -- do not run.
#
# This filtered hg38_fixed.gtf, a 2022 download whose provenance was lost: there is
# no record of which UCSC table or RepeatMasker run it came from, so the annotation
# could not be regenerated from a named source. Stage 1 now starts from a UCSC Table
# Browser dump of the hg38 rmsk table instead. Kept only to document how the old
# hg38_L1_Alu_SVA_canonical.gtf (archived in backup_hg38_fixed_provenance/) was made.
# The selection rules below are unchanged in the replacement, just keyed off table
# columns rather than the packed attribute string.
#
# Filter the RepeatMasker hg38 GTF down to canonical L1, Alu and SVA elements.
#
# The attribute field packs subfamily/class/family into one pipe-delimited string:
#   gene_id |3|subfamily_id|L1MC5a|class_id|LINE|family_id|L1|
# Note that "subfamily_id" contains "family_id" as a substring, so an unanchored
# regex on "family_id" silently matches the subfamily slot. Split on "|" and key
# off exact field names instead.
#
# "Canonical" excludes family members that are not the element proper:
#   Alu family also contains ancestral monomers FLAM_A/FLAM_C/FRAM/FAM
#   L1  family also contains half-L1s HAL1* and a mis-filed X9_LINE
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IN="${1:-$PROJ/hg38_fixed.gtf}"
OUT="${2:-$PROJ/hg38_L1_Alu_SVA_canonical.gtf}"

awk -F'\t' '
/^#/ { print; next }
{
  # NB: "sub" is an awk builtin, so the subfamily variable is named sfam
  n = split($9, a, "|"); cls = ""; fam = ""; sfam = "";
  for (i = 1; i < n; i++) {
    if      (a[i] == "class_id")     cls  = a[i+1];
    else if (a[i] == "family_id")    fam  = a[i+1];
    else if (a[i] == "subfamily_id") sfam = a[i+1];
  }
  if ((cls == "SINE"       && fam == "Alu" && sfam ~ /^Alu/) ||
      (cls == "LINE"       && fam == "L1"  && sfam ~ /^L1/)  ||
      (cls == "Retroposon" && fam == "SVA"))
    print
}' "$IN" > "$OUT"

echo "wrote $OUT ($(wc -l < "$OUT") elements)"
