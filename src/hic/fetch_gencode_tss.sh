#!/usr/bin/env bash
# GENCODE hg38 -> TSS BED, for eigenvector phasing and as the promoter anchor set.
#
#   bash src/hic/fetch_gencode_tss.sh [RELEASE] [OUTDIR]
#
# The project carries no gene annotation: hg38_L1_Alu_SVA_canonical.gtf is TE-only
# by construction and the archived hg38_fixed.gtf is 5.5M repeat_element records
# with no gene models. So TSS have to come from outside.
#
# LC_ALL=C is not optional here. Under a comma-decimal locale awk parses GTF
# coordinates and any scientific notation wrongly -- the same trap CLAUDE.md
# records for the TSV sanity checks. Set it before awk touches anything.
#
# Note also that `sub` is an awk builtin and cannot be used as a variable name;
# it fails with a syntax error and, mid-pipeline, yields an empty output file.
# Nothing below needs it, but that is why the field variables are named as they are.
set -euo pipefail
export LC_ALL=C

REL="${1:-47}"
OUT="${2:-data/hic}"
URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${REL}/gencode.v${REL}.basic.annotation.gtf.gz"
GTF="${OUT}/gencode.v${REL}.basic.annotation.gtf.gz"
BED="${OUT}/tss.bed"
BED_PC="${OUT}/tss_protein_coding.bed"

mkdir -p "$OUT"

if [[ -s "$GTF" ]]; then
  echo "[fetch_gencode_tss] reusing $GTF"
else
  echo "[fetch_gencode_tss] downloading GENCODE v${REL} ..."
  curl -fL --retry 3 --continue-at - -o "${GTF}.part" "$URL"
  mv "${GTF}.part" "$GTF"
fi

# transcript lines -> one TSS per transcript, BED 0-based half-open.
# TSS is the 5' end: start-1 on +, end-1 on -.
echo "[fetch_gencode_tss] extracting TSS ..."
zcat "$GTF" | awk -F'\t' -v OFS='\t' '
  $0 ~ /^#/          { next }
  $3 != "transcript" { next }
  {
    gname = ""; gtype = ""
    n = split($9, attrs, ";")
    for (k = 1; k <= n; k++) {
      if (attrs[k] ~ /gene_name/)      { gsub(/.*gene_name "|"$/,  "", attrs[k]); gname = attrs[k] }
      else if (attrs[k] ~ /gene_type/) { gsub(/.*gene_type "|"$/,  "", attrs[k]); gtype = attrs[k] }
    }
    pos = ($7 == "+") ? $4 - 1 : $5 - 1
    if (pos < 0) pos = 0
    print $1, pos, pos + 1, gname, gtype, $7
  }' | sort -k1,1 -k2,2n -u > "$BED"

awk -F'\t' '$5 == "protein_coding"' "$BED" > "$BED_PC"

echo "[fetch_gencode_tss] wrote:"
printf '  %-40s %s TSS\n' "$BED"    "$(wc -l < "$BED")"
printf '  %-40s %s TSS\n' "$BED_PC" "$(wc -l < "$BED_PC")"
echo
echo "Use $BED_PC for --tss unless you specifically want lncRNA/pseudogene promoters:"
echo "  the phasing track and the promoter anchor set are both cleaner without them."
