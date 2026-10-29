#!/usr/bin/env bash
# EM quantification of TE subfamilies with TEtranscripts' TEcount, sense and
# antisense, from the all-alignments BAMs of 12_align_multi.sh.
#
# TEcount only counts a TE in the orientation its GTF gives, so -- as with the
# two SAFs of make_saf.py -- antisense comes from a second TE GTF with the strand
# flipped, counted with the SAME --stranded setting. Both GTFs are built here
# from data/te_copy_key.tsv in TEtranscripts' format (gene_id = subfamily,
# transcript_id = copy id, family_id, class_id); the original
# hg38_L1_Alu_SVA_canonical.gtf is not needed.
#
# --stranded reverse corresponds to featureCounts -s 2 (the verdict from
# 04_strandedness.sh for every library). Override per sample with
# EMSTRAND_<sample>=forward|reverse.
#
# Needs the 'em' pixi env (TEcount requires python 3.10):
#   cd deps_chrrna && pixi install -e em
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAMDIR="${1:-$ROOT/work/bam_multi}"
OUT="${2:-$ROOT/work/counts_em}"
GENE_GTF="${GENE_GTF:-/mnt/raid1/rafail/aluaca/ref/gencode.v47.primary_assembly.annotation.gtf.gz}"
TECOUNT="${TECOUNT:-$ROOT/deps_chrrna/.pixi/envs/em/bin/TEcount}"
EMSTRAND_DEFAULT="${EMSTRAND_DEFAULT:-reverse}"
JOBS="${JOBS:-6}"
EMDIR="$ROOT/data/em"

[ -x "$TECOUNT" ] || { echo "missing $TECOUNT -- cd deps_chrrna && pixi install -e em" >&2; exit 1; }
mkdir -p "$OUT" "$EMDIR"
export LC_ALL=C

if [ ! -s "$EMDIR/te_antisense.gtf" ]; then
    echo "== building TEtranscripts GTFs from te_copy_key.tsv"
    awk -F'\t' -v OFS='\t' -v dir="$EMDIR" 'NR>1{
        cls = ($3=="Alu") ? "SINE" : ($3=="L1") ? "LINE" : "Retroposon"
        a = sprintf("gene_id \"%s\"; transcript_id \"%s\"; family_id \"%s\"; class_id \"%s\";", $2, $1, $3, cls)
        flip = ($7=="+") ? "-" : "+"
        print $4, "rmsk", "exon", $5, $6, ".", $7,   ".", a > dir"/te_sense.gtf"
        print $4, "rmsk", "exon", $5, $6, ".", flip, ".", a > dir"/te_antisense.gtf"
    }' "$ROOT/data/te_copy_key.tsv"
fi
if [ ! -s "$EMDIR/genes.gtf" ]; then
    zcat -f "$GENE_GTF" > "$EMDIR/genes.gtf"
fi

run_one() {
    local sample="$1" orient="$2"
    local var="EMSTRAND_${sample}"
    local strand="${!var:-$EMSTRAND_DEFAULT}"
    local prefix="$OUT/${sample}.${orient}"
    if [ -s "${prefix}.cntTable" ]; then
        echo "== $sample $orient already counted, skipping"
        return
    fi
    echo "== $sample $orient (--stranded $strand)"
    "$TECOUNT" -b "$BAMDIR/${sample}.Aligned.out.bam" \
        --GTF "$EMDIR/genes.gtf" --TE "$EMDIR/te_${orient}.gtf" \
        --mode multi --stranded "$strand" --format BAM \
        --outdir "$OUT" --project "${sample}.${orient}" > "${prefix}.log" 2>&1
}
export -f run_one
export OUT EMDIR TECOUNT EMSTRAND_DEFAULT BAMDIR
for v in $(compgen -e | grep '^EMSTRAND_'); do export "$v"; done

tail -n +2 "$ROOT/config/runs.tsv" | cut -f1 | grep -v '^$' | while read -r s; do
    [ -s "$BAMDIR/${s}.Aligned.out.bam" ] || { echo "missing $BAMDIR/${s}.Aligned.out.bam" >&2; exit 1; }
    printf '%s sense\n%s antisense\n' "$s" "$s"
done | xargs -P "$JOBS" -L 1 bash -c 'run_one "$0" "$1"'

echo
echo "TEcount tables in $OUT; next: python src/14_em_matrix.py $OUT results/em"
