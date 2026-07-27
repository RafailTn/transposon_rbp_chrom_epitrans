#!/usr/bin/env bash
# Adapter-, quality- and length-trim the FASTQs before alignment.
#
# THIS STEP IS NOT OPTIONAL, for a reason specific to this pipeline. 03_align.sh
# runs STAR with --alignEndsType EndToEnd, which is the authors' setting and
# forbids soft-clipping. An untrimmed adapter therefore cannot be clipped off the
# end of a read -- the read simply fails to align, and does so silently. The loss
# is worst for short inserts, which for a fragmented chromatin RNA library is a
# large fraction. Skipping this does not error; it quietly deletes data.
#
# The authors used Cutadapt v1.9.1. This uses whatever modern Cutadapt the env
# provides (see deps_chrrna/pixi.toml): 1.9.1 predates Python 3 support and
# cannot be installed alongside the rest of the toolchain. Adapter removal is
# not where the versions differ meaningfully -- the alignment parameters are,
# and STAR is pinned to the authors' exact 2.7.2b.
#
# Adapters are the standard Illumina TruSeq pair, matching the NEBNext Ultra II
# Directional and Takara kits used for these libraries.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FQ="${1:-$PWD/fastq}"
OUT="${2:-$PWD/fastq_trimmed}"
THREADS="${THREADS:-8}"

A1="${ADAPTER_R1:-AGATCGGAAGAGCACACGTCTGAACTCCAGTCA}"
A2="${ADAPTER_R2:-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT}"
MINLEN="${MINLEN:-25}"
QUAL="${QUAL:-20}"

command -v cutadapt >/dev/null || { echo "missing cutadapt" >&2; exit 1; }
mkdir -p "$OUT"

tail -n +2 "$HERE/runs.tsv" | while IFS=$'\t' read -r sample fraction rep run srx layout; do
    [ -z "${sample:-}" ] && continue
    o1="$OUT/${sample}_1.fastq.gz"
    o2="$OUT/${sample}_2.fastq.gz"
    if [ -s "$o2" ]; then
        echo "== $sample already trimmed, skipping"
        continue
    fi
    r1="$FQ/${sample}_1.fastq.gz"
    r2="$FQ/${sample}_2.fastq.gz"
    for f in "$r1" "$r2"; do
        [ -s "$f" ] || { echo "missing $f -- run 00_fetch.sh first" >&2; exit 1; }
    done
    echo "== trimming $sample"
    cutadapt -j "$THREADS" \
        -a "$A1" -A "$A2" \
        -q "$QUAL" -m "$MINLEN" \
        --pair-filter=any \
        -o "$o1" -p "$o2" \
        "$r1" "$r2" > "$OUT/${sample}.cutadapt.log"
    grep -E 'Total read pairs|Pairs written|Reads with adapters' \
        "$OUT/${sample}.cutadapt.log" | sed 's/^/    /'
done

echo
echo "trimmed reads in $OUT"
echo "A 'Pairs written' fraction below ~80% is worth investigating before aligning."
