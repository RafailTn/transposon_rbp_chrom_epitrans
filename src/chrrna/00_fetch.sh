#!/usr/bin/env bash
# Fetch the 6 K562 runs from SRA as paired gzipped FASTQ.
#
# GSE296212 deposits only .txt.gz count tables, and those tables contain L1 loci
# and gene symbols but ZERO Alu and ZERO SVA rows -- the authors ran featureCounts
# against an L1-only SAF. There is no bigWig for any of the 80 samples either, so
# there is no coverage shortcut. Re-quantifying from raw reads is the only way to
# get Alu/SVA out of this dataset.
#
# ~17 GB of SRA, expanding to roughly 90-110 GB of FASTQ. Run on the server.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-$PWD/fastq}"
THREADS="${THREADS:-8}"

for t in prefetch fasterq-dump; do
    command -v "$t" >/dev/null || { echo "missing $t (sra-tools)" >&2; exit 1; }
done

mkdir -p "$OUT"
cd "$OUT"

tail -n +2 "$HERE/runs.tsv" | while IFS=$'\t' read -r sample fraction rep run srx layout; do
    [ -z "${run:-}" ] && continue
    if [ -s "${sample}_2.fastq.gz" ]; then
        echo "== $sample ($run) already present, skipping"
        continue
    fi
    echo "== $sample ($run)"
    prefetch --max-size 100G "$run"
    # --split-3 so a run that is unexpectedly single-end fails loudly here rather
    # than silently producing a half-empty pair downstream.
    fasterq-dump --split-3 --threads "$THREADS" --progress "$run"
    if [ ! -s "${run}_2.fastq" ]; then
        echo "ERROR: $run did not yield a read-2 file; runs.tsv says PAIRED" >&2
        exit 1
    fi
    mv "${run}_1.fastq" "${sample}_1.fastq"
    mv "${run}_2.fastq" "${sample}_2.fastq"
    rm -f "${run}.fastq"
    pigz -p "$THREADS" "${sample}_1.fastq" "${sample}_2.fastq" 2>/dev/null \
        || gzip "${sample}_1.fastq" "${sample}_2.fastq"
    rm -rf "$run"
done

echo
echo "fetched into $OUT:"
ls -la "$OUT"/*.fastq.gz
