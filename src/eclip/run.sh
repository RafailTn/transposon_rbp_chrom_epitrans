#!/usr/bin/env bash
# eCLIP run: observed crosslink peaks, K562 + HepG2, from the two
# eCLIP_Peaks_K562_HepG2_* folders (BED6, hg38, 2 replicates per RBP).
#
# The dataset-specific work is entirely in prep_replicates.py: intersect
# replicates, and rename into the <RBP>_<CELL>.bed convention that every script's
# basename[:-4].rsplit("_", 1) parser requires. Everything after that is the
# generic TE-enrichment pipeline and makes no assumption about where the peaks
# came from -- except the two flagged below, which are sized to THIS peak set.
#
# Skip the prep step (it is deterministic and takes ~1 min) with --no-prep.
#
# Usage: bash src/eclip/run.sh [--no-prep] [--pmax N]
set -euo pipefail

PROJ="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY="$PROJ/deps/.pixi/envs/default/bin/python"
if [ ! -x "$PY" ]; then
    echo "missing interpreter: $PY -- run 'cd $PROJ/deps && pixi install'" >&2
    exit 1
fi
DATA="$PROJ/data/eclip"
OUT="$PROJ/results/eclip"
export LC_ALL=C          # comma-decimal locales make awk misparse these TSVs

PREP=1
if [ "${1:-}" = "--no-prep" ]; then PREP=0; shift; fi

# The BH floor rule (1/(P_MAX+1) * candidates-in-a-cell < 0.05) scales with the
# candidate count, which is a property of the peak set, not of the method. This one
# screens 1,037 candidates in HepG2, so it needs P_MAX > 20,740 -- far above
# enrich_permutation.py's own default. Defaulting it here rather than documenting a
# flag: getting this wrong does not error, it silently leaves 56% of hits tied on
# the p-floor and therefore unrankable. Re-derive it if the peak set changes.
PERM_ARGS=("$@")
case " $* " in *" --pmax "*) ;; *) PERM_ARGS+=(--pmax 25000) ;; esac

mkdir -p "$OUT"
# Stage 1-2 (rmsk_to_canonical_gtf.sh, build_te_index.py) are cached artefacts:
# the TE annotation does not depend on the peaks, and rebuilding it per run would
# only risk the cached index and the GTF drifting apart.
[ "$PREP" = 1 ] && "$PY" "$PROJ/src/eclip/prep_replicates.py" "$PROJ/data" "$DATA"

"$PY" "$PROJ/src/enrich_stranded.py"        "$DATA" "$OUT"
"$PY" "$PROJ/src/enrich_permutation.py"     "$DATA" "$OUT" "${PERM_ARGS[@]}"
"$PY" "$PROJ/src/report_loci_per_cell.py"   "$DATA" "$OUT"
"$PY" "$PROJ/src/plot_enrichment.py"        "$OUT"
"$PY" "$PROJ/src/plot_across_cell_lines.py" "$OUT"
echo "eCLIP run complete -> $OUT"
