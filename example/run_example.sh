#!/usr/bin/env bash
# Test a Betanfer install end to end on a small synthetic dataset.
#
# Usage:
#   bash example/run_example.sh [OUTPUT_DIR]
#
# Optional environment variables:
#   BETANFER    command used to run betanfer   (default: betanfer)
#               e.g. BETANFER="python -m betanfer" to test a source checkout
#   PYTHON      Python used for the helper scripts (default: python)
#   JOBS        parallel `betanfer index` processes (default: number of CPUs)
#   INFER_ARGS  extra arguments for `betanfer infer`, e.g. "--cpu"
#
# Takes about 1-5 minutes, mostly spent starting Python for 150 index jobs.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-betanfer_example}"
BETANFER="${BETANFER:-betanfer}"
PYTHON="${PYTHON:-python}"
JOBS="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
INFER_ARGS="${INFER_ARGS:-}"

mkdir -p "$OUT"
echo "Output directory: $OUT"

echo
echo "[1/4] Generating synthetic summary statistics..."
"$PYTHON" "$HERE/make_example_data.py" --out "$OUT"

echo
echo "[2/4] Indexing 150 trait files with 'betanfer index' ($JOBS in parallel)..."
mkdir -p "$OUT/indexed"
if ! find "$OUT/sumstats" -name '*.txt' -print0 \
    | xargs -0 -P "$JOBS" -I{} $BETANFER index --file {} --output-dir "$OUT/indexed" \
    > "$OUT/index.log" 2>&1; then
    echo "ERROR: indexing failed. See $OUT/index.log"
    exit 1
fi
n_tbi=$(find "$OUT/indexed" -name '*.gz.tbi' | wc -l | tr -d ' ')
echo "  Indexed files: $n_tbi"
if [ "$n_tbi" -ne 150 ]; then
    echo "ERROR: expected 150 indexed files. See $OUT/index.log"
    exit 1
fi

COMMON=(--dir "$OUT/indexed" --chrom 22 --start 20000000 --end 20200000
        --pval-threshold 5e-8 --min-traits 100 --plink2-format)

run_infer() {
    local label="$1" min_pair="$2"
    local prefix="$OUT/ld_minpair${min_pair}"
    # shellcheck disable=SC2086
    if ! $BETANFER infer "${COMMON[@]}" --min-pair-traits "$min_pair" \
        --output "$prefix" $INFER_ARGS > "$prefix.log" 2>&1; then
        echo "ERROR: betanfer infer failed. See $prefix.log"
        exit 1
    fi
    grep -E "Retained .* traits|Device  |set to NaN" "$prefix.log" \
        | sed -E 's/^.* - INFO - +/  /'
    echo "  Checking against the true LD matrix:"
    "$PYTHON" "$HERE/check_example.py" --inferred "$prefix" --data "$OUT" \
        --min-pair-traits "$min_pair"
}

status=0
echo
echo "[3/4] 'betanfer infer' with --min-pair-traits 2 (every pair reported)..."
run_infer "all pairs" 2 || status=1

echo
echo "[4/4] 'betanfer infer' with --min-pair-traits 100 (rare-variant pairs set to NaN)..."
run_infer "filtered" 100 || status=1

echo
if [ "$status" -eq 0 ]; then
    echo "SUCCESS: Betanfer is installed and working."
else
    echo "FAILED: see the checks above and the logs in $OUT/"
fi
exit "$status"
