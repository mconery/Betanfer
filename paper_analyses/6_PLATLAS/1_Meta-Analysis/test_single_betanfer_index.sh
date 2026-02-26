#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# test_single_betanfer_index.sh
#
# Debug PBS job to run betanfer index on a single small file.  Output and
# error are written directly to the PBS log/err files so the full betanfer
# diagnostic output is visible without any parsl indirection.
#
# Dynamically selects the smallest unindexed .out file in the EUR results
# directory so the script always tests on a file that still needs indexing.
#
# Submit with:  qsub test_single_betanfer_index.sh
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=0:30:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q debug
#PBS -N betanfer_index_test
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/betanfer_index_test.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/betanfer_index_test.err

source $HOME/.bashrc
module load frameworks

# Pick the smallest unindexed .out file in the EUR results directory
results_base=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results
test_file=$(
  find "${results_base}/EUR" -maxdepth 1 -name "*.out" \
       ! -name "*.log.out" ! -name "*.gc.out" ! -name "*.err.out" \
  | while read f; do
      [ -f "${f}.gz.tbi" ] && continue
      stat --format="%s %n" "$f"
    done \
  | sort -n | head -1 | awk '{print $2}'
)

if [ -z "${test_file}" ]; then
    echo "ERROR: No unindexed EUR .out file found.  Nothing to test."
    exit 1
fi

echo "=== betanfer index test started: $(date) ==="
echo "Node: $(hostname)"
echo "Python: $(which python)"
echo "betanfer: $(which betanfer 2>/dev/null || echo 'NOT FOUND IN PATH')"
echo "betanfer version: $(python -m betanfer --version 2>&1)"
echo "Test file: ${test_file}"
echo "Test file exists: $(test -f ${test_file} && echo yes || echo no)"
echo "Test file size: $(du -h ${test_file} | cut -f1)"
echo "Test file header: $(head -1 ${test_file})"
echo "Output dir writable: $(test -w $(dirname ${test_file}) && echo yes || echo no)"
echo "==="

python -m betanfer index \
    --file ${test_file} \
    --variant-col rs_number \
    --variant-id-format CHR:POS:REF:ALT \
    --force

exit_code=$?
echo "=== Test finished: $(date) | exit code: ${exit_code} ==="
exit ${exit_code}
