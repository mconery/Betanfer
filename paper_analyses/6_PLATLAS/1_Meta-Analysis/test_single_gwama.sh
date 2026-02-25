#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# test_single_gwama.sh
#
# Debug PBS job to test run_single_gwama.py on a single small trait/population
# combination (Phe_288_3.EUR, 1 input file ~990 MB compressed).
#
# Submit with:  qsub test_single_gwama.sh
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=0:30:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N gwama_test
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwama_test.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwama_test.err

source $HOME/.bashrc
conda activate finemap_old

script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/1_Meta-Analysis

echo "=== GWAMA single-job test started: $(date) ==="
echo "Node: $(hostname)"
echo "GWAMA binary: $(which GWAMA)"

python ${script_dir}/run_single_gwama.py \
    --trait_pop Phe_288_3.EUR

exit_code=$?
echo "=== Test finished: $(date) | exit code: ${exit_code} ==="
exit ${exit_code}
