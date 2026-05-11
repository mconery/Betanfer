#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_diagnose_pair.sh
#
# PBS submission script for diagnose_pair.py.
# Fill in VARIANT1, VARIANT2, LOCUS, and POPULATION from the worst_pair.tsv
# output produced by run_find_worst_pair.sh before submitting.
#
# Submit with:   qsub run_diagnose_pair.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=12
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N diagnose_pair
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/diagnose_pair.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/diagnose_pair.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap_old

# ── Pair to diagnose (fill in from worst_pair.tsv) ────────────────────────────
VARIANT1=FILL_IN
VARIANT2=FILL_IN
LOCUS=FILL_IN      # e.g. chr1_26750001_27500000
POPULATION=FILL_IN # EUR or AFR

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check
output=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545/worst_pair_diagnosis.tsv

# ── Run ───────────────────────────────────────────────────────────────────────
echo "NOTE: Diagnosing pair ${VARIANT1} vs ${VARIANT2} at ${LOCUS} (${POPULATION})..."
python3 ${script_dir}/diagnose_pair.py \
    -1 ${VARIANT1} \
    -2 ${VARIANT2} \
    -l ${LOCUS} \
    -p ${POPULATION} \
    -n 12 \
    -o ${output}

if [ $? -ne 0 ]; then
    echo "ERROR: diagnose_pair.py failed."
    exit 1
fi

echo "NOTE: Per-trait beta table written to ${output}"
