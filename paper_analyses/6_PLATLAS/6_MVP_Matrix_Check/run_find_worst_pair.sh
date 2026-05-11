#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_find_worst_pair.sh
#
# PBS submission script for find_worst_pair.py.
# Streams the 50 GB mvp_matrix_comparison.tsv to find the common variant pair
# (MAF > 0.05 in the relevant population) with the greatest |betanfer_r - mvp_r|.
#
# Submit with:   qsub run_find_worst_pair.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=1
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N find_worst_pair
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/find_worst_pair.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/find_worst_pair.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap_old

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check
ritm_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545
output=${ritm_dir}/worst_pair.tsv

# ── Run ───────────────────────────────────────────────────────────────────────
echo "NOTE: Finding worst common-variant pair (MAF > 0.05)..."
python3 ${script_dir}/find_worst_pair.py \
    -c ${ritm_dir}/mvp_matrix_comparison.tsv \
    -m ${ritm_dir}/maf_lookup.tsv \
    -t 0.05 \
    -o ${output}

if [ $? -ne 0 ]; then
    echo "ERROR: find_worst_pair.py failed."
    exit 1
fi

echo "NOTE: Worst pair written to ${output}"
