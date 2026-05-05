#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_mvp_locus_plots.sh
#
# PBS submission script for generating MVP locus comparison plots.
#
# Randomly selects 100 loci per population (AFR, EUR) and produces a
# 3-panel plot per locus comparing original SuSiE credible sets against
# new SuSiE and RSparsePro results.
#
# Submit with:   qsub run_mvp_locus_plots.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=10:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q debug-scaling
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_locus_plots.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_locus_plots.err

# ── Environment ───────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap_old

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/5_MVP_Fine-Mapping
output_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/locus_plots

mkdir -p ${output_dir}/AFR ${output_dir}/EUR

# ── Run plotting ──────────────────────────────────────────────────────────────
echo "NOTE: Starting MVP locus plot generation..."
python ${script_dir}/plot_mvp_loci.py \
    --output-dir ${output_dir} \
    --workers    32

if [ $? -ne 0 ]; then
    echo "ERROR: plot_mvp_loci.py failed."
    exit 1
fi

echo "NOTE: Plots written to ${output_dir}"
