#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_plot_mvp_matrix_comparison.sh
#
# PBS submission script for plot_mvp_matrix_comparison.R.
# Produces two TIFFs comparing Betanfer vs raw MVP LD matrices:
#   scatter.betanfer_vs_mvp.tiff
#   MAE.by_genotype_bucket.tiff
#
# Submit with:   qsub run_plot_mvp_matrix_comparison.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=1
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N plot_mvp_matrix_comparison
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/plot_mvp_matrix_comparison.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/plot_mvp_matrix_comparison.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/r_plotting

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check

# ── Run plotting script ────────────────────────────────────────────────────────
echo "NOTE: Running plot_mvp_matrix_comparison.R..."
Rscript ${script_dir}/plot_mvp_matrix_comparison.R

if [ $? -ne 0 ]; then
    echo "ERROR: plot_mvp_matrix_comparison.R failed."
    exit 1
fi

echo "NOTE: Plots written successfully."
