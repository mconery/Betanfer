#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_plot_mvp_ld_error_scatter.sh
#
# PBS submission script for plot_mvp_ld_error_scatter.R.
# Produces mvp_ld_error_scatter.{png,tiff}: betanfer beta correlation vs
# real MVP genotype correlation, one point per sampled variant pair
# (Figure 4B).
#
# Submit with:   qsub run_plot_mvp_ld_error_scatter.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=1
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N plot_mvp_ld_error_scatter
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/plot_mvp_ld_error_scatter.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/plot_mvp_ld_error_scatter.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/r_plotting

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check

# ── Run plotting script ────────────────────────────────────────────────────────
echo "NOTE: Running plot_mvp_ld_error_scatter.R..."
Rscript ${script_dir}/plot_mvp_ld_error_scatter.R tiff

if [ $? -ne 0 ]; then
    echo "ERROR: plot_mvp_ld_error_scatter.R failed."
    exit 1
fi

echo "NOTE: Plot written successfully."
