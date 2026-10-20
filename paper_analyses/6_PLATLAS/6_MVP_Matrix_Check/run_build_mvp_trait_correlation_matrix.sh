#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_build_mvp_trait_correlation_matrix.sh
#
# PBS submission script for build_mvp_trait_correlation_matrix.py.
# Pulls a genome-wide random sample of common-variant BETA/SE values from
# every MVP EUR trait, producing the matrix that
# collate_effective_trait_counts.py and decompose_beta_vs_z_ld_error.py use
# for the Supplementary Note on beta- vs z-based LD inference.
# A comparable earlier run over 1,443 traits with 48 workers took about
# 6 minutes; the 1-hour debug-queue walltime below leaves ample margin.
#
# Submit with:   qsub run_build_mvp_trait_correlation_matrix.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=48
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N build_mvp_trait_correlation_matrix
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/build_mvp_trait_correlation_matrix.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/build_mvp_trait_correlation_matrix.err

# ── Environment ────────────────────────────────────────────────────────────────
# needs pysam, which is not installed in the finemap / finemap_old environments
source $HOME/.bashrc
conda activate /home/jianhuawang/miniforge3/envs/betainfer

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check
out_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/effective_traits

# ── Run ───────────────────────────────────────────────────────────────────────
echo "NOTE: Building the MVP trait correlation matrix..."
python3 ${script_dir}/build_mvp_trait_correlation_matrix.py \
    -n 48 \
    -o ${out_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: build_mvp_trait_correlation_matrix.py failed."
    exit 1
fi

echo "NOTE: mvp_trait_correlation_matrix.npz written to ${out_dir}"
