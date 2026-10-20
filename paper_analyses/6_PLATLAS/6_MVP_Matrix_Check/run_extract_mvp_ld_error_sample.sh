#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_extract_mvp_ld_error_sample.sh
#
# PBS submission script for extract_mvp_ld_error_sample.py.
# Streams mvp_matrix_comparison.tsv (produced by compare_matrices.py) and
# writes a downsampled pair table plus exact per-population moments of the
# betanfer-vs-MVP correlation error, for the Figure 4B scatter plot.
#
# Submit with:   qsub run_extract_mvp_ld_error_sample.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=16
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N extract_mvp_ld_error_sample
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/extract_mvp_ld_error_sample.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/extract_mvp_ld_error_sample.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap_old

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check
out_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545

# ── Run ───────────────────────────────────────────────────────────────────────
echo "NOTE: Extracting the sampled betanfer-vs-MVP LD pairs..."
python3 ${script_dir}/extract_mvp_ld_error_sample.py \
    -o ${out_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: extract_mvp_ld_error_sample.py failed."
    exit 1
fi

echo "NOTE: mvp_ld_error_sample.tsv.gz and mvp_ld_error_moments.tsv written to ${out_dir}"
