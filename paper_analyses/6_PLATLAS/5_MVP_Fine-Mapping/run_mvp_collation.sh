#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_mvp_collation.sh
#
# PBS submission script for collating original vs. new MVP fine-mapping results.
#
# Compares variants in credible sets from the original analysis (Data_S1.xlsx)
# against new SuSiE and RSparsePro results produced with updated LD matrices.
# Only AFR and EUR populations are compared (AMR and EAS were not re-run).
#
# Outputs (written to mvp_results/):
#   mvp_collation.signal_recovery.tsv  -- signal-level recovery check
#   mvp_collation.pip_comparison.tsv   -- variant-level PIP comparison
#
# Submit with:   qsub run_mvp_collation.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A gpu_hack
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q debug
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_collation.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_collation.err

# ── Environment ────────────────────────────────────────────────────────────────
# finemap_old provides pandas and openpyxl for reading the Excel input.
# If openpyxl is missing: conda install -c conda-forge openpyxl
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap_old

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/5_MVP_Fine-Mapping
results_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/MVP
excel_file=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/Data_S1.xlsx
output_prefix=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/mvp_collation

# ── Run collation ──────────────────────────────────────────────────────────────
echo "NOTE: Starting MVP fine-mapping collation..."
python ${script_dir}/collate_mvp_finemap_results.py \
    --excel-file    ${excel_file} \
    --results-dir   ${results_dir} \
    --output-prefix ${output_prefix}

if [ $? -ne 0 ]; then
    echo "ERROR: collate_mvp_finemap_results.py failed."
    exit 1
fi

echo "NOTE: Collation complete. Outputs written to:"
echo "  ${output_prefix}.signal_recovery.tsv"
echo "  ${output_prefix}.pip_comparison.tsv"
