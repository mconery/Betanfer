#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_effective_trait_pipeline.sh
#
# PBS submission script chaining the three Supplementary Note steps that run
# on the already-built mvp_trait_correlation_matrix.npz:
#   collate_effective_trait_counts.py  -> mvp_teff_summary.tsv, mvp_teff_vs_t.tsv,
#                                          mvp_trait_scale.tsv
#   decompose_beta_vs_z_ld_error.py    -> mvp_teff_decomposition.tsv
#   plot_mvp_teff_vs_t.py              -> mvp_teff_vs_t.{pdf,png}
# All three are cheap linear-algebra passes over a 1,443-trait matrix, so
# they are chained in one small job rather than three separate submissions.
#
# Run build_mvp_trait_correlation_matrix.py first if
# mvp_trait_correlation_matrix.npz does not already exist.
#
# Submit with:   qsub run_effective_trait_pipeline.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=4
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N effective_trait_pipeline
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/effective_trait_pipeline.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/effective_trait_pipeline.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /home/jianhuawang/miniforge3/envs/betainfer

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check
out_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/effective_traits

# ── Run ───────────────────────────────────────────────────────────────────────
echo "NOTE: Collating effective trait counts..."
python3 ${script_dir}/collate_effective_trait_counts.py -o ${out_dir}
if [ $? -ne 0 ]; then
    echo "ERROR: collate_effective_trait_counts.py failed."
    exit 1
fi

echo "NOTE: Decomposing the beta-vs-z LD error..."
python3 ${script_dir}/decompose_beta_vs_z_ld_error.py -o ${out_dir}
if [ $? -ne 0 ]; then
    echo "ERROR: decompose_beta_vs_z_ld_error.py failed."
    exit 1
fi

echo "NOTE: Plotting T_eff vs number of traits..."
python3 ${script_dir}/plot_mvp_teff_vs_t.py -i ${out_dir}
if [ $? -ne 0 ]; then
    echo "ERROR: plot_mvp_teff_vs_t.py failed."
    exit 1
fi

echo "NOTE: Effective-trait pipeline outputs written to ${out_dir}"
