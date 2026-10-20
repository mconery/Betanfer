#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_collate_mvp_traits_per_locus.sh
#
# PBS submission script for collate_mvp_traits_per_locus.py.
# For each locus in sampled_loci_750kb.csv and each population, counts how
# many MVP traits have a genome-wide-significant hit in the locus (and are
# therefore excluded from betanfer's LD inference there). Feeds the
# per-trait-count coloring in the Figure 4C plot.
#
# Submit with:   qsub run_collate_mvp_traits_per_locus.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=48
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N collate_mvp_traits_per_locus
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/collate_mvp_traits_per_locus.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/collate_mvp_traits_per_locus.err

# ── Environment ────────────────────────────────────────────────────────────────
# needs pysam, which is not installed in the finemap / finemap_old environments
source $HOME/.bashrc
conda activate /home/jianhuawang/miniforge3/envs/betainfer

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check
out_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545

# ── Run ───────────────────────────────────────────────────────────────────────
echo "NOTE: Counting genome-wide-significant traits per locus..."
python3 ${script_dir}/collate_mvp_traits_per_locus.py \
    -n 48 \
    -o ${out_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: collate_mvp_traits_per_locus.py failed."
    exit 1
fi

echo "NOTE: mvp_traits_per_locus.tsv written to ${out_dir}"
