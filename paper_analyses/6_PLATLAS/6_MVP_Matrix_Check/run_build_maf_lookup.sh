#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_build_maf_lookup.sh
#
# PBS submission script for build_maf_lookup.py.
# Builds an rsID -> MAF lookup TSV from MVP hg19 summary statistics for
# EUR and AFR populations.  Output is used by plot_mvp_matrix_comparison.R
# to color the scatter plot by minimum minor allele frequency.
#
# Submit with:   qsub run_build_maf_lookup.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=1
#PBS -q debug
#PBS -M mconery@anl.gov
#PBS -N build_maf_lookup
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/build_maf_lookup.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/build_maf_lookup.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap_old

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check
sumstats_base=/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19
output=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545/maf_lookup.tsv

# ── Run lookup builder ────────────────────────────────────────────────────────
echo "NOTE: Building rsID -> MAF lookup..."
python3 ${script_dir}/build_maf_lookup.py \
    -r ${script_dir}/rsid_pos_lookup.json \
    -e ${sumstats_base}/EUR/MVP_R4.1000G_AGR.A1C_Max_INT.EUR.GIA.dbGaP.sumstats.hg19.txt \
    -a ${sumstats_base}/AFR/MVP_R4.1000G_AGR.A1C_Max_INT.AFR.GIA.dbGaP.sumstats.hg19.txt \
    -o ${output}

if [ $? -ne 0 ]; then
    echo "ERROR: build_maf_lookup.py failed."
    exit 1
fi

echo "NOTE: MAF lookup written to ${output}"
