#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_mvp_fine_mapping.sh
#
# PBS submission script for the PLATLAS MVP fine-mapping pipeline.
#
# Runs RSparsePro and SuSiE fine-mapping on MVP AFR and EUR loci defined in
# the master locus files produced by step 2 (run_locus_defining.sh).
#
# Pipeline steps:
#   Step 1: create_finemap_commands.py
#     Reads sig_loci.MVP_AFR.master.csv and sig_loci.MVP_EUR.master.csv,
#     validates that the corresponding sumstats (.gz) and LD matrix (.vcor1)
#     files exist, checks for completed sentinels, and writes a sorted command
#     file (finemap_commands.txt) containing bash invocations of
#     run_rsparsepro.sh and run_susie.sh for each pending locus.
#
#   Step 2: finemap_parsl_allocator.py
#     Uses a LocalProvider over the pre-allocated PBS nodes to distribute the
#     fine-mapping commands (4 workers/node = 64 concurrent jobs with 16 nodes).
#     Already-completed sentinels are skipped, so resubmitting this job is safe.
#
# Inputs:
#   Locus lists:  /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/
#                   sig_loci.MVP_{AFR,EUR}.master.csv
#   Sumstats:     /flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19/{AFR,EUR}/
#                   MVP_R4.1000G_AGR.{TRAIT}.{POP}.GIA.dbGaP.sumstats.hg19.txt.gz
#   LD matrices:  /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices/
#                   {TRAIT}.MVP_{POP}/chr{CHR}_{START}_{END}.MVP.{POP}.vcor1
#
# Outputs:
#   /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/MVP/
#     rsparsepro/{AFR,EUR}/{TRAIT}.chr{CHR}_{START}_{END}.rsparsepro.txt
#     rsparsepro/{AFR,EUR}/timing/...
#     susie/{AFR,EUR}/{TRAIT}.chr{CHR}_{START}_{END}.susie.txt
#     susie/{AFR,EUR}/timing/...
#
# Scale note:
#   ~22,866 loci × 2 populations × 2 tools ≈ 91,000 total jobs.
#   With 16 nodes × 4 workers = 64 concurrent jobs, expect multiple walltime
#   submissions.  The allocator skips completed sentinels, so resubmission is safe.
#   Increase the node count (select=N) to reduce the number of submissions needed.
#
# Submit with:   qsub run_mvp_fine_mapping.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A gpu_hack
#PBS -l walltime=72:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q capacity
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_fine_mapping.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_fine_mapping.err

# ── Environment ────────────────────────────────────────────────────────────────
# The finemap conda env provides parsl (for the allocator/interchange on the
# head node) and all fine-mapping tools (tabix, Python/RSparsePro, Rscript/SuSiE).
# worker_init in parsl_config.py activates the same env on worker nodes so that
# parsl manager and interchange run identical Python versions.
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/5_MVP_Fine-Mapping
locus_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci
sumstats_dir=/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19
ld_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices
output_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/MVP
log_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs
commands_file=${script_dir}/finemap_commands.txt

# ── Create output base directory ───────────────────────────────────────────────
mkdir -p ${output_dir}

# ── Change to log directory so parsl internal logs land there too ─────────────
cd ${log_dir}

# ── Step 1: Generate fine-mapping command list ─────────────────────────────────
# Runs on every submission — skips already-completed outputs so only pending
# tasks are written to the commands file (fast on re-runs because completed
# loci are detected with a single os.path.exists() and skipped immediately).
echo "NOTE: Generating fine-mapping command list..."
python ${script_dir}/create_finemap_commands.py \
    --master_afr  ${locus_dir}/sig_loci.MVP_AFR.master.csv \
    --master_eur  ${locus_dir}/sig_loci.MVP_EUR.master.csv \
    --sumstats_dir ${sumstats_dir} \
    --ld_dir      ${ld_dir} \
    --output_dir  ${output_dir} \
    --script_dir  ${script_dir} \
    --output_file ${commands_file}

if [ $? -ne 0 ]; then
    echo "ERROR: create_finemap_commands.py failed.  Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${commands_file}"
echo "NOTE: Total commands: $(wc -l < ${commands_file})"

# ── Step 2: Run the parsl allocator ───────────────────────────────────────────
# The allocator uses LocalProvider (defined in parsl_config.py) to spread
# 4 workers/node across the pre-allocated nodes.
# Already-completed sentinels are skipped; resubmission is safe.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/finemap_parsl_allocator.py \
    --command_file ${commands_file} \
    --log_dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: finemap_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: finemap_parsl_allocator.py completed successfully."
