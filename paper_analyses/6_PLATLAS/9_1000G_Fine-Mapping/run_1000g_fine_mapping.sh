#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_1000g_fine_mapping.sh
#
# PBS submission script for the PLATLAS 1000G fine-mapping pipeline.
#
# Runs RSparsePro and SuSiE fine-mapping on MVP AFR and EUR loci defined in
# the master locus files, using 1000G LD matrices instead of the Betanfer-
# derived matrices.  The 1000G LD matrices must have been produced by step 8
# (run_1000g_matrix_making.sh) before submitting this job.
#
# Pipeline steps:
#   Step 1: create_1000g_finemap_commands.py
#     Reads sig_loci.MVP_AFR.master.csv and sig_loci.MVP_EUR.master.csv,
#     validates that the corresponding sumstats (.gz) and 1000G LD matrix
#     (.unphased.vcor1 + .unphased.vcor1.vars) files exist, checks for
#     completed sentinels, and writes a sorted command file containing bash
#     invocations of run_rsparsepro_1000g.sh and run_susie_1000g.sh for each
#     pending locus.
#
#   Step 2: finemap_1000g_parsl_allocator.py
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
#                   {TRAIT}.MVP_{POP}/chr{CHR}_{START}_{END}.1000G.{POP}.unphased.vcor1
#
# Outputs:
#   /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/MVP_1000G/
#     rsparsepro/{AFR,EUR}/{TRAIT}.chr{CHR}_{START}_{END}.rsparsepro.txt
#     rsparsepro/{AFR,EUR}/timing/...
#     susie/{AFR,EUR}/{TRAIT}.chr{CHR}_{START}_{END}.susie.txt
#     susie/{AFR,EUR}/timing/...
#
# Scale note:
#   ~22,866 loci x 2 populations x 2 tools = ~91,000 total jobs.
#   With 16 nodes x 4 workers = 64 concurrent jobs, expect multiple walltime
#   submissions.  The allocator skips completed sentinels, so resubmission is safe.
#   Increase the node count (select=N) to reduce the number of submissions needed.
#
# Submit with:   qsub run_1000g_fine_mapping.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=5:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q debug-scaling
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/1000g_fine_mapping.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/1000g_fine_mapping.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/9_1000G_Fine-Mapping
locus_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci
sumstats_dir=/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19
ld_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices
output_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/MVP_1000G
log_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs
commands_file=${script_dir}/finemap_1000g_commands.txt

# ── Create output base directory ───────────────────────────────────────────────
mkdir -p ${output_dir}

# ── Change to log directory so parsl internal logs land there too ─────────────
cd ${log_dir}

# ── Step 1: Generate fine-mapping command list ─────────────────────────────────
echo "NOTE: Generating 1000G fine-mapping command list..."
python ${script_dir}/create_1000g_finemap_commands.py \
    --master_afr  ${locus_dir}/sig_loci.MVP_AFR.master.csv \
    --master_eur  ${locus_dir}/sig_loci.MVP_EUR.master.csv \
    --sumstats_dir ${sumstats_dir} \
    --ld_dir      ${ld_dir} \
    --output_dir  ${output_dir} \
    --script_dir  ${script_dir} \
    --output_file ${commands_file}

if [ $? -ne 0 ]; then
    echo "ERROR: create_1000g_finemap_commands.py failed.  Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${commands_file}"
echo "NOTE: Total commands: $(wc -l < ${commands_file})"

# ── Step 2: Run the parsl allocator ───────────────────────────────────────────
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/finemap_1000g_parsl_allocator.py \
    --command_file ${commands_file} \
    --log_dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: finemap_1000g_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: finemap_1000g_parsl_allocator.py completed successfully."
