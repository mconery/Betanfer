#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_gwas_rsparsepro.sh
#
# PBS submission script for the PLATLAS GWAS RSparsePro fine-mapping pipeline.
#
# Runs RSparsePro on meta-analysis loci for the EUR, AFR, and ALL populations
# defined in the master locus files from step 2 (run_locus_defining.sh).
# Unlike the SuSiEx/MultiSuSiE pipeline, RSparsePro takes the population-level
# meta-analysis summary statistics and the merged weighted-average LD matrices
# produced by step 4 (run_ld_merge.sh).
#
# Pipeline steps:
#   Step 1: create_rsparsepro_commands.py
#     Reads sig_loci.{EUR,AFR,ALL}.master.csv and, for each pending locus,
#     validates that the meta-analysis sumstats (.out.gz) and merged LD matrix
#     (.merged.merged.vcor1) files exist and emits one RSparsePro command,
#     sorted by locus size ascending.  Completed sentinels (.rsparsepro.txt
#     output files) are skipped.
#
#   Step 2: rsparsepro_parsl_allocator.py
#     Uses a LocalProvider over the pre-allocated PBS nodes to distribute the
#     RSparsePro commands (4 workers/node).  Completed sentinels are skipped,
#     so resubmitting this job is safe.
#
# Inputs:
#   Locus lists:   /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/
#                    sig_loci.{EUR,AFR,ALL}.master.csv
#   Meta-analysis: /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results/
#                    {EUR,AFR,ALL}/{TRAIT}.out.gz  (bgzip+tabix indexed)
#   LD matrices:   /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices/
#                    {TRAIT}.{ANALYSIS_POP}/chr{CHR}_{START}_{END}.merged.merged.vcor1
#
# Outputs:
#   /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/GWAS/
#     rsparsepro/{EUR,AFR,ALL}/{TRAIT}.chr{CHR}_{START}_{END}.rsparsepro.txt
#     rsparsepro/{EUR,AFR,ALL}/timing/...
#
# Submit with:   qsub run_gwas_rsparsepro.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=24:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q capacity
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwas_rsparsepro.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwas_rsparsepro.err

# ── Environment ────────────────────────────────────────────────────────────────
# The finemap conda env provides parsl (head-node interchange), pandas, numpy,
# and scipy.  Each run_rsparsepro.py command also uses the finemap Python path.
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/7_GWAS_Fine-Mapping
locus_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci
meta_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results
ld_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices
output_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/GWAS
log_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs
commands_file=${script_dir}/rsparsepro.commands.txt

# ── Create output base directory ───────────────────────────────────────────────
mkdir -p ${output_dir}

# ── Change to log directory so parsl internal logs land there too ─────────────
cd ${log_dir}

# ── Step 1: Generate RSparsePro command list ───────────────────────────────────
# Runs on every submission — skips already-completed outputs, so pending
# tasks are the only ones written to the commands file (fast on re-runs).
echo "NOTE: Generating RSparsePro command list..."
python ${script_dir}/create_rsparsepro_commands.py \
    --locus-dir   ${locus_dir} \
    --meta-dir    ${meta_dir} \
    --ld-dir      ${ld_dir} \
    --output-dir  ${output_dir} \
    --script-dir  ${script_dir} \
    --output-file ${commands_file}

if [ $? -ne 0 ]; then
    echo "ERROR: create_rsparsepro_commands.py failed.  Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${commands_file}"
echo "NOTE: Total commands: $(wc -l < ${commands_file})"

# ── Step 2: Run the parsl allocator ───────────────────────────────────────────
# The allocator uses LocalProvider (defined in parsl_config.py) to spread
# 4 workers/node across the pre-allocated nodes.
# Completed sentinels are skipped; resubmission is safe.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/rsparsepro_parsl_allocator.py \
    --command-file ${commands_file} \
    --log-dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: rsparsepro_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: rsparsepro_parsl_allocator.py completed successfully."
