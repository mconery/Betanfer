#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_gwas_fine_mapping.sh
#
# PBS submission script for the PLATLAS GWAS meta-analysis fine-mapping pipeline.
#
# Runs SuSiEx and MultiSuSiE fine-mapping on meta-analysis loci for the EUR,
# AFR, and ALL populations defined in the master locus files produced by step 2
# (run_locus_defining.sh).  Inputs are the per-cohort summary statistics and
# per-cohort LD matrices produced by steps 1–4 of the pipeline.
#
# Pipeline steps:
#   Step 1: create_finemap_commands.py
#     Reads sig_loci.{EUR,AFR,ALL}.master.csv, the linker file, and
#     sample_sizes.tsv.  For each pending locus it validates that the
#     per-cohort sumstats (.txt.gz) and LD matrix (.vcor1) files exist and
#     emits one SuSiEx and one MultiSuSiE command, sorted by locus size
#     ascending.  Completed sentinels (output .txt files) are skipped.
#
#   Step 2: finemap_parsl_allocator.py
#     Uses a LocalProvider over the pre-allocated PBS nodes to distribute the
#     fine-mapping commands (4 workers/node).  Completed sentinels are skipped,
#     so resubmitting this job is safe.
#
# Inputs:
#   Locus lists:     /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/
#                      sig_loci.{EUR,AFR,ALL}.master.csv
#   Linker:          /flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv
#   Sumstats:        /flare/GeomicVar/multi-biobank-meta/sumstats/{COHORT}/{COHORT_POP}/
#                      *.txt.gz  (bgzip+tabix indexed)
#   LD matrices:     /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices/
#                      {TRAIT}.{ANALYSIS_POP}/chr{CHR}_{START}_{END}.{COHORT}.{COHORT_POP}.vcor1
#   Sample sizes:    /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices/
#                      sample_sizes.tsv
#
# Outputs:
#   /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/GWAS/
#     susiex/{EUR,AFR,ALL}/{TRAIT}.chr{CHR}_{START}_{END}.susiex.txt
#     susiex/{EUR,AFR,ALL}/timing/...
#     multisusie/{EUR,AFR,ALL}/{TRAIT}.chr{CHR}_{START}_{END}.multisusie.txt
#     multisusie/{EUR,AFR,ALL}/timing/...
#
# Scale note:
#   ~46,000 loci x 3 populations x 2 methods ≈ up to ~276,000 commands.
#   In practice many loci are unique to a single population.  Increase the
#   node count (select=N) to reduce the number of submissions needed.
#   The allocator skips completed sentinels so resubmission is always safe.
#
# Submit with:   qsub run_gwas_fine_mapping.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A gpu_hack
#PBS -l walltime=72:00:00
#PBS -l filesystems=home:flare
#PBS -l select=16:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q capacity
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwas_fine_mapping.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwas_fine_mapping.err

# ── Environment ────────────────────────────────────────────────────────────────
# The finemap conda env provides parsl (head-node interchange) and pandas
# (for create_finemap_commands.py).  Each command in the command file uses
# its own conda env's Python via full path:
#   SuSiEx:     /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/python
#   MultiSuSiE: /flare/GeomicVar/mconery/tools/miniconda3/envs/MultiSuSiE/bin/python
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/7_GWAS_Fine-Mapping
locus_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci
ld_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices
sumstats_base=/flare/GeomicVar/multi-biobank-meta/sumstats
linker_file=/flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv
sample_sizes=${ld_dir}/sample_sizes.tsv
output_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/GWAS
log_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs
commands_file=${script_dir}/finemap.commands.txt

# ── Create output base directory ───────────────────────────────────────────────
mkdir -p ${output_dir}

# ── Change to log directory so parsl internal logs land there too ─────────────
cd ${log_dir}

# ── Step 1: Generate fine-mapping command list ─────────────────────────────────
# Runs on every submission — skips already-completed outputs, so pending
# tasks are the only ones written to the commands file (fast on re-runs).
echo "NOTE: Generating fine-mapping command list..."
python ${script_dir}/create_finemap_commands.py \
    --linker-file   ${linker_file} \
    --locus-dir     ${locus_dir} \
    --ld-dir        ${ld_dir} \
    --sumstats-base ${sumstats_base} \
    --sample-sizes  ${sample_sizes} \
    --output-dir    ${output_dir} \
    --script-dir    ${script_dir} \
    --output-file   ${commands_file}

if [ $? -ne 0 ]; then
    echo "ERROR: create_finemap_commands.py failed.  Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${commands_file}"
echo "NOTE: Total commands: $(wc -l < ${commands_file})"

# ── Step 2: Run the parsl allocator ───────────────────────────────────────────
# The allocator uses LocalProvider (defined in parsl_config.py) to spread
# 4 workers/node across the pre-allocated nodes.
# Completed sentinels are skipped; resubmission is safe.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/finemap_parsl_allocator.py \
    --command-file ${commands_file} \
    --log-dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: finemap_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: finemap_parsl_allocator.py completed successfully."
