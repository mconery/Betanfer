#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_intersection_rsparsepro.sh
#
# PBS submission script for the trait-aware intersection LD merge and
# RSparsePro fine-mapping pipeline.
#
# Replaces, for RSparsePro only, the union-merged LD from step 4
# (4_LD_Matrix_Merging), which fills missing cells with NaN and averages
# only over the cohorts with coverage. This pipeline instead restricts each
# locus to the SNPs present in every contributing per-cohort LD matrix and
# in every stratum's and the meta-analysis's sumstats, and takes the
# sample-size-weighted LD average on that strict intersection.
#
# Pipeline steps:
#   Step 1: create_intersection_commands.py
#     Reads sig_loci.{EUR,AFR,ALL}.master.csv, the PLATLAS linker, and the
#     per-cohort sample sizes, and for each locus with complete inputs emits
#     one merge command and one RSparsePro command, sorted by total LD input
#     size descending.
#
#   Step 2: intersection_parsl_allocator.py on the merge commands
#     Uses a LocalProvider over the pre-allocated PBS nodes to run
#     merge_intersection_ld.py for every locus (4 workers/node). Completed
#     sentinels are skipped, so resubmitting this job is safe.
#
#   Step 3: intersection_parsl_allocator.py on the RSparsePro commands
#     Runs run_rsparsepro_intersection.py on the merged LD from step 2.
#     Must run after step 2 completes, since it reads step 2's output.
#
# Inputs:
#   Linker:        /flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv
#   Locus lists:   /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/
#                    sig_loci.{EUR,AFR,ALL}.master.csv
#   Sample sizes:  /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices/sample_sizes.tsv
#   Per-cohort LD: /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices/
#                    {TRAIT}.{POP}/{LOCUS}.{COHORT}.{COHORT_POP}.vcor1
#   Meta-analysis: /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results/
#                    {POP}/{TRAIT}.out.gz
#
# Outputs:
#   /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/intersection_finemap/
#     ld/{TRAIT}.{POP}/{LOCUS}.trait-intersection-maf005.vcor1[.vars]
#     rsparsepro/{EUR,AFR,ALL}/{TRAIT}.{LOCUS}.rsparsepro.txt
#
# Submit with:   qsub run_intersection_rsparsepro.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=24:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q capacity
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/intersection_finemap/script_logs/intersection_rsparsepro.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/intersection_finemap/script_logs/intersection_rsparsepro.err

set -euo pipefail

# ── Environment ────────────────────────────────────────────────────────────────
# The finemap conda env provides parsl (head-node interchange), pandas, and numpy.
# Each merge/RSparsePro command also embeds the finemap Python path directly.
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/10_Intersection_LD_RSparsePro
output_root=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/intersection_finemap
log_dir=${output_root}/script_logs
manifest=${script_dir}/intersection.manifest.tsv
merge_commands=${script_dir}/merge.commands.txt
rsparsepro_commands=${script_dir}/rsparsepro.commands.txt

# ── Create output + log directories ────────────────────────────────────────────
mkdir -p ${output_root} ${log_dir}
cd ${log_dir}

################################################
## Step 1: Generate merge and RSparsePro command lists
################################################
# Runs on every submission; only loci with complete linker/LD/sumstats
# inputs are included, sorted largest-LD-input first.
echo "NOTE: Generating intersection command lists..."
python ${script_dir}/create_intersection_commands.py \
    --output-root         ${output_root} \
    --manifest            ${manifest} \
    --merge-commands      ${merge_commands} \
    --rsparsepro-commands ${rsparsepro_commands}

echo "NOTE: Manifest written to ${manifest}"
echo "NOTE: Merge commands: $(wc -l < ${merge_commands})"
echo "NOTE: RSparsePro commands: $(wc -l < ${rsparsepro_commands})"

################################################
## Step 2: Run the merge commands via parsl
################################################
# Must complete before Step 3: RSparsePro reads this step's LD output.
echo "NOTE: Launching parsl allocator on merge commands..."
python ${script_dir}/intersection_parsl_allocator.py \
    --command-file ${merge_commands} \
    --log-dir      ${log_dir}

echo "NOTE: Merge step complete."

################################################
## Step 3: Run the RSparsePro commands via parsl
################################################
echo "NOTE: Launching parsl allocator on RSparsePro commands..."
python ${script_dir}/intersection_parsl_allocator.py \
    --command-file ${rsparsepro_commands} \
    --log-dir      ${log_dir}

echo "NOTE: run_intersection_rsparsepro.sh completed successfully."
