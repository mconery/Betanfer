#!/bin/bash
# run_ld_merge.sh
#
# PBS outer script for the PLATLAS LD matrix merging pipeline.
#
# Requests Aurora CPU nodes (no GPUs needed) and runs three steps:
#   1. extract_sample_sizes.py  — one-time extraction of max(N) from sumstats
#                                 (no-op if sample_sizes.tsv already exists)
#   2. create_merge_commands.py — completeness + checkpoint checks, writes
#                                 commands file for jobs still to run
#   3. merge_parsl_allocator.py — distributes merge jobs across CPU nodes
#                                 via parsl
#
# Current settings: 2-node DEBUG run for initial testing.
# For production, change:
#   -l select=2  →  -l select=8   (or more; see note below)
#   -q debug-scaling  →  -q capacity
#   walltime=1:00:00  →  walltime=24:00:00
#
# Node count guidance (production):
#   4 workers/node × N_nodes = concurrent merge jobs.
#   Estimate: ~20,000 ALL/AFR/EUR merge jobs total (depends on loci).
#   At ~5 min/job:  20,000 / (4 × 8 nodes) × 5 min ≈ 104 hours with 8 nodes
#   At ~5 min/job:  20,000 / (4 × 32 nodes) × 5 min ≈ 26 hours with 32 nodes
#   Resubmit the job after walltime to continue; completed merges are skipped.
#
#PBS -N ld_merge
#PBS -l select=16:ncpus=208
#PBS -l filesystems=home:flare
#PBS -l walltime=3:00:00
#PBS -q capacity
#PBS -A gpu_hack
#PBS -j oe
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/ld_merge.log

SCRIPT_DIR="/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/4_LD_Matrix_Merging"
LOG_DIR="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"
CMDS_FILE="${SCRIPT_DIR}/ld_merge.commands.txt"

# ── Environment setup ────────────────────────────────────────────────────────
source ~/.bashrc
conda activate /flare/GeomicVar/mconery/tools/parsl_gpu

# ── Step 1: Extract sample sizes (checkpoints if file already exists) ────────
echo "NOTE: Extracting sample sizes from sumstats..."
python "${SCRIPT_DIR}/extract_sample_sizes.py"

# ── Step 2: Generate merge command file ──────────────────────────────────────
echo "NOTE: Generating merge command list..."
python "${SCRIPT_DIR}/create_merge_commands.py" \
    --command_file "${CMDS_FILE}"

# ── Step 3: Distribute merge jobs via parsl ──────────────────────────────────
echo "NOTE: Launching parsl allocator..."
python "${SCRIPT_DIR}/merge_parsl_allocator.py" \
    --command_file "${CMDS_FILE}" \
    --log_dir "${LOG_DIR}"
