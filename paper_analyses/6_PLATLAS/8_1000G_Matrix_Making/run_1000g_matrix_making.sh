#!/bin/bash
# run_1000g_matrix_making.sh
#
# PBS outer script for the 1000G LD matrix-making pipeline.
#
# Requests Aurora CPU nodes and runs:
#   1. create_1000g_commands.py  — build per-population keep files and generate
#                                  plink2 --r-unphased commands
#   2. matrix_1000g_parsl_allocator.py — distribute jobs across nodes via parsl
#
# Each plink2 job gets a full node (1 worker/node) with all 208 CPUs available
# via --threads 208, since plink2 is natively multi-threaded.
#
# Current settings: 16-node run.
# For production, change:
#   -l select=16  →  -l select=64   (or more; see note below)
#   walltime=12:00:00  →  walltime=24:00:00
#
# Node count guidance (production):
#   ~45,732 total (locus × population) commands (22,866 MVP_AFR + 22,866 MVP_EUR).
#   1 job/node; plink2 on a 1 Mb region with 208 threads ≈ 1–3 min/locus.
#   At ~2 min/locus:  45,732 / 64 nodes × 2 min ≈ 24 hours  ✓
#
#PBS -N 1000g_matrix_making
#PBS -l select=256:ncpus=208
#PBS -l filesystems=home:flare
#PBS -l walltime=1:00:00
#PBS -q debug-scaling
#PBS -A GeomicVar
#PBS -j oe
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/1000g_matrix_making.log

SCRIPT_DIR="/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/8_1000G_Matrix_Making"
LOG_DIR="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"
CMDS_FILE="${SCRIPT_DIR}/1000g_matrix_making.commands.txt"

# ── Environment setup ───────────────────────────────────────────────────────────
source ~/.bashrc
conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap

# ── Step 1: Generate keep files and command file ────────────────────────────────
echo "NOTE: Generating 1000G keep files and matrix-making command list..."
python "${SCRIPT_DIR}/create_1000g_commands.py" \
    --command_file "${CMDS_FILE}" \
    --log_dir "${LOG_DIR}"

# ── Step 2: Distribute plink2 jobs via parsl ────────────────────────────────────
echo "NOTE: Launching parsl allocator..."
python "${SCRIPT_DIR}/matrix_1000g_parsl_allocator.py" \
    --command_file "${CMDS_FILE}" \
    --log_dir "${LOG_DIR}"
