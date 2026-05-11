#!/bin/bash
# run_matrix_making.sh
#
# PBS outer script for the PLATLAS LD matrix-making pipeline.
#
# Requests Aurora GPU nodes and runs:
#   1. create_matrix_commands.py  — build file lists + generate betanfer infer commands
#   2. matrix_parsl_allocator.py  — distribute jobs across GPU tiles via parsl
#
# Current settings: 2-node DEBUG run for initial testing.
# For production, change:
#   -l select=2  →  -l select=64   (or more; see note below)
#   -q debug     →  -q capacity
#   walltime=1:00:00  →  walltime=24:00:00
#
# Node count guidance (production):
#   ~162,000 (locus × cohort+pop) commands total.
#   12 GPU tiles/node → 12 concurrent betanfer runs/node.
#   At ~5 min/locus:  162,000 / (12 × 64 nodes) × 5 min ≈ 18 hours  ✓
#
#PBS -N matrix_making
#PBS -l select=1:ncpus=208:ngpus=6
#PBS -l filesystems=home:flare
#PBS -l walltime=1:00:00
#PBS -q debug
#PBS -A gpu_hack
#PBS -j oe
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/matrix_making.log

SCRIPT_DIR="/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/3_Matrix_Making"
BETANFER_DIR="/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer"
LOG_DIR="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"
CMDS_FILE="${SCRIPT_DIR}/matrix_making.commands.txt"

# ── Environment setup ───────────────────────────────────────────────────────────
source ~/.bashrc
conda activate /flare/GeomicVar/mconery/tools/parsl_gpu

# ── Step 1: Generate file lists and command file ────────────────────────────────
echo "NOTE: Generating matrix-making command list..."
python "${SCRIPT_DIR}/create_matrix_commands.py" \
    --command_file "${CMDS_FILE}" \
    --log_dir "${LOG_DIR}"

# ── Step 2: Distribute betanfer infer jobs via parsl ───────────────────────────
echo "NOTE: Launching parsl allocator..."
python "${SCRIPT_DIR}/matrix_parsl_allocator.py" \
    --command_file "${CMDS_FILE}" \
    --log_dir "${LOG_DIR}"
