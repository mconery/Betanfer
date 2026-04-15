#!/bin/bash
# run_restricted_matrix_making.sh
#
# PBS outer script for re-running the 181 large-locus betanfer infer jobs that
# could not complete in the standard matrix-making pass.
#
# Runs:
#   1. create_restricted_commands.py  — identify failing jobs, build per-locus
#                                       variant lists from meta-analysis output,
#                                       and generate --variant-file commands
#   2. matrix_parsl_allocator.py      — distribute the restricted jobs across
#                                       GPU tiles via parsl (reused as-is)
#
#PBS -N restricted_matrix_making
#PBS -l select=1:ncpus=208:ngpus=6
#PBS -l filesystems=home:flare
#PBS -l walltime=2:00:00
#PBS -q capacity
#PBS -A GeomicVar
#PBS -j oe
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/restricted_matrix_making.log

SCRIPT_DIR="/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/3_Matrix_Making"
LOG_DIR="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"
MAIN_CMDS_FILE="${SCRIPT_DIR}/matrix_making.commands.txt"
RESTRICTED_CMDS_FILE="${SCRIPT_DIR}/matrix_making.restricted.commands.txt"

# ── Environment setup ────────────────────────────────────────────────────────────
source ~/.bashrc
conda activate /flare/GeomicVar/mconery/tools/parsl_gpu

# ── Step 1: Build variant lists and generate restricted command file ─────────────
echo "NOTE: Generating restricted command list..."
python "${SCRIPT_DIR}/create_restricted_commands.py" \
    --input_commands "${MAIN_CMDS_FILE}" \
    --command_file   "${RESTRICTED_CMDS_FILE}" \
    --log_dir        "${LOG_DIR}"

# ── Step 2: Distribute restricted betanfer infer jobs via parsl ──────────────────
echo "NOTE: Launching parsl allocator for restricted commands..."
python "${SCRIPT_DIR}/matrix_parsl_allocator.py" \
    --command_file "${RESTRICTED_CMDS_FILE}" \
    --log_dir      "${LOG_DIR}"
