#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# run_mvp_recompute_ld.sh
#
# PBS script to recompute LD matrices for the 175 MVP loci that failed
# fine-mapping due to either:
#   1. CHR:POS:REF:ALT .vars format (restricted-mode betanfer, pre-chunk-size)
#   2. Truncated .vcor1 files (GPU OOM during large-locus computation)
#
# Pipeline steps:
#   Step 1: create_mvp_recompute_ld_commands.py
#     Reads finemap_commands.txt and writes regular-mode betanfer infer commands
#     for the 175 failing loci (--chunk-size 15000, --plink2-format).
#
#   Step 2: Delete existing corrupted/restricted LD files
#     Parses --output prefixes from the commands file and removes .vcor1 and
#     .vcor1.vars for each locus so betanfer can write fresh output.
#     (matrix_parsl_allocator.py skips loci where .vcor1.vars already exists.)
#
#   Step 3: matrix_parsl_allocator.py
#     Distributes the 175 betanfer infer jobs across the allocated GPU nodes
#     (6 workers/node via parsl GPU-tile binding).
#
# After this job completes:
#   1. Re-run create_finemap_commands.py to rebuild finemap_commands.txt
#   2. Submit run_mvp_fine_mapping.sh (parsl allocator step only) to complete
#      fine-mapping for these 175 loci.
#
# Node sizing: 175 loci / (8 nodes × 6 workers) = ~4 rounds at ~5 min/locus
#              ≈ 20 minutes compute time; 4 hours walltime provides ample margin.
#
# Submit with:   qsub run_mvp_recompute_ld.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=2:00:00
#PBS -l filesystems=home:flare
#PBS -l select=16:ncpus=208:ngpus=6
#PBS -l place=scatter:excl
#PBS -q capacity
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_recompute_ld.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_recompute_ld.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate /flare/GeomicVar/mconery/tools/parsl_gpu

# ── Key paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/3_Matrix_Making
FINEMAP_CMDS=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/5_MVP_Fine-Mapping/finemap_commands.txt
CMDS_FILE=${SCRIPT_DIR}/mvp_recompute_ld.commands.txt
LOG_DIR=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs

# ── Change to log directory so parsl internal logs land there too ──────────────
cd ${LOG_DIR}

# ── Step 1: Generate betanfer infer commands for the 175 failing loci ──────────
echo "NOTE: Generating recompute commands..."
python ${SCRIPT_DIR}/create_mvp_recompute_ld_commands.py \
    --finemap_commands ${FINEMAP_CMDS} \
    --output_file      ${CMDS_FILE}

if [ $? -ne 0 ]; then
    echo "ERROR: create_mvp_recompute_ld_commands.py failed. Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${CMDS_FILE}"
echo "NOTE: Total commands: $(wc -l < ${CMDS_FILE})"

# ── Step 2: Delete corrupted/restricted LD files so betanfer can overwrite ─────
# Parses --output prefix from each command and removes .vcor1 and .vcor1.vars.
# matrix_parsl_allocator.py skips loci where .vcor1.vars already exists, so
# deletion is required before the allocator runs.
echo "NOTE: Removing existing LD files for the 175 loci..."
python -c "
import os, sys
removed = 0
with open(sys.argv[1]) as fh:
    for line in fh:
        parts = line.strip().split()
        for i, p in enumerate(parts):
            if p == '--output' and i + 1 < len(parts):
                prefix = parts[i + 1]
                for ext in ('.vcor1', '.vcor1.vars'):
                    path = prefix + ext
                    if os.path.exists(path):
                        os.remove(path)
                        removed += 1
print(f'NOTE: Removed {removed} LD files')
" ${CMDS_FILE}

# ── Step 3: Run betanfer infer via parsl ───────────────────────────────────────
echo "NOTE: Launching parsl allocator for LD recomputation..."
python ${SCRIPT_DIR}/matrix_parsl_allocator.py \
    --command_file ${CMDS_FILE} \
    --log_dir      ${LOG_DIR}

if [ $? -ne 0 ]; then
    echo "ERROR: matrix_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: LD recomputation complete."
echo "NOTE: Next steps:"
echo "      1. Re-run create_finemap_commands.py to generate fresh fine-mapping commands"
echo "      2. Submit run_mvp_fine_mapping.sh to complete fine-mapping for these 175 loci"
