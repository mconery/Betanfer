#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_gwama_analysis.sh
#
# PBS submission script for the GWAMA fixed-effects meta-analysis pipeline.
#
# This script requests ONE compute node to:
#   1. Generate the per-job command list via create_gwama_commands.py.
#   2. Launch the parsl allocator (gwama_parsl_allocator.py), which then
#      submits its own PBS worker jobs (8 nodes, 4 workers/node = 32 concurrent
#      GWAMA runs) to Aurora's capacity queue.
#
# Temp decompressed files go to TEMP_DIR (separate from the original sumstats
# files and their .tbi indices, which are never modified).
#
# The allocator skips trait/population combinations whose .out output file
# already exists, so re-submitting this script is safe and will only run
# remaining jobs.
#
# Submit with:   qsub run_gwama_analysis.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=48:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q capacity
#PBS -M mconery@anl.gov
#PBS -N gwama_meta_analysis
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwama_meta_analysis.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwama_meta_analysis.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
conda activate finemap_old

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/1_Meta-Analysis
linker_file=/flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv
sumstats_base=/flare/GeomicVar/multi-biobank-meta/sumstats
output_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results
temp_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/gwama_temp
log_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs
commands_file=${script_dir}/gwama_meta_analysis.commands.txt

# ── Create temp directory if it doesn't exist ─────────────────────────────────
mkdir -p ${temp_dir}

# ── Change to log directory so parsl internal logs land there too ─────────────
cd ${log_dir}

# ── Step 1: Generate the commands file ────────────────────────────────────────
echo "NOTE: Generating GWAMA command list..."
python ${script_dir}/create_gwama_commands.py \
    --linker_file   ${linker_file}   \
    --sumstats_base ${sumstats_base} \
    --output_file   ${commands_file} \
    --output_dir    ${output_dir}    \
    --temp_dir      ${temp_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: create_gwama_commands.py failed.  Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${commands_file}"
echo "NOTE: Total commands: $(wc -l < ${commands_file})"

# ── Step 2: Run the parsl allocator ──────────────────────────────────────────
# The allocator submits its own PBS worker jobs (defined in parsl_config.py:
#   8 nodes/block × 4 workers/node × 1 max block = 32 concurrent GWAMA runs).
# It skips any trait/population combination whose .out file already exists.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/gwama_parsl_allocator.py \
    --command_file ${commands_file} \
    --log_dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: gwama_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: gwama_parsl_allocator.py completed successfully."
