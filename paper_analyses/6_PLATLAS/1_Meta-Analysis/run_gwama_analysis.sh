#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_gwama_analysis.sh
#
# PBS submission script for the GWAMA fixed-effects meta-analysis pipeline.
#
# Requests 128 nodes.  The parsl allocator (gwama_parsl_allocator.py) uses a
# LocalProvider to distribute work across these pre-allocated nodes directly
# (4 workers/node = 512 concurrent GWAMA runs) — no secondary PBS jobs are
# submitted.
#
# Temp decompressed files go to TEMP_DIR (separate from the original sumstats
# files and their .tbi indices, which are never modified).
#
# The allocator skips trait/population combinations whose .done sentinel file
# already exists, so re-submitting this script is safe and will only run
# remaining jobs.
#
# Submit with:   qsub run_gwama_analysis.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=128:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q debug-scaling
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
# The allocator uses LocalProvider (defined in parsl_config.py) to spread
# 4 workers/node across the 128 pre-allocated nodes = 512 concurrent GWAMA runs.
# It skips any trait/population combination whose .done sentinel already exists.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/gwama_parsl_allocator.py \
    --command_file ${commands_file} \
    --log_dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: gwama_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: gwama_parsl_allocator.py completed successfully."
