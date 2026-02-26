#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_betanfer_indexing.sh
#
# PBS submission script for the betanfer index pipeline.
#
# Requests 32 nodes.  The parsl allocator (betanfer_index_parsl_allocator.py)
# uses a LocalProvider to distribute work across these pre-allocated nodes
# (8 workers/node = 256 concurrent betanfer index runs) — no secondary PBS
# jobs are submitted.
#
# betanfer index reads each .out file into pandas, sorts by chromosome and
# position, writes a BGZF-compressed .gz file, and creates a tabix .tbi
# index.  Output is written to the same directory as the input file.
#
# The allocator skips files whose .gz.tbi index already exists, so
# re-submitting this script is safe and will only run remaining jobs.
#
# Submit with:   qsub run_betanfer_indexing.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=32:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q debug-scaling
#PBS -N betanfer_indexing
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/betanfer_indexing.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/betanfer_indexing.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
module load frameworks

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/1_Meta-Analysis
results_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results
log_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs
commands_file=${script_dir}/betanfer_index.commands.txt

# ── Change to log directory so parsl internal logs land there too ─────────────
cd ${log_dir}

# ── Step 1: Generate the commands file ────────────────────────────────────────
echo "NOTE: Generating betanfer index command list..."
python ${script_dir}/create_index_commands.py \
    --results_dir ${results_dir} \
    --output_file ${commands_file}

if [ $? -ne 0 ]; then
    echo "ERROR: create_index_commands.py failed.  Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${commands_file}"
echo "NOTE: Total commands: $(wc -l < ${commands_file})"

# ── Step 2: Run the parsl allocator ───────────────────────────────────────────
# The allocator uses LocalProvider (defined in parsl_config_betanfer.py) to
# spread 8 workers/node across the 32 pre-allocated nodes = 256 concurrent
# betanfer index runs.  It skips any file whose .gz.tbi already exists.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/betanfer_index_parsl_allocator.py \
    --command_file ${commands_file} \
    --log_dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: betanfer_index_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: betanfer_index_parsl_allocator.py completed successfully."
