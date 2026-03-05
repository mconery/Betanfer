#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_mvp_hg19_indexing.sh
#
# PBS submission script for betanfer indexing of the MVP hg19 summary
# statistics produced by run_mvp_hg19_conversion.sh.
#
# Requests 32 nodes.  The parsl allocator (mvp_hg19_index_parsl_allocator.py)
# uses a LocalProvider to distribute work across these pre-allocated nodes
# (4 workers/node = 128 concurrent betanfer index runs) — no secondary PBS
# jobs are submitted.
#
# betanfer index reads each .sumstats.hg19.txt file into pandas, sorts by
# chromosome and position, writes a BGZF-compressed .gz file, and creates
# a tabix .tbi index.  Output is written to the same directory as the input.
#
# The allocator skips files whose .gz.tbi index already exists, so
# re-submitting this script is safe and will only run remaining jobs.
#
# Run AFTER run_mvp_hg19_conversion.sh has completed successfully.
#
# Submit with:   qsub run_mvp_hg19_indexing.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=2:00:00
#PBS -l filesystems=home:flare
#PBS -l select=32:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q capacity
#PBS -N mvp_hg19_indexing
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_hg19_indexing.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_hg19_indexing.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
module load frameworks

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/0_MVP_hg19_Conversion
hg19_dir=/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19
log_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs
commands_file=${script_dir}/mvp_hg19_index.commands.txt

# ── Change to log directory so parsl internal logs land there too ─────────────
cd ${log_dir}

# ── Step 1: Generate the commands file ────────────────────────────────────────
echo "NOTE: Generating betanfer index command list for MVP hg19 sumstats..."
python ${script_dir}/create_mvp_hg19_index_commands.py \
    --hg19_dir    ${hg19_dir}    \
    --output_file ${commands_file}

if [ $? -ne 0 ]; then
    echo "ERROR: create_mvp_hg19_index_commands.py failed.  Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${commands_file}"
echo "NOTE: Total commands: $(wc -l < ${commands_file})"

# ── Step 2: Run the parsl allocator ───────────────────────────────────────────
# The allocator uses LocalProvider (defined in parsl_config_mvp_hg19_index.py)
# to spread 4 workers/node across the 32 pre-allocated nodes = 128 concurrent
# betanfer index runs.  It skips any file whose .gz.tbi already exists.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/mvp_hg19_index_parsl_allocator.py \
    --command_file ${commands_file} \
    --log_dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: mvp_hg19_index_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: mvp_hg19_index_parsl_allocator.py completed successfully."
