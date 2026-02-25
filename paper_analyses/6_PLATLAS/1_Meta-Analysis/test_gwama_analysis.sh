#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# test_gwama_analysis.sh
#
# Debug PBS job to test the full GWAMA pipeline (LocalProvider parsl approach)
# on the debug-scaling queue with a 1-hour walltime.
#
# Requests 8 nodes — same as the production script — so that the LocalProvider
# configuration (4 workers/node = 32 concurrent GWAMA runs) is exercised under
# real conditions.
#
# Submit with:  qsub test_gwama_analysis.sh
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=8:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q debug-scaling
#PBS -M mconery@anl.gov
#PBS -N gwama_analysis_test
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwama_analysis_test.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/gwama_analysis_test.err

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

# ── Diagnostics ────────────────────────────────────────────────────────────────
echo "=== GWAMA analysis test started: $(date) ==="
echo "Head node: $(hostname)"
echo "PBS_NODEFILE: ${PBS_NODEFILE}"
echo "Allocated nodes ($(sort -u ${PBS_NODEFILE} | wc -l)):"
sort -u ${PBS_NODEFILE}

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
# 4 workers/node across the 8 pre-allocated nodes = 32 concurrent GWAMA runs.
# It skips any trait/population combination whose .out file already exists.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/gwama_parsl_allocator.py \
    --command_file ${commands_file} \
    --log_dir      ${log_dir}

exit_code=$?
echo "=== Test finished: $(date) | exit code: ${exit_code} ==="
exit ${exit_code}
