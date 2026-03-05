#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_mvp_hg19_conversion.sh
#
# PBS submission script for converting MVP summary statistics from GRCh38
# to hg19 coordinates.
#
# Requests 32 nodes.  The parsl allocator (mvp_hg19_parsl_allocator.py) uses
# a LocalProvider to distribute work across these pre-allocated nodes
# (1 worker/node = 32 concurrent conversions) — no secondary PBS jobs are
# submitted.
#
# Each conversion job loads the full annotation file into memory to map
# GRCh38 positions to hg19 positions and RSIDs.  1 worker per node is used
# to avoid memory pressure from concurrent annotation loads.
#
# Output .sumstats.hg19.txt files are written to:
#   /flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19/{AFR,EUR}/
#
# The allocator skips files whose output .sumstats.hg19.txt already exists,
# so re-submitting this script is safe and will only run remaining jobs.
#
# Submit with:   qsub run_mvp_hg19_conversion.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=4:00:00
#PBS -l filesystems=home:flare
#PBS -l select=32:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q capacity
#PBS -N mvp_hg19_conversion
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_hg19_conversion.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/mvp_hg19_conversion.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
module load frameworks

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/0_MVP_hg19_Conversion
sumstats_dir=/flare/GeomicVar/multi-biobank-meta/sumstats/MVP
output_dir=/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19
annotation_file=/flare/GeomicVar/multi-biobank-meta/Annotation_gwPheWAS_dbGaP_b38.noAF_R2.txt.gz
log_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs
commands_file=${script_dir}/mvp_hg19.commands.txt

# ── Change to log directory so parsl internal logs land there too ─────────────
cd ${log_dir}

# ── Step 1: Generate the commands file ────────────────────────────────────────
echo "NOTE: Generating MVP hg19 conversion command list..."
python ${script_dir}/create_mvp_hg19_commands.py \
    --sumstats_dir    ${sumstats_dir}    \
    --output_dir      ${output_dir}      \
    --annotation_file ${annotation_file} \
    --output_file     ${commands_file}

if [ $? -ne 0 ]; then
    echo "ERROR: create_mvp_hg19_commands.py failed.  Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${commands_file}"
echo "NOTE: Total commands: $(wc -l < ${commands_file})"

# ── Step 2: Run the parsl allocator ───────────────────────────────────────────
# The allocator uses LocalProvider (defined in parsl_config_mvp_hg19.py) to
# spread 1 worker/node across the 32 pre-allocated nodes = 32 concurrent
# conversion runs.  It skips any file whose output already exists.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/mvp_hg19_parsl_allocator.py \
    --command_file ${commands_file} \
    --log_dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: mvp_hg19_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: mvp_hg19_parsl_allocator.py completed successfully."
