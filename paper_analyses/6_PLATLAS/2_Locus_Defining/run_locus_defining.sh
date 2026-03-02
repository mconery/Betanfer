#!/bin/bash -l
# ─────────────────────────────────────────────────────────────────────────────
# run_locus_defining.sh
#
# PBS submission script for the PLATLAS locus-defining pipeline.
#
# Requests 32 nodes.  The parsl allocator (locus_parsl_allocator.py) uses a
# LocalProvider to distribute work across these pre-allocated nodes directly
# (4 workers/node = 128 concurrent locus-defining runs) — no secondary PBS
# jobs are submitted.
#
# Covers five population inputs:
#   Meta-analysis results (GWAMA output):   EUR, AFR, ALL  [tile-based discovery]
#   MVP single-cohort:                      AFR, EUR       [predefined from S11]
#
# MVP master locus files are written directly from Supplementary_Table-S11.txt
# by create_mvp_loci_from_s11.py (Step 0), bypassing per-trait computation.
# Meta-analysis loci use the tile-based algorithm via define_sig_loci_for_trait.py.
#
# The allocator skips any trait/population whose .sig_loci.csv sentinel file
# already exists, so re-submitting this script is safe and will only process
# remaining traits.  After all per-trait jobs finish, the allocator combines
# results into per-population master CSV files in the output directory.
#
# Submit with:   qsub run_locus_defining.sh
# Monitor with:  qstat -u mconery
# ─────────────────────────────────────────────────────────────────────────────

#PBS -A GeomicVar
#PBS -l walltime=1:00:00
#PBS -l filesystems=home:flare
#PBS -l select=1:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q debug-scaling
#PBS -M mconery@anl.gov
#PBS -N locus_defining
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/locus_defining.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/locus_defining.err

# ── Environment ────────────────────────────────────────────────────────────────
source $HOME/.bashrc
module load frameworks

# ── Key paths ─────────────────────────────────────────────────────────────────
script_dir=/lus/flare/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/6_PLATLAS/2_Locus_Defining
meta_results_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results
output_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci
s11_file=${output_dir}/Supplementary_Table-S11.txt
log_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs
commands_file=${script_dir}/locus_defining.commands.txt

# ── Create output directory if it doesn't exist ───────────────────────────────
mkdir -p ${output_dir}

# ── Change to log directory so parsl internal logs land there too ─────────────
cd ${log_dir}

# ── Step 0: Write MVP master locus files from Supplementary Table S11 ─────────
echo "NOTE: Creating MVP master locus files from Supplementary Table S11..."
python ${script_dir}/create_mvp_loci_from_s11.py \
    --s11_file   ${s11_file} \
    --output_dir ${output_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: create_mvp_loci_from_s11.py failed.  Aborting."
    exit 1
fi

# ── Step 1: Generate the commands file (meta-analysis populations only) ────────
echo "NOTE: Generating locus-defining command list..."
python ${script_dir}/create_locus_commands.py \
    --meta_results_dir ${meta_results_dir} \
    --output_dir       ${output_dir}       \
    --output_file      ${commands_file}

if [ $? -ne 0 ]; then
    echo "ERROR: create_locus_commands.py failed.  Aborting."
    exit 1
fi

echo "NOTE: Commands written to ${commands_file}"
echo "NOTE: Total commands: $(wc -l < ${commands_file})"

# ── Step 2: Run the parsl allocator ──────────────────────────────────────────
# The allocator uses LocalProvider (defined in parsl_config.py) to spread
# 4 workers/node across the 32 pre-allocated nodes = 128 concurrent runs.
# It skips any trait/population whose .sig_loci.csv sentinel already exists.
# After all per-trait futures resolve, it combines results into master files.
echo "NOTE: Launching parsl allocator..."
python ${script_dir}/locus_parsl_allocator.py \
    --command_file ${commands_file} \
    --log_dir      ${log_dir}

if [ $? -ne 0 ]; then
    echo "ERROR: locus_parsl_allocator.py exited with an error."
    exit 1
fi

echo "NOTE: locus_parsl_allocator.py completed successfully."
