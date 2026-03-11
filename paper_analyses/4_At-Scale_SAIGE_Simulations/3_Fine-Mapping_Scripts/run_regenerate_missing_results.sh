#!/bin/bash

#SBATCH -A bif154
#SBATCH -J Regenerate_Missing
#SBATCH -o run_regenerate_missing_results.log
#SBATCH -N 50
#SBATCH -p batch
#SBATCH -t 4:00:00

# NOTE: Adjust -N (node count) based on the number of missing replicates.
# As a rough guide:
#   < 1,000 missing   →  5–10 nodes
#   1,000–10,000      →  20–50 nodes
#   > 10,000          →  50–120 nodes

################################################################################################################
################################################### Load Modules ###############################################
################################################################################################################

module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

PYTHON_SCRIPT="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/regenerate_missing_results.py"

MISSING_LIST="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/summary_results/missing_results_files.txt"

RESULTS_BASE="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"
CMD_FILES_DIR="${RESULTS_BASE}/command_files"
FLUX_OUTPUT_DIR="${RESULTS_BASE}/flux_output"

# Conda environments
CARMA_ENV="/lustre/orion/bif154/proj-shared/mconery/tools/environments/carma"
SUSIE_ENV="/lustre/orion/bif154/proj-shared/mconery/tools/environments/susie"
RSPARSEPRO_ENV="/lustre/orion/bif154/proj-shared/mconery/tools/environments/rsparsepro"
MULTISUSIE_ENV="/lustre/orion/bif154/proj-shared/mconery/tools/environments/multisusie"
# MESuSiE uses the carma environment; SuSiEx uses the default miniforge3 base env.

NUM_NODES=${SLURM_NNODES:-50}

mkdir -p "$CMD_FILES_DIR"
mkdir -p "$FLUX_OUTPUT_DIR"

################################################################################################################
############################################# Validate Setup ###################################################
################################################################################################################

echo "=================================================="
echo "REGENERATE MISSING FINE-MAPPING RESULTS"
echo "=================================================="
echo "Start time: $(date)"
echo ""
echo "Missing list : $MISSING_LIST"
echo "Command dir  : $CMD_FILES_DIR"
echo "Flux logs    : $FLUX_OUTPUT_DIR"
echo "Nodes        : $NUM_NODES"
echo ""

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

if [ ! -f "$MISSING_LIST" ]; then
    echo "ERROR: Missing list file not found: $MISSING_LIST"
    exit 1
fi

NUM_MISSING=$(wc -l < "$MISSING_LIST")
echo "Missing replicates to process: $NUM_MISSING"
echo ""

################################################################################################################
################################## Step 1: Generate Per-Method Command Files ##################################
################################################################################################################

echo "=================================================="
echo "STEP 1: Generating per-method command files"
echo "=================================================="

python "$PYTHON_SCRIPT" \
    --missing-list "$MISSING_LIST" \
    --output-dir   "$CMD_FILES_DIR" \
    --populations  "AFR,EUR" \
    --plink-conv   "CHR#.plink.maf"

if [ $? -ne 0 ]; then
    echo "ERROR: Command generation script failed."
    exit 1
fi

echo ""
echo "Command files written to: $CMD_FILES_DIR"
echo ""

################################################################################################################
####################################### Helper: run one method via Flux #######################################
################################################################################################################

# Usage: run_method METHOD CMD_FILE ENV_PATH
#   ENV_PATH may be empty string for methods that use the base miniforge3 env.
run_method() {
    local method="$1"
    local cmd_file="$2"
    local env_path="$3"

    echo "=================================================="
    echo "Running method: ${method^^}"
    echo "Command file  : $cmd_file"
    echo "Start time    : $(date)"
    echo "=================================================="

    if [ ! -f "$cmd_file" ] || [ ! -s "$cmd_file" ]; then
        echo "  No commands for $method (file missing or empty). Skipping."
        echo ""
        return 0
    fi

    local num_cmds
    num_cmds=$(wc -l < "$cmd_file")
    echo "  Commands to run: $num_cmds"

    # Activate the method-specific conda environment (if any)
    if [ -n "$env_path" ]; then
        echo "  Activating environment: $env_path"
        source activate "$env_path"
    else
        echo "  Using default miniforge3 environment (no activation)"
    fi

    # Launch Flux and submit all commands for this method
    srun -N ${NUM_NODES} -n ${NUM_NODES} -c 56 flux start bash -c "
        echo 'Flux started for ${method} at \$(date)'
        flux resource list
        echo ''

        mapfile -t cmds < '${cmd_file}'
        total_cmds=\${#cmds[@]}
        echo \"Submitting \${total_cmds} ${method} jobs to Flux...\"

        for ((i=0; i<\${#cmds[@]}; i++)); do
            flux submit -N1 -n1 -c2 \
                --output='${FLUX_OUTPUT_DIR}/${method}_\${i}.log' \
                \${cmds[\$i]}
        done

        echo \"All ${method} jobs submitted. Waiting for completion...\"
        flux queue drain
        echo \"${method} complete at \$(date)\"
    "

    local exit_code=$?

    # Deactivate environment if one was activated
    if [ -n "$env_path" ]; then
        conda deactivate
    fi

    if [ $exit_code -ne 0 ]; then
        echo "ERROR: Flux execution failed for $method (exit code $exit_code)"
        echo "Check logs in: $FLUX_OUTPUT_DIR"
        return $exit_code
    fi

    echo "  $method completed successfully at $(date)"
    echo ""
    return 0
}

################################################################################################################
#################### Step 2: Execute Methods Sequentially (each needs its own environment) ####################
################################################################################################################

echo "=================================================="
echo "STEP 2: Executing per-method command files"
echo "(Methods run sequentially; commands within each method run in parallel)"
echo "=================================================="
echo ""

# ── CARMA ──────────────────────────────────────────────────────────────────────
run_method "carma" \
    "${CMD_FILES_DIR}/carma_regenerate_commands.txt" \
    "$CARMA_ENV"
if [ $? -ne 0 ]; then exit 1; fi

# ── SuSiE ──────────────────────────────────────────────────────────────────────
run_method "susie" \
    "${CMD_FILES_DIR}/susie_regenerate_commands.txt" \
    "$SUSIE_ENV"
if [ $? -ne 0 ]; then exit 1; fi

# ── RSparsePro ─────────────────────────────────────────────────────────────────
run_method "rsparsepro" \
    "${CMD_FILES_DIR}/rsparsepro_regenerate_commands.txt" \
    "$RSPARSEPRO_ENV"
if [ $? -ne 0 ]; then exit 1; fi

# ── MESuSiE (uses CARMA environment) ───────────────────────────────────────────
run_method "mesusie" \
    "${CMD_FILES_DIR}/mesusie_regenerate_commands.txt" \
    "$CARMA_ENV"
if [ $? -ne 0 ]; then exit 1; fi

# ── SuSiEx (no special environment — uses base miniforge3) ─────────────────────
run_method "susiex" \
    "${CMD_FILES_DIR}/susiex_regenerate_commands.txt" \
    ""
if [ $? -ne 0 ]; then exit 1; fi

# ── MultiSuSiE ─────────────────────────────────────────────────────────────────
run_method "multisusie" \
    "${CMD_FILES_DIR}/multisusie_regenerate_commands.txt" \
    "$MULTISUSIE_ENV"
if [ $? -ne 0 ]; then exit 1; fi

################################################################################################################
################################################ Summary #######################################################
################################################################################################################

echo "=================================================="
echo "REGENERATION COMPLETE"
echo "=================================================="
echo "End time: $(date)"
echo ""
echo "Command files in : $CMD_FILES_DIR"
echo "Flux logs in     : $FLUX_OUTPUT_DIR"
echo ""
echo "To verify results, re-run find_missing_results.py and compare the"
echo "missing count to the value before this job."
echo "=================================================="
