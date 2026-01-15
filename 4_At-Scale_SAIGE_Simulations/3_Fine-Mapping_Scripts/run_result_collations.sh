#!/bin/bash

#SBATCH -A bif154
#SBATCH -J Collate_finemap_parallel
#SBATCH -o run_result_collations.log
#SBATCH -N 20
#SBATCH -p batch
#SBATCH -t 2:00:00
#SBATCH -q debug

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

# Define base paths
CAUSAL_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants"
RESULTS_BASE_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"
OUTPUT_BASE_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/summary_results"

# Fine-mapping collator script
COLLATE_SCRIPT="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/fine_mapping_collator.py"

# Flux output directory
FLUX_OUTPUT_DIR="${OUTPUT_BASE_DIR}/flux_logs"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_BASE_DIR"
mkdir -p "$FLUX_OUTPUT_DIR"

################################################################################################################
############################################# Validate Setup ###################################################
################################################################################################################

echo "=================================================="
echo "VALIDATING SETUP"
echo "=================================================="

if [ ! -d "$CAUSAL_DIR" ]; then
    echo "ERROR: Causal variants directory does not exist: $CAUSAL_DIR"
    exit 1
fi

if [ ! -d "$RESULTS_BASE_DIR" ]; then
    echo "ERROR: Results base directory does not exist: $RESULTS_BASE_DIR"
    exit 1
fi

if [ ! -f "$COLLATE_SCRIPT" ]; then
    echo "ERROR: Collation script does not exist: $COLLATE_SCRIPT"
    exit 1
fi

echo "✓ Causal variants directory: $CAUSAL_DIR"
echo "✓ Results base directory: $RESULTS_BASE_DIR"
echo "✓ Output base directory: $OUTPUT_BASE_DIR"
echo "✓ Collation script: $COLLATE_SCRIPT"
echo "✓ Flux logs directory: $FLUX_OUTPUT_DIR"

################################################################################################################
########################### Find All Result Directories (fixed* and random*) ##################################
################################################################################################################

echo ""
echo "=================================================="
echo "DISCOVERING RESULT DIRECTORIES"
echo "=================================================="

# Create array to hold all result directories
declare -a RESULT_DIRS

# Find all directories starting with "fixed" or "random"
cd "$RESULTS_BASE_DIR" || exit 1

for dir in fixed* random* ; do
    if [ -d "$dir" ]; then
        RESULT_DIRS+=("$dir")
        echo "Found: $dir"
    fi
done

num_dirs=${#RESULT_DIRS[@]}

if [ $num_dirs -eq 0 ]; then
    echo "ERROR: No directories found starting with 'fixed' or 'random'"
    exit 1
fi

echo ""
echo "Found $num_dirs result directories to process"

################################################################################################################
################################### Create Collation Scripts for Flux ########################################
################################################################################################################

echo ""
echo "=================================================="
echo "CREATING COLLATION SCRIPTS FOR FLUX SUBMISSION"
echo "=================================================="

# Create temporary directory for collation scripts
SCRIPT_DIR="${FLUX_OUTPUT_DIR}/collation_scripts"
mkdir -p "$SCRIPT_DIR"

# Create individual collation scripts
script_count=0

for result_dir in "${RESULT_DIRS[@]}"; do
    script_count=$((script_count + 1))
    
    # Define paths for this directory
    RESULTS_DIR="${RESULTS_BASE_DIR}/${result_dir}"
    OUTPUT_FILE="${OUTPUT_BASE_DIR}/${result_dir}.results_file.txt"
    
    # Create a script for this collation job
    COLLATION_SCRIPT="${SCRIPT_DIR}/collate_${result_dir}.sh"
    
    cat > "$COLLATION_SCRIPT" << 'SCRIPT_EOF'
#!/bin/bash
set -e

# Get environment variables passed in
CAUSAL_DIR="$1"
RESULTS_DIR="$2"
OUTPUT_FILE="$3"
COLLATE_SCRIPT="$4"
DIR_NAME="$5"

# Load modules
module load miniforge3

echo "=========================================="
echo "Processing: $DIR_NAME"
echo "=========================================="
echo "Start time: $(date)"
echo "Input directory: $RESULTS_DIR"
echo "Output file: $OUTPUT_FILE"
echo ""

# First validate the directory structure
echo "Validating directory structure..."
python "$COLLATE_SCRIPT" \
    --causal-dir "$CAUSAL_DIR" \
    --results-dir "$RESULTS_DIR" \
    --validate-only

if [ $? -ne 0 ]; then
    echo "ERROR: Directory validation failed for $DIR_NAME"
    exit 1
fi

# Run the collation script
echo "Starting fine-mapping results collation for $DIR_NAME..."
python "$COLLATE_SCRIPT" \
    --causal-dir "$CAUSAL_DIR" \
    --results-dir "$RESULTS_DIR" \
    --output "$OUTPUT_FILE" \
    --verbose

if [ $? -eq 0 ]; then
    echo "✓ Collation completed successfully for $DIR_NAME!"
    echo "Output file: $OUTPUT_FILE"
    echo "File size: $(du -h "$OUTPUT_FILE" | cut -f1)"
    echo "Number of lines: $(wc -l < "$OUTPUT_FILE")"
else
    echo "ERROR: Collation failed for $DIR_NAME"
    exit 1
fi

echo "Job completed at: $(date)"
SCRIPT_EOF
    
    chmod +x "$COLLATION_SCRIPT"
    
    echo "Created: collate_${result_dir}.sh"
done

echo ""
echo "Created $script_count collation scripts"

################################################################################################################
################################### Submit Jobs to Flux #######################################################
################################################################################################################

echo ""
echo "=================================================="
echo "SUBMITTING COLLATION JOBS TO FLUX"
echo "=================================================="
echo "Start time: $(date)"
echo ""

# Get number of SLURM nodes available
NUM_NODES=${SLURM_NNODES:-1}

echo "Available SLURM nodes: $NUM_NODES"
echo "Collation directories to process: $num_dirs"
echo ""

# Write directory list to temporary file for passing to Flux
TEMP_DIR_FILE=$(mktemp)
for result_dir in "${RESULT_DIRS[@]}"; do
    echo "$result_dir" >> "$TEMP_DIR_FILE"
done

# Start Flux and submit all jobs
srun -N ${NUM_NODES} -n ${NUM_NODES} -c 56 flux start bash << 'FLUX_EOF'

# Display Flux resource allocation
echo "Flux instance started at $(date)"
echo "Available resources:"
flux resource list
echo ""

# Submit all collation jobs to Flux
echo "Submitting collation jobs to Flux queue..."
job_count=0

# Read directories from the temp file passed via environment
while IFS= read -r result_dir; do
    # Define paths
    RESULTS_DIR="${RESULTS_BASE_DIR}/${result_dir}"
    OUTPUT_FILE="${OUTPUT_BASE_DIR}/${result_dir}.results_file.txt"
    COLLATION_SCRIPT="${SCRIPT_DIR}/collate_${result_dir}.sh"
    
    # Submit job to Flux - one task per node
    flux submit -N1 -n1 -c2 \
        --output="${FLUX_OUTPUT_DIR}/collate_${result_dir}.log" \
        bash "${COLLATION_SCRIPT}" \
        "${CAUSAL_DIR}" \
        "${RESULTS_DIR}" \
        "${OUTPUT_FILE}" \
        "${COLLATE_SCRIPT}" \
        "${result_dir}"
    
    ((job_count++))
    
    # Print progress every 10 jobs
    if [ $((job_count % 10)) -eq 0 ]; then
        echo " Submitted ${job_count} jobs"
    fi
done < "${TEMP_DIR_FILE}"

echo "All collation jobs submitted to Flux queue at $(date)"
echo ""

# Wait for all jobs to complete
echo "Waiting for all collations to complete..."
flux queue drain

echo ""
echo "All collation jobs completed at $(date)"
echo "Final Flux resource status:"
flux resource list
FLUX_EOF

# Cleanup temp file
rm -f "$TEMP_DIR_FILE"

################################################################################################################
################################## Verify Completion and Summary ##############################################
################################################################################################################

echo ""
echo "=================================================="
echo "COLLATION COMPLETE - VERIFYING RESULTS"
echo "=================================================="

# Count successful output files
successful_count=0
failed_dirs=()

for result_dir in "${RESULT_DIRS[@]}"; do
    OUTPUT_FILE="${OUTPUT_BASE_DIR}/${result_dir}.results_file.txt"
    
    if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
        successful_count=$((successful_count + 1))
        file_size=$(du -h "$OUTPUT_FILE" | cut -f1)
        num_lines=$(wc -l < "$OUTPUT_FILE")
        echo "✓ $result_dir → $file_size ($num_lines lines)"
    else
        failed_dirs+=("$result_dir")
        echo "✗ $result_dir → MISSING OR EMPTY"
    fi
done

echo ""
echo "=================================================="
echo "SUMMARY"
echo "=================================================="
echo "Total directories: $num_dirs"
echo "Successfully processed: $successful_count"
echo "Failed: $((num_dirs - successful_count))"
echo ""

if [ ${#failed_dirs[@]} -gt 0 ]; then
    echo "Failed directories:"
    for dir in "${failed_dirs[@]}"; do
        echo "  - $dir (check ${FLUX_OUTPUT_DIR}/collate_${dir}.log)"
    done
else
    echo "✓ All collations completed successfully!"
fi

echo ""
echo "Output directory: $OUTPUT_BASE_DIR"
echo "Flux logs directory: $FLUX_OUTPUT_DIR"
echo "Completion time: $(date)"
echo "=================================================="
