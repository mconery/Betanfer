#!/bin/bash

#SBATCH -A bif154
#SBATCH -J gwama_meta_analysis
#SBATCH -o run_gwama_meta_analysis_%j.log
#SBATCH -N 40
#SBATCH -p batch
#SBATCH -t 2:00:00

##################################################################################################################
################################################### Set Needed Variables #########################################
##################################################################################################################

# Comma-separated population names
POP_NAMES="AFR,EUR"
#Enter meta-analysis types
ANALYSIS_TYPE=random

# Comma-separated population directories matching these names
POP_DIRS="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/AFR,/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/EUR"

# Output directory for meta-analysis
OUTPUT_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"

# Command files directory
command_files_dir="${OUTPUT_DIR}/command_files"
mkdir -p "$command_files_dir"

# Path to the updated generate_meta_commands.py script
CMD_SCRIPT="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/generate_meta_commands.py"

# Python splitter for batch creation
python_splitter="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/split_commands_parallel.py"

# Batch configuration
commands_per_batch=100

# Master commands file
master_commands_file="${command_files_dir}/master_meta_analysis_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Number of nodes allocated by SLURM
num_nodes=${SLURM_NNODES}

##################################################################################################################
############################################## Load Modules and Setup ###########################################
##################################################################################################################

module load python/3.10-miniforge3
module load flux

mkdir -p ${OUTPUT_DIR}

##################################################################################################################
################################################### Generate Commands ############################################
##################################################################################################################

echo "Generating meta-analysis commands with populations: ${POP_NAMES}"

python ${CMD_SCRIPT} \
    --pop_names ${POP_NAMES} \
    --pop_dirs ${POP_DIRS} \
    --output_dir ${OUTPUT_DIR} \
    --commands_file ${master_commands_file} \
    --analysis_type ${ANALYSIS_TYPE}

if [ ! -f ${master_commands_file} ] || [ ! -s ${master_commands_file} ]; then
    echo "Error: No commands generated or commands file is empty"
    exit 1
fi

NUM_COMMANDS=$(wc -l $master_commands_file)
echo "Generated ${NUM_COMMANDS} meta-analysis commands"

##################################################################################################################
################################## Create batched command files ################################################
##################################################################################################################

echo "Creating batched command files using parallel Python splitter..."

# Clean up any existing meta analysis batch files
rm -f "${command_files_dir}"/meta_analysis_batch_*.sh

# Use Python splitter for fast parallel batch creation
python "${python_splitter}" \
    -i "${master_commands_file}" \
    -o "${command_files_dir}" \
    -b ${commands_per_batch} \
    -p "meta_analysis_batch_" \
    -s ".sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Python batch splitter failed"
    exit 1
fi

# Verify batch creation
num_batch_files=$(ls "${command_files_dir}" | grep -c "^meta_analysis_batch_" || true)
echo "Verified: Created ${num_batch_files} batch files in ${command_files_dir}"
echo ""

if [ ${num_batch_files} -eq 0 ]; then
    echo "ERROR: No batch files were created"
    exit 1
fi

##################################################################################################################
####################################### Execute Commands in Parallel #########################################
##################################################################################################################

echo "Starting flux-based execution of batched GWAMA meta-analysis commands at $(date)"

# Create flux output directory
flux_output_dir="${OUTPUT_DIR}/flux_output"
mkdir -p "$flux_output_dir"

# Start Flux and submit batch jobs
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start bash -c "
    # Display Flux resource allocation
    echo \"Flux instance started at \$(date)\"
    echo \"Available resources:\"
    flux resource list
    echo \"\"
    
    # Submit all batch files to Flux
    echo \"Submitting ${num_batch_files} batches to Flux queue...\"
    batch_count=0
    for batch in ${command_files_dir}/meta_analysis_batch_*.sh; do
        batch_name=\$(basename \${batch} .sh)
        flux submit -N1 -n1 -c4 \
            --output=${flux_output_dir}/\${batch_name}.log \
            bash \${batch}
        
        ((batch_count++))
        
        # Print progress every 10 batches
        if [ \$((batch_count % 10)) -eq 0 ]; then
            echo \"  Submitted \${batch_count}/${num_batch_files} batches\"
        fi
    done
    
    echo \"All ${num_batch_files} batches submitted to Flux queue at \$(date)\"
    echo \"\"
    
    # Wait for all jobs to complete
    echo \"Waiting for all batches to complete...\"
    flux queue drain
    
    echo \"\"
    echo \"All batches completed at \$(date)\"
    
    # Display final resource status
    echo \"Final Flux resource status:\"
    flux resource list
"

##################################################################################################################
################################## Post-processing and cleanup ##################################################
##################################################################################################################

echo ""
echo "Post-processing results..."

# Count successful batch completions by checking for output files
successful_batches=$(ls -1 "${flux_output_dir}"/meta_analysis_batch_*.log 2>/dev/null | wc -l)
echo "Successfully completed batches: ${successful_batches}/${num_batch_files}"

# Check for any failed batches
if [ ${successful_batches} -ne ${num_batch_files} ]; then
    echo "WARNING: Some batches may have failed. Check logs in ${flux_output_dir}"
fi

echo ""
echo "=================================================="
echo "GWAMA META-ANALYSIS WORKFLOW COMPLETED"
echo "=================================================="
echo "End time: $(date)"
echo "Total meta-analysis commands: ${NUM_COMMANDS}"
echo "Total batches: ${num_batch_files}"
echo "Batch files retained in: ${command_files_dir}"
echo "Flux logs available in: ${flux_output_dir}"
echo "Meta-analysis results in: ${OUTPUT_DIR}"
echo "=================================================="
