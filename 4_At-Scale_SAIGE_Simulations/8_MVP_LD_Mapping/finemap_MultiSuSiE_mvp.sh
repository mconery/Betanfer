#!/bin/bash

#SBATCH -A bif154_mde
#SBATCH -J MultiSuSiE_finemap_mvp
#SBATCH -o finemap_MultiSuSiE_mvp.log
#SBATCH -N 40
#SBATCH -p batch-spi
#SBATCH -t 2:00:00

################################################################################################################
################################################### Load Modules ###############################################
################################################################################################################

module load Core/25.03
module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set script locations
map_script="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/MultiSuSiE/finemap_MultiSuSiE.py"
command_script="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/MultiSuSiE/generate_multisusie_commands.py"
python_splitter="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Meta_Fine-Mapping/5_LD_Inference/split_commands_parallel.py"

# Set environment variables
ENV_NAME=multisusie
TAR_FILE=${ENV_NAME}.tar.gz
TAR_DIR=/sw/frontier/mde/bif154_mde/tools/environments

# Set directories
sumstats_dir="/sw/frontier/mde/bif154_mde/betanfer_project/step2_results"
matrix_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/LD_matrices"
pheno_dir="/sw/frontier/mde/bif154_mde/betanfer_project/simulated_phenotypes"

# Base output directory for fine-mapping results
base_finemap_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/fine-mapping_results"

# Command files directory
command_files_dir="${base_finemap_dir}/command_files"
mkdir -p "$command_files_dir"

# Master commands file
master_commands_file="${command_files_dir}/master_multisusie_realistic_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Populations
populations="AFR,EUR"

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Define heritabilities
heritabilities=("5e-5" "3e-4")
heritability_values=(0.00005 0.0003)

################################################################################################################
################################## Generate Commands for All Combinations ####################################
################################################################################################################

echo "=================================================="
echo "STARTING MULTISUSIE FINE-MAPPING PIPELINE (All Configurations)"
echo "=================================================="
echo "Start time: $(date)"
echo "Per-SNP heritabilities: ${heritabilities[@]}"
echo "Master commands file: $master_commands_file"
echo "Batch size: ${commands_per_batch} commands per batch"
echo "=================================================="

# Counters
total_command_sets=0
total_commands=0

# Loop over heritabilities
for i in "${!heritabilities[@]}"; do
    h_label=${heritabilities[$i]}
    h_value=${heritability_values[$i]}
    
    echo ""
    echo "========================================="
    echo "Processing heritability: $h_label"
    echo "========================================="

        echo ""
        echo "Generating commands for heritability: $h_label"

        # Define matrix types for this configuration
	matrix_types="mvp_external"        

        # Define output directory and command file
        if [[ $h_label == "5e-5" ]]; then
            h_dir="fixed_homo_5e-5"
        else
            h_dir="fixed_homo_3e-4"
        fi

        out_dir="${base_finemap_dir}/${h_dir}/multisusie_results"
        command_file="${command_files_dir}/multisusie_realistic_${h_label}_commands.txt"

        mkdir -p "$out_dir"

        # Generate commands for this configuration
        echo "  Running command generation script..."
        python "$command_script" \
            --wrapper-script "$map_script" \
            --gwas-dir "$sumstats_dir" \
            --ld-dir "$matrix_dir" \
            --pheno-dir "$pheno_dir" \
            --output-dir "$out_dir" \
            --commands-file "$command_file" \
            --matrix-types "$matrix_types" \
            --populations "$populations" \
            --skip-existing \
            --heritability "$h_value"

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to generate commands for heritability: $h_label"
            exit 1
        fi

        if [[ ! -f "$command_file" ]]; then
            echo "ERROR: Commands file is missing or empty: $command_file"
            exit 1
        fi

        num_commands=$(wc -l < "$command_file")
        echo "  Generated $num_commands commands"

        ((total_commands += num_commands))
        ((total_command_sets++))

        # Append to master commands file
        grep -v "^#" "$command_file" >> "$master_commands_file"

done

################################################################################################################
################################## Summary and Final Status ###################################################
################################################################################################################

echo ""
echo "=================================================="
echo "COMMAND FILE GENERATION COMPLETE"
echo "=================================================="
echo "End time: $(date)"
echo "Total command sets generated: $total_command_sets"
echo "Total commands in master file: $(wc -l < "$master_commands_file")"
echo "Master commands file: $master_commands_file"
echo "Command directory: $command_files_dir"
echo "=================================================="

# Verify master commands file
if [[ ! -f "$master_commands_file" ]] || [[ ! -s "$master_commands_file" ]]; then
    echo "ERROR: Master commands file is missing or empty!"
    exit 1
fi

################################################################################################################
################################## Create batched command files ################################################
################################################################################################################

echo ""
echo "Creating batched command files using parallel Python splitter..."

# Clean up any existing multisusie batch files
rm -f "${command_files_dir}"/multisusie_realistic_finemap_batch_*.sh

#Calculate number of batches
num_batches=$(( num_nodes * 28 ))
# Use Python splitter for fast parallel batch creation
python "${python_splitter}" \
    -i "${master_commands_file}" \
    -o "${command_files_dir}" \
    -n ${num_batches} \
    -p "multisusie_realistic_finemap_batch_" \
    -s ".sh" \

if [ $? -ne 0 ]; then
    echo "ERROR: Python batch splitter failed"
    exit 1
fi

# Verify batch creation
num_batch_files=$(ls "${command_files_dir}" | grep -c "multisusie_realistic_finemap_batch_")
echo "Verified: Created ${num_batch_files} batch files in ${command_files_dir}"

if [ ${num_batch_files} -eq 0 ]; then
    echo "ERROR: No batch files were created"
    exit 1
fi

################################################################################################################
############################################### Load Conda Environment #########################################
################################################################################################################

##### START OF SBCAST AND CONDA-UNPACK #####
#Set variables
NNODES=${SLURM_NNODES}
# Move a copy of the env to the NVMe on each node
echo "copying ${ENV_NAME} to each node in the job"
sbcast -pf ${TAR_DIR}/${TAR_FILE} /mnt/bb/${USER}/${TAR_FILE}
if [ ! "$?" == "0" ]; then
    # CHECK EXIT CODE. When SBCAST fails, it may leave partial files on the compute nodes, and if you continue to launch srun,
    # your application may pick up partially complete shared library files, which would give you confusing errors.
    echo "SBCAST failed!"
    exit 1
fi

# Untar the environment file (only need 1 task per node to do this)
srun -N ${NNODES} --ntasks-per-node 1 mkdir -p /mnt/bb/${USER}/${ENV_NAME}
echo "untaring ${ENV_NAME}"
srun -N ${NNODES} --ntasks-per-node 1 tar -xzf /mnt/bb/${USER}/${TAR_FILE} -C /mnt/bb/${USER}/${ENV_NAME}
##### END OF SBCAST AND CONDA-UNPACK #####

################################################################################################################
####################################### Execute Commands in Parallel #########################################
################################################################################################################

echo ""
echo "Starting flux-based execution of batched MultiSuSiE fine-mapping commands at $(date)"

# Create flux output directory
flux_output_dir="${base_finemap_dir}/flux_output"
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
    for batch in ${command_files_dir}/multisusie_realistic_finemap_batch_*.sh; do
        batch_name=\$(basename \${batch} .sh)
        flux submit -N1 -n1 -c2 \
            --output=${flux_output_dir}/\${batch_name}.log \
            bash \${batch}
        
        ((batch_count++))
        
        # Print progress every 20 batches
        if [ \$((batch_count % 20)) -eq 0 ]; then
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

################################################################################################################
################################## Post-processing and cleanup ##################################################
################################################################################################################

echo ""
echo "Post-processing results..."

# Count successful batch completions by checking for output files
successful_batches=$(ls -1 "${flux_output_dir}"/multisusie_realistic_finemap_batch_*.log 2>/dev/null | wc -l)
echo "Successfully completed batches: ${successful_batches}/${num_batch_files}"

# Check for any failed batches
if [ ${successful_batches} -ne ${num_batch_files} ]; then
    echo "WARNING: Some batches may have failed. Check logs in ${flux_output_dir}"
fi

echo ""
echo "MultiSuSiE fine-mapping workflow completed at $(date)"
echo "Batch files retained in: ${command_files_dir}"
echo "Flux logs available in: ${flux_output_dir}"
echo "Fine-mapping results in: ${base_finemap_dir}"
