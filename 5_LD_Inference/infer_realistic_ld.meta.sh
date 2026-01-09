#!/bin/bash

#SBATCH -A bif154
#SBATCH -J infer_realistic_ld_meta
#SBATCH -o infer_realistic_ld.meta.log
#SBATCH -N 64
#SBATCH -p extended
#SBATCH -t 12:00:00

################################################################################################################
###################################### Load Modules and Environment ############################################
################################################################################################################

module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux

source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/infersparsepro

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set script locations
matrix_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/infer_matrix_with_pval_filter.py"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/create_realistic_beta_inference_commands.py"
python_splitter="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/split_commands_parallel.py"

# Base directories
base_ld_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# EUR population directories
eur_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_realistic_EUR_sumstats"
eur_variant_dir="${base_ld_dir}/raw_in-sample/EUR"

# AFR population directories
afr_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_realistic_AFR_sumstats"
afr_variant_dir="${base_ld_dir}/raw_in-sample/AFR"

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Define p-value thresholds and buffer sizes
pval_thresholds=(5e-6 5e-8) #Need to check 0 later
buffer_sizes=(0 100000 250000 500000)

# Directory for all command files
commands_dir="${base_ld_dir}/command_files"
mkdir -p "$commands_dir"

# Batch configuration
commands_per_batch=20

# Master commands file
master_commands_file="${commands_dir}/master_realistic_meta_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Track total commands generated
total_commands=0
total_command_files=0

################################################################################################################
################################## Generate Commands for All Combinations ####################################
################################################################################################################

echo "================================================"
echo "STARTING REALISTIC TRAIT LD MATRIX INFERENCE PIPELINE"
echo "================================================"
echo "Start time: $(date)"
echo "P-value thresholds: ${pval_thresholds[@]}"
echo "Buffer sizes (bp): ${buffer_sizes[@]}"
echo "Command files will be written to: $commands_dir"
echo "Master commands file: $master_commands_file"
echo "Batch size: ${commands_per_batch} commands per batch"
echo "Number of nodes: $num_nodes"
echo "Indexed sumstats: $indexed_dir"
echo "================================================"

# Loop over p-value thresholds
for pval_thresh in "${pval_thresholds[@]}"; do

    echo ""
    echo "====================================="
    echo "Processing p-value threshold: $pval_thresh"
    echo "====================================="

    # Loop over buffer sizes
    for buffer_size in "${buffer_sizes[@]}"; do

        echo ""
        echo "Generating commands for p-value: $pval_thresh, buffer: ${buffer_size}bp"

        # Format buffer size for output directory (convert to kb)
        if [[ $buffer_size -eq 0 ]]; then
            buffer_kb="0kb"
        else
            buffer_kb="$((buffer_size / 1000))kb"
        fi

        # Format p-value for output directory
        if [[ $pval_thresh == "0" ]]; then
            pval_str="0"
        else
            pval_str="$pval_thresh"
        fi

        # Set output directory following naming convention
        output_dir="${base_ld_dir}/beta-inferred_realistic_${pval_str}_${buffer_kb}/Meta"
        commands_file="${commands_dir}/realistic_pval${pval_str}_buffer${buffer_kb}_commands.txt"

        mkdir -p "$output_dir"

        # Generate commands for this configuration
        echo "Creating commands file: $commands_file"
        echo "Output directory: $output_dir"

        python "$command_script" \
            --variant_dir "$variant_dir" \
            --output_dir "$output_dir" \
            --commands_file "$commands_file" \
            --matrix_script "$matrix_script" \
            --indexed_dir "$indexed_dir" \
            --pval_threshold "$pval_thresh" \
            --buffer "$buffer_size" \
            --pval_col p-value \
            --variant_col '#rs_number' \
            --beta_col beta \
            --device cuda

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to generate commands for p-value: $pval_thresh, buffer: ${buffer_size}bp"
            exit 1
        fi

        if [[ ! -f "$commands_file" ]] || [[ ! -s "$commands_file" ]]; then
            echo "ERROR: Commands file is missing or empty: $commands_file"
            exit 1
        fi

        num_commands=$(wc -l < "$commands_file")
        echo "Generated $num_commands commands"

        ((total_commands += num_commands))
        ((total_command_files++))

        # Add commands file reference to master file
        cat "$commands_file" >> "$master_commands_file"

    done

done

################################################################################################################
################################## Create batched command files ################################################
################################################################################################################

echo ""
echo "Creating batched command files using parallel Python splitter..."

# Clean up any existing meta_infer batch files
rm -f "${commands_dir}"/meta_realistic_batch_*.sh

# Use Python splitter for fast parallel batch creation
python "${python_splitter}" \
    -i "${master_commands_file}" \
    -o "${commands_dir}" \
    -b ${commands_per_batch} \
    -p "meta_realistic_batch_" \
    -s ".sh" 

if [ $? -ne 0 ]; then
    echo "ERROR: Python batch splitter failed"
    exit 1
fi

# Verify batch creation
num_batch_files=$(ls "${commands_dir}" | grep -c "meta_realistic_batch")
echo "Verified: Created ${num_batch_files} batch files in ${commands_dir}"

################################################################################################################
####################################### Execute Commands in Parallel #########################################
################################################################################################################

echo ""
echo "Starting flux-based execution of batched meta-analysis commands at $(date)"

# Create flux output directory
mkdir -p "${base_ld_dir}/flux_output"

# Start Flux and submit batch jobs
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start --broker-opts=-Sstatedir=/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/jobstate bash -c "
    # Display Flux resource allocation
    echo \"Flux instance started at \$(date)\"
    echo \"Available resources:\"
    flux resource list
    echo \"\"
    
    # Submit all batch files to Flux
    echo \"Submitting ${num_batch_files} batches to Flux queue...\"
    batch_count=0
    for batch in ${commands_dir}/meta_realistic_batch_*.sh; do
        batch_name=\$(basename \${batch} .sh)
        flux submit -N1 -n1 -c7 --gpus-per-task=1 \
            --output=${base_ld_dir}/flux_output/\${batch_name}.log \
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

################################################################################################################
################################## Post-processing and cleanup ##################################################
################################################################################################################

echo ""
echo "Post-processing results..."

# Count successful batch completions by checking for output files
successful_batches=$(ls -1 "${base_ld_dir}"/flux_output/meta_realistic_batch_*.log 2>/dev/null | wc -l)
echo "Successfully completed batches: ${successful_batches}/${num_batch_files}"

# Check for any failed batches (optional - you can add error checking logic here)
if [ ${successful_batches} -ne ${num_batch_files} ]; then
    echo "WARNING: Some batches may have failed. Check logs in ${base_ld_dir}/flux_output/"
fi

echo ""
echo "Meta-analysis LD matrix inference workflow completed at $(date)"
echo "Batch files retained in: ${commands_dir}"
echo "Flux logs available in: ${base_ld_dir}/flux_output"
