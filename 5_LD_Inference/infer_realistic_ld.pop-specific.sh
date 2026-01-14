#!/bin/bash

#SBATCH -A bif154
#SBATCH -J infer_realistic_ld_pop-specific
#SBATCH -o infer_realistic_ld.pop-specific.log
#SBATCH -N 64
#SBATCH -p batch
#SBATCH -t 2:00:00

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
pval_thresholds=(5e-8 5e-6 5e-4) #Need to check 0 later
buffer_sizes=(0 100000 250000 500000)

# Directory for all command files
commands_dir="${base_ld_dir}/command_files"
mkdir -p "$commands_dir"

# Batch configuration
commands_per_batch=20

# Master commands file
master_commands_file="${commands_dir}/master_realistic_pop_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Track total commands generated
total_commands=0
total_command_files=0

################################################################################################################
################################## Generate Commands for All Combinations ####################################
################################################################################################################

echo "================================================"
echo "STARTING MULTI-TRAIT LD MATRIX INFERENCE PIPELINE"
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


        # EUR population
        echo "  → EUR population..."
        eur_output_dir="${base_ld_dir}/beta-inferred_realistic_${pval_str}_${buffer_kb}/EUR"
        eur_commands_file="${commands_dir}/realistic_pval${pval_str}_buffer${buffer_kb}_eur_commands.txt"

        mkdir -p "$eur_output_dir"

        python "$command_script"             --variant_dir "$eur_variant_dir"             --output_dir "$eur_output_dir"             --commands_file "$eur_commands_file"             --matrix_script "$matrix_script"             --indexed_dir "$eur_scratch_dir"             --pval_threshold "$pval_thresh"             --buffer "$buffer_size"             --pval_col "p.value"             --variant_col MarkerID             --beta_col BETA             --device cuda

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to generate EUR commands for p-value: $pval_thresh, buffer: ${buffer_size}bp"
            exit 1
        fi

        if [[ ! -f "$eur_commands_file" ]] || [[ ! -s "$eur_commands_file" ]]; then
            echo "ERROR: EUR commands file is missing or empty: $eur_commands_file"
            exit 1
        fi

        eur_num_commands=$(wc -l < "$eur_commands_file")
        echo "    Generated $eur_num_commands EUR commands"
        ((total_commands += eur_num_commands))
        ((total_command_files++))

        # Append to master commands file
        grep -v "^#" "$eur_commands_file" >> "$master_commands_file"

        # AFR population
        echo "  → AFR population..."
        afr_output_dir="${base_ld_dir}/beta-inferred_realistic_${pval_str}_${buffer_kb}/AFR"
        afr_commands_file="${commands_dir}/realistic_pval${pval_str}_buffer${buffer_kb}_afr_commands.txt"

        mkdir -p "$afr_output_dir"

        python "$command_script"             --variant_dir "$afr_variant_dir"             --output_dir "$afr_output_dir"             --commands_file "$afr_commands_file"             --matrix_script "$matrix_script"             --indexed_dir "$afr_scratch_dir"             --pval_threshold "$pval_thresh"             --buffer "$buffer_size"             --pval_col "p.value"             --variant_col MarkerID             --beta_col BETA             --device cuda

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to generate AFR commands for p-value: $pval_thresh, buffer: ${buffer_size}bp"
            exit 1
        fi

        if [[ ! -f "$afr_commands_file" ]] || [[ ! -s "$afr_commands_file" ]]; then
            echo "ERROR: AFR commands file is missing or empty: $afr_commands_file"
            exit 1
        fi

        afr_num_commands=$(wc -l < "$afr_commands_file")
        echo "    Generated $afr_num_commands AFR commands"
        ((total_commands += afr_num_commands))
        ((total_command_files++))

        # Append to master commands file
        grep -v "^#" "$afr_commands_file" >> "$master_commands_file"

    done
done

################################################################################################################
################################## Summary and Final Status ###################################################
################################################################################################################

echo ""
echo "================================================"
echo "COMMAND FILE GENERATION COMPLETE"
echo "================================================"
echo "End time: $(date)"
echo "Total command files generated: $total_command_files"
echo "Total commands in master file: $(wc -l < "$master_commands_file")"
echo "Master commands file: $master_commands_file"
echo "Command directory: $commands_dir"
echo "================================================"

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

# Clean up any existing pop_infer batch files
rm -f "${commands_dir}"/pop_realistic_batch_*.sh

# Use Python splitter for fast parallel batch creation
python "${python_splitter}" \
    -i "${master_commands_file}" \
    -o "${commands_dir}" \
    -b ${commands_per_batch} \
    -p "pop_realistic_batch_" \
    -s ".sh" 

if [ $? -ne 0 ]; then
    echo "ERROR: Python batch splitter failed"
    exit 1
fi

# Verify batch creation
num_batch_files=$(ls "${commands_dir}" | grep -c "^pop_realistic_batch")
echo "Verified: Created ${num_batch_files} batch files in ${commands_dir}"

if [ ${num_batch_files} -eq 0 ]; then
    echo "ERROR: No batch files were created"
    exit 1
fi

################################################################################################################
####################################### Execute Commands in Parallel #########################################
################################################################################################################

echo ""
echo "Starting flux-based execution of batched population-specific commands at $(date)"

# Create flux output directory
mkdir -p "${base_ld_dir}/flux_output"
mkdir -p "${base_ld_dir}/jobstate/real_pop"
rm "${base_ld_dir}/jobstate/real_pop/*"

# Start Flux and submit batch jobs
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start --broker-opts=-Sstatedir=${base_ld_dir}/jobstate/real_pop bash -c "
    # Display Flux resource allocation
    echo \"Flux instance started at \$(date)\"
    echo \"Available resources:\"
    flux resource list
    echo \"\"
    
    # Submit all batch files to Flux
    echo \"Submitting ${num_batch_files} batches to Flux queue...\"
    batch_count=0
    for batch in ${commands_dir}/pop_realistic_batch_*.sh; do
        batch_name=\$(basename \${batch} .sh)
        flux submit -N1 -n1 -c56 --gpus-per-task=1 \
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
successful_batches=$(ls -1 "${base_ld_dir}"/flux_output/pop_realistic_batch_*.log 2>/dev/null | wc -l)
echo "Successfully completed batches: ${successful_batches}/${num_batch_files}"

# Check for any failed batches (optional - you can add error checking logic here)
if [ ${successful_batches} -ne ${num_batch_files} ]; then
    echo "WARNING: Some batches may have failed. Check logs in ${base_ld_dir}/flux_output/"
fi

echo ""
echo "Population-specific LD matrix inference workflow completed at $(date)"
echo "Batch files retained in: ${commands_dir}"
echo "Flux logs available in: ${base_ld_dir}/flux_output"
