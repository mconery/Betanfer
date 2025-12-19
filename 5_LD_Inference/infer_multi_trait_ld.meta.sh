#!/bin/bash

#SBATCH -A bif154
#SBATCH -J infer_multi_trait_ld_meta
#SBATCH -o infer_multi_trait_ld.meta.log
#SBATCH -N 1882
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
matrix_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/infer_matrix_from_null_traits.py"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/create_beta_inference_commands.py"
python_splitter="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/split_commands_parallel.py"

# Base directories
base_ld_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# Meta-analysis specific directories
fixed_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats_fixed"
random_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats_random"

# Variant directory (same for all meta-analyses)
variant_dir="${base_ld_dir}/raw_in-sample/weighted"

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Define meta-analysis types, trait types and number of traits
meta_types=("fixed" "random")
trait_types=("all" "continuous" "binary")
num_traits_list=(50 100 200 300 400 500 750 1000 1250 1500 1750 2000)

# Directory for all command files
commands_dir="${base_ld_dir}/command_files"
mkdir -p "$commands_dir"

# Batch configuration
commands_per_batch=5

# Master commands file
master_commands_file="${commands_dir}/master_meta_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Track total commands generated
total_commands=0
total_command_files=0

################################################################################################################
################################## Generate Commands for All Combinations ####################################
################################################################################################################

echo "================================================"
echo "STARTING META-ANALYSIS LD MATRIX INFERENCE PIPELINE"
echo "================================================"
echo "Start time: $(date)"
echo "Meta-analysis types: ${meta_types[@]}"
echo "Trait types: ${trait_types[@]}"
echo "Number of traits: ${num_traits_list[@]}"
echo "Command files will be written to: $commands_dir"
echo "Master commands file: $master_commands_file"
echo "Batch size: ${commands_per_batch} commands per batch"
echo "================================================"

# Loop over meta-analysis types
for meta_type in "${meta_types[@]}"; do
    echo ""
    echo "====================================="
    echo "Processing meta-analysis type: $meta_type"
    echo "====================================="

    # Set indexed directory based on meta-analysis type
    if [[ $meta_type == "fixed" ]]; then
        indexed_dir="$fixed_scratch_dir"
    else
        indexed_dir="$random_scratch_dir"
    fi

    # Loop over trait types
    for trait_type in "${trait_types[@]}"; do
        echo ""
        echo "Processing trait type: $trait_type"

        # Loop over number of traits
        for n_traits in "${num_traits_list[@]}"; do
            echo ""
            echo "Generating commands for $n_traits traits, meta-type: $meta_type, trait-type: $trait_type"

            # Set output directory
            output_dir="${base_ld_dir}/beta-inferred_${n_traits}_${meta_type}_${trait_type}/Meta"
            commands_file="${commands_dir}/${meta_type}_${n_traits}_${trait_type}_meta_commands.txt"
            mkdir -p "$output_dir"

            # Generate commands for this configuration
            python "$command_script" \
                --variant_dir "$variant_dir" \
                --output_dir "$output_dir" \
                --commands_file "$commands_file" \
                --matrix_script "$matrix_script" \
                --indexed_dir "$indexed_dir" \
                --n_traits "$n_traits" \
                --trait_type "$trait_type" \
                --variant_col '"#rs_number"' \
                --beta_col beta \
                --device cuda

            if [[ $? -ne 0 ]]; then
                echo "ERROR: Failed to generate commands for $n_traits traits, meta-type: $meta_type, trait-type: $trait_type"
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

            # Append to master commands file
            grep -v "^#" "$commands_file" >> "$master_commands_file"

        done
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

# Clean up any existing meta_infer batch files
rm -f "${commands_dir}"/meta_infer_batch_*.sh

# Use Python splitter for fast parallel batch creation
python "${python_splitter}" \
    -i "${master_commands_file}" \
    -o "${commands_dir}" \
    -b ${commands_per_batch} \
    -p "meta_infer_batch_" \
    -s ".sh" 

if [ $? -ne 0 ]; then
    echo "ERROR: Python batch splitter failed"
    exit 1
fi

# Verify batch creation
num_batch_files=$(ls "${commands_dir}" | grep -c "meta_infer_batch")
echo "Verified: Created ${num_batch_files} batch files in ${commands_dir}"

################################################################################################################
####################################### Execute Commands in Parallel #########################################
################################################################################################################

echo ""
echo "Starting flux-based execution of batched meta-analysis commands at $(date)"

# Create flux output directory
mkdir -p "${base_ld_dir}/flux_output"

# Start Flux and submit batch jobs
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start bash -c "
    # Display Flux resource allocation
    echo \"Flux instance started at \$(date)\"
    echo \"Available resources:\"
    flux resource list
    echo \"\"
    
    # Submit all batch files to Flux
    echo \"Submitting ${num_batch_files} batches to Flux queue...\"
    batch_count=0
    for batch in ${commands_dir}/meta_infer_batch_*.sh; do
        batch_name=\$(basename \${batch} .sh)
        flux submit -N1 -n1 -c56 --gpus-per-task=8 \
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
successful_batches=$(ls -1 "${base_ld_dir}"/flux_output/meta_infer_batch_*.log 2>/dev/null | wc -l)
echo "Successfully completed batches: ${successful_batches}/${num_batch_files}"

# Check for any failed batches (optional - you can add error checking logic here)
if [ ${successful_batches} -ne ${num_batch_files} ]; then
    echo "WARNING: Some batches may have failed. Check logs in ${base_ld_dir}/flux_output/"
fi

echo ""
echo "Meta-analysis LD matrix inference workflow completed at $(date)"
echo "Batch files retained in: ${commands_dir}"
echo "Flux logs available in: ${base_ld_dir}/flux_output"
