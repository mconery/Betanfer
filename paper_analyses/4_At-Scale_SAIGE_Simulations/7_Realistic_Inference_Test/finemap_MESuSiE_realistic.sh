#!/bin/bash

#SBATCH -A bif154
#SBATCH -J MESuSiE_finemap_realistic
#SBATCH -o finemap_MESuSiE_realistic.log
#SBATCH -N 40
#SBATCH -p batch
#SBATCH -t 2:00:00

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

# Set script locations
map_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/MESuSiE/finemap_MESuSiE.R"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/MESuSiE/generate_mesusie_commands.py"
python_splitter="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/5_LD_Inference/split_commands_parallel.py"

# Set environment location
map_enviro="/lustre/orion/bif154/proj-shared/mconery/tools/environments/carma"

# Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"

# Base output directory for fine-mapping results
base_finemap_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"

# Command files directory
command_files_dir="${base_finemap_dir}/command_files"
mkdir -p "$command_files_dir"

# Master commands file
master_commands_file="${command_files_dir}/master_mesusie_realistic_commands.txt"
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
echo "STARTING MESUSIE FINE-MAPPING PIPELINE (All Configurations)"
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
	matrix_types="beta-inferred_realistic_0_0kb,beta-inferred_realistic_5e-4_0kb,beta-inferred_realistic_5e-4_100kb,beta-inferred_realistic_5e-4_250kb,beta-inferred_realistic_5e-4_500kb,beta-inferred_realistic_5e-6_0kb,beta-inferred_realistic_5e-6_100kb,beta-inferred_realistic_5e-6_250kb,beta-inferred_realistic_5e-6_500kb,beta-inferred_realistic_5e-8_0kb,beta-inferred_realistic_5e-8_100kb,beta-inferred_realistic_5e-8_250kb,beta-inferred_realistic_5e-8_500kb"        

        # Define output directory and command file
        if [[ $h_label == "5e-5" ]]; then
            h_dir="fixed_homo_5e-5"
        else
            h_dir="fixed_homo_3e-4"
        fi

        out_dir="${base_finemap_dir}/${h_dir}/mesusie_results"
        command_file="${command_files_dir}/mesusie_realistic_${h_label}_commands.txt"

        mkdir -p "$out_dir"

        # Generate commands for this configuration
        echo "  Running command generation script..."
        python "$command_script" \
            --wrapper-script "$map_script" \
            --gwas-dir "$sumstats_dir" \
            --ld-dir "$matrix_dir" \
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

# Clean up any existing mesusie batch files
rm -f "${command_files_dir}"/mesusie_realistic_finemap_batch_*.sh

#Calculate number of batches
num_batches=$(( num_nodes * 28 ))
# Use Python splitter for fast parallel batch creation
python "${python_splitter}" \
    -i "${master_commands_file}" \
    -o "${command_files_dir}" \
    -n ${num_batches} \
    -p "mesusie_realistic_finemap_batch_" \
    -s ".sh" \

if [ $? -ne 0 ]; then
    echo "ERROR: Python batch splitter failed"
    exit 1
fi

# Verify batch creation
num_batch_files=$(ls "${command_files_dir}" | grep -c "mesusie_realistic_finemap_batch_")
echo "Verified: Created ${num_batch_files} batch files in ${command_files_dir}"

if [ ${num_batch_files} -eq 0 ]; then
    echo "ERROR: No batch files were created"
    exit 1
fi

################################################################################################################
####################################### Execute Commands in Parallel #########################################
################################################################################################################

echo ""
echo "Starting flux-based execution of batched MESuSiE fine-mapping commands at $(date)"

# Create flux output directory
flux_output_dir="${base_finemap_dir}/flux_output"
mkdir -p "$flux_output_dir"

# Activate environment before launching flux
source activate $map_enviro

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
    for batch in ${command_files_dir}/mesusie_realistic_finemap_batch_*.sh; do
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
successful_batches=$(ls -1 "${flux_output_dir}"/mesusie_realistic_finemap_batch_*.log 2>/dev/null | wc -l)
echo "Successfully completed batches: ${successful_batches}/${num_batch_files}"

# Check for any failed batches
if [ ${successful_batches} -ne ${num_batch_files} ]; then
    echo "WARNING: Some batches may have failed. Check logs in ${flux_output_dir}"
fi

echo ""
echo "MESuSiE fine-mapping workflow completed at $(date)"
echo "Batch files retained in: ${command_files_dir}"
echo "Flux logs available in: ${flux_output_dir}"
echo "Fine-mapping results in: ${base_finemap_dir}"
