#!/bin/bash
#SBATCH -A bif154
#SBATCH -J merge_ld_matrices
#SBATCH -o merge_ld_matrices.log
#SBATCH -N 1882
#SBATCH -p batch
#SBATCH -t 40:00

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
command_generator_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/generate_merge_commands.py"
matrix_combiner_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/combine_ld_matrices.py"
# Set directory locations
input_dirs=(
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/PC-corrected"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_50_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_100_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_200_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_300_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_400_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_500_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_750_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1000_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1250_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1500_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1750_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_2000_random_all"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_50_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_100_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_200_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_300_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_400_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_500_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_750_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1000_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1250_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1500_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1750_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_2000_random_binary"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_50_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_100_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_200_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_300_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_400_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_500_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_750_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1000_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1250_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1500_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_1750_random_continuous"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_2000_random_continuous"
)

# Output location for commands and flux output
output_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
commands_file="${output_dir}/ld_merge_commands.txt"

# Batch configuration
batch_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/command_files"
commands_per_batch=25

# Sample sizes (both AFR and EUR = 400,000)
afr_sample_size=400000
eur_sample_size=400000

# Set number of nodes
num_nodes=${SLURM_NNODES}

################################################################################################################
################################## Create output directories and prep commands ##################################
################################################################################################################

#echo "Starting LD matrix merge workflow at $(date)"
#echo "Input directories: ${input_dirs[@]}"
#echo "Output directory: ${output_dir}"
#echo "Using ${num_nodes} nodes"
#
## Ensure weighted subdirectories exist in input directories
#for input_dir in "${input_dirs[@]}"; do
#    if [ -d "${input_dir}" ]; then
#        echo "Checking directory: ${input_dir}"
#        mkdir -p "${input_dir}/weighted"
#
#        # List contents for verification
#        echo "  AFR matrices: $(find ${input_dir}/AFR -name '*.unphased.vcor1' 2>/dev/null | wc -l)"
#        echo "  EUR matrices: $(find ${input_dir}/EUR -name '*.unphased.vcor1' 2>/dev/null | wc -l)"
#    else
#        echo "WARNING: Input directory does not exist: ${input_dir}"
#    fi
#done
#
#echo ""

################################################################################################################
################################## Generate merge commands ###################################################
################################################################################################################

#echo "Generating merge commands..."
#python "${command_generator_script}" \
#    --input-dirs "${input_dirs[@]}" \
#    --output-commands "${commands_file}" \
#    --afr-sample-size ${afr_sample_size} \
#    --eur-sample-size ${eur_sample_size} \
#    --combiner-script "${matrix_combiner_script}"
#
## Check if command generation was successful
#if [ ! -f "${commands_file}" ]; then
#    echo "ERROR: Command generation failed. Commands file not found: ${commands_file}"
#    exit 1
#fi
#
## Count the number of commands
#num_commands=$(grep -v '^#' "${commands_file}" | grep -v '^$' | wc -l)
#echo "Generated ${num_commands} merge commands"
#
#if [ ${num_commands} -eq 0 ]; then
#    echo "WARNING: No commands were generated. Check input directories and file structure."
#    exit 1
#fi

################################################################################################################
################################## Create batched command files ################################################
################################################################################################################

#echo ""
#echo "Creating batched command files..."
#
## Create batch directory
#mkdir -p "${batch_dir}"
#
## Clean up any existing batch files
#rm -f "${batch_dir}"/batch_*.sh
#
## Split commands into batch files (skip comments and empty lines)
#grep -v '^#' "${commands_file}" | grep -v '^$' | \
#    split -l ${commands_per_batch} -d --additional-suffix=.sh - "${batch_dir}/batch_"
#
# Verify batch creation
num_batch_files=$(ls -1 "${batch_dir}"/batch_*.sh 2>/dev/null | wc -l)
echo "Created ${num_batch_files} batch files in ${batch_dir}"

if [ ${num_batch_files} -eq 0 ]; then
    echo "ERROR: No batch files were created"
    exit 1
fi

################################################################################################################
################################## Execute merge commands using Flux #########################################
################################################################################################################

echo ""
echo "Starting flux-based execution of batched merge commands at $(date)"

# Create flux output directory
mkdir -p "${output_dir}/flux_output"

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
    for batch in ${batch_dir}/batch_*.sh; do
        batch_name=\$(basename \${batch} .sh)
        flux submit -N1 -n1 -c56 \
            --output=${output_dir}/flux_output/\${batch_name}.log \
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
