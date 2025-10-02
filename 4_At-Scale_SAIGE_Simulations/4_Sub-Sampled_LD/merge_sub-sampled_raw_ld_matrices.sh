#!/bin/bash
#SBATCH -A bif154
#SBATCH -J merge_ld_matrices
#SBATCH -o merge_sub-sampled_raw_ld_matrices.log
#SBATCH -N 20
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
command_generator_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/generate_merge_commands.py"
matrix_combiner_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/combine_ld_matrices.py"
# Set directory locations
input_dirs=(
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/sub-sample_1000"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/sub-sample_5000"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/sub-sample_10000"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/sub-sample_50000"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/sub-sample_100000"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/sub-sample_200000"
    "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/sub-sample_300000"
)

# Output location for commands and flux output
output_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
commands_file="${output_dir}/ld_merge_commands.txt"

# Sample sizes (both AFR and EUR = 400,000) #Actual sample size doesn't matter only the ratio does for weighting
afr_sample_size=400000
eur_sample_size=400000

# Set number of nodes
num_nodes=${SLURM_NNODES}

################################################################################################################
################################## Create output directories and prep commands ##################################
################################################################################################################

echo "Starting LD matrix merge workflow at $(date)"
echo "Input directories: ${input_dirs[@]}"
echo "Output directory: ${output_dir}"
echo "Using ${num_nodes} nodes"

# Ensure weighted subdirectories exist in input directories
for input_dir in "${input_dirs[@]}"; do
    if [ -d "${input_dir}" ]; then
        echo "Checking directory: ${input_dir}"
        mkdir -p "${input_dir}/weighted"

        # List contents for verification
        echo "  AFR matrices: $(find ${input_dir}/AFR -name '*.unphased.vcor1' 2>/dev/null | wc -l)"
        echo "  EUR matrices: $(find ${input_dir}/EUR -name '*.unphased.vcor1' 2>/dev/null | wc -l)"
    else
        echo "WARNING: Input directory does not exist: ${input_dir}"
    fi
done

################################################################################################################
################################## Generate merge commands ###################################################
################################################################################################################

echo "Generating merge commands..."
python "${command_generator_script}" \
    --input-dirs "${input_dirs[@]}" \
    --output-commands "${commands_file}" \
    --afr-sample-size ${afr_sample_size} \
    --eur-sample-size ${eur_sample_size} \
    --combiner-script "${matrix_combiner_script}" \

# Check if command generation was successful
if [ ! -f "${commands_file}" ]; then
    echo "ERROR: Command generation failed. Commands file not found: ${commands_file}"
    exit 1
fi

# Count the number of commands
num_commands=$(grep -v '^#' "${commands_file}" | grep -v '^$' | wc -l)
echo "Generated ${num_commands} merge commands"

if [ ${num_commands} -eq 0 ]; then
    echo "WARNING: No commands were generated. Check input directories and file structure."
    exit 1
fi

################################################################################################################
################################## Execute merge commands using Flux #########################################
################################################################################################################

echo "Starting flux-based execution of merge commands at $(date)"

# Read commands into array (skip comments and empty lines)
mapfile -t cmds < ${commands_file}
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start "flux resource list; mapfile -t cmds < ${commands_file}; for ((i=0; i<\${#cmds[@]}; i++)); do flux submit -N1 -n1 -c56 --output=flux_output/merge_in-sample_\$i.txt \${cmds[\$i]}; done; flux queue drain;"

