#!/bin/bash

#SBATCH -A bif154
#SBATCH -J infer_2000-trait_ld
#SBATCH -o infer_2000-trait_ld.log
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
matrix_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/infer_matrix_from_null_traits.py"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/create_beta_inference_commands.py"

# Set directory locations
scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats_fixed"
variant_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/weighted"
output_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_2000/Fixed_Meta"
commands_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_2000/fixed_meta_matrix_commands.txt"

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

################################################################################################################
################################## Prep Commands and launch LD matrix inference ##############################
################################################################################################################

echo "================================================"
echo "STARTING LD MATRIX INFERENCE PIPELINE"
echo "================================================"
echo "Start time: $(date)"
echo "Input variant directory: $variant_dir"
echo "Output directory: $output_dir"
echo "Indexed sumstats directory: $scratch_dir"
echo "Number of nodes: $num_nodes"
echo "================================================"

# Generate commands file using Python script
echo "Generating commands file..."
python $command_script \
    --variant_dir "$variant_dir" \
    --output_dir "$output_dir" \
    --commands_file "$commands_file" \
    --matrix_script "$matrix_script" \
    --indexed_dir "$scratch_dir" \
    --n_traits 2000 \
    --variant_col '#rs_number' \
    --beta_col beta \
    --device cuda

# Check if command generation was successful
if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to generate commands file!"
    exit 1
fi

# Verify commands file exists and has content
if [[ ! -f "$commands_file" ]] || [[ ! -s "$commands_file" ]]; then
    echo "ERROR: Commands file is missing or empty: $commands_file"
    exit 1
fi

num_commands=$(wc -l < "$commands_file")
echo "Generated $num_commands commands"

################################################################################################################
################################## Execute Commands in Parallel ###############################################
################################################################################################################

echo "Starting parallel execution..."

# Read commands into array  
mapfile -t cmds < "$commands_file"

# Launch flux and execute commands in parallel
echo "Launching flux with $num_nodes nodes..."
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start "flux resource list; mapfile -t cmds < ${commands_file}; for ((i=0; i<\${#cmds[@]}; i++)); do flux submit -N1 -n1 -c7 --gpus-per-task=1 --output=flux_output/beta-infer_2000_fixed_\$i.txt \${cmds[\$i]}; done; flux queue drain;"

echo "All LD matrix inference jobs submitted and completed."