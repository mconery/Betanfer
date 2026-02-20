#!/bin/bash

#SBATCH -A bif154
#SBATCH -J infer_2000-trait_ld_pop-specific
#SBATCH -o infer_2000-trait_ld.pop-specific.log
#SBATCH -N 30
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
matrix_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/5_LD_Inference/infer_matrix_from_null_traits.py"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/5_LD_Inference/create_beta_inference_commands.py"

# Set directory locations for EUR population
eur_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_EUR_sumstats"
eur_variant_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/EUR"
eur_output_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_2000_random_all/EUR"
eur_commands_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_2000/EUR/eur_matrix_commands.txt"

# Set directory locations for AFR population
afr_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_AFR_sumstats"
afr_variant_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/AFR"
afr_output_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_2000_random_all/AFR"
afr_commands_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/beta-inferred_2000_random_all/AFR/afr_matrix_commands.txt"

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Create output directories if they don't exist
mkdir -p "$eur_output_dir"
mkdir -p "$afr_output_dir"
mkdir -p flux_output

################################################################################################################
################################## Prep Commands for EUR population #########################################
################################################################################################################

echo "================================================"
echo "STARTING LD MATRIX INFERENCE PIPELINE - EUR"
echo "================================================"
echo "Start time: $(date)"
echo "EUR Input variant directory: $eur_variant_dir"
echo "EUR Output directory: $eur_output_dir"
echo "EUR Indexed sumstats directory: $eur_scratch_dir"
echo "Number of nodes: $num_nodes"
echo "================================================"

# Generate EUR commands file using Python script
echo "Generating EUR commands file..."
python $command_script \
--variant_dir "$eur_variant_dir" \
--output_dir "$eur_output_dir" \
--commands_file "$eur_commands_file" \
--matrix_script "$matrix_script" \
--indexed_dir "$eur_scratch_dir" \
--n_traits 2000 \
--variant_col MarkerID \
--beta_col BETA \
--device cuda

# Check if EUR command generation was successful
if [[ $? -ne 0 ]]; then
echo "ERROR: Failed to generate EUR commands file!"
exit 1
fi

# Verify EUR commands file exists and has content
if [[ ! -f "$eur_commands_file" ]] || [[ ! -s "$eur_commands_file" ]]; then
echo "ERROR: EUR commands file is missing or empty: $eur_commands_file"
exit 1
fi

eur_num_commands=$(wc -l < "$eur_commands_file")
echo "Generated $eur_num_commands EUR commands"

################################################################################################################
################################## Prep Commands for AFR population #########################################
################################################################################################################

echo "================================================"
echo "STARTING LD MATRIX INFERENCE PIPELINE - AFR"
echo "================================================"
echo "AFR Input variant directory: $afr_variant_dir"
echo "AFR Output directory: $afr_output_dir"
echo "AFR Indexed sumstats directory: $afr_scratch_dir"
echo "================================================"

# Generate AFR commands file using Python script
echo "Generating AFR commands file..."
python $command_script \
--variant_dir "$afr_variant_dir" \
--output_dir "$afr_output_dir" \
--commands_file "$afr_commands_file" \
--matrix_script "$matrix_script" \
--indexed_dir "$afr_scratch_dir" \
--n_traits 2000 \
--variant_col MarkerID \
--beta_col BETA \
--device cuda

# Check if AFR command generation was successful
if [[ $? -ne 0 ]]; then
echo "ERROR: Failed to generate AFR commands file!"
exit 1
fi

# Verify AFR commands file exists and has content
if [[ ! -f "$afr_commands_file" ]] || [[ ! -s "$afr_commands_file" ]]; then
echo "ERROR: AFR commands file is missing or empty: $afr_commands_file"
exit 1
fi

afr_num_commands=$(wc -l < "$afr_commands_file")
echo "Generated $afr_num_commands AFR commands"

################################################################################################################
################################## Execute EUR Commands in Parallel ########################################
################################################################################################################

echo "Starting EUR parallel execution..."

# Read EUR commands into array
mapfile -t eur_cmds < "$eur_commands_file"

# Launch flux and execute EUR commands in parallel
echo "Launching flux for EUR with $num_nodes nodes..."
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start "flux resource list; mapfile -t eur_cmds < ${eur_commands_file}; for ((i=0; i<\${#eur_cmds[@]}; i++)); do flux submit -N1 -n1 -c7 --gpus-per-task=1 --output=flux_output/beta-infer_2000_eur_\$i.txt \${eur_cmds[\$i]}; done; flux queue drain;"

echo "All EUR LD matrix inference jobs submitted and completed."

################################################################################################################
################################## Execute AFR Commands in Parallel ########################################
################################################################################################################

echo "Starting AFR parallel execution..."

# Read AFR commands into array
mapfile -t afr_cmds < "$afr_commands_file"

# Launch flux and execute AFR commands in parallel
echo "Launching flux for AFR with $num_nodes nodes..."
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start "flux resource list; mapfile -t afr_cmds < ${afr_commands_file}; for ((i=0; i<\${#afr_cmds[@]}; i++)); do flux submit -N1 -n1 -c7 --gpus-per-task=1 --output=flux_output/beta-infer_2000_afr_\$i.txt \${afr_cmds[\$i]}; done; flux queue drain;"

echo "All AFR LD matrix inference jobs submitted and completed."

################################################################################################################
######################################## Final Summary ######################################################
################################################################################################################

echo "================================================"
echo "LD MATRIX INFERENCE PIPELINE SUMMARY"
echo "================================================"
echo "Completion time: $(date)"
echo "EUR commands processed: $eur_num_commands"
echo "AFR commands processed: $afr_num_commands"
echo "Total commands processed: $((eur_num_commands + afr_num_commands))"
echo ""
echo "Output directories:"
echo "  EUR LD matrices: $eur_output_dir"
echo "  AFR LD matrices: $afr_output_dir"
echo ""
echo "Command files:"
echo "  EUR commands: $eur_commands_file"
echo "  AFR commands: $afr_commands_file"
echo ""
echo "Flux output logs: flux_output/"
echo "================================================"
echo "LD matrix inference completed for both populations."