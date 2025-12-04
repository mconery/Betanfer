#!/bin/bash

#SBATCH -A bif154
#SBATCH -J infer_multi_trait_ld_pop-specific
#SBATCH -o infer_multi_trait_ld.pop-specific.log
#SBATCH -N 10
#SBATCH -p extended
#SBATCH -t 24:00:00

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

# Base directories
base_ld_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# EUR population directories
eur_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_EUR_sumstats"
eur_variant_dir="${base_ld_dir}/raw_in-sample/EUR"

# AFR population directories
afr_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_AFR_sumstats"
afr_variant_dir="${base_ld_dir}/raw_in-sample/AFR"

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Define trait types and number of traits
trait_types=("all" "continuous" "binary")
num_traits_list=(50 100 200 300 400 500 750 1000)

# Directory for all command files
commands_dir="${base_ld_dir}/command_files"
mkdir -p "$commands_dir"

# Master commands file
master_commands_file="${commands_dir}/master_commands.txt"
> "$master_commands_file"  # Clear the file if it exists

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
echo "Trait types: ${trait_types[@]}"
echo "Number of traits: ${num_traits_list[@]}"
echo "Command files will be written to: $commands_dir"
echo "Master commands file: $master_commands_file"
echo "================================================"

# Loop over trait types
for trait_type in "${trait_types[@]}"; do
    echo ""
    echo "===================================="
    echo "Processing trait type: $trait_type"
    echo "===================================="

    # Loop over number of traits
    for n_traits in "${num_traits_list[@]}"; do
        echo ""
        echo "Generating commands for $n_traits traits, type: $trait_type"

        # EUR population
        echo "  → EUR population..."
        eur_output_dir="${base_ld_dir}/beta-inferred_${n_traits}_random_${trait_type}/EUR"
        eur_commands_file="${commands_dir}/eur_${n_traits}_${trait_type}_commands.txt"

        mkdir -p "$eur_output_dir"

        python "$command_script"             --variant_dir "$eur_variant_dir"             --output_dir "$eur_output_dir"             --commands_file "$eur_commands_file"             --matrix_script "$matrix_script"             --indexed_dir "$eur_scratch_dir"             --n_traits "$n_traits"             --trait_type "$trait_type"             --variant_col MarkerID             --beta_col BETA             --device cuda

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to generate EUR commands for $n_traits traits, type: $trait_type"
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
        afr_output_dir="${base_ld_dir}/beta-inferred_${n_traits}_random_${trait_type}/AFR"
        afr_commands_file="${commands_dir}/afr_${n_traits}_${trait_type}_commands.txt"

        mkdir -p "$afr_output_dir"

        python "$command_script"             --variant_dir "$afr_variant_dir"             --output_dir "$afr_output_dir"             --commands_file "$afr_commands_file"             --matrix_script "$matrix_script"             --indexed_dir "$afr_scratch_dir"             --n_traits "$n_traits"             --trait_type "$trait_type"             --variant_col MarkerID             --beta_col BETA             --device cuda

        if [[ $? -ne 0 ]]; then
            echo "ERROR: Failed to generate AFR commands for $n_traits traits, type: $trait_type"
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
####################################### Execute Commands in Parallel ########################################
################################################################################################################

echo "Starting AFR parallel execution..."

# Read AFR commands into array
mapfile -t master_cmds < "$master_commands_file"

# Launch flux and execute all commands in parallel
echo "Launching flux with $num_nodes nodes..."
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start "flux resource list; mapfile -t master_cmds < ${master_commands_file}; for ((i=0; i<\${#master_cmds[@]}; i++)); do flux submit -N1 -n1 -c7 --gpus-per-task=1 --output=flux_output/beta-infer_\$i.txt \${master_cmds[\$i]}; done; flux queue drain;"

echo "All LD matrix inference jobs submitted and completed."
