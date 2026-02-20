#!/bin/bash

#SBATCH -A bif154
#SBATCH -J make_matrix
#SBATCH -o make_raw_in-sample_matrices.log
#SBATCH -N 10
#SBATCH -p batch
#SBATCH -t 1:00:00

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set input directories and files
var_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants/master_causal_list.txt"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample"
plink_dirs="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/AFR,/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/EUR"

# Script directories and files
script_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/scripts"
script_list_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/script_list.txt"

# Set script locations
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/make_raw_in-sample_matrix_command_writer.py"
dedup_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/remove_matrix_duplicates.py"

# Set populations
populations=AFR,EUR

# Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
############################################## Make Basic LD Matrices ##########################################
################################################################################################################

# Load modules and activate conda directory
module load flux
module load miniforge3
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno

# Create output directory for flux logs
mkdir -p flux_output

# Generate individual shell scripts
echo "Generating individual shell scripts..."
python $command_script \
    --script_dir $script_dir \
    --script_list_file $script_list_file \
    --matrix_dir $matrix_dir \
    --plink_dirs $plink_dirs \
    --populations $populations \
    --variant_lists_file $var_file \
    --variant-file-suffix ".variant_set"

# Check if script generation was successful
if [ ! -f "$script_list_file" ]; then
    echo "ERROR: Script list file not created. Exiting."
    exit 1
fi

# Count number of scripts
num_scripts=$(wc -l < "$script_list_file")
echo "Generated $num_scripts individual shell scripts"

# Use flux to run matrix creation
echo "Starting flux and submitting jobs..."
mapfile -t script_paths < "${script_list_file}"

srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start \
"flux resource list; mapfile -t script_paths < ${script_list_file}; for ((i=0; i<\${#script_paths[@]}; i++)); do flux submit -N1 -n1 -c4 --output=flux_output/make_raw_matrix_\$i.txt \${script_paths[\$i]}; done; flux queue drain;"

# Clean-up log and nosex files
echo "Cleaning up temporary files..."
rm -f $matrix_dir/*/*.nosex
rm -f $matrix_dir/*/*.log
rm -f $matrix_dir/scripts/*
rm -f $matrix_dir/*/*temp_dedup*

echo 'SUCCESS: All raw in-sample matrices created.'
