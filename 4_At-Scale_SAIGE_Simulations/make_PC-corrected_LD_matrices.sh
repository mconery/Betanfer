#!/bin/bash

#SBATCH -A bif154
#SBATCH -J make_PC_LD_matrix
#SBATCH -o make_PC-corrected_LD_matrices.log
#SBATCH -N 20
#SBATCH -p batch
#SBATCH -t 1:00:00

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set input directories and files
var_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants/master_causal_list.txt"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/PC-corrected"
plink_dirs="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/AFR,/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/EUR"
extract_cmd_out_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/PC-corrected/extract_genotypes.cmds.txt"
ld_cmd_out_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/PC-corrected/calculate_ld.cmds.txt"

# Set PC file locations for each population
pc_files="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/AFR/merged.plink.maf.ld.pcs.eigenvec,/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/EUR/merged.plink.maf.ld.pcs.eigenvec"

# Set script locations
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/make_PC-corrected_LD_command_writer.py"
residualize_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/residualize_and_compute_ld.py"

# Set populations
populations=AFR,EUR

# Set number of PCs to regress out (top 10 PCs as mentioned in your GWAS)
num_pcs=10

# Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
###################################### Make Full Matrix-Making Commands ########################################
################################################################################################################

#Load modules and activate conda directory
module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/infersparsepro

echo "Starting PC-residualized LD matrix generation..."
echo "Using $num_pcs principal components for residualization"
echo "PC files: $pc_files"

# Create output directory if it doesn't exist
mkdir -p "$matrix_dir"

# Make command file with PC residualization
python $command_script \
    --extract_cmd_out $extract_cmd_out_file \
    --ld_cmd_out $ld_cmd_out_file \
    --matrix_dir $matrix_dir \
    --plink_dirs $plink_dirs \
    --populations $populations \
    --variant_lists_file $var_file \
    --pc_files $pc_files \
    --residualize_script $residualize_script \
    --num_pcs $num_pcs \
    --variant-file-suffix ".variant_set"

echo "Extraction command file created: $extract_cmd_out_file"
echo "Number of command sets: $(grep -c "plink2" $extract_cmd_out_file)"
echo "LD command file created: $ld_cmd_out_file"
echo "Number of command sets: $(grep -c "python" $ld_cmd_out_file)"

################################################################################################################
################################################# Extract Genotypes ############################################
################################################################################################################

#Use flux to extract genotypes
mapfile -t cmds < $extract_cmd_out_file
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start "flux resource list; mapfile -t cmds < ${extract_cmd_out_file}; for ((i=0; i<${#cmds[@]}; i++)); do flux submit -N1 -n1 -c56 --output=flux_output/extract_pc_geno_\$i.txt \${cmds[\$i]}; done; flux queue drain;" 
echo 'NOTE: All genotypes extracted.'


################################################################################################################
############################################ Make Residualized Matrices ########################################
################################################################################################################

#Use flux to extract genotypes
mapfile -t cmds < $ld_cmd_out_file
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start "flux resource list; mapfile -t cmds < ${ld_cmd_out_file}; for ((i=0; i<${#cmds[@]}; i++)); do flux submit -N1 -n1 -c7 --gpus-per-task=1 --output=flux_output/make_pc_matrix_\$i.txt \${cmds[\$i]}; done; flux queue drain;" 
echo 'NOTE: All genotypes extracted.'

#Clean up temp
rm -rf $matrix_dir/*/temp

