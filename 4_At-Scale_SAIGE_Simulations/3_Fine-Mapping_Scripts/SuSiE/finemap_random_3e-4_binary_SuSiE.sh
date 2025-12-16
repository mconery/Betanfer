#!/bin/bash

#SBATCH -A bif154
#SBATCH -J SuSiE_finemap
#SBATCH -o finemap_SuSiE.log
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
map_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/SuSiE/finemap_SuSiE.R"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/SuSiE/generate_susie_commands.py"

#Set environment location
map_enviro="/lustre/orion/bif154/proj-shared/mconery/tools/environments/susie"

#Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
out_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results/random_homo_3e-4/susie_results"

#Set command file
command_file=$out_dir/command_files/"susie_binary_commands.txt"

# Set number of nodes
num_nodes=${SLURM_NNODES}

################################################################################################################
############################################## Prep Mapping Commands ###########################################
################################################################################################################

python $command_script --wrapper-script $map_script --gwas-dir $sumstats_dir --ld-dir $matrix_dir --output-dir $out_dir --commands-file $command_file --matrix-types raw_in-sample,PC-corrected,beta-inferred_50_random_binary,beta-inferred_100_random_binary,beta-inferred_200_random_binary,beta-inferred_300_random_binary,beta-inferred_400_random_binary,beta-inferred_500_random_binary,beta-inferred_750_random_binary,beta-inferred_1000_random_binary,beta-inferred_1250_random_binary,beta-inferred_1500_random_binary,beta-inferred_1750_random_binary,beta-inferred_2000_random_binary --skip-existing --heritability 0.0003 --analysis-type random

################################################################################################################
################################## Execute Mapping Commands in Parallel ########################################
################################################################################################################

# Read mapping commands into array
mapfile -t cmds < "$command_file"

# Activate environment
source activate $map_enviro

# Launch flux and execute SuSiE commands in parallel
echo "Launching flux for SuSiE with $num_nodes nodes..."
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start "flux resource list; mapfile -t cmds < ${command_file}; for ((i=0; i<\${#cmds[@]}; i++)); do flux submit -N1 -n1 -c2 --output=flux_output/SuSiE_\$i.txt \${cmds[\$i]}; done; flux queue drain;"

echo "All SuSiE mapping jobs submitted and completed."