#!/bin/bash

#SBATCH -A bif154
#SBATCH -J MESuSiE_finemap
#SBATCH -o finemap_MESuSiE.log
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
map_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/MESuSiE/finemap_MESuSiE.R"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/MESuSiE/generate_mesusie_commands.py"

#Set environment location
map_enviro="/lustre/orion/bif154/proj-shared/mconery/tools/environments/carma"

#Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
geno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
out_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results/fixed_homo_5e-5/mesusie_results"

#Set command file
command_file=$out_dir/command_files/"mesusie_binary_commands.txt"

#Set populations
populations="AFR,EUR"

# Set number of nodes
num_nodes=${SLURM_NNODES}

################################################################################################################
############################################## Prep Mapping Commands ###########################################
################################################################################################################

python $command_script --wrapper-script $map_script --gwas-dir $sumstats_dir --ld-dir $matrix_dir --output-dir $out_dir --commands-file $command_file --matrix-types raw_in-sample,PC-corrected,beta-inferred_50_fixed_binary,beta-inferred_100_fixed_binary,beta-inferred_200_fixed_binary,beta-inferred_300_fixed_binary,beta-inferred_400_fixed_binary,beta-inferred_500_fixed_binary,beta-inferred_750_fixed_binary,beta-inferred_1000_fixed_binary,beta-inferred_1250_fixed_binary,beta-inferred_1500_fixed_binary,beta-inferred_1750_fixed_binary,beta-inferred_2000_fixed_binary --populations $populations --skip-existing --heritability 0.00005 

################################################################################################################
################################## Execute Mapping Commands in Parallel ########################################
################################################################################################################

# Read mapping commands into array
mapfile -t cmds < "$command_file"

# Activate environment
source activate $map_enviro

# Launch flux and execute MESuSiE commands in parallel
echo "Launching flux for MESuSiE with $num_nodes nodes..."
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start "flux resource list; mapfile -t cmds < ${command_file}; for ((i=0; i<\${#cmds[@]}; i++)); do flux submit -N1 -n1 -c2 --output=flux_output/MESuSiE_\$i.txt \${cmds[\$i]}; done; flux queue drain;"

echo "All MESuSiE mapping jobs submitted and completed."