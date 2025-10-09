#!/bin/bash

#SBATCH -A bif154
#SBATCH -J MultiSuSiE_finemap
#SBATCH -o finemap_MultiSuSiE.log
#SBATCH -N 40
#SBATCH -p batch
#SBATCH -t 1:00:00

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
map_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/finemap_MultiSuSiE.py"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/generate_multisusie_commands.py"

#Set environment location
map_enviro="/lustre/orion/bif154/proj-shared/mconery/tools/environments/multisusie"

#Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"
out_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results/random_homo_5e-5/multisusie_results"

#Set command file
command_file=$out_dir/"multisusie_commands.txt"

#Set populations
populations="AFR,EUR"

# Set number of nodes
num_nodes=${SLURM_NNODES}

################################################################################################################
############################################## Prep Mapping Commands ###########################################
################################################################################################################

python $command_script --wrapper-script $map_script --gwas-dir $sumstats_dir --ld-dir $matrix_dir --pheno-dir $pheno_dir --output-dir $out_dir --commands-file $command_file --matrix-types raw_in-sample,PC-corrected,beta-inferred_2000 --populations $populations --skip-existing --heritability 0.00005

################################################################################################################
################################## Execute Mapping Commands in Parallel ########################################
################################################################################################################

# Read mapping commands into array
mapfile -t cmds < "$command_file"

# Activate environment
source activate $map_enviro

# Launch flux and execute MultiSuSiE commands in parallel
echo "Launching flux for MultiSuSiE with $num_nodes nodes..."
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start "flux resource list; mapfile -t cmds < ${command_file}; for ((i=0; i<\${#cmds[@]}; i++)); do flux submit -N1 -n1 -c2 --output=flux_output/MultiSuSiE_\$i.txt \${cmds[\$i]}; done; flux queue drain;"

echo "All MultiSuSiE mapping jobs submitted and completed."
