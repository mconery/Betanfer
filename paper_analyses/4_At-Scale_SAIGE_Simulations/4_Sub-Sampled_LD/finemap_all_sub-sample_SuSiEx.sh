#!/bin/bash

#SBATCH -A bif154
#SBATCH -J SuSiEx_sub-sample
#SBATCH -o finemap_all_sub-sample_SuSiEx.log
#SBATCH -N 10
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
map_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/SuSiEx/finemap_SuSiEx.py"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/SuSiEx/generate_susiex_commands.py"

#Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample"
geno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
out_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results/random_sub-sample_LD_3e-4/susiex_results"

#Set command file
command_file=$out_dir/"susiex_commands.txt"

#Set populations
populations="AFR,EUR"

# Set number of nodes
num_nodes=${SLURM_NNODES}

################################################################################################################
############################################## Prep Mapping Commands ###########################################
################################################################################################################

python $command_script --wrapper-script $map_script --gwas-dir $sumstats_dir --ld-dir $matrix_dir --geno-dir $geno_dir --output-dir $out_dir --commands-file $command_file --matrix-types sub-sample_1000,sub-sample_5000,sub-sample_10000,sub-sample_50000,sub-sample_100000,sub-sample_200000,sub-sample_300000 --populations $populations --skip-existing --heritability 0.0003

################################################################################################################
################################## Execute Mapping Commands in Parallel ########################################
################################################################################################################

# Read mapping commands into array
mapfile -t cmds < "$command_file"

# Launch flux and execute SuSiEx commands in parallel
echo "Launching flux for SuSiEx with $num_nodes nodes..."
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start "flux resource list; mapfile -t cmds < ${command_file}; for ((i=0; i<\${#cmds[@]}; i++)); do flux submit -N1 -n1 -c2 --output=flux_output/SuSiEx_\$i.txt \${cmds[\$i]}; done; flux queue drain;"

echo "All SuSiEx mapping jobs submitted and completed."
