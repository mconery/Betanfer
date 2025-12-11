#!/bin/bash
#SBATCH -A bif154
#SBATCH -J index_gwama_mpi
#SBATCH -o index_all_null_gwama_random.log
#SBATCH -N 40
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

#Set script locations
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/create_gwama_index_commands.py"
indexer_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/index_single_gwama_file.py"

#Set other file location
scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats_random"
meta_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"
commands_file=$scratch_dir/gwama_index_commands.txt

#Set number of nodes
num_nodes=${SLURM_NNODES}

################################################################################################################
################################## Prep Commands and launch indexing ###########################################
################################################################################################################

#Prep commands file
python $command_script --scratch_dir $scratch_dir --pattern $meta_dir/'*/Pheno*.meta.out' --commands_file $commands_file --filename_regex "^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.random_meta\\.out$" --indexer_script $indexer_script

#Execute commands
mapfile -t cmds < $commands_file
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start "flux resource list; mapfile -t cmds < ${commands_file}; for ((i=0; i<${#cmds[@]}; i++)); do flux submit -N1 -n1 -c2 --output=flux_output/index_gwama_\$i.txt \${cmds[\$i]}; done; flux queue drain;" 
echo 'NOTE: All files indexed.'
