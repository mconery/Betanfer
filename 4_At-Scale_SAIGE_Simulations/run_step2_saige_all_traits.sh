#!/bin/bash
#SBATCH -A bif154
#SBATCH -J saige_step2
#SBATCH -o run_step2_saige_all_traits.log
#SBATCH -N 809
#SBATCH -q debug
#SBATCH -t 2:00:00
#SBATCH -C nvme

##################################################################################################################
################################################### Set Needed Variables #########################################
##################################################################################################################
#Set script locations
make_cmd_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/make_step2_commands.py"
execute_cmd_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/parsl_run_1_command_per_node.py"
#Set command file
cmd_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/step2_commands.txt" #This isn't actually used to ensure flux works (need to modify below)

##################################################################################################################
################################################### Prepare Command File #########################################
##################################################################################################################

#Reset environment
source deactivate
#Activate needed modules
module load python/3.10-miniforge3 r/4.4.1 amd/6.3.1 rocm/6.3.1 hwloc/2.11.1-gpu flux
#Call command making script
#python $make_cmd_script --causal-var-suffix ".variant_set.txt" #Commenting out to only run null phenotypes for now

##################################################################################################################
############################################## Load Modules and Conda Env ########################################
##################################################################################################################

##### START OF SBCAST AND CONDA-UNPACK #####
#Set variables
ENV_NAME=RSAIGE_1.3.3_amd6.3.1_gpu
TAR_FILE=${ENV_NAME}.tar.gz
TAR_DIR=/lustre/orion/bif154/proj-shared/arodriguez/tools/conda_envs
NNODES=${SLURM_NNODES}
# Move a copy of the env to the NVMe on each node
echo "copying ${ENV_NAME} to each node in the job"
sbcast -pf ${TAR_DIR}/${TAR_FILE} /mnt/bb/${USER}/${TAR_FILE}
if [ ! "$?" == "0" ]; then
    # CHECK EXIT CODE. When SBCAST fails, it may leave partial files on the compute nodes, and if you continue to launch srun,
    # your application may pick up partially complete shared library files, which would give you confusing errors.
    echo "SBCAST failed!"
    exit 1
fi

# Untar the environment file (only need 1 task per node to do this)
srun -N ${NNODES} --ntasks-per-node 1 mkdir -p /mnt/bb/${USER}/${ENV_NAME}
echo "untaring ${ENV_NAME}"
srun -N ${NNODES} --ntasks-per-node 1 tar -xzf /mnt/bb/${USER}/${TAR_FILE} -C /mnt/bb/${USER}/${ENV_NAME}

# Unpack the env
source activate /mnt/bb/${USER}/${ENV_NAME}
srun -N ${NNODES} --ntasks-per-node 1 conda-unpack
##### END OF SBCAST AND CONDA-UNPACK #####

##################################################################################################################
################################################## Execute SAIGE Step 2 ##########################################
##################################################################################################################

#Read in the commands
mapfile -t cmds < /lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/step2_commands.txt
#Start flux and Run the flux script
srun -N $NNODES -n $NNODES -c 56 --gpus-per-node=8 flux start "flux resource list; mapfile -t cmds < /lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/step2_commands.txt; for ((i=0; i<${#cmds[@]}; i++)); do flux submit -N1 -n1 -c56 --output=flux_output/steptwo_\$i.txt \${cmds[\$i]}; done; flux queue drain;"