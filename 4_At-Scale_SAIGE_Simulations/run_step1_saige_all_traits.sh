#!/bin/bash
#SBATCH -A bif154
#SBATCH -J saige_step1
#SBATCH -o run_step1_saige_all_traits.log
#SBATCH -N 1882
#SBATCH -p batch
#SBATCH -t 3:00:00
#SBATCH -C nvme

##################################################################################################################
################################################### Set Needed Variables #########################################
##################################################################################################################
#Input directories
plink_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"
out_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step1_results"

#Set up SAIGE-GPU variables
R_LIB=/lustre/orion/bif154/proj-shared/arodriguez/tools/conda_envs/RSAIGE_1.3.3_amd_gpu/lib/R/library
path_to_saige="/lustre/orion/bif154/proj-shared/arodriguez/tools/SAIGE-GPU/src/SAIGE"
saige_step1_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/run_step1_saige_single_trait.sh"

#Set population names
populations=("AFR" "EUR")
trait_types=("binary" "continuous")

#Set number of gpus needed per phenotype
gpus=3

##################################################################################################################
############################################## Load Modules and Conda Env ########################################
##################################################################################################################
#Activate needed modules
module load python/3.10-miniforge3
module load r/4.4.1
module load amd/6.4.0

##### START OF SBCAST AND CONDA-UNPACK #####
#Set variables
ENV_NAME=RSAIGE_1.3.3_amd_gpu
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
#################################################### Run SAIGE Step 1 ############################################
##################################################################################################################

#Get total number of gpus
total_gpus=$(($NNODES * 8))
#Get maximum number of parallel jobs
max_parallel=$(($total_gpus / $gpus))

cmds=()
#Loop over the populations, create the mapping commands for all phenotypes, and launch them
for pop in ${populations[@]}; do
	plink_file_pre=$plink_dir/$pop/merged.plink.maf.ld
	plink_file_post=/mnt/bb/${USER}/$pop
	#sbcast genotype files to local scratch
	sbcast -pf $plink_file_pre.bim $plink_file_post.bim
	sbcast -pf $plink_file_pre.fam $plink_file_post.fam
	sbcast -pf $plink_file_pre.bed $plink_file_post.bed
	#Change into output directory
	pop_out_dir=$out_dir/$pop
	mkdir -p $pop_out_dir
	cd $pop_out_dir
	touch $pop_out_dir/timing.csv  # Create timing file
	#Create commands to process step 1 for all phenotypes of the population
	#Loop over the trait types
	for trait_type in ${trait_types[@]}; do
		#Make a list of the phenotype files for the given trait type and population
		ls $pheno_dir/$pop/Pheno*.$trait_type*txt > $pheno_dir/$pop/$trait_type.phenotype_files.txt
		#Loop over the files and create the commands
		while read file; do
			temp=$(basename $file)
			prefix=${temp%.txt}
			#Check if output files exist and add command if they do not
			if [ ! -f $pop_out_dir/$prefix.step1_output.rda ] || [ ! -f $pop_out_dir/$prefix.step1_output.varianceRatio.txt ]; then
				cmds+=("$saige_step1_script $plink_file_post $file $prefix $pop_out_dir $path_to_saige $trait_type")
			fi
		done < $pheno_dir/$pop/$trait_type.phenotype_files.txt
	done
done

#Run all the commands for the populations
job_size=${#cmds[@]}
for (( i=0; i<${#cmds[@]}; i++ )); do
	j=$(( $i + 1 ))
	srun -N1 -n $gpus --gpus-per-task 1 --gpu-bind closest --exclusive -C nvme --job-name=$j ${cmds[$i]} &
	if [ $(($j % $max_parallel)) == 0 ] || [ $j == $job_size ]; then 
	    	wait
        fi
done
