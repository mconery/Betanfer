#!/bin/bash
#SBATCH -A bif154
#SBATCH -J saige_step2
#SBATCH -o run_step2_saige_all_traits.log
#SBATCH -N 1
#SBATCH -q debug
#SBATCH -t 1:00:00
#SBATCH -C nvme

##################################################################################################################
################################################### Set Needed Variables #########################################
##################################################################################################################
#Input directories
plink_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"
out_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/null_saige_results"

#Set up SAIGE-GPU variables
R_LIB=/lustre/orion/bif154/proj-shared/arodriguez/tools/conda_envs/RSAIGE_1.3.3_amd_gpu/lib/R/library
path_to_saige="/lustre/orion/bif154/proj-shared/arodriguez/tools/SAIGE-GPU/src/SAIGE"
saige_step1_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/run_step1_saige_single_trait.sh"

#Set population names
populations=("AFR" "EUR")
trait_types=("binary" "continuous")

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
################################################## Execute SAIGE Step 2 ##########################################
##################################################################################################################

#Testing code
pop=AFR
trait_type=binary
chromo=1
plink_file=$plink_dir/$pop/chr"$chromo".plink.maf.ld
#Change into output directory
cd $out_dir
pop_out_dir=$out_dir/$pop

#Create manifest file
manifest_file=$out_dir/$pop/$trait_type.manifest.txt
>$manifest_file
ls $out_dir/$pop/*.step1_output.rda | while read file; do
	temp=${file%.step1_output.rda}
	pheno=$(basename $temp)
	echo $file $out_dir/$pop/$pheno.step1_output.varianceRatio.txt $out_dir/$pop/$pheno.step2_output.txt | awk -F " " '{print $1"\t"$2"\t"$3}' >> $manifest_file
done

Rscript $path_to_saige/extdata/step2_SPAtests.R -bgenFile=$plink_file.bgen --bgenFileIndex=$plink_file.bgen.bgi \
--manifestFile $manifest_file --chrom=$chromo --AlleleOrder=ref-first --minMAF=0 --minMAC=0.5 --LOCO=TRUE --is_fastTest=TRUE --is_Firth_beta=TRUE
