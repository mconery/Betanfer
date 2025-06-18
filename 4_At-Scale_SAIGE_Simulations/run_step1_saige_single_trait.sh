#!/bin/bash
#Extract command line parameters
plink_file=$1
pheno_file=$2
phenotype=$3
out_dir=$4

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
srun -N ${NNODES} --ntasks-per-node 1 mkdir /mnt/bb/${USER}/${ENV_NAME}
echo "untaring ${ENV_NAME}"
srun -N ${NNODES} --ntasks-per-node 1 tar -xzf /mnt/bb/${USER}/${TAR_FILE} -C /mnt/bb/${USER}/${ENV_NAME}

# Unpack the env
source activate /mnt/bb/${USER}/${ENV_NAME}
srun -N ${NNODES} --ntasks-per-node 1 conda-unpack
##### END OF SBCAST AND CONDA-UNPACK #####


#Set up SAIGE-GPU variables
R_LIB=/lustre/orion/bif154/proj-shared/arodriguez/tools/conda_envs/RSAIGE_1.3.3_amd_gpu/lib/R/library
path_to_saige="/lustre/orion/bif154/proj-shared/arodriguez/tools/SAIGE-GPU/src/SAIGE"
#Change into output directory
cd $out_dir

#sbcast genotype files to local scratch
sbcast -pf $plink_file.bim /mnt/bb/${USER}/temp.bim
sbcast -pf $plink_file.fam /mnt/bb/${USER}/temp.fam
sbcast -pf $plink_file.bed /mnt/bb/${USER}/temp.bed

#Run step 1
Rscript $path_to_saige/extdata/step1_fitNULLGLMM.R --plinkFile=/mnt/bb/${USER}/temp --phenoFile=$pheno_file --invNormalize=FALSE --phenoCol=$phenotype --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 --sampleIDColinphenoFile=IID --traitType=binary --outputPrefix=./$phenotype.step1_output --minMAFforGRM 0.01 --LOCO T  --IsOverwriteVarianceRatioFile=TRUE --nThreads=1