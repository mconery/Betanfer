#!/bin/bash
#SBATCH -A bif154
#SBATCH -J index_gwama_mpi
#SBATCH -o index_all_null_gwama.log
#SBATCH -N 8
#SBATCH -p batch
#SBATCH -t 1:00:00                    

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

mpi_index_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/index_gwama_files_mpi.py"
gwama_pattern="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results/*/Pheno*.meta.out"
gwama_regex="^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.meta\\.out$"

srun -N $SLURM_NNODES -n $SLURM_NNODES -c 60 "python $mpi_index_script \
  --scratch_dir /lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats \
  --pattern $gwama_pattern \
  --filename_regex $gwama_regex \
  --local_workers 60"