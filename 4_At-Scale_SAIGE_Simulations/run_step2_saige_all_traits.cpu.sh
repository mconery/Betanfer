#!/bin/bash
#SBATCH -A bif154
#SBATCH -J saige_step2_cpu
#SBATCH -o run_step2_saige_all_traits.cpu.log
#SBATCH -N 10
#SBATCH -p extended
#SBATCH -t 24:00:00
#SBATCH -C nvme

#Source bashrc and activate the conda environment
source ~/.bashrc
conda activate parsl_bash

#Run the script
saige_step2_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/run_step2_saige_all_traits.cpu.py"
python $saige_step2_script