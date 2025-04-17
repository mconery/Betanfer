#!/bin/bash -l
#PBS -A GeomicVar
#PBS -l walltime=3:00:00
#PBS -l filesystems=flare
#PBS -l select=32
#PBS -q prod
#PBS -M mconery@anl.gov
#PBS -N onekg_test 
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/multinode_1000G_simulations.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/multinode_1000G_simulations.err

#Set up conda environment
source ~/.bashrc
conda activate finemap_old

#Set key directories and locations
onekg_script=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/2_Correlation_Inference/OnekG_Simulation_Wrapper.py
plink_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/plink_files
out_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/simulation_results
keep_dir=/flare/GeomicVar/mconery/resources/1000G

#Change directory into directory where script logs are stored
cd $(dirname $onekg_script)
cd ../../script_logs

#Launch wrapper script
mpiexec -n 32 python $onekg_script --plink_dir $plink_dir --out_dir $out_dir --keep_dir $keep_dir