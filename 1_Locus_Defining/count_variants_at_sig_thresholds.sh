#!/bin/bash -l
#PBS -N countvarsigthresh
#PBS -l select=1:ncpus=64
#PBS -l walltime=24:00:00
#PBS -q preemptable
#PBS -A GeomicVar

source ~/.bashrc
conda activate finemap
python /grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/1_Locus_Defining/count_variants_at_sig_thresholds.py -l /grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.master.csv -o /grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.plus_var_counts.csv -p 1 -m 0