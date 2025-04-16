#!/bin/bash -l
#PBS -N countvarsigthresh
#PBS -l select=4:ncpus=256:mpiprocs=256
#PBS -l walltime=1:00:00
#PBS -q debug-scaling
#PBS -A GeomicVar
#PBS -l filesystems=home:flare

source ~/.bashrc
conda activate finemap_old
mpiexec -n 256 python -m mpi4py.futures /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/1_Locus_Defining/count_variants_at_sig_thresholds.py -l /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.ALL.master.csv -g /flare/GeomicVar/multi-biobank-meta/metaanalysis_for_export/ALL -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.ALL.plus_var_counts.csv -p 1 -m 0