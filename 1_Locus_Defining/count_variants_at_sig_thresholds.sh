#!/bin/bash -l
#PBS -N countvarsigthresh
#PBS -l select=2:ncpus=416
#PBS -l walltime=1:00:00
#PBS -q debug
#PBS -A GeomicVar
#PBS -l filesystems=home:flare

source ~/.bashrc
conda activate finemap_old
python -m mpi4py.futures /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/1_Locus_Defining/count_variants_at_sig_thresholds.py -l /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.EAS.master.csv -g /flare/GeomicVar/multi-biobank-meta/metaanalysis_for_export/EAS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.EAS.plus_var_counts.csv -p 1 -m 0