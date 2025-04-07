#!/bin/bash -l
#PBS -N countvarsig
#PBS -l select=1:ncpus=64
#PBS -l walltime=24:00:00
#PBS -q preemptable
#PBS -A GeomicVar
#PBS -l filesystems=home:grand

source ~/.bashrc
conda activate finemap
python /grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/1_Locus_Defining/get_variants_in_sig_loci.py -l /grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.master.csv -o /grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci/sig_loci.trait_variant_identities.csv