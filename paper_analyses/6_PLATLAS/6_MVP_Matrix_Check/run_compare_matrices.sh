#!/bin/bash

#SBATCH -A bif154_mde
#SBATCH -J mvp_matrix_compare
#SBATCH -o mvp_matrix_compare.log
#SBATCH -N 1
#SBATCH -p batch-spi
#SBATCH -q debug
#SBATCH -t 4:00:00

#################################################################### Load Needed Modules #################################################################

module load Core/25.03
module load miniforge3

################################################################# Set files and locations #################################################################

# Base directories
base_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project"
script_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check"

# Input: betanfer matrices transferred from Aurora
betanfer_dir="/sw/frontier/mde/bif154_mde/betanfer_project/LD_matrices/mvp_platlas_betanfer"

# Input: raw plink-based MVP matrices produced by make_matrices.py
mvp_dir="${base_dir}/LD_matrices/mvp_platlas_loci"

# Locus list
locus_csv="${script_dir}/sampled_loci_750kb.csv"

# rsID position lookup (built on Aurora by build_rsid_lookup.py, then transferred)
lookup_file="/sw/frontier/mde/bif154_mde/betanfer_project/LD_matrices/rsid_pos_lookup.tsv"

# Output
output_file="${base_dir}/LD_matrices/mvp_matrix_comparison.tsv.gz"

# Populations
populations="AFR,EUR"

echo "Starting Betanfer vs raw MVP LD matrix comparison"
echo "Betanfer matrices: ${betanfer_dir}"
echo "Raw MVP matrices:  ${mvp_dir}"
echo "rsID lookup:       ${lookup_file}"
echo "Output:            ${output_file}"

###########################################################################################################################################################

python3 "${script_dir}/compare_matrices.py" \
    -b "${betanfer_dir}" \
    -m "${mvp_dir}" \
    -l "${locus_csv}" \
    -k "${lookup_file}" \
    -o "${output_file}" \
    -a "${populations}"

if [ $? -ne 0 ]; then
    echo "ERROR: compare_matrices.py failed"
    exit 1
fi

echo "Comparison complete. Output: ${output_file}"
