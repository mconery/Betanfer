#!/bin/bash

#SBATCH -A bif154_mde
#SBATCH -J mvp_subset_plot
#SBATCH -o mvp_subset_plot.log
#SBATCH -N 1
#SBATCH -p batch-spi
#SBATCH -q debug
#SBATCH -t 2:00:00

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

# Input: rsID position lookup built on Aurora and transferred
lookup_file="/sw/frontier/mde/bif154_mde/betanfer_project/LD_matrices/rsid_pos_lookup.json"

# Locus list
locus_csv="${script_dir}/sampled_loci_750kb.csv"

# Output: shared-variant subset matrices
subset_dir="${base_dir}/LD_matrices/mvp_platlas_subset"

# Output: side-by-side heatplots
plot_dir="${base_dir}/mvp_matrix_plots"

# Populations
populations="AFR,EUR"

# Number of parallel workers
n_workers=20

echo "Starting MVP subset matrix extraction and heatplot generation"
echo "Betanfer matrices: ${betanfer_dir}"
echo "Raw MVP matrices:  ${mvp_dir}"
echo "rsID lookup:       ${lookup_file}"
echo "Subset output:     ${subset_dir}"
echo "Plot output:       ${plot_dir}"

###########################################################################################################################################################
# Step 1: Extract shared-variant subset matrices
###########################################################################################################################################################

echo ""
echo "Step 1: Extracting shared-variant subset matrices..."

python3 "${script_dir}/subset_matrices.py" \
    -b "${betanfer_dir}" \
    -m "${mvp_dir}" \
    -l "${locus_csv}" \
    -k "${lookup_file}" \
    -o "${subset_dir}" \
    -a "${populations}" \
    -n "${n_workers}"

if [ $? -ne 0 ]; then
    echo "ERROR: subset_matrices.py failed"
    exit 1
fi

echo "Subset matrices written to: ${subset_dir}"

###########################################################################################################################################################
# Step 2: Generate side-by-side heatplots
###########################################################################################################################################################

echo ""
echo "Step 2: Generating heatplots..."

python3 "${script_dir}/plot_matrices.py" \
    -i "${subset_dir}" \
    -o "${plot_dir}" \
    -l "${locus_csv}" \
    -a "${populations}" \
    -n "${n_workers}"

if [ $? -ne 0 ]; then
    echo "ERROR: plot_matrices.py failed"
    exit 1
fi

echo "Plots written to: ${plot_dir}"
echo ""
echo "Pipeline complete."
