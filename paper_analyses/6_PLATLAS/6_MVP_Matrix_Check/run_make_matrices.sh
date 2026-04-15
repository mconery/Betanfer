#!/bin/bash

#SBATCH -A bif154_mde
#SBATCH -J mvp_platlas_matrices
#SBATCH -o mvp_platlas_matrices.log
#SBATCH -N 10
#SBATCH -p batch-spi
#SBATCH -q debug
#SBATCH -t 2:00:00

#################################################################### Load Needed Modules #################################################################

module load Core/25.03
module load r/4.4.1
module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux

################################################################# Set files and locations #################################################################

# Base directories
base_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project"
script_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Betanfer/paper_analyses/6_PLATLAS/6_MVP_Matrix_Check"

# Input files
locus_csv="${script_dir}/sampled_loci_750kb.csv"
registry_file="/gpfs/arx2/bif154_mde/proj-shared/RITM0309091/pgen/hg19.chunk_register.csv"
individual_dir="/gpfs/arx2/bif154_mde/proj-shared/MVP_meta-analysis/real_genotypes"

# Plink files and binaries
plink_dir="/gpfs/arx2/bif154_mde/proj-shared/RITM0309091/pgen"
plink_bin_dir="/gpfs/arx2/bif154_mde/proj-shared/tools/bin"

# Output locations
matrix_dir="${base_dir}/LD_matrices/mvp_platlas_loci"
command_file="${matrix_dir}/make_matrices.master.txt"
command_files_dir="${base_dir}/LD_matrices/command_files"

# Populations to process (comma-separated)
populations="AFR,EUR"

# Create output directories if they don't exist
mkdir -p "${matrix_dir}"
mkdir -p "${command_files_dir}"

echo "Starting MVP PLATLAS locus LD matrix generation"
echo "Locus CSV:      ${locus_csv}"
echo "Matrix dir:     ${matrix_dir}"
echo "Populations:    ${populations}"

###########################################################################################################################################################

# Step 1: Generate the per-job command file
python3 "${script_dir}/make_matrices.py" \
    -l "${locus_csv}" \
    -o "${matrix_dir}" \
    -a "${populations}" \
    -i "${individual_dir}" \
    -r "${registry_file}" \
    -c "${command_file}" \
    -p "${plink_dir}" \
    -b "${plink_bin_dir}"

if [ $? -ne 0 ]; then
    echo "ERROR: make_matrices.py failed"
    exit 1
fi

echo "Command file generated: ${command_file}"
echo "Total commands: $(wc -l < ${command_file})"

################################################################################################################
################################## Create batched command files ################################################
################################################################################################################

echo ""
echo "Creating batched command files..."

# Clean up any existing batch files for this pipeline
rm -f "${command_files_dir}"/mvp_platlas_matrix_*.sh

# Split commands — one command per batch file (each runs as a separate flux job)
python_splitter="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Betanfer/paper_analyses/5_LD_Inference/split_commands_parallel.py"

python "${python_splitter}" \
    -i "${command_file}" \
    -o "${command_files_dir}" \
    -b 1 \
    -p "mvp_platlas_matrix_" \
    -s ".sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Python batch splitter failed"
    exit 1
fi

num_batch_files=$(ls "${command_files_dir}" | grep -c "^mvp_platlas_matrix_" || true)
echo "Verified: Created ${num_batch_files} batch files in ${command_files_dir}"

if [ ${num_batch_files} -eq 0 ]; then
    echo "ERROR: No batch files were created"
    exit 1
fi

################################################################################################################
####################################### Execute Commands in Parallel ###########################################
################################################################################################################

NNODES=${SLURM_NNODES}

mkdir -p flux_output

srun -N $NNODES -n $NNODES -c 56 flux start \
    "flux resource list; \
     batch_count=0; \
     for batch in ${command_files_dir}/mvp_platlas_matrix_*.sh; do \
         batch_name=\$(basename \${batch} .sh); \
         flux submit -N1 -n1 -c56 \
             --output=flux_output/\${batch_name}.log \
             bash \${batch}; \
         ((batch_count++)); \
         if [ \$((batch_count % 20)) -eq 0 ]; then \
             echo \"  Submitted \${batch_count}/${num_batch_files} batches\"; \
         fi; \
     done; \
     flux queue drain;"

echo "All jobs completed"
echo "LD matrices located in: ${matrix_dir}/AFR/ and ${matrix_dir}/EUR/"
