#!/bin/bash

#SBATCH -A bif154_mde
#SBATCH -J make_matrices
#SBATCH -o make_matrices.log
#SBATCH -N 40
#SBATCH -p batch-spi
#SBATCH -q debug
#SBATCH -t 2:00:00

module load miniforge3

################################################################# Set files and locations #################################################################

# Base directories
base_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project"
script_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/8_MVP_LD_Mapping"
var_dir="/sw/frontier/mde/bif154_mde/betanfer_project/LD_matrices/AFR"

# Script locations
wrapper_script="${script_dir}/make_matrices.py"
execution_script="${script_dir}/execute_command.py"

# Files and directories
plink_dir="/gpfs/arx2/bif154_mde/proj-shared/RITM039091/pgen"
matrix_dir="${base_dir}/LD_matrices"
individual_dir="/gpfs/arx2/bif154_mde/proj-shared/MVP_meta-analysis/real_genotypes"
variant_file="${var_dir}/Pheno1.1_causal.vars"  # Update this path to your variant file
registry_file="/gpfs/arx2/bif154_mde/proj-shared/RITM039091/pgen/hg19.chunk_register.csv"

# Plink binary directory
plink_bin_dir="/gpfs/arx2/bif154_mde/proj-shared/tools/bin"

# Output prefix for matrices
output_prefix="finemapping"

# Populations to process (comma-separated)
populations="AFR,EUR"

# Create matrix directory if it doesn't exist
mkdir -p "${matrix_dir}"
mkdir -p "${matrix_dir}/logs"

echo "Starting LD matrix generation"
echo "Matrix directory: ${matrix_dir}"
echo "Populations: ${populations}"
echo "Variant file: ${variant_file}"

###########################################################################################################################################################

# Call python wrapper script with given parameters
# This will generate make_matrices.AFR.txt, make_matrices.EUR.txt, and make_matrices.master.txt
python3 "${wrapper_script}" \
    -v "${variant_file}" \
    -o "${matrix_dir}" \
    -a "${populations}" \
    -i "${individual_dir}" \
    -r "${registry_file}" \
    -c "${matrix_dir}/make_matrices.txt" \
    -p "${plink_dir}" \
    -b "${plink_bin_dir}" \
    -n "${output_prefix}"

if [ $? -ne 0 ]; then
    echo "ERROR: Python wrapper script failed"
    exit 1
fi

echo "Command files generated successfully"

# Read the master command file and execute commands in parallel
master_command_file="${matrix_dir}/make_matrices.master.txt"

if [ ! -f "${master_command_file}" ]; then
    echo "ERROR: Master command file not found: ${master_command_file}"
    exit 1
fi

# Count the number of commands
num_commands=$(wc -l < "${master_command_file}")
echo "JOB SIZE: ${num_commands}"

# Array to hold commands
declare -a cmd_array=()

# Read commands from file
while IFS= read -r line; do
    if [ ! -z "$line" ]; then
        cmd_array+=("$line")
    fi
done < "${master_command_file}"

echo "Read ${#cmd_array[@]} commands from ${master_command_file}"

# Execute commands in parallel using srun
# Adjust -n1 and -c32 based on your resource allocation
for (( i=0; i<${#cmd_array[@]}; i++ )); do
    # Get a simple identifier for the job name
    job_id=$(printf "%04d" $i)
    
    echo "Submitting job ${job_id}/${#cmd_array[@]}"
    
    srun -n1 -c8 --exclusive --job-name="${job_id}" bash -c "${cmd_array[$i]}" &
    
    # Wait after every 480 tasks to prevent overwhelming the scheduler
    # Adjust this number based on your resource allocation
    if [ $((i+1)) -ge 480 ]; then
        echo "Waiting for batch to complete..."
        wait
    fi
done

echo "Waiting for all remaining jobs to complete..."
wait

echo "All jobs completed"
echo "LD matrices should be located in: ${matrix_dir}/AFR/ and ${matrix_dir}/EUR/"
