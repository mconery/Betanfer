#!/bin/bash

#SBATCH -A bif154_mde
#SBATCH -J make_matrices
#SBATCH -o make_matrices.log
#SBATCH -N 300
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
script_dir="/gpfs/arx2/bif154_mde/proj-shared/betanfer_project/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/8_MVP_LD_Mapping"
var_dir="/sw/frontier/mde/bif154_mde/betanfer_project/LD_matrices"

# Script locations
wrapper_script="${script_dir}/make_matrices.py"
execution_script="${script_dir}/execute_command.py"

# Files and directories
plink_dir="/gpfs/arx2/bif154_mde/proj-shared/RITM0309091/pgen"
matrix_dir="${base_dir}/LD_matrices/mvp_external"
individual_dir="/gpfs/arx2/bif154_mde/proj-shared/MVP_meta-analysis/real_genotypes"
registry_file="/gpfs/arx2/bif154_mde/proj-shared/RITM0309091/pgen/hg19.chunk_register.csv"
command_file="${matrix_dir}/make_matrices.master.txt"

# Plink binary directory
plink_bin_dir="/gpfs/arx2/bif154_mde/proj-shared/tools/bin"

# Populations to process (comma-separated)
populations="AFR,EUR"

# Create matrix directory if it doesn't exist
mkdir -p "${matrix_dir}"

echo "Starting LD matrix generation"
echo "Matrix directory: ${matrix_dir}"
echo "Populations: ${populations}"
echo "Variant directory: ${var_dir}"

###########################################################################################################################################################

# Call python wrapper script with given parameters
# This will generate make_matrices.AFR.txt, make_matrices.EUR.txt, and make_matrices.master.txt
python3 "${wrapper_script}" \
    -v "${var_dir}" \
    -o "${matrix_dir}" \
    -a "${populations}" \
    -i "${individual_dir}" \
    -r "${registry_file}" \
    -c "${command_file}" \
    -p "${plink_dir}" \
    -b "${plink_bin_dir}"

if [ $? -ne 0 ]; then
    echo "ERROR: Python wrapper script failed"
    exit 1
fi

echo "Command files generated successfully"

#Read in the commands
mapfile -t cmds < ${command_file}
#Start flux and Run the flux script
NNODES=${SLURM_NNODES}
srun -N $NNODES -n $NNODES -c 56 flux start "flux resource list; mapfile -t cmds < ${command_file}; for ((i=0; i<${#cmds[@]}; i++)); do flux submit -N1 -n1 -c56 --output=flux_output/make_matrix_\$i.txt \${cmds[\$i]}; done; flux queue drain;"

echo "All jobs completed"
echo "LD matrices should be located in: ${matrix_dir}/AFR/ and ${matrix_dir}/EUR/"
