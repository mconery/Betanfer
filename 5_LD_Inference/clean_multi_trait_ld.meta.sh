#!/bin/bash

#SBATCH -A bif154
#SBATCH -J clean_multi_trait_ld_meta
#SBATCH -o clean_multi_trait_ld.meta.log
#SBATCH -N 1
#SBATCH -p batch
#SBATCH -t 2:00:00
#SBATCH -q debug

################################################################################################################
###################################### Load Modules and Environment ############################################
################################################################################################################

module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux

source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/infersparsepro

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set script locations
matrix_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/infer_matrix_from_null_traits.py"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/create_beta_inference_commands.py"
cleanup_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/cleanup_junk_ld_files.py"
python_splitter="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/split_commands_parallel.py"

# Base directories
base_ld_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# Meta-analysis specific directories
fixed_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats_fixed"
random_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats_random"

# Variant directory (same for all meta-analyses)
variant_dir="${base_ld_dir}/raw_in-sample/weighted"

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Define meta-analysis types, trait types and number of traits
meta_types=("fixed" "random")
trait_types=("all" "continuous" "binary")
#trait_types=("all")
num_traits_list=(50 100 200 300 400 500 750 1000 1250 1500 1750 2000)
#num_traits_list=(2500 3000 3500 4000)

# Directory for all command files
commands_dir="${base_ld_dir}/command_files"
mkdir -p "$commands_dir"

# Batch configuration
commands_per_batch=2

# Master commands file
master_commands_file="${commands_dir}/master_meta_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Track total commands generated
total_commands=0
total_command_files=0

# Cleanup log file
cleanup_log_file="${commands_dir}/cleanup_junk_files.meta.log"
> "$cleanup_log_file" # Clear the file if it exists

################################################################################################################
################################## Generate Commands for All Combinations ####################################
################################################################################################################

# Create a temporary directory for job scripts
temp_job_dir="${commands_dir}/.temp-job-scripts"
mkdir -p "${temp_job_dir}"

# Create a temporary manifest file to track all job combinations
job_manifest="${temp_job_dir}/job-manifest.txt"
> "${job_manifest}"

job_count=0

# Loop over meta-analysis types
for meta_type in "${meta_types[@]}"; do

echo ""
echo "====================================="
echo "Processing meta-analysis type: $meta_type"
echo "====================================="

# Set indexed directory based on meta-analysis type
if [[ $meta_type == "fixed" ]]; then
indexed_dir="$fixed_scratch_dir"
else
indexed_dir="$random_scratch_dir"
fi

# Loop over trait types
for trait_type in "${trait_types[@]}"; do

echo ""
echo "Processing trait type: $trait_type"

# Loop over number of traits
for n_traits in "${num_traits_list[@]}"; do

echo ""
echo "Generating commands for $n_traits traits, meta-type: $meta_type, trait-type: $trait_type"

# Set output directory
output_dir="${base_ld_dir}/beta-inferred_${n_traits}_${meta_type}_${trait_type}/Meta"
commands_file="${commands_dir}/${meta_type}_${n_traits}_${trait_type}_meta_commands.txt"
kdir -p "$output_dir"

# Create a job script for this combination
job_script="${temp_job_dir}/gen-job-${meta_type}-${n_traits}-${trait_type}.sh"

# Add to manifest for tracking
echo "${job_script}|${commands_file}|${output_dir}|${meta_type}|${n_traits}|${trait_type}" >> "${job_manifest}"

((job_count++))

done

done

done

echo "Prepared ${job_count} command generation jobs"
echo ""

################################################################################################################
################################### Clean Junk Files from Output Directories ##################################
################################################################################################################

echo ""
echo "================================================"
echo "CLEANING JUNK FILES FROM OUTPUT DIRECTORIES"
echo "================================================"
echo "Start time: $(date)"
echo ""

# Loop through the manifest to clean each output directory
while IFS='|' read -r job_script commands_file output_dir meta_type n_traits trait_type; do
  
  echo ""
  echo "Cleaning: ${meta_type}_${n_traits}_${trait_type}"
  
  # Run cleanup script in safe mode (no deletion, just reporting)
  python "${cleanup_script}" \
    --variant_dir "${variant_dir}" \
    --output_dir "${output_dir}" \
    --n_traits "${n_traits}" \
    --delete-mode >> "$cleanup_log_file" 2>&1
  
  if [ $? -ne 0 ]; then
    echo "WARNING: Cleanup scan failed for ${output_dir}" >&2
    echo "Check ${cleanup_log_file} for details"
  else
    echo "✓ Cleanup scan completed for ${meta_type}_${n_traits}_${trait_type}"
  fi

done < "${job_manifest}"

echo ""
echo "Cleanup scan completed at: $(date)"
echo "Review cleanup report: ${cleanup_log_file}"
echo ""
