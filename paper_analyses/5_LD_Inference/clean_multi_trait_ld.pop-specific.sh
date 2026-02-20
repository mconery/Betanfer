#!/bin/bash

#SBATCH -A bif154
#SBATCH -J clean_multi_trait_ld_pop-specific
#SBATCH -o clean_multi_trait_ld.pop-specific.log
#SBATCH -N 10
#SBATCH -p extended
#SBATCH -t 6:00:00

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
matrix_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/5_LD_Inference/infer_matrix_from_null_traits.py"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/5_LD_Inference/create_beta_inference_commands.py"
cleanup_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/5_LD_Inference/cleanup_junk_ld_files.py"
python_splitter="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/5_LD_Inference/split_commands_parallel.py"

# Base directories
base_ld_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# EUR population directories
eur_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_EUR_sumstats"
eur_variant_dir="${base_ld_dir}/raw_in-sample/EUR"

# AFR population directories
afr_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_AFR_sumstats"
afr_variant_dir="${base_ld_dir}/raw_in-sample/AFR"

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Define trait types and number of traits
#trait_types=("all" "continuous" "binary")
trait_types=("all")
#num_traits_list=(50 100 200 300 400 500 750 1000 1250 1500 1750 2000)
num_traits_list=(2500 3000 3500 4000)

# Directory for all command files
commands_dir="${base_ld_dir}/command_files"
mkdir -p "$commands_dir"

# Batch configuration
commands_per_batch=1

# Master commands file
master_commands_file="${commands_dir}/master_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Track total commands generated
total_commands=0
total_command_files=0

# Cleanup log file
cleanup_log_file="${commands_dir}/cleanup_junk_files.pop-specific.log"
> "$cleanup_log_file" # Clear the file if it exists

################################################################################################################
################################## Generate Commands for All Combinations ####################################
################################################################################################################

# Create a temporary directory for job scripts
temp_job_dir="${commands_dir}/.temp-job-scripts"
mkdir -p "${temp_job_dir}"

# Create a temporary manifest file to track all job combinations
job_manifest="${temp_job_dir}/job-manifest.pop-specific.txt"
> "${job_manifest}"

# Loop over trait types
for trait_type in "${trait_types[@]}"; do
    echo ""
    echo "===================================="
    echo "Processing trait type: $trait_type"
    echo "===================================="

    # Loop over number of traits
    for n_traits in "${num_traits_list[@]}"; do
        echo ""
        echo "Cleaning junk for $n_traits traits, type: $trait_type"

        # EUR population
        echo "  → EUR population..."
        eur_output_dir="${base_ld_dir}/beta-inferred_${n_traits}_random_${trait_type}/EUR"
        eur_commands_file="${commands_dir}/eur_${n_traits}_${trait_type}_commands.txt"
        # Add to manifest for tracking
        echo "${eur_commands_file}|${eur_output_dir}|${eur_variant_dir}|EUR|${n_traits}|${trait_type}" >> "${job_manifest}"
	       
	
        # AFR population
        echo "  → AFR population..."
        afr_output_dir="${base_ld_dir}/beta-inferred_${n_traits}_random_${trait_type}/AFR"
        afr_commands_file="${commands_dir}/afr_${n_traits}_${trait_type}_commands.txt"
        # Add to manifest for tracking
        echo "${afr_commands_file}|${afr_output_dir}|${afr_variant_dir}|AFR|${n_traits}|${trait_type}" >> "${job_manifest}"


    done
done

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
while IFS='|' read -r commands_file output_dir variant_dir pop n_traits trait_type; do
  
  echo ""
  echo "Cleaning: ${pop}_${n_traits}_${trait_type}"
  
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
    echo "✓ Cleanup scan completed for ${pop}_${n_traits}_${trait_type}"
  fi

done < "${job_manifest}"

echo ""
echo "Cleanup scan completed at: $(date)"
echo "Review cleanup report: ${cleanup_log_file}"
echo ""

################################################################################################################
################################### Clean Junk for Raw and PC-Corrected In-Sample ##############################
################################################################################################################

for type in raw_in-sample PC-corrected; do for pop in AFR EUR; do for suffix in bin bim frq; do rm "${base_ld_dir}/$type/$pop/*.$suffix"; done; done; done	
