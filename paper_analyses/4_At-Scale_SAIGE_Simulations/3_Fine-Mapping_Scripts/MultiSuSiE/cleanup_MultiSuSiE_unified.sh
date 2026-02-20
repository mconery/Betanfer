#!/bin/bash

#SBATCH -A bif154
#SBATCH -J MultiSuSiE_cleanup_unified
#SBATCH -o cleanup_MultiSuSiE_all.log
#SBATCH -N 1
#SBATCH -p batch
#SBATCH -t 2:00:00

################################################################################################################
################################################### Load Modules ###############################################
################################################################################################################

module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set script locations
map_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/MultiSuSiE/finemap_MultiSuSiE.py"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/MultiSuSiE/generate_multisusie_commands.py"
cleanup_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/MultiSuSiE/cleanup_multisusie_junk_files.py"
python_splitter="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/5_LD_Inference/split_commands_parallel.py"

# Set environment location
map_enviro="/lustre/orion/bif154/proj-shared/mconery/tools/environments/multisusie"

# Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"

# Base output directory for fine-mapping results
base_finemap_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"

# Command files directory
command_files_dir="${base_finemap_dir}/command_files"
mkdir -p "$command_files_dir"

# Cleanup log file
cleanup_log_file="${command_files_dir}/cleanup_multisusie_junk_files.log"
> "$cleanup_log_file" # Clear the file if it exists

# Batch configuration
commands_per_batch=50

# Master commands file
master_commands_file="${command_files_dir}/master_multisusie_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Populations
populations="AFR,EUR"

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Define heritabilities and trait types
heritabilities=("5e-5" "3e-4")
heritability_values=(0.00005 0.0003)
trait_types=("all" "continuous" "binary")

################################################################################################################
######################################## Clean Junk Files from Output Directories #############################
################################################################################################################

echo ""
echo "=================================================="
echo "CLEANING JUNK FILES FROM FINE-MAPPING OUTPUT DIRECTORIES"
echo "=================================================="
echo "Start time: $(date)"
echo "Running in SAFE MODE (no files will be deleted)"
echo ""

# Loop through each heritability and trait type combination
cleanup_count=0

for i in "${!heritabilities[@]}"; do
  h_label=${heritabilities[$i]}
  h_value=${heritability_values[$i]}
  
  if [[ $h_label == "5e-5" ]]; then
    h_dir="fixed_homo_5e-5"
  else
    h_dir="fixed_homo_3e-4"
  fi
  
  for trait_type in "${trait_types[@]}"; do
    out_dir="${base_finemap_dir}/${h_dir}/multisusie_results"
    
    echo ""
    echo "Scanning: ${trait_type}_${h_label}"
    
    # Run cleanup script in safe mode (no deletion, just reporting)
    python "${cleanup_script}" \
      --gwas-dir "${sumstats_dir}" \
      --pheno-dir "${pheno_dir}" \
      --ld-dir "${matrix_dir}" \
      --output-dir "${out_dir}" \
      --populations "${populations}" \
      --heritability "${h_value}" \
      --delete-mode >> "$cleanup_log_file" 2>&1
    
    if [ $? -ne 0 ]; then
      echo "WARNING: Cleanup scan failed for ${out_dir}" >&2
      echo "Check ${cleanup_log_file} for details"
    else
      echo "✓ Cleanup scan completed for ${trait_type}_${h_label}"
    fi
    
    ((cleanup_count++))
  done
done

echo ""
echo "Cleanup scan completed: $cleanup_count output directory(ies) scanned"
echo "Cleanup scan completed at: $(date)"
echo "Review cleanup report: ${cleanup_log_file}"
echo ""
