#!/bin/bash

#SBATCH -A bif154
#SBATCH -J CARMA_cleanup_unified
#SBATCH -o cleanup_CARMA_unified.log
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
map_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/CARMA/finemap_CARMA.R"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/CARMA/generate_carma_commands.py"
cleanup_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/CARMA/cleanup_carma_junk_files.py"
python_splitter="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/5_LD_Inference/split_commands_parallel.py"
node_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/flux_node_worker.sh"

# Set environment location
ENV_NAME=carma
TAR_FILE=${ENV_NAME}.tar.gz
TAR_DIR="/lustre/orion/bif154/proj-shared/mconery/tools/environments"

# Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# Base output directory for fine-mapping results
base_finemap_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"

# Command files directory
command_files_dir="${base_finemap_dir}/command_files"
batch_dir="${command_files_dir}/node_batches"
mkdir -p "$command_files_dir"
mkdir -p $batch_dir

# Cleanup log file
cleanup_log_file="${command_files_dir}/cleanup_carma_junk_files.log"
> "$cleanup_log_file" # Clear the file if it exists

# Master commands file
master_commands_file="${command_files_dir}/master_carma_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Set number of nodes
num_nodes=${SLURM_NNODES}
NNODES=${SLURM_NNODES}

# Define heritabilities, trait types, and analysis types
heritabilities=("5e-5" "3e-4")
heritability_values=(0.00005 0.0003)
trait_types=("all" "continuous" "binary")
analysis_types=("fixed" "random")

# Total configurations: 2 heritabilities × 3 trait types × 2 analysis types = 12

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

# Loop through each heritability, trait type, and analysis type combination
cleanup_count=0

for h_idx in "${!heritabilities[@]}"; do
  h_label=${heritabilities[$h_idx]}
  h_value=${heritability_values[$h_idx]}
  
  for trait_type in "${trait_types[@]}"; do
    
    for analysis_type in "${analysis_types[@]}"; do
      
      h_dir="${analysis_type}_homo_${h_label}"
      out_dir="${base_finemap_dir}/${h_dir}/carma_results"
      
      echo ""
      echo "Scanning: ${analysis_type}_${trait_type}_${h_label}"
      
      # Run cleanup script in safe mode (no deletion, just reporting)
      python "${cleanup_script}" \
        --gwas-dir "${sumstats_dir}" \
        --ld-dir "${matrix_dir}" \
        --output-dir "${out_dir}" \
        --analysis-type "${analysis_type}" \
        --heritability "${h_value}" \
        --delete-mode >> "$cleanup_log_file" 2>&1
      
      if [ $? -ne 0 ]; then
        echo "WARNING: Cleanup scan failed for ${out_dir}" >&2
        echo "Check ${cleanup_log_file} for details"
      else
        echo "✓ Cleanup scan completed for ${analysis_type}_${trait_type}_${h_label}"
      fi
      
      ((cleanup_count++))
    done
  done
done

echo ""
echo "Cleanup scan completed: $cleanup_count output directory(ies) scanned"
echo "Cleanup scan completed at: $(date)"
echo "Review cleanup report: ${cleanup_log_file}"
echo ""
