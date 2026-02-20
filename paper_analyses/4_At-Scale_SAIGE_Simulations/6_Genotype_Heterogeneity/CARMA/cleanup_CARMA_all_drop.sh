#!/bin/bash

#SBATCH -A bif154
#SBATCH -J CARMA_cleanup_drop
#SBATCH -o cleanup_CARMA_all_drop.log
#SBATCH -N 1
#SBATCH -p batch
#SBATCH -t 0:30:00

################################################################################################################
################################################### Load Modules ###############################################
################################################################################################################

module load miniforge3

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set script location
cleanup_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/5_Effect-Size_Homogeneity/CARMA/cleanup_carma_meta_junk_files.py"

# Set environment location
map_enviro="/lustre/orion/bif154/proj-shared/mconery/tools/environments/carma"

# Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# Base output directory for fine-mapping results
base_finemap_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"

# Define heritabilities, trait types, and analysis types
heritabilities=("5e-5" "3e-4")
heritability_values=(0.00005 0.0003)
trait_types=("all" "continuous" "binary")
analysis_types=("fixed" "random")
match_types=("drop_causal_EUR" "drop_causal_AFR")

# Total cleanup jobs: 2 heritabilities × 2 analysis types × 2 match types = 8
# (Note: trait_types don't affect output filenames, so we clean all together)

################################################################################################################
############################ Create and Execute Cleanup Jobs in Parallel ######################################
################################################################################################################

echo "=================================================="
echo "STARTING CARMA FINE-MAPPING JUNK FILE CLEANUP (Matched Effect-Size Analysis)"
echo "=================================================="
echo "Start time: $(date)"
echo "Per-SNP heritabilities: ${heritabilities[@]}"
echo "Analysis types: ${analysis_types[@]}"
echo "Match types: ${match_types[@]}"
echo "Cleanup mode: SAFE MODE (files will NOT be deleted)"
echo "=================================================="
echo ""

# Create jobs directory for parallel execution
cleanup_jobs_dir="${base_finemap_dir}/cleanup_jobs_carma"
mkdir -p "$cleanup_jobs_dir"
rm -f "$cleanup_jobs_dir"/cleanup_job_*.sh

echo "Creating cleanup jobs..."
echo ""

# Counter for cleanup jobs
total_cleanup_jobs=0

# Activate environment before launching cleanup jobs
source activate $map_enviro

# Create individual cleanup job scripts for parallel execution
for h_idx in "${!heritabilities[@]}"; do
  h_label=${heritabilities[$h_idx]}
  h_value=${heritability_values[$h_idx]}
  
  for analysis_type in "${analysis_types[@]}"; do
    
    for match_type in "${match_types[@]}"; do
      
      # Define output directory for this configuration
      h_dir="${analysis_type}_${match_type}_${h_label}"
      out_dir="${base_finemap_dir}/${h_dir}/carma_results"
      
      # Create job script
      job_script="${cleanup_jobs_dir}/cleanup_job_${total_cleanup_jobs}.sh"
      
      cat > "$job_script" << 'JOBEOF'
#!/bin/bash

set -e

CLEANUP_SCRIPT="PLACEHOLDER_CLEANUP_SCRIPT"
GWAS_DIR="PLACEHOLDER_GWAS_DIR"
MATRIX_DIR="PLACEHOLDER_MATRIX_DIR"
OUTPUT_DIR="PLACEHOLDER_OUTPUT_DIR"
H_VALUE="PLACEHOLDER_H_VALUE"
ANALYSIS_TYPE="PLACEHOLDER_ANALYSIS_TYPE"
MATCH_TYPE="PLACEHOLDER_MATCH_TYPE"

echo "Scanning for junk files in: $OUTPUT_DIR"

python "$CLEANUP_SCRIPT" \
  --gwas-dir "$GWAS_DIR" \
  --ld-dir "$MATRIX_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --heritability "$H_VALUE" \
  --analysis-type "$ANALYSIS_TYPE" \
  --filter-keyword "$MATCH_TYPE" \
  --delete-mode

if [ $? -ne 0 ]; then
  echo "ERROR: Cleanup scan failed for $OUTPUT_DIR"
  exit 1
fi

echo "✓ Cleanup scan completed for ${ANALYSIS_TYPE}_${MATCH_TYPE}_${H_VALUE}"

JOBEOF

      # Replace placeholders
      sed -i "s|PLACEHOLDER_CLEANUP_SCRIPT|$cleanup_script|g" "$job_script"
      sed -i "s|PLACEHOLDER_GWAS_DIR|$sumstats_dir|g" "$job_script"
      sed -i "s|PLACEHOLDER_MATRIX_DIR|$matrix_dir|g" "$job_script"
      sed -i "s|PLACEHOLDER_OUTPUT_DIR|$out_dir|g" "$job_script"
      sed -i "s|PLACEHOLDER_H_VALUE|$h_value|g" "$job_script"
      sed -i "s|PLACEHOLDER_ANALYSIS_TYPE|$analysis_type|g" "$job_script"
      sed -i "s|PLACEHOLDER_MATCH_TYPE|$match_type|g" "$job_script"
      
      chmod +x "$job_script"
      
      ((total_cleanup_jobs++))
      
    done
    
  done
  
done

echo "Created $total_cleanup_jobs cleanup jobs"
echo ""

echo "=================================================="
echo "EXECUTING CLEANUP JOBS"
echo "=================================================="
echo ""

# Initialize job counter
cleanup_job_count=0

# Run cleanup jobs sequentially in safe mode
for job_script in "${cleanup_jobs_dir}"/cleanup_job_*.sh; do
  
  job_name=$(basename "$job_script" .sh)
  
  echo "[$(($cleanup_job_count + 1))/${total_cleanup_jobs}] Running $job_name..."
  
  bash "$job_script"
  
  if [ $? -ne 0 ]; then
    echo "ERROR: $job_name failed"
    exit 1
  fi
  
  ((cleanup_job_count++))
  
done

echo ""
echo "=================================================="
echo "CLEANUP SCAN COMPLETE"
echo "=================================================="
echo "Cleanup jobs completed: $cleanup_job_count/${total_cleanup_jobs}"
echo "End time: $(date)"
echo ""

# Cleanup temporary job scripts
rm -rf "$cleanup_jobs_dir"
