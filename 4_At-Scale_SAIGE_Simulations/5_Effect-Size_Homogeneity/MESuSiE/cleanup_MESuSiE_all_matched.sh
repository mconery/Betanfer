#!/bin/bash

#SBATCH -A bif154
#SBATCH -J MESuSiE_cleanup_matched
#SBATCH -o cleanup_MESuSiE_all_matched.log
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
cleanup_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/5_Effect-Size_Homogeneity/MESuSiE/cleanup_mesusie_meta_junk_files.py"

# Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# Base output directory for fine-mapping results
base_finemap_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"

# Define heritabilities and match types for MESuSiE matched analysis
heritabilities=("5e-5" "3e-4")
heritability_values=(0.00005 0.0003)
match_types=("matched_to_EUR" "matched_to_AFR")

# Total cleanup jobs: 2 heritabilities × 2 match types = 4
# Note: MESuSiE uses analysis_type-based subdirectories within each heritability/match_type directory

################################################################################################################
############################ Create and Execute Cleanup Jobs in Parallel ######################################
################################################################################################################

echo "=================================================="
echo "STARTING MESUSIE FINE-MAPPING JUNK FILE CLEANUP (Matched Effect-Size Analysis)"
echo "=================================================="
echo "Start time: $(date)"
echo "Per-SNP heritabilities: ${heritabilities[@]}"
echo "Match types: ${match_types[@]}"
echo "Cleanup mode: SAFE MODE (files will NOT be deleted)"
echo "=================================================="
echo ""

# Create jobs directory for parallel execution
cleanup_jobs_dir="${base_finemap_dir}/cleanup_jobs_mesusie"
mkdir -p "$cleanup_jobs_dir"
rm -f "$cleanup_jobs_dir"/cleanup_job_*.sh

echo "Creating cleanup jobs..."
echo ""

# Counter for cleanup jobs
total_cleanup_jobs=0

# Create individual cleanup job scripts for parallel execution
for h_idx in "${!heritabilities[@]}"; do
  h_label=${heritabilities[$h_idx]}
  h_value=${heritability_values[$h_idx]}
  
  for match_type in "${match_types[@]}"; do
    
    # For MESuSiE, output goes to analysis_type-based subdirectories
    h_dir="fixed_${match_type}_${h_label}"
    out_dir="${base_finemap_dir}/${h_dir}/mesusie_results"
    
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
MATCH_TYPE="PLACEHOLDER_MATCH_TYPE"

echo "Scanning for junk files in: $OUTPUT_DIR"
echo "  Heritability: $H_VALUE"
echo "  Match type: $MATCH_TYPE"

python "$CLEANUP_SCRIPT" \
  --gwas-dir "$GWAS_DIR" \
  --ld-dir "$MATRIX_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --heritability "$H_VALUE" \
  --filter-keyword "$MATCH_TYPE"

if [ $? -ne 0 ]; then
  echo "ERROR: Cleanup scan failed for $OUTPUT_DIR"
  exit 1
fi

echo "✓ Cleanup scan completed for ${MATCH_TYPE} (h²=${H_VALUE})"

JOBEOF

      # Replace placeholders
      sed -i "s|PLACEHOLDER_CLEANUP_SCRIPT|$cleanup_script|g" "$job_script"
      sed -i "s|PLACEHOLDER_GWAS_DIR|$sumstats_dir|g" "$job_script"
      sed -i "s|PLACEHOLDER_MATRIX_DIR|$matrix_dir|g" "$job_script"
      sed -i "s|PLACEHOLDER_OUTPUT_DIR|$out_dir|g" "$job_script"
      sed -i "s|PLACEHOLDER_H_VALUE|$h_value|g" "$job_script"
      sed -i "s|PLACEHOLDER_MATCH_TYPE|$match_type|g" "$job_script"
      
      chmod +x "$job_script"
      
      ((total_cleanup_jobs++))
      
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
