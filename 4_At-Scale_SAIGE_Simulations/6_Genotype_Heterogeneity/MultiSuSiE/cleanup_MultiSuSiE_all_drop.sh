#!/bin/bash

#SBATCH -A bif154
#SBATCH -J MultiSuSiE_cleanup_drop
#SBATCH -o cleanup_MultiSuSiE_all_drop.log
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
cleanup_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/6_Genotype_Heterogeneity/MultiSuSiE/cleanup_multisusie_drop_causal_junk_files.py"

# Set directories - NOTE: MultiSuSiE uses step2_results (population subdirectories), NOT meta_analysis_results
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# Base output directory for fine-mapping results
base_finemap_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"

# Define heritabilities and central populations for MultiSuSiE drop_causal analysis
heritabilities=("5e-5" "3e-4")
heritability_values=(0.00005 0.0003)

# Central populations (the ones with drop_causal files)
central_pops=("EUR" "AFR")

# Map central pop to drop_causal match_type for directory naming
declare -A match_type_map
match_type_map["EUR"]="drop_causal_EUR"
match_type_map["AFR"]="drop_causal_AFR"

# Total cleanup jobs: 2 heritabilities × 2 central_pops = 4

# Note: MultiSuSiE uses flat output structure (base_finemap_dir/{h_dir}/multisusie_results)
# where h_dir = fixed_drop_causal_{central_pop}_{heritability}

################################################################################################################
############################ Create and Execute Cleanup Jobs in Parallel ######################################
################################################################################################################

echo "=================================================="
echo "STARTING MULTISUSIE FINE-MAPPING JUNK FILE CLEANUP (Drop Causal Analysis)"
echo "=================================================="
echo "Start time: $(date)"
echo "Per-SNP heritabilities: ${heritabilities[@]}"
echo "Central populations: ${central_pops[@]}"
echo "Cleanup mode: SAFE MODE (files will NOT be deleted)"
echo "=================================================="
echo ""

# Create jobs directory for parallel execution
cleanup_jobs_dir="${base_finemap_dir}/cleanup_jobs_multisusie_drop"
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
    
    for central_pop in "${central_pops[@]}"; do
        # Get the match_type for directory naming
        match_type=${match_type_map[$central_pop]}
        
        # Output directory: fixed_{match_type}_{h_label}/multisusie_results
        h_dir="fixed_${match_type}_${h_label}"
        out_dir="${base_finemap_dir}/${h_dir}/multisusie_results"
        
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
CENTRAL_POP="PLACEHOLDER_CENTRAL_POP"

echo "Scanning for junk files in: $OUTPUT_DIR"
echo " Heritability: $H_VALUE"
echo " Central population (drop_causal): $CENTRAL_POP"

python3 "$CLEANUP_SCRIPT" \
    --gwas-dir "$GWAS_DIR" \
    --ld-dir "$MATRIX_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --central-pop "$CENTRAL_POP" \
    --heritability "$H_VALUE" 

if [ $? -ne 0 ]; then
    echo "ERROR: Cleanup scan failed for $OUTPUT_DIR"
    exit 1
fi

echo "✓ Cleanup scan completed for ${CENTRAL_POP} drop_causal (h²=${H_VALUE})"

JOBEOF

        # Replace placeholders
        sed -i "s|PLACEHOLDER_CLEANUP_SCRIPT|$cleanup_script|g" "$job_script"
        sed -i "s|PLACEHOLDER_GWAS_DIR|$sumstats_dir|g" "$job_script"
        sed -i "s|PLACEHOLDER_MATRIX_DIR|$matrix_dir|g" "$job_script"
        sed -i "s|PLACEHOLDER_OUTPUT_DIR|$out_dir|g" "$job_script"
        sed -i "s|PLACEHOLDER_H_VALUE|$h_value|g" "$job_script"
        sed -i "s|PLACEHOLDER_CENTRAL_POP|$central_pop|g" "$job_script"
        
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

# Run cleanup jobs sequentially in safe mode (no deletion by default)
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
echo "To delete junk files, run with --delete-mode flag:"
echo "  python3 cleanup_multisusie_drop_causal_junk_files.py --central-pop EUR --delete-mode"
echo ""

# Cleanup temporary job scripts
rm -rf "$cleanup_jobs_dir"
