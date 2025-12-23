#!/bin/bash

#SBATCH -A bif154
#SBATCH -J RSparsePro_finemap_all
#SBATCH -o finemap_RSparsePro_all.log
#SBATCH -N 14
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
map_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/RSparsePro/finemap_RSparsePro.sh"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/RSparsePro/generate_rsparsepro_commands.py"
python_splitter="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/split_commands_parallel.py"

# Set environment location
map_enviro="/lustre/orion/bif154/proj-shared/mconery/tools/environments/rsparsepro"

# Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# Base output directory for fine-mapping results
base_finemap_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"

# Command files directory
command_files_dir="${base_finemap_dir}/command_files"
mkdir -p "$command_files_dir"

# Batch configuration
commands_per_batch=20

# Master commands file
master_commands_file="${command_files_dir}/master_rsparsepro_commands.txt"
> "$master_commands_file" # Clear the file if it exists

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Define heritabilities, trait types, and analysis types
heritabilities=("5e-5" "3e-4")
heritability_values=(0.00005 0.0003)
trait_types=("all" "continuous" "binary")
analysis_types=("fixed" "random")

# Total configurations: 2 heritabilities × 3 trait types × 2 analysis types = 12

################################################################################################################
################################## Generate Commands in Parallel Using Flux ##################################
################################################################################################################

echo "=================================================="
echo "STARTING RSPARSEPRO FINE-MAPPING PIPELINE (All Configurations)"
echo "=================================================="
echo "Start time: $(date)"
echo "Per-SNP heritabilities: ${heritabilities[@]}"
echo "Trait types: ${trait_types[@]}"
echo "Analysis types: ${analysis_types[@]}"
echo "Total configurations: $((${#heritabilities[@]} * ${#trait_types[@]} * ${#analysis_types[@]}))"
echo "Master commands file: $master_commands_file"
echo "Batch size: ${commands_per_batch} commands per batch"
echo "=================================================="

# Create jobs directory for parallel execution
jobs_dir="${command_files_dir}/command_gen_jobs"
mkdir -p "$jobs_dir"
rm -f "$jobs_dir"/*.sh

echo ""

# Create individual job scripts for parallel execution
total_jobs=0
for h_idx in "${!heritabilities[@]}"; do
    h_label=${heritabilities[$h_idx]}
    h_value=${heritability_values[$h_idx]}
    
    for trait_type in "${trait_types[@]}"; do
        for analysis_type in "${analysis_types[@]}"; do
            # Define matrix types for this configuration
            if [[ $trait_type == "all" ]]; then
                matrix_types="raw_in-sample,PC-corrected,beta-inferred_50_${analysis_type}_all,beta-inferred_100_${analysis_type}_all,beta-inferred_200_${analysis_type}_all,beta-inferred_300_${analysis_type}_all,beta-inferred_400_${analysis_type}_all,beta-inferred_500_${analysis_type}_all,beta-inferred_750_${analysis_type}_all,beta-inferred_1000_${analysis_type}_all,beta-inferred_1250_${analysis_type}_all,beta-inferred_1500_${analysis_type}_all,beta-inferred_1750_${analysis_type}_all,beta-inferred_2000_${analysis_type}_all"
            elif [[ $trait_type == "continuous" ]]; then
                matrix_types="raw_in-sample,PC-corrected,beta-inferred_50_${analysis_type}_continuous,beta-inferred_100_${analysis_type}_continuous,beta-inferred_200_${analysis_type}_continuous,beta-inferred_300_${analysis_type}_continuous,beta-inferred_400_${analysis_type}_continuous,beta-inferred_500_${analysis_type}_continuous,beta-inferred_750_${analysis_type}_continuous,beta-inferred_1000_${analysis_type}_continuous,beta-inferred_1250_${analysis_type}_continuous,beta-inferred_1500_${analysis_type}_continuous,beta-inferred_1750_${analysis_type}_continuous,beta-inferred_2000_${analysis_type}_continuous"
            else  # binary
                matrix_types="raw_in-sample,PC-corrected,beta-inferred_50_${analysis_type}_binary,beta-inferred_100_${analysis_type}_binary,beta-inferred_200_${analysis_type}_binary,beta-inferred_300_${analysis_type}_binary,beta-inferred_400_${analysis_type}_binary,beta-inferred_500_${analysis_type}_binary,beta-inferred_750_${analysis_type}_binary,beta-inferred_1000_${analysis_type}_binary,beta-inferred_1250_${analysis_type}_binary,beta-inferred_1500_${analysis_type}_binary,beta-inferred_1750_${analysis_type}_binary,beta-inferred_2000_${analysis_type}_binary"
            fi

            # Define output directory
            h_dir="${analysis_type}_homo_${h_label}"
            out_dir="${base_finemap_dir}/${h_dir}/rsparsepro_results"
            command_file="${command_files_dir}/rsparsepro_${analysis_type}_${trait_type}_${h_label}_commands.txt"

            # Create job script
            job_script="${jobs_dir}/job_${total_jobs}.sh"
            cat > "$job_script" << 'JOBEOF'
#!/bin/bash
set -e

OUT_DIR="PLACEHOLDER_OUT_DIR"
CMD_FILE="PLACEHOLDER_CMD_FILE"
MASTER_CMD_FILE="PLACEHOLDER_MASTER_CMD_FILE"

mkdir -p "$OUT_DIR"

echo "Generating commands for PLACEHOLDER_ANALYSIS_TYPE_PLACEHOLDER_TRAIT_TYPE_PLACEHOLDER_H_LABEL..."

python PLACEHOLDER_CMD_SCRIPT \
    --wrapper-script PLACEHOLDER_MAP_SCRIPT \
    --gwas-dir PLACEHOLDER_SUMSTATS_DIR \
    --ld-dir PLACEHOLDER_MATRIX_DIR \
    --output-dir "$OUT_DIR" \
    --commands-file "$CMD_FILE" \
    --matrix-types PLACEHOLDER_MATRIX_TYPES \
    --skip-existing \
    --heritability PLACEHOLDER_H_VALUE \
    --analysis-type PLACEHOLDER_ANALYSIS_TYPE

if [[ ! -f "$CMD_FILE" ]] || [[ ! -s "$CMD_FILE" ]]; then
    echo "ERROR: Commands file is missing or empty: $CMD_FILE"
    exit 1
fi

num_commands=$(wc -l < "$CMD_FILE")
echo "✓ Generated $num_commands commands for PLACEHOLDER_ANALYSIS_TYPE_PLACEHOLDER_TRAIT_TYPE_PLACEHOLDER_H_LABEL"

# Append to master commands file
grep -v "^#" "$CMD_FILE" >> "$MASTER_CMD_FILE"
JOBEOF

            # Replace placeholders
            sed -i "s|PLACEHOLDER_OUT_DIR|$out_dir|g" "$job_script"
            sed -i "s|PLACEHOLDER_CMD_FILE|$command_file|g" "$job_script"
            sed -i "s|PLACEHOLDER_MASTER_CMD_FILE|$master_commands_file|g" "$job_script"
            sed -i "s|PLACEHOLDER_CMD_SCRIPT|$command_script|g" "$job_script"
            sed -i "s|PLACEHOLDER_MAP_SCRIPT|$map_script|g" "$job_script"
            sed -i "s|PLACEHOLDER_SUMSTATS_DIR|$sumstats_dir|g" "$job_script"
            sed -i "s|PLACEHOLDER_MATRIX_DIR|$matrix_dir|g" "$job_script"
            sed -i "s|PLACEHOLDER_MATRIX_TYPES|$matrix_types|g" "$job_script"
            sed -i "s|PLACEHOLDER_H_VALUE|$h_value|g" "$job_script"
            sed -i "s|PLACEHOLDER_ANALYSIS_TYPE|$analysis_type|g" "$job_script"
            sed -i "s|PLACEHOLDER_TRAIT_TYPE|$trait_type|g" "$job_script"
            sed -i "s|PLACEHOLDER_H_LABEL|$h_label|g" "$job_script"

            chmod +x "$job_script"
            ((total_jobs++))
        done
    done
done

echo "Created $total_jobs command generation jobs"
echo ""
echo "Executing command generation jobs in parallel using Flux..."
echo "Start time: $(date)"

# Use Flux to parallelize command generation across nodes
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start bash -c "
    echo \"Flux instance started for parallel command generation at \$(date)\"
    echo \"Available resources:\"
    flux resource list
    echo \"\"
    
    echo \"Submitting ${total_jobs} command generation jobs...\"
    job_count=0
    for job_script in ${jobs_dir}/job_*.sh; do
        flux submit -N1 -n1 -c56 \
            bash \${job_script}
        
        ((job_count++))
        
        # Print progress every 2 jobs
        if [ \$((job_count % 2)) -eq 0 ]; then
            echo \"  Submitted \${job_count}/${total_jobs} command generation jobs\"
        fi
    done
    
    echo \"All ${total_jobs} command generation jobs submitted at \$(date)\"
    echo \"\"
    
    # Wait for all jobs to complete
    echo \"Waiting for all command generation jobs to complete...\"
    flux queue drain
    
    echo \"All command generation jobs completed at \$(date)\"
"

if [ $? -ne 0 ]; then
    echo "ERROR: Command generation failed"
    exit 1
fi

echo "Command generation completed at $(date)"

################################################################################################################
################################## Summary and Final Status ###################################################
################################################################################################################

echo ""
echo "=================================================="
echo "COMMAND FILE GENERATION COMPLETE"
echo "=================================================="
echo "Total configurations processed: $total_jobs"
echo "Total commands in master file: $(wc -l < "$master_commands_file")"
echo "Master commands file: $master_commands_file"
echo "Command directory: $command_files_dir"
echo "=================================================="

# Verify master commands file
if [[ ! -f "$master_commands_file" ]] || [[ ! -s "$master_commands_file" ]]; then
    echo "ERROR: Master commands file is missing or empty!"
    exit 1
fi

################################################################################################################
################################## Create batched command files ################################################
################################################################################################################

echo ""
echo "Creating batched command files using parallel Python splitter..."

# Clean up any existing rsparsepro batch files
rm -f "${command_files_dir}"/rsparsepro_finemap_batch_*.sh

# Use Python splitter for fast parallel batch creation
python "${python_splitter}" \
    -i "${master_commands_file}" \
    -o "${command_files_dir}" \
    -b ${commands_per_batch} \
    -p "rsparsepro_finemap_batch_" \
    -s ".sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Python batch splitter failed"
    exit 1
fi

# Verify batch creation
num_batch_files=$(ls "${command_files_dir}" | grep -c "^rsparsepro_finemap_batch_" || true)
echo "Verified: Created ${num_batch_files} batch files in ${command_files_dir}"

if [ ${num_batch_files} -eq 0 ]; then
    echo "ERROR: No batch files were created"
    exit 1
fi

################################################################################################################
####################################### Execute Commands in Parallel #########################################
################################################################################################################

echo ""
echo "Starting flux-based execution of batched RSparsePro fine-mapping commands at $(date)"

# Create flux output directory
flux_output_dir="${base_finemap_dir}/flux_output"
mkdir -p "$flux_output_dir"

# Activate environment
source activate $map_enviro

# Start Flux and submit batch jobs
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start bash -c "
    # Display Flux resource allocation
    echo \"Flux instance started at \$(date)\"
    echo \"Available resources:\"
    flux resource list
    echo \"\"
    
    # Submit all batch files to Flux
    echo \"Submitting ${num_batch_files} batches to Flux queue...\"
    batch_count=0
    for batch in ${command_files_dir}/rsparsepro_finemap_batch_*.sh; do
        batch_name=\$(basename \${batch} .sh)
        flux submit -N1 -n1 -c2 \
            --output=${flux_output_dir}/\${batch_name}.log \
            bash \${batch}
        
        ((batch_count++))
        
        # Print progress every 20 batches
        if [ \$((batch_count % 20)) -eq 0 ]; then
            echo \"  Submitted \${batch_count}/${num_batch_files} batches\"
        fi
    done
    
    echo \"All ${num_batch_files} batches submitted to Flux queue at \$(date)\"
    echo \"\"
    
    # Wait for all jobs to complete
    echo \"Waiting for all batches to complete...\"
    flux queue drain
    
    echo \"\"
    echo \"All batches completed at \$(date)\"
    
    # Display final resource status
    echo \"Final Flux resource status:\"
    flux resource list
"

################################################################################################################
################################## Post-processing and cleanup ##################################################
################################################################################################################

echo ""
echo "Post-processing results..."

# Count successful batch completions by checking for output files
successful_batches=$(ls -1 "${flux_output_dir}"/rsparsepro_finemap_batch_*.log 2>/dev/null | wc -l)
echo "Successfully completed batches: ${successful_batches}/${num_batch_files}"

# Check for any failed batches
if [ ${successful_batches} -ne ${num_batch_files} ]; then
    echo "WARNING: Some batches may have failed. Check logs in ${flux_output_dir}"
fi

# Cleanup temporary job scripts
rm -rf "$jobs_dir"

echo ""
echo "RSparsePro fine-mapping workflow completed at $(date)"
echo "Batch files retained in: ${command_files_dir}"
echo "Flux logs available in: ${flux_output_dir}"
echo "Fine-mapping results in: ${base_finemap_dir}"
