#!/bin/bash

#SBATCH -A bif154
#SBATCH -J SuSiE_finemap_distributed
#SBATCH -o finemap_SuSiE_distributed.log
#SBATCH -N 1
#SBATCH -p batch
#SBATCH -t 2:00:00

################################################################################
################## DISTRIBUTED SUSIE FINE-MAPPING PIPELINE #####################
#################### Using Lustre-based Work Queue System ######################
################################################################################

set -e  # Exit on error

echo "=========================================="
echo "DISTRIBUTED SUSIE FINE-MAPPING PIPELINE"
echo "=========================================="
echo "Start time: $(date)"

################################################################################################################
###################################### Load Modules ##########################################################
################################################################################################################

module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux

echo "Modules loaded successfully"

################################################################################################################
########################## Define directories and key file/script locations #################################
################################################################################################################

# Set script locations
map_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/SuSiE/finemap_SuSiE.R"
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/SuSiE/generate_susie_commands.py"

# Scripts for distributed workflow (copy paths to your locations)
node_manager_script="./node_manager.py"
distribute_script="./generate_and_distribute_commands.py"
monitor_script="./monitor_queue.py"

# Set environment location
map_enviro="/lustre/orion/bif154/proj-shared/mconery/tools/environments/susie"

# Set directories
sumstats_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"

# Base output directory for fine-mapping results
base_finemap_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"

# Directory for command files and work queue
command_files_dir="${base_finemap_dir}/command_files"
work_queue_dir="${base_finemap_dir}/work_queue"

mkdir -p "$command_files_dir"

# Configuration
commands_per_batch=41
num_nodes=${SLURM_NNODES:-1000}  # Default to 1000 if not in SLURM
workers_per_node=28

# Define heritabilities, trait types, and analysis types
heritabilities=("5e-5" "3e-4")
heritability_values=(0.00005 0.0003)
trait_types=("all" "continuous" "binary")
analysis_types=("fixed" "random")

# Total configurations: 2 heritabilities × 3 trait types × 2 analysis types = 12

echo "Configuration:"
echo "  Base finemap dir: $base_finemap_dir"
echo "  Work queue dir: $work_queue_dir"
echo "  Num nodes: $num_nodes"
echo "  Workers per node: $workers_per_node"
echo "  Commands per batch: $commands_per_batch"
echo ""

################################################################################################################
####################################### Verify Scripts and Paths #############################################
################################################################################################################

echo "Verifying required scripts and paths..."

for script in "$map_script" "$command_script" "$node_manager_script" "$distribute_script" "$monitor_script"; do
    if [[ ! -f "$script" ]]; then
        echo "ERROR: Required script not found: $script"
        exit 1
    fi
done

for dir in "$sumstats_dir" "$matrix_dir"; do
    if [[ ! -d "$dir" ]]; then
        echo "ERROR: Required directory not found: $dir"
        exit 1
    fi
done

echo "✓ All scripts and paths verified"
echo ""

################################################################################################################
################################## PHASE 1: Generate All Commands ###########################################
################################################################################################################

echo "=========================================="
echo "PHASE 1: GENERATE ALL SUSIE COMMANDS"
echo "=========================================="
echo "Start time: $(date)"

# Activate conda environment
source activate $map_enviro

# Create jobs directory for parallel command generation
jobs_dir="${command_files_dir}/command_gen_jobs"
mkdir -p "$jobs_dir"
rm -f "$jobs_dir"/*.sh

master_commands_file="${command_files_dir}/master_susie_commands.txt"
> "$master_commands_file"  # Clear if exists

total_jobs=0

# Create individual job scripts for command generation
for h_idx in "${!heritabilities[@]}"; do
    h_label=${heritabilities[$h_idx]}
    h_value=${heritability_values[$h_idx]}
    
    for trait_type in "${trait_types[@]}"; do
        for analysis_type in "${analysis_types[@]}"; do
            # Define matrix types for this configuration
            if [[ $trait_type == "all" ]]; then
                matrix_types="raw_in-sample,PC-corrected,beta-inferred_50_${analysis_type}_all,beta-inferred_100_${analysis_type}_all,beta-inferred_200_${analysis_type}_all,beta-inferred_300_${analysis_type}_all,beta-inferred_400_${analysis_type}_all,beta-inferred_500_${analysis_type}_all,beta-inferred_750_${analysis_type}_all,beta-inferred_1000_${analysis_type}_all,beta-inferred_1250_${analysis_type}_all,beta-inferred_1500_${analysis_type}_all,beta-inferred_1750_${analysis_type}_all,beta-inferred_2000_${analysis_type}_all,beta-inferred_2500_${analysis_type}_all,beta-inferred_3000_${analysis_type}_all,beta-inferred_3500_${analysis_type}_all,beta-inferred_4000_${analysis_type}_all"
            elif [[ $trait_type == "continuous" ]]; then
                matrix_types="raw_in-sample,PC-corrected,beta-inferred_50_${analysis_type}_continuous,beta-inferred_100_${analysis_type}_continuous,beta-inferred_200_${analysis_type}_continuous,beta-inferred_300_${analysis_type}_continuous,beta-inferred_400_${analysis_type}_continuous,beta-inferred_500_${analysis_type}_continuous,beta-inferred_750_${analysis_type}_continuous,beta-inferred_1000_${analysis_type}_continuous,beta-inferred_1250_${analysis_type}_continuous,beta-inferred_1500_${analysis_type}_continuous,beta-inferred_1750_${analysis_type}_continuous,beta-inferred_2000_${analysis_type}_continuous"
            else # binary
                matrix_types="raw_in-sample,PC-corrected,beta-inferred_50_${analysis_type}_binary,beta-inferred_100_${analysis_type}_binary,beta-inferred_200_${analysis_type}_binary,beta-inferred_300_${analysis_type}_binary,beta-inferred_400_${analysis_type}_binary,beta-inferred_500_${analysis_type}_binary,beta-inferred_750_${analysis_type}_binary,beta-inferred_1000_${analysis_type}_binary,beta-inferred_1250_${analysis_type}_binary,beta-inferred_1500_${analysis_type}_binary,beta-inferred_1750_${analysis_type}_binary,beta-inferred_2000_${analysis_type}_binary"
            fi
            
            # Define output directory
            h_dir="${analysis_type}_homo_${h_label}"
            out_dir="${base_finemap_dir}/${h_dir}/susie_results"
            command_file="${command_files_dir}/susie_${analysis_type}_${trait_type}_${h_label}_commands.txt"
            
            # Create job script for this configuration
            job_script="${jobs_dir}/job_${total_jobs}.sh"
            
            cat > "$job_script" << 'JOBEOF'
#!/bin/bash
set -e
OUT_DIR="PLACEHOLDER_OUT_DIR"
CMD_FILE="PLACEHOLDER_CMD_FILE"

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
JOBEOF
            
            # Replace placeholders
            sed -i "s|PLACEHOLDER_OUT_DIR|$out_dir|g" "$job_script"
            sed -i "s|PLACEHOLDER_CMD_FILE|$command_file|g" "$job_script"
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

# Execute jobs in parallel using Flux
echo "Executing command generation jobs in parallel..."
srun -N 1 -n 1 -c 56 flux start bash -c "
    echo \"Flux instance started for parallel command generation at \$(date)\"
    
    job_count=0
    for job_script in ${jobs_dir}/job_*.sh; do
        job_name=\$(basename \${job_script} .sh)
        flux submit -N1 -n1 -c56 bash \${job_script}
        ((job_count++))
        
        if [ \$((job_count % 2)) -eq 0 ]; then
            echo \" Submitted \${job_count}/${total_jobs} command generation jobs\"
        fi
    done
    
    echo \"All ${total_jobs} command generation jobs submitted at \$(date)\"
    
    # Wait for all jobs to complete
    echo \"Waiting for all command generation jobs to complete...\"
    flux queue drain
    
    echo \"All command generation jobs completed at \$(date)\"
"

if [ $? -ne 0 ]; then
    echo "ERROR: Command generation failed"
    exit 1
fi

echo "✓ Command generation completed"
echo ""

################################################################################################################
################################ PHASE 2: Merge Command Files #################################################
################################################################################################################

echo "=========================================="
echo "PHASE 2: MERGE COMMAND FILES"
echo "=========================================="

echo "Merging individual command files into master file..."

# Loop through all individual command files and append to master
ls "${command_files_dir}" | grep "^susie_" | grep "_commands.txt$" | while read command_file; do
    if [[ -f "${command_files_dir}/$command_file" ]] && [[ -s "${command_files_dir}/$command_file" ]]; then
        num_lines=$(wc -l < "${command_files_dir}/$command_file")
        grep -v "^#" "${command_files_dir}/$command_file" >> "$master_commands_file"
        echo "✓ Merged $num_lines commands from $(basename "$command_file")"
    fi
done

if [[ ! -f "$master_commands_file" ]] || [[ ! -s "$master_commands_file" ]]; then
    echo "ERROR: Master commands file is missing or empty!"
    exit 1
fi

total_commands=$(wc -l < "$master_commands_file")
echo "✓ Total commands merged: $total_commands"
echo ""

################################################################################################################
############################ PHASE 3: Create Distributed Work Queue ###########################################
################################################################################################################

echo "=========================================="
echo "PHASE 3: CREATE DISTRIBUTED WORK QUEUE"
echo "=========================================="

python "$distribute_script" \
    --all-commands "$master_commands_file" \
    --queue-dir "$work_queue_dir" \
    --num-nodes "$num_nodes" \
    --commands-per-batch "$commands_per_batch" \
    --node-manager "$node_manager_script" \
    --workers-per-node "$workers_per_node"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to create distributed work queue"
    exit 1
fi

echo ""

################################################################################################################
############################ PHASE 4: Monitor Queue Progress #################################################
################################################################################################################

echo "=========================================="
echo "PHASE 4: MONITOR QUEUE PROGRESS"
echo "=========================================="
echo ""
echo "Monitoring queue progress. This will run until all work is complete."
echo "Alternatively, you can monitor in another terminal with:"
echo "  python $monitor_script --queue-dir $work_queue_dir"
echo ""

# Monitor queue progress
python "$monitor_script" \
    --queue-dir "$work_queue_dir" \
    --interval 10 \
    --no-curses

echo ""

################################################################################################################
################################## COMPLETION AND SUMMARY ###################################################
################################################################################################################

echo "=========================================="
echo "FINE-MAPPING PIPELINE COMPLETE"
echo "=========================================="
echo "End time: $(date)"
echo ""

# Count results
echo "Results summary:"
total_output_files=$(find "$base_finemap_dir" -name "*.txt" -type f | wc -l)
echo "  Total output files: $total_output_files"

echo ""
echo "Output directory: $base_finemap_dir"
echo "Work queue directory: $work_queue_dir"
echo "Logs directory: $work_queue_dir/logs"
echo ""

echo "✓ Pipeline completed successfully!"
