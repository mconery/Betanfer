#!/bin/bash

#SBATCH -A bif154
#SBATCH -J Collate_finemap
#SBATCH -o run_result_collations.log
#SBATCH -N 1
#SBATCH -p batch
#SBATCH -t 1:00:00

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

# Define paths - MODIFY THESE TO MATCH YOUR ACTUAL PATHS
CAUSAL_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants"
RESULTS_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"
OUTPUT_FILE="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/summary_results/homogeneity_results_file.txt"

#Fine-mapping collator script
collate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/fine_mapping_collator.py"

################################################################################################################
############################################# Execute Results Collation ########################################
################################################################################################################

# First, validate the directory structure
echo "Validating directory structure..."
python $collate_script \
    --causal-dir "$CAUSAL_DIR" \
    --results-dir "$RESULTS_DIR" \
    --validate-only

if [ $? -ne 0 ]; then
    echo "Directory validation failed. Exiting."
    exit 1
fi

# Run the collation script
echo "Starting fine-mapping results collation..."
python $collate_script \
    --causal-dir "$CAUSAL_DIR" \
    --results-dir "$RESULTS_DIR" \
    --output "$OUTPUT_FILE" \
    --verbose

if [ $? -eq 0 ]; then
    echo "Collation completed successfully!"
    echo "Output file: $OUTPUT_FILE"
    echo "File size: $(du -h $OUTPUT_FILE | cut -f1)"
    echo "Number of lines: $(wc -l < $OUTPUT_FILE)"
else
    echo "Collation failed!"
    exit 1
fi

echo "Job completed at: $(date)"