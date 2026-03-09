#!/bin/bash

#SBATCH -A bif154
#SBATCH -J Find_Missing_Results
#SBATCH -o run_find_missing_results.log
#SBATCH -N 1
#SBATCH -q debug
#SBATCH -t 2:00:00


################################################################################################################
################################################### Load Modules ###############################################
################################################################################################################

module load miniforge3

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

RESULTS_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"
SUMMARY_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/summary_results"
OUTPUT_FILE="${SUMMARY_DIR}/missing_results_files.txt"

PYTHON_SCRIPT="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/find_missing_results.py"

################################################################################################################
############################################# Validate Setup ###################################################
################################################################################################################

echo "=================================================="
echo "FIND MISSING FINE-MAPPING RESULTS"
echo "=================================================="
echo "Start time: $(date)"
echo ""

if [ ! -d "$RESULTS_DIR" ]; then
    echo "ERROR: Results directory does not exist: $RESULTS_DIR"
    exit 1
fi

if [ ! -d "$SUMMARY_DIR" ]; then
    echo "ERROR: Summary directory does not exist: $SUMMARY_DIR"
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "ERROR: Python script does not exist: $PYTHON_SCRIPT"
    exit 1
fi

echo "Results directory : $RESULTS_DIR"
echo "Summary directory : $SUMMARY_DIR"
echo "Output file       : $OUTPUT_FILE"
echo "Python script     : $PYTHON_SCRIPT"
echo ""

################################################################################################################
################################################# Run Script ###################################################
################################################################################################################

python "$PYTHON_SCRIPT" \
    --results-dir "$RESULTS_DIR" \
    --summary-dir "$SUMMARY_DIR" \
    --output "$OUTPUT_FILE"

EXIT_CODE=$?

echo ""
echo "=================================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "Script completed successfully."
    if [ -f "$OUTPUT_FILE" ]; then
        NUM_MISSING=$(wc -l < "$OUTPUT_FILE")
        echo "Missing files listed: $NUM_MISSING"
        echo "Output file: $OUTPUT_FILE"
    fi
else
    echo "ERROR: Script exited with code $EXIT_CODE"
fi
echo "End time: $(date)"
echo "=================================================="

exit $EXIT_CODE
