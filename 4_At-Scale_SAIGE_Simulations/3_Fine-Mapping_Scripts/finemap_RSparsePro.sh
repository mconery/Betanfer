#!/bin/bash

# Single-locus RSparsePro fine-mapping script
# Modified from the original map_rsparsepro.sh to process a single locus
# Usage: ./finemapping_rsparsepro.sh [gwas_file] [ld_file] [output_file] [varcol] [zscore_col] [beta_col] [se_col]

################################################################################################################
########################## Parse command line arguments ######################################################
################################################################################################################

# Check if minimum required arguments are provided
if [ $# -lt 4 ]; then
    echo "Usage: $0 <gwas_file> <ld_file> <output_file> [varcol] [zscore_col] [beta_col] [se_col]"
    echo ""
    echo "Required arguments:"
    echo "  gwas_file     - Path to GWAS summary statistics file"
    echo "  ld_file       - Path to LD matrix file (.vcor1 format)"
    echo "  output_file   - Full path for the final output file"
    echo ""
    echo "Optional arguments (provide either zscore_col OR both beta_col and se_col):"
    echo "  varcol        - Column name for variant IDs (default: rsid)"
    echo "  zscore_col    - Column name for Z-scores (if available)"
    echo "  beta_col      - Column name for beta coefficients (default: beta)"
    echo "  se_col        - Column name for standard errors (default: se)"
    echo ""
    echo "Examples:"
    echo "  $0 locus1 gwas.tsv matrix.vcor1 results/locus1.rsparsepro.txt"
    echo "  $0 locus1 gwas.tsv matrix.vcor1 results/locus1.rsparsepro.txt SNP Z"
    echo "  $0 locus1 gwas.tsv matrix.vcor1 results/locus1.rsparsepro.txt SNP "" beta se"
    exit 1
fi

# Required arguments
GWAS_FILE="$1"
LD_FILE="$2"
OUTPUT_FILE="$3"

# Optional arguments with defaults (matching typical finemap parameter patterns)
VARCOL_NAME="${4:-rsid}"
ZSCORE_COL="${5:-}"
BETA_COL="${6:-beta}"
SE_COL="${7:-se}"

################################################################################################################
########################## Define script locations ############################################################
################################################################################################################

# Set script locations (adjust paths as needed)
Z_SCRIPT="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/3_Fine-Mapping_Robustness_Testing/prep_z-score_files.py"
RSPARSEPRO_SCRIPT="/lustre/orion/bif154/proj-shared/mconery/tools/RSparsePro_LD/src/rsparsepro_ld.py"

################################################################################################################
########################## Validate inputs and setup ##########################################################
################################################################################################################

# Check if input files exist
if [ ! -f "$GWAS_FILE" ]; then
    echo "ERROR: GWAS file not found: $GWAS_FILE"
    exit 1
fi

if [ ! -f "$LD_FILE" ]; then
    echo "ERROR: LD file not found: $LD_FILE"
    exit 1
fi

# Extract output directory from output file path
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Create z-score_files subdirectory
Z_SCORE_DIR="${OUTPUT_DIR}/z-score_files"
mkdir -p "$Z_SCORE_DIR"

# Generate variant list file from LD matrix
VARS_FILE="${LD_FILE}.vars"
if [ ! -f "$VARS_FILE" ]; then
    echo "ERROR: Variant list file not found: $VARS_FILE"
    echo "Please ensure the .vars file exists alongside the LD matrix file"
    exit 1
fi

# Define z-scores file in the subdirectory
Z_SCORES_FILE="${Z_SCORE_DIR}/$(basename $OUTPUT_FILE)"

################################################################################################################
########################## Step 1: Generate Z-scores file ####################################################
################################################################################################################

echo "=== Step 1: Preparing Z-scores file ==="
echo "GWAS file: $GWAS_FILE"
echo "Variants file: $VARS_FILE"
echo "Z-scores output: $Z_SCORES_FILE"

# Build the z-score preparation command
if [ -n "$ZSCORE_COL" ]; then
    # Use existing Z-score column
    Z_CMD="python $Z_SCRIPT -v $VARS_FILE -g $GWAS_FILE -o $Z_SCORES_FILE --varcol=$VARCOL_NAME --zscore=$ZSCORE_COL"
else
    # Calculate Z-scores from beta and SE
    Z_CMD="python $Z_SCRIPT -v $VARS_FILE -g $GWAS_FILE -o $Z_SCORES_FILE --varcol=$VARCOL_NAME --beta=$BETA_COL --se=$SE_COL"
fi

echo "Command: $Z_CMD"

# Execute z-score preparation
if [ ! -f "$Z_SCORES_FILE" ]; then
    eval $Z_CMD
    if [ $? -ne 0 ]; then
        echo "ERROR: Z-score preparation failed"
        exit 1
    fi
    echo "Z-scores file created successfully"
else
    echo "Z-scores file already exists, skipping generation"
fi

################################################################################################################
########################## Step 2: Run RSparsePro fine-mapping ###############################################
################################################################################################################

echo ""
echo "=== Step 2: Running RSparsePro fine-mapping ==="
echo "Z-scores file: $Z_SCORES_FILE"
echo "LD matrix file: $LD_FILE"
echo "Final output file: $OUTPUT_FILE"

# Extract the base name for RSparsePro save parameter (without extension)
OUTPUT_BASE=$(basename "$OUTPUT_FILE" .txt)
OUTPUT_BASE=$(basename "$OUTPUT_BASE" .rsparsepro)

# Build the RSparsePro command
RSPARSEPRO_CMD="python $RSPARSEPRO_SCRIPT --z $Z_SCORES_FILE --ld $LD_FILE --save ${OUTPUT_DIR}/${OUTPUT_BASE}"

echo "Command: $RSPARSEPRO_CMD"

# Execute RSparsePro fine-mapping
if [ ! -f "$OUTPUT_FILE" ]; then
    eval $RSPARSEPRO_CMD
    if [ $? -ne 0 ]; then
        echo "ERROR: RSparsePro fine-mapping failed"
        exit 1
    fi
    echo "RSparsePro fine-mapping completed successfully"
else
    echo "RSparsePro output file already exists, skipping analysis"
fi
