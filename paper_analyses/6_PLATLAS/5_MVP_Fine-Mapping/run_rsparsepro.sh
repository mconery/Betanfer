#!/bin/bash

# run_rsparsepro.sh
#
# Per-locus RSparsePro fine-mapping wrapper for MVP population-specific analyses.
#
# 1. Tabix-extracts the locus region from the bgzip-compressed MVP summary
#    statistics into a temporary file.
# 2. Prepares a Z-scores file via prep_z-score_files.py (filtering to variants
#    present in the LD matrix .vars file).
# 3. Runs rsparsepro_ld.py to perform fine-mapping.
# 4. Renames/reformats the output (RSID → VARIANT_ID header fix).
# 5. Writes a timing file alongside the output.
# 6. Cleans up all temporary files on exit.
#
# Usage:
#   bash run_rsparsepro.sh <sumstats_gz> <ld_file> <output_file> <CHR> <START> <END>
#
# Arguments:
#   sumstats_gz  Path to bgzip-compressed MVP sumstats (.txt.gz with .tbi index)
#   ld_file      Path to LD correlation matrix (.vcor1); .vars file must exist alongside
#   output_file  Full path for the final RSparsePro output (.rsparsepro.txt)
#   CHR          Chromosome (integer, e.g. 10)
#   START        Locus start position
#   END          Locus end position
#
# ***Run in finemap conda environment***

set -euo pipefail

SUMSTATS_GZ="$1"
LD_FILE="$2"
OUTPUT_FILE="$3"
CHR="$4"
START="$5"
END="$6"

################################################
## Tool locations (Aurora / flare paths)
################################################
Z_SCRIPT="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/3_Fine-Mapping_Robustness_Testing/prep_z-score_files.py"
RSPARSEPRO_SCRIPT="/flare/GeomicVar/mconery/tools/RSparsePro_LD/src/rsparsepro_ld.py"
ALIGN_SCRIPT="$(dirname "$0")/align_zscore_and_ld.py"
# Full path to tabix: conda activate is unreliable in non-interactive parsl worker shells
TABIX="/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/tabix"

################################################
## Skip if output already complete
################################################
if [ -f "$OUTPUT_FILE" ]; then
    echo "NOTE: Output already exists, skipping: $OUTPUT_FILE"
    exit 0
fi

################################################
## Validate inputs
################################################
if [ ! -f "$SUMSTATS_GZ" ]; then
    echo "ERROR: sumstats file not found: $SUMSTATS_GZ"
    exit 1
fi
if [ ! -f "$LD_FILE" ]; then
    echo "ERROR: LD matrix file not found: $LD_FILE"
    exit 1
fi
if [ ! -f "${LD_FILE}.vars" ]; then
    echo "ERROR: .vars file not found: ${LD_FILE}.vars"
    exit 1
fi

################################################
## Setup directories and temp files
################################################
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
BASE_NAME=$(basename "$OUTPUT_FILE" .txt)   # strips one .txt
BASE_NAME=$(basename "$BASE_NAME" .rsparsepro)  # strips .rsparsepro → TRAIT.chrN_S_E
TIMING_DIR="${OUTPUT_DIR}/timing"

mkdir -p "$OUTPUT_DIR" "$TIMING_DIR"

# Unique temp files (cleaned up on exit via trap)
TEMP_SUMSTATS=$(mktemp /tmp/mvp_rsp_sumstats_XXXXXX.txt)
Z_FILE=$(mktemp /tmp/mvp_rsp_zscore_XXXXXX.tsv)
ALIGNED_Z=$(mktemp /tmp/mvp_rsp_aligned_z_XXXXXX.tsv)
ALIGNED_LD=$(mktemp /tmp/mvp_rsp_aligned_ld_XXXXXX.txt)
trap "rm -f '$TEMP_SUMSTATS' '$Z_FILE' '$ALIGNED_Z' '$ALIGNED_LD'" EXIT

################################################
## Step 1: Tabix-extract locus region
################################################
echo "NOTE: Extracting locus ${CHR}:${START}-${END} from $(basename $SUMSTATS_GZ)"
$TABIX -h "$SUMSTATS_GZ" "${CHR}:${START}-${END}" > "$TEMP_SUMSTATS"

if [ ! -s "$TEMP_SUMSTATS" ]; then
    echo "ERROR: tabix extraction produced an empty file for ${CHR}:${START}-${END}"
    exit 1
fi

################################################
## Step 2: Prepare Z-scores file
################################################
echo "NOTE: Preparing Z-scores file (varcol=RSID, beta=BETA, se=SE)"
python "$Z_SCRIPT" \
    -v "${LD_FILE}.vars" \
    -g "$TEMP_SUMSTATS" \
    -o "$Z_FILE" \
    --varcol=RSID \
    --beta=BETA \
    --se=SE

if [ ! -s "$Z_FILE" ]; then
    echo "ERROR: Z-scores file is empty — no variants matched between sumstats and .vars"
    exit 1
fi

################################################
## Step 2b: Align z-score order and subset LD matrix
################################################
# betanfer builds LD matrices with variants in lexicographic RSID order, but
# prep_z-score_files.py outputs in SE-ascending order.  align_zscore_and_ld.py
# reorders the z-score file to match the .vars order and subsets the LD matrix
# to the intersection of variants, so rsparsepro_ld.py receives matching inputs.
echo "NOTE: Aligning z-score order and subsetting LD matrix"
python "$ALIGN_SCRIPT" \
    --zscore     "$Z_FILE" \
    --vars       "${LD_FILE}.vars" \
    --ld         "$LD_FILE" \
    --out_zscore "$ALIGNED_Z" \
    --out_ld     "$ALIGNED_LD"

if [ ! -s "$ALIGNED_Z" ] || [ ! -s "$ALIGNED_LD" ]; then
    echo "ERROR: align_zscore_and_ld.py produced empty output"
    exit 1
fi

################################################
## Step 3: Run RSparsePro
################################################
# rsparsepro_ld.py --save BASE produces BASE.rsparsepro.txt and BASE.rsparsepro.log
OUTPUT_BASE="${OUTPUT_DIR}/${BASE_NAME}"

echo "NOTE: Running RSparsePro"
START_TIME=$(date +%s%N)
python "$RSPARSEPRO_SCRIPT" --z "$ALIGNED_Z" --ld "$ALIGNED_LD" --save "$OUTPUT_BASE"
END_TIME=$(date +%s%N)

ELAPSED_NS=$(( END_TIME - START_TIME ))
ELAPSED_SEC=$(echo "scale=3; $ELAPSED_NS / 1000000000" | bc)

################################################
## Step 4: Rename output and fix header
################################################
# rsparsepro produces {OUTPUT_BASE}.rsparsepro.txt; rename to OUTPUT_FILE
RSP_RAW="${OUTPUT_BASE}.rsparsepro.txt"
if [ ! -f "$RSP_RAW" ]; then
    echo "ERROR: Expected RSparsePro output not found: $RSP_RAW"
    exit 1
fi

mv "$RSP_RAW" "$OUTPUT_FILE"
# Rename variant ID column to VARIANT_ID (tool outputs RSID)
sed -i 's/^RSID\t/VARIANT_ID\t/' "$OUTPUT_FILE"

# Remove the log file
rm -f "${OUTPUT_BASE}.rsparsepro.log"

################################################
## Step 5: Write timing file
################################################
TIMING_FILE="${TIMING_DIR}/${BASE_NAME}.rsparsepro.txt"
printf "rsparsepro_time_seconds\n%s\n" "$ELAPSED_SEC" > "$TIMING_FILE"

echo "SUCCESS: RSparsePro completed in ${ELAPSED_SEC}s → $OUTPUT_FILE"
