#!/bin/bash

# run_rsparsepro_1000g.sh
#
# Per-locus RSparsePro fine-mapping wrapper for MVP loci using 1000G LD matrices.
#
# 1. Tabix-extracts the locus region from the bgzip-compressed MVP summary
#    statistics into a temporary file.
# 2. Calls harmonize_1000g.py to match variants between the MVP sumstats
#    (ID column: CHR:POS:REF:ALT) and the 1000G LD matrix (.vars file:
#    CHR:POS:REF:ALT), handling allele flips by negating z-scores.  Outputs
#    an aligned z-score file and an aligned LD submatrix.
# 3. Runs rsparsepro_ld.py on the aligned inputs.
# 4. Renames/reformats the output (RSID → VARIANT_ID header fix).
# 5. Writes a timing file alongside the output.
# 6. Cleans up all temporary files on exit.
#
# Usage:
#   bash run_rsparsepro_1000g.sh <sumstats_gz> <ld_file> <output_file> <CHR> <START> <END>
#
# Arguments:
#   sumstats_gz  Path to bgzip-compressed MVP sumstats (.txt.gz with .tbi index)
#   ld_file      Path to 1000G LD correlation matrix (.unphased.vcor1);
#                .unphased.vcor1.vars must exist alongside
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
RSPARSEPRO_SCRIPT="/flare/GeomicVar/mconery/tools/RSparsePro_LD/src/rsparsepro_ld.py"
HARMONIZE_SCRIPT="$(dirname "$0")/harmonize_1000g.py"
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
BASE_NAME=$(basename "$OUTPUT_FILE" .txt)
BASE_NAME=$(basename "$BASE_NAME" .rsparsepro)  # → TRAIT.chrN_S_E
TIMING_DIR="${OUTPUT_DIR}/timing"

mkdir -p "$OUTPUT_DIR" "$TIMING_DIR"

TEMP_SUMSTATS=$(mktemp /tmp/1000g_rsp_sumstats_XXXXXX.txt)
TEMP_ZSCORE=$(mktemp /tmp/1000g_rsp_zscore_XXXXXX.tsv)
TEMP_LD=$(mktemp /tmp/1000g_rsp_ld_XXXXXX.txt)
trap "rm -f '$TEMP_SUMSTATS' '$TEMP_ZSCORE' '$TEMP_LD'" EXIT

################################################
## Step 1: Tabix-extract locus region (all columns)
################################################
echo "NOTE: Extracting locus ${CHR}:${START}-${END} from $(basename $SUMSTATS_GZ)"
$TABIX -h "$SUMSTATS_GZ" "${CHR}:${START}-${END}" > "$TEMP_SUMSTATS"

if [ ! -s "$TEMP_SUMSTATS" ]; then
    echo "ERROR: tabix extraction produced an empty file for ${CHR}:${START}-${END}"
    exit 1
fi

################################################
## Step 2: Harmonize sumstats with 1000G LD matrix
## Matches MVP sumstats ID (CHR:POS:REF:ALT) to .vars IDs, handles allele
## flips by negating z-scores, and subsets the LD matrix to matched variants.
################################################
echo "NOTE: Harmonizing sumstats with 1000G LD matrix"
python "$HARMONIZE_SCRIPT" \
    --sumstats "$TEMP_SUMSTATS" \
    --vars     "${LD_FILE}.vars" \
    --ld       "$LD_FILE" \
    --out_zscore "$TEMP_ZSCORE" \
    --out_ld_rsp "$TEMP_LD"

if [ ! -s "$TEMP_ZSCORE" ] || [ ! -s "$TEMP_LD" ]; then
    echo "ERROR: harmonize_1000g.py produced empty output"
    exit 1
fi

################################################
## Step 3: Run RSparsePro
## rsparsepro_ld.py --save BASE produces BASE.rsparsepro.txt and BASE.rsparsepro.log
################################################
OUTPUT_BASE="${OUTPUT_DIR}/${BASE_NAME}"

echo "NOTE: Running RSparsePro"
START_TIME=$(date +%s%N)
python "$RSPARSEPRO_SCRIPT" --z "$TEMP_ZSCORE" --ld "$TEMP_LD" --save "$OUTPUT_BASE"
END_TIME=$(date +%s%N)

ELAPSED_NS=$(( END_TIME - START_TIME ))
ELAPSED_SEC=$(echo "scale=3; $ELAPSED_NS / 1000000000" | bc)

################################################
## Step 4: Fix header in output file
## rsparsepro_ld.py --save BASE writes directly to BASE.rsparsepro.txt,
## which equals OUTPUT_FILE, so no rename is needed.
################################################
if [ ! -f "$OUTPUT_FILE" ]; then
    echo "ERROR: Expected RSparsePro output not found: $OUTPUT_FILE"
    exit 1
fi

# Rename variant ID column header (tool outputs RSID; standardise to VARIANT_ID)
sed -i 's/^RSID\t/VARIANT_ID\t/' "$OUTPUT_FILE"

# Remove the log file
rm -f "${OUTPUT_BASE}.rsparsepro.log"

################################################
## Step 5: Write timing file
################################################
TIMING_FILE="${TIMING_DIR}/${BASE_NAME}.rsparsepro.txt"
printf "rsparsepro_time_seconds\n%s\n" "$ELAPSED_SEC" > "$TIMING_FILE"

echo "SUCCESS: RSparsePro completed in ${ELAPSED_SEC}s → $OUTPUT_FILE"
