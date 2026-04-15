#!/bin/bash

# run_susie.sh
#
# Per-locus SuSiE fine-mapping wrapper for MVP population-specific analyses.
#
# 1. Tabix-extracts the locus region from the bgzip-compressed MVP summary
#    statistics into a temporary file (avoids loading the full 2-5 GB file
#    into R memory).
# 2. Calls finemap_SuSiE.R on the extracted locus, passing the MVP column
#    names (RSID, BETA, SE, N).  The R script filters to LD matrix variants,
#    runs susie_rss(), and writes output + timing files.
# 3. Cleans up the temporary file on exit.
#
# SuSiE output columns: VARIANT_ID  Z  PIP  converged  cs
# Timing file is written by finemap_SuSiE.R to {output_dir}/timing/{basename}.
#
# Usage:
#   bash run_susie.sh <sumstats_gz> <ld_file> <output_file> <CHR> <START> <END>
#
# Arguments:
#   sumstats_gz  Path to bgzip-compressed MVP sumstats (.txt.gz with .tbi index)
#   ld_file      Path to LD correlation matrix (.vcor1); .vars file must exist alongside
#   output_file  Full path for the final SuSiE output (.susie.txt)
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
## Tool location
################################################
SUSIE_SCRIPT="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/SuSiE/finemap_SuSiE.R"

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
## Setup directories and temp file
################################################
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
TIMING_DIR="${OUTPUT_DIR}/timing"

mkdir -p "$OUTPUT_DIR" "$TIMING_DIR"

# Unique temp file (cleaned up on exit via trap)
TEMP_SUMSTATS=$(mktemp /tmp/mvp_susie_sumstats_XXXXXX.txt)
trap "rm -f '$TEMP_SUMSTATS'" EXIT

################################################
## Step 1: Tabix-extract locus region, then select the 4 columns SuSiE needs.
## Two issues with passing the raw tabix output to R:
##   1. The header starts with '#ID'; R's read.table() treats '#' as a comment
##      character by default and silently discards the entire header line.
##   2. Several columns (e.g. EFFECT, P_HETERO) are blank for many rows;
##      R's default sep="" (whitespace) collapses consecutive tabs, making rows
##      appear shorter than the header and causing scan() to error.
## Piping through awk strips the '#', builds a column-index map, and outputs
## only RSID/BETA/SE/N — a clean 4-column file R can parse without issue.
################################################
echo "NOTE: Extracting locus ${CHR}:${START}-${END} from $(basename $SUMSTATS_GZ)"
tabix -h "$SUMSTATS_GZ" "${CHR}:${START}-${END}" | awk '
    BEGIN { OFS="\t" }
    NR==1 {
        sub(/^#/, "", $1)
        for (i=1; i<=NF; i++) idx[$i] = i
        print "RSID", "BETA", "SE", "N"
        next
    }
    { print $idx["RSID"], $idx["BETA"], $idx["SE"], $idx["N"] }
' > "$TEMP_SUMSTATS"

if [ ! -s "$TEMP_SUMSTATS" ]; then
    echo "ERROR: tabix extraction produced an empty file for ${CHR}:${START}-${END}"
    exit 1
fi

################################################
## Step 2: Run SuSiE
## finemap_SuSiE.R args: sumstats ld_file output_file varcol beta_col se_col n_col
## The R script writes the timing file to {output_dir}/timing/{basename(output_file)}
################################################
echo "NOTE: Running SuSiE (varcol=RSID, beta=BETA, se=SE, n=N)"
Rscript "$SUSIE_SCRIPT" \
    "$TEMP_SUMSTATS" \
    "$LD_FILE" \
    "$OUTPUT_FILE" \
    RSID \
    BETA \
    SE \
    N

echo "SUCCESS: SuSiE completed → $OUTPUT_FILE"
