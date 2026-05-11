#!/bin/bash

# run_susie_1000g.sh
#
# Per-locus SuSiE fine-mapping wrapper for MVP loci using 1000G LD matrices.
#
# 1. Tabix-extracts the locus region from the bgzip-compressed MVP summary
#    statistics into a temporary file.
# 2. Calls harmonize_1000g.py to match variants between the MVP sumstats
#    (ID column: CHR:POS:REF:ALT) and the 1000G LD matrix (.vars file:
#    CHR:POS:REF:ALT), handling allele flips by negating effect sizes.
#    Outputs an aligned sumstats file (VARIANT_ID, BETA, SE, N) and an
#    aligned LD submatrix (tab-delimited) alongside a .vars file.
# 3. Calls finemap_SuSiE.R on the pre-aligned inputs so it does not need to
#    load the full N×N LD matrix or perform its own variant filtering.
# 4. Cleans up all temporary files on exit.
#
# SuSiE output columns: VARIANT_ID  Z  PIP  converged  cs
# Timing file is written by finemap_SuSiE.R to {output_dir}/timing/{basename}.
#
# Usage:
#   bash run_susie_1000g.sh <sumstats_gz> <ld_file> <output_file> <CHR> <START> <END>
#
# Arguments:
#   sumstats_gz  Path to bgzip-compressed MVP sumstats (.txt.gz with .tbi index)
#   ld_file      Path to 1000G LD correlation matrix (.unphased.vcor1);
#                .unphased.vcor1.vars must exist alongside
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
## Tool locations (Aurora / flare paths)
################################################
SUSIE_SCRIPT="/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/SuSiE/finemap_SuSiE.R"
HARMONIZE_SCRIPT="$(dirname "$0")/harmonize_1000g.py"
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
TIMING_DIR="${OUTPUT_DIR}/timing"

mkdir -p "$OUTPUT_DIR" "$TIMING_DIR"

TEMP_SUMSTATS=$(mktemp /tmp/1000g_susie_sumstats_XXXXXX.txt)
TEMP_ALIGNED_STATS=$(mktemp /tmp/1000g_susie_aligned_stats_XXXXXX.tsv)
TEMP_LD=$(mktemp /tmp/1000g_susie_ld_XXXXXX.txt)
# harmonize_1000g.py writes TEMP_LD.vars automatically; include in cleanup
trap "rm -f '$TEMP_SUMSTATS' '$TEMP_ALIGNED_STATS' '$TEMP_LD' '${TEMP_LD}.vars'" EXIT

################################################
## Step 1: Tabix-extract locus region (all columns)
## harmonize_1000g.py needs the full row including ID, BETA, SE, N.
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
## flips by negating BETA, and writes a pre-aligned LD submatrix (tab-delimited)
## plus a .vars file alongside it (so finemap_SuSiE.R can locate it).
################################################
echo "NOTE: Harmonizing sumstats with 1000G LD matrix"
python "$HARMONIZE_SCRIPT" \
    --sumstats     "$TEMP_SUMSTATS" \
    --vars         "${LD_FILE}.vars" \
    --ld           "$LD_FILE" \
    --out_sumstats "$TEMP_ALIGNED_STATS" \
    --out_ld_susie "$TEMP_LD"

if [ ! -s "$TEMP_ALIGNED_STATS" ] || [ ! -s "$TEMP_LD" ]; then
    echo "ERROR: harmonize_1000g.py produced empty output"
    exit 1
fi

################################################
## Step 3: Run SuSiE
## finemap_SuSiE.R args: sumstats ld_file output_file varcol beta_col se_col n_col
## Reads ${TEMP_LD}.vars (written by harmonize_1000g.py) for variant order.
## Writes timing file to {output_dir}/timing/{basename(output_file)}.
################################################
echo "NOTE: Running SuSiE (varcol=VARIANT_ID, beta=BETA, se=SE, n=N)"
Rscript "$SUSIE_SCRIPT" \
    "$TEMP_ALIGNED_STATS" \
    "$TEMP_LD" \
    "$OUTPUT_FILE" \
    VARIANT_ID \
    BETA \
    SE \
    N

echo "SUCCESS: SuSiE completed → $OUTPUT_FILE"
