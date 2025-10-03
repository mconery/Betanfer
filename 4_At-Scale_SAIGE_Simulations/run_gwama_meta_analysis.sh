#!/bin/bash

#SBATCH -A bif154
#SBATCH -J gwama_meta_analysis
#SBATCH -o gwama_meta_analysis_%j.log
#SBATCH -N 20
#SBATCH -p batch
#SBATCH -t 2:00:00

##################################################################################################################
################################################### Set Needed Variables #########################################
##################################################################################################################

# Comma-separated population names
POP_NAMES="AFR,EUR"

# Comma-separated population directories matching these names
POP_DIRS="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/AFR,/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/EUR"

# Output directory for meta-analysis
OUTPUT_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"

# Command file to write GWAMA commands
COMMANDS_FILE="${OUTPUT_DIR}/meta_analysis_commands.txt"

# Path to the updated generate_meta_commands.py script
CMD_SCRIPT="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/generate_meta_commands.py"

# Number of nodes allocated by SLURM
NNODES=${SLURM_NNODES}

#Set analysis type
analysis_type=fixed

##################################################################################################################
############################################## Load Modules and Setup ###########################################
##################################################################################################################

module load python/3.10-miniforge3
module load flux

mkdir -p ${OUTPUT_DIR}

##################################################################################################################
################################################### Generate Commands ############################################
##################################################################################################################

echo "Generating meta-analysis commands with populations: ${POP_NAMES}"

python ${CMD_SCRIPT} \
    --pop_names ${POP_NAMES} \
    --pop_dirs ${POP_DIRS} \
    --output_dir ${OUTPUT_DIR} \
    --commands_file ${COMMANDS_FILE}
    --analysis_type $analysis_type

if [ ! -f ${COMMANDS_FILE} ] || [ ! -s ${COMMANDS_FILE} ]; then
    echo "Error: No commands generated or commands file is empty"
    exit 1
fi

NUM_COMMANDS=$(wc -l < ${COMMANDS_FILE})
echo "Generated ${NUM_COMMANDS} meta-analysis commands"

##################################################################################################################
################################################### Run Meta-Analysis #############################################
##################################################################################################################

echo "Starting meta-analysis using Flux..."
echo "Submitting ${NUM_COMMANDS} meta-analysis jobs to Flux..."
mapfile -t cmds < ${COMMANDS_FILE}
srun -N ${NNODES} -n ${NNODES} -c 56 flux start \
    "flux resource list; \
     mapfile -t cmds < ${COMMANDS_FILE}; \
     for ((i=0; i<${#cmds[@]}; i++)); do \
         flux submit -N1 -n1 -c4 --output=flux_output/meta_output_\$i.txt \${cmds[\$i]}; \
     done; \
     flux queue drain;" 
echo 'SUCCESS: All meta-analysis jobs completed.'
