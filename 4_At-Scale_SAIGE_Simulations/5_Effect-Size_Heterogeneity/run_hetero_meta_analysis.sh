#!/bin/bash

#SBATCH -A bif154
#SBATCH -J hetero_meta
#SBATCH -o run_hetero_meta_analysis.log
#SBATCH -N 40
#SBATCH -p batch
#SBATCH -t 2:00:00

##################################################################################################################
################################################### Set Needed Variables #########################################
##################################################################################################################

# Comma-separated population names
POP_NAMES="AFR,EUR"

# Comma-separated population directories matching these names
POP_DIRS="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/AFR,/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/EUR"

# Combinations of comma-separated per-variant heritabilities for each population
HERITABILITY_SETS=("0.00005,0.0003" "0.0003,0.00005" "0,0.0003" "0.0003,0" "0,0.00005" "0.00005,0")

# Output directory for meta-analysis
OUTPUT_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"

# Command file to write GWAMA commands
COMMANDS_DIR="${OUTPUT_DIR}/command_files"
MASTER_COMMANDS_FILE=$COMMANDS_DIR/master.commands_file.txt

# Path to the updated generate_meta_commands.py script
CMD_SCRIPT="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/5_Effect-Size_Heterogeneity/generate_hetero_meta_commands.py"

# Number of nodes allocated by SLURM
NNODES=${SLURM_NNODES}

# Set analysis type
analysis_type=random

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

#Loop over heritability sets and make command files
for HERITABILITIES in ${HERITABILITY_SETS[@]}; do
	heritabilities=$(echo $HERITABILITIES | sed 's/,/_/g')
	COMMANDS_FILE=$COMMANDS_DIR/heritability_"$heritabilities".commands_file.txt
	python ${CMD_SCRIPT} \
	--pop_names ${POP_NAMES} \
	--pop_dirs ${POP_DIRS} \
	--heritabilities ${HERITABILITIES} \
	--output_dir ${OUTPUT_DIR} \
	--commands_file ${COMMANDS_FILE} \
	--analysis_type $analysis_type
done

#Combine all commands files
cat $COMMANDS_DIR/heritability_*.commands_file.txt > $MASTER_COMMANDS_FILE

if [ ! -f ${MASTER_COMMANDS_FILE} ] || [ ! -s ${MASTER_COMMANDS_FILE} ]; then
echo "Error: No commands generated or commands file is empty"
exit 1
fi

NUM_COMMANDS=$(wc -l < ${MASTER_COMMANDS_FILE})
echo "Generated ${NUM_COMMANDS} meta-analysis commands"

##################################################################################################################
################################################### Run Meta-Analysis #############################################
##################################################################################################################

echo "Starting meta-analysis using Flux..."
echo "Submitting ${NUM_COMMANDS} meta-analysis jobs to Flux..."
mapfile -t cmds < ${MASTER_COMMANDS_FILE}
srun -N ${NNODES} -n ${NNODES} -c 56 flux start \
    "flux resource list; \
     mapfile -t cmds < ${MASTER_COMMANDS_FILE}; \
     for ((i=0; i<${#cmds[@]}; i++)); do \
         flux submit -N1 -n1 -c4 --output=flux_output/meta_output_\$i.txt \${cmds[\$i]}; \
     done; \
     flux queue drain;" 
echo 'SUCCESS: All meta-analysis jobs completed.'