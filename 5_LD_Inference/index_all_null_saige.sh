#!/bin/bash

#SBATCH -A bif154
#SBATCH -J index_saige_mpi
#SBATCH -o index_all_null_saige.log
#SBATCH -N 100
#SBATCH -q debug
#SBATCH -t 1:00:00

################################################################################################################
###################################### Load Modules and Environment ############################################
################################################################################################################

module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/infersparsepro

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set script locations
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/create_saige_index_commands.py"
indexer_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/5_LD_Inference/index_single_saige_file.py"

# Set input and output directories
eur_input_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/EUR"
afr_input_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/AFR"
eur_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_EUR_sumstats"
afr_scratch_dir="/lustre/orion/bif154/scratch/mconery/indexed_AFR_sumstats"

# Set command files
eur_commands_file=$eur_scratch_dir/eur_saige_index_commands.txt
afr_commands_file=$afr_scratch_dir/afr_saige_index_commands.txt

# Set number of nodes
num_nodes=${SLURM_NNODES}

# Create output directories if they don't exist
mkdir -p $eur_scratch_dir
mkdir -p $afr_scratch_dir
mkdir -p flux_output

################################################################################################################
##################################### Prep Commands for EUR population #########################################
################################################################################################################

echo "Generating indexing commands for EUR population..."
python $command_script \
    --scratch_dir $eur_scratch_dir \
    --pattern $eur_input_dir/"Pheno*.txt" \
    --commands_file $eur_commands_file \
    --filename_regex "^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.step2\\.txt$" \
    --indexer_script $indexer_script \
    --chrom_col CHR \
    --pos_col POS

# Check if EUR commands were generated
if [ ! -f "$eur_commands_file" ] || [ ! -s "$eur_commands_file" ]; then
    echo "ERROR: No EUR commands generated or file is empty: $eur_commands_file"
    exit 1
fi

eur_cmd_count=$(wc -l < $eur_commands_file)
echo "Generated $eur_cmd_count commands for EUR population"

################################################################################################################
##################################### Prep Commands for AFR population #########################################
################################################################################################################

echo "Generating indexing commands for AFR population..."
python $command_script \
    --scratch_dir $afr_scratch_dir \
    --pattern $afr_input_dir/"Pheno*.txt" \
    --commands_file $afr_commands_file \
    --filename_regex "^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.step2\\.txt$" \
    --indexer_script $indexer_script \
    --chrom_col CHR \
    --pos_col POS

# Check if AFR commands were generated
if [ ! -f "$afr_commands_file" ] || [ ! -s "$afr_commands_file" ]; then
    echo "ERROR: No AFR commands generated or file is empty: $afr_commands_file"
    exit 1
fi

afr_cmd_count=$(wc -l < $afr_commands_file)
echo "Generated $afr_cmd_count commands for AFR population"

################################################################################################################
################################## Execute EUR population indexing #############################################
################################################################################################################

echo "Starting EUR population indexing with $num_nodes nodes..."
echo "Processing $eur_cmd_count EUR files..."

mapfile -t eur_cmds < $eur_commands_file

srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start "flux resource list; mapfile -t eur_cmds < ${eur_commands_file}; for ((i=0; i<\${#eur_cmds[@]}; i++)); do flux submit -N1 -n1 -c2 --output=flux_output/index_eur_\$i.txt \${eur_cmds[\$i]}; done; flux queue drain;"

echo "EUR population indexing completed."

################################################################################################################
################################## Execute AFR population indexing #############################################
################################################################################################################

echo "Starting AFR population indexing with $num_nodes nodes..."
echo "Processing $afr_cmd_count AFR files..."

mapfile -t afr_cmds < $afr_commands_file

srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start "flux resource list; mapfile -t afr_cmds < ${afr_commands_file}; for ((i=0; i<\${#afr_cmds[@]}; i++)); do flux submit -N1 -n1 -c2 --output=flux_output/index_afr_\$i.txt \${afr_cmds[\$i]}; done; flux queue drain;"

echo "AFR population indexing completed."

