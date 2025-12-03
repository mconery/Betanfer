#!/bin/bash
#SBATCH -A bif154
#SBATCH -J filter_gwas
#SBATCH -o filter_gwas.log
#SBATCH -N 4
#SBATCH -p batch
#SBATCH -t 2:00:00
#SBATCH -C nvme

# ============================================================================
# CONFIGURATION - MODIFY AS NEEDED
# ============================================================================

# Path to the Python scripts
PYTHON_SCRIPT="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/6_Genotype_Heterogeneity/drop_variants.py"
COMMAND_GENERATOR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/6_Genotype_Heterogeneity/generate_filter_commands.py"

# Path to causal variants directory
CAUSAL_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants"

# Path to summary statistics root directory
RESULTS_DIR="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"

# Column names (should match your data)
MARKER_COL="MarkerID"
CHR_COL="CHR"
POS_COL="POS"

# Number of nodes
num_nodes="$SLURM_NNODES"

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Load modules
echo "Loading modules..."
module load miniforge3
module load amd/6.3.1
module load rocm/6.3.1
module load hwloc/2.11.1-gpu
module load flux

echo "Starting GWAS summary statistics filtering pipeline (flux-parallelized)"
echo "Command generator script: $COMMAND_GENERATOR"
echo "Python script: $PYTHON_SCRIPT"
echo "Causal variants directory: $CAUSAL_DIR"
echo "Results directory: $RESULTS_DIR"
    
# Validate paths
if [[ ! -f "$COMMAND_GENERATOR" ]]; then
	echo "Command generator script not found: $COMMAND_GENERATOR"
        exit 1
fi
    
if [[ ! -f "$PYTHON_SCRIPT" ]]; then
        echo "Python script not found: $PYTHON_SCRIPT"
        exit 1
fi
    
if [[ ! -d "$CAUSAL_DIR" ]]; then
        echo "Causal variants directory not found: $CAUSAL_DIR"
        exit 1
fi
    
if [[ ! -d "$RESULTS_DIR" ]]; then
        echo "Results directory not found: $RESULTS_DIR"
        exit 1
fi
  
# Command file
cmd_file="$RESULTS_DIR/drop_causal_commands.txt"
    
echo "Generating filter commands..."
    
# Run command generator
if ! python "$COMMAND_GENERATOR" \
        -r "$RESULTS_DIR" \
        -c "$CAUSAL_DIR" \
        -p "$PYTHON_SCRIPT" \
        -o "$cmd_file" \
        -m "$MARKER_COL" \
        --chr-col "$CHR_COL" \
        --pos-col "$POS_COL" \
        -v; then
	echo "Failed to generate filter commands"
        exit 1
fi
    
# Check if command file has content
if [[ ! -f "$cmd_file" ]] || [[ ! -s "$cmd_file" ]]; then
        echo "Command file is empty or not created"
        exit 1
fi
    
# Count commands
num_commands=$(wc -l < "$cmd_file")
echo "Generated $num_commands commands"
    
if [[ $num_commands -eq 0 ]]; then
        echo "No new commands to process. Exiting."
        exit 0
fi
    
echo "Starting flux scheduler..."
echo "Submitting $num_commands filter jobs to flux..."
    
# Launch flux scheduler with command execution
# Read commands into array and submit to flux
mapfile -t cmds < $cmd_file
srun -N "${num_nodes}" -n "${num_nodes}" -c 56 flux start "mapfile -t cmds < $cmd_file; flux resource list; for ((i=0; i<\${#cmds[@]}; i++)); do flux submit -N1 -n1 -c1 --output=flux_output/filter_\$i.txt \${cmds[\${i}]}; done; flux queue drain"
    
echo "Filtering completed"
echo "Check flux output files in: flux_output"

