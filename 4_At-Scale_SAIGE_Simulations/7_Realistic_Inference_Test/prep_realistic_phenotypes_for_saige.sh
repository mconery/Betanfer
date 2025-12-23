#!/bin/bash

#SBATCH -A bif154
#SBATCH -J prep_realistic_pheno
#SBATCH -o prep_realistic_phenotypes_for_saige.log
#SBATCH -N 1
#SBATCH -p batch
#SBATCH -t 4:00:00

# Load plink tools and conda environment
module load miniforge3
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set directories
geno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"

# Set file path
variant_file="$pheno_dir/causal_variants/master_variants_list.txt"

# Set script location
simulate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/7_Realistic_Inference_Test/simulate_phenotypes_many_causal.py"

# Get number of nodes
num_nodes="$SLURM_NNODES"

# Phenotype simulation parameters
N_PHENOTYPES=2500
N_CAUSAL=500
H2_SNP_MIN=5e-5
H2_SNP_MAX=1.5e-4
SEED=42

################################################################################################################
############################### Simulate Realistic Phenotypes ################################################
################################################################################################################

echo "Starting phenotype simulations for AFR and EUR populations..."
echo "Parameters: N_PHENOTYPES=$N_PHENOTYPES, N_CAUSAL=$N_CAUSAL, H2_SNP_MIN=$H2_SNP_MIN, H2_SNP_MAX=$H2_SNP_MAX"

# Create empty command array
cmds=()

# Loop over the populations
for pop in AFR EUR; do
    echo "Preparing simulation for population: $pop"
    
    # Create output directory if it doesn't exist
    mkdir -p "$pheno_dir/$pop"
    cd "$pheno_dir/$pop"
    
    # Get PC file path based on population
    pc_file="$geno_dir/$pop/merged.plink.maf.ld.pcs.eigenvec"
    
    # Intermediate output file (all phenotypes in one file)
    intermediate_output="$pheno_dir/$pop/merged.realistic.continuous.txt"
    
    # Add simulation command
    cmds+=("python $simulate_script \
        --plink $geno_dir/$pop/merged.plink.maf.ld \
        --output $intermediate_output \
        --variants $variant_file \
        --n-phenotypes $N_PHENOTYPES \
        --n-causal $N_CAUSAL \
        --h2-snp-min $H2_SNP_MIN \
        --h2-snp-max $H2_SNP_MAX \
        --trait-type continuous \
        --pc-file $pc_file \
        --n-jobs 128 \
        --seed $SEED \
        --verbose")
done

# Launch commands
job_size=${#cmds[@]}
echo "JOB SIZE: $job_size"

for (( i=0; i<${#cmds[@]}; i++ )); do
    j=$(( $i + 1 ))
    srun -n1 -N1 --ntasks-per-node=1 --job-name=$j bash -c "${cmds[$i]}" &
    
    if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then
        wait
    fi
done

wait

echo "Phenotype simulations completed."

################################################################################################################
#################################### Split Phenotypes into Individual Files ##################################
################################################################################################################

echo "Splitting merged phenotype files into individual phenotype files..."

# Create commands to split the phenotypes
cmds=()

# Loop over populations and extract individual phenotypes
for pop in AFR EUR; do
    cd "$pheno_dir/$pop"
    merged_file="$pheno_dir/$pop/merged.realistic.continuous.txt"
    
    # Loop over phenotype indices (1 to N_PHENOTYPES)
    for i in $(seq 1 $N_PHENOTYPES); do
        # Column index: phenotype columns start at column 13 (after #FID, IID, PC1-PC10)
        # PHENO_0 is column 13, PHENO_1 is column 14, etc.
        pheno_col=$(( $i + 12 ))
        
        output_file="$pheno_dir/$pop/Pheno$i.realistic_trait.continuous.txt"
        
        # Only create command if output doesn't already exist
        if [ ! -f "$output_file" ]; then
            cmds+=("cut -f 1-12,$pheno_col $merged_file > $output_file")
        fi
    done
done

# Calculate number of possible jobs
job_size=${#cmds[@]}
num_parallel_jobs=$((($SLURM_CPUS_ON_NODE / 4) * num_nodes))

echo "JOB SIZE (file splitting): $job_size"
echo "Parallel jobs: $num_parallel_jobs"

# Launch commands
for (( i=0; i<${#cmds[@]}; i++ )); do
    j=$(( $i + 1 ))
    srun -n1 -c1 -N1 --exclusive --job-name=$j bash -c "${cmds[$i]}" &
    
    if [ $(($j % $num_parallel_jobs)) == 0 ] || [ $j == $job_size ]; then
        wait
    fi
done

wait

echo "File splitting completed."

################################################################################################################
##################### Replace Phenotype Column Names in Individual Files ######################################
################################################################################################################

echo "Renaming phenotype columns in individual phenotype files..."

# Create commands to replace phenotype column names
cmds=()

# Loop over populations and rename phenotype columns
for pop in AFR EUR; do
    cd "$pheno_dir/$pop"
    
    # Loop over phenotype indices
    for i in $(seq 1 $N_PHENOTYPES); do
        output_file="$pheno_dir/$pop/Pheno$i.realistic_trait.continuous.txt"
        
        # Replace the column name PHENO_X with PHENO
        cmds+=("sed -i 's/PHENO_[0-9]*\$/PHENO/' $output_file")
    done
done

# Calculate number of possible jobs
job_size=${#cmds[@]}
num_parallel_jobs=$((($SLURM_CPUS_ON_NODE / 4) * num_nodes))

echo "JOB SIZE (renaming): $job_size"

# Launch commands
for (( i=0; i<${#cmds[@]}; i++ )); do
    j=$(( $i + 1 ))
    srun -n1 -c1 -N1 --exclusive --job-name=$j bash -c "${cmds[$i]}" &
    
    if [ $(($j % $num_parallel_jobs)) == 0 ] || [ $j == $job_size ]; then
        wait
    fi
done

wait

echo "Phenotype renaming completed."

################################################################################################################
##################################### Summary and Verification ##############################################
################################################################################################################

echo ""
echo "=========================================="
echo "Phenotype Generation Summary"
echo "=========================================="

for pop in AFR EUR; do
    echo ""
    echo "Population: $pop"
    echo "  Total phenotype files: $(ls $pheno_dir/$pop/Pheno*.realistic_trait.continuous.txt 2>/dev/null | wc -l)"
    
    # Check first phenotype file
    first_file=$(ls $pheno_dir/$pop/Pheno1.realistic_trait.continuous.txt 2>/dev/null)
    if [ -f "$first_file" ]; then
        echo "  First phenotype dimensions: $(wc -l < $first_file) samples"
        echo "  Sample columns: $(head -1 $first_file | cut -d$'\t' -f1-12)"
    fi
    
    # Check metadata
    metadata_file="$pheno_dir/$pop/merged.realistic.continuous_metadata.txt"
    if [ -f "$metadata_file" ]; then
        echo "  Metadata file: $metadata_file"
    fi
done

echo ""
echo "=========================================="
echo "Pipeline Complete!"
echo "=========================================="
