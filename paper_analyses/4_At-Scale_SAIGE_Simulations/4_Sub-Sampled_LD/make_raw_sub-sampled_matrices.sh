#!/bin/bash

#SBATCH -A bif154
#SBATCH -J make_matrix
#SBATCH -o make_raw_sub-sampled_matrices.log
#SBATCH -N 20
#SBATCH -p batch
#SBATCH -t 1:00:00

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################

# Set input directories and files
var_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants/master_causal_list.txt"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample"
plink_dirs="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/AFR,/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/EUR"

# Script directories and files
script_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/scripts"

# Set script locations
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/make_raw_in-sample_matrix_command_writer.py"
dedup_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/remove_matrix_duplicates.py"
sample_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/4_Sub-Sampled_LD/sample_individuals.py"

# Set populations
populations=AFR,EUR

#Set sample sizes
sizes=1000,5000,10000,50000,100000,200000,300000

# Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
############################################ Sub-Sample the Population #########################################
################################################################################################################

# Load modules and activate conda directory
module load flux
module load miniforge3
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno

# Conver populations and sizes to arrays
IFS=',' read -ra pop_arr <<< "$populations"
IFS=',' read -ra plink_arr <<< "$plink_dirs"
IFS=',' read -ra size_arr <<< "$sizes"
#Run sampling script 
for ((i=0; i<${#pop_arr[@]}; i++)); do for size in ${size_arr[@]}; do
	python $sample_script ${plink_arr[$i]}/chr1.plink.maf.fam $size ${plink_arr[$i]}/sub-sample_"$size".list.txt
done; done

################################################################################################################
############################################## Make Basic LD Matrices ##########################################
################################################################################################################

# Create output directory for flux logs
mkdir -p flux_output

#Create output directories for 
for size in ${size_arr[@]}; do 
	mkdir -p $matrix_dir/sub-sample_"$size"
	for pop in ${pop_arr[@]}; do
		mkdir -p $matrix_dir/sub-sample_"$size"/$pop
done; done

# Generate individual shell scripts
echo "Generating individual shell scripts..."
for size in ${size_arr[@]}; do 
	script_list_file=$script_dir/script_list.sub-sample_"$size".txt
	keep_files=$(echo "$plink_dirs" | sed 's|[^,]\+|&/sub-sample_'"$size"'.list.txt|g')
	python $command_script \
    	--script_dir $script_dir \
    	--script_list_file $script_list_file \
    	--matrix_dir $matrix_dir/sub-sample_"$size" \
    	--plink_dirs $plink_dirs \
    	--populations $populations \
    	--variant_lists_file $var_file \
    	--variant-file-suffix ".variant_set" \
	--keep_individuals $keep_files \
    	--dedup_script $dedup_script
done

#Merge all the script files together
cat $script_dir/script_list.sub-sample_*.txt > $script_dir/all.script_list.sub-sample.txt

# Count number of scripts
num_scripts=$(wc -l < "$script_dir/all.script_list.sub-sample.txt")
echo "Generated $num_scripts individual shell scripts"

# Use flux to run matrix creation
echo "Starting flux and submitting jobs..."
mapfile -t script_paths < $script_dir/all.script_list.sub-sample.txt

srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start \
"flux resource list; mapfile -t script_paths < $script_dir/all.script_list.sub-sample.txt; for ((i=0; i<\${#script_paths[@]}; i++)); do flux submit -N1 -n1 -c4 --output=flux_output/make_sampled_raw_matrix_\$i.txt \${script_paths[\$i]}; done; flux queue drain;"

# Clean-up log and nosex files
for size in ${size_arr[@]}; do 
	echo "Cleaning up temporary files for $size sub-sample..."
	rm -f $matrix_dir/sub-sample_"$size"/*/*.nosex
	rm -f $matrix_dir/sub-sample_"$size"/*/*.log
	rm -f $matrix_dir/sub-sample_"$size"/scripts/*
	rm -f $matrix_dir/sub-sample_"$size"/*/*temp_dedup*
done

echo 'SUCCESS: All raw in-sample matrices created.'
