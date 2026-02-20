#!/bin/bash
#SBATCH -A bif154
#SBATCH -J make_matrix
#SBATCH -o make_all_matrices.log
#SBATCH -N 10
#SBATCH -p batch
#SBATCH -t 1:00:00

#Activate conda directory
source ~/.bashrc
conda activate rsparsepro

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set directories
gwas_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_gwas"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/ld_matrices"

#Set general file locations
locus_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/temp_files/make_matrix.loci.txt"
stddev_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/ld_matrices/correlation_error_standard_deviations.csv"
noise_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/3_Fine-Mapping_Robustness_Testing/make_noisy_matrices.py"

#Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
############################################## Make Basic LD Matrices ##########################################
################################################################################################################
#Change directory into ld matrices
cd $matrix_dir

#Send loci to file
ls $gwas_dir | grep "_causal_per_locus_" > $locus_file
#Loop over the directories and make the raw ld matrices
cmds=()
while read directory; do
	plink_file=$gwas_dir/$directory/$directory"_out_filt_sim"
	chr=$(echo "$directory" | sed 's/.*chr\([0-9]\+\).*/\1/')
	start=$(echo "$directory" | grep -oE '[0-9]+' | tail -2 | head -1)
	end=$(echo "$directory" | grep -oE '[0-9]+' | tail -1)
	if [ ! -f $matrix_dir/$directory".unphased.vcor1" ]; then
		cmds+=("plink2 --bfile $plink_file --chr $chr --from-bp $start --to-bp $end --r-unphased square --out $matrix_dir/$directory")
	fi
done < $locus_file

#Launch commands
job_size=${#cmds[@]}
echo "JOB SIZE: $job_size"
for (( i=0; i<${#cmds[@]}; i++ )); do
	j=$(( $i + 1 ))
	srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
	if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then 
	    wait
        fi
done
wait

################################################################################################################
############################################## Make Noisy LD Matrices ##########################################
################################################################################################################
#Get list of trait counts from 1000G simulations
num_trait_levels=$(cut -d "," -f 1 $stddev_file | awk 'NR>1' | sort -u -V | wc -l)

#Loop over the directories and make the raw ld matrices
cmds=()
shopt -s nullglob
while read directory; do
	matrix_file=$matrix_dir/$directory".unphased.vcor1"
	#Check if noisy matrices need to be made for the current locus by counting number of trait files
	files=($matrix_dir/$directory.*_traits.unphased.vcor1)
	count_trait_matrix_files=${#files[@]}
	#Create command if needed for current locus
	if [ $num_trait_levels != $count_trait_matrix_files ]; then cmds+=("python $noise_script --matrix=$matrix_file --stddev_file=$stddev_file --outdir=$matrix_dir"); fi
done < $locus_file

#Calculate number of possible jobs
job_size=${#cmds[@]}
num_parallel_jobs=$((($SLURM_CPUS_ON_NODE / 2) * num_nodes))
echo "JOB SIZE: $job_size"
#Launch commands
for (( i=0; i<${#cmds[@]}; i++ )); do
	j=$(( $i + 1 ))
	srun -n1 -c1 --exclusive --job-name=$j ${cmds[$i]} &
	if [ $(($j % $num_parallel_jobs)) == 0 ] || [ $j == $job_size ]; then 
	    wait
        fi
done
wait
