#!/bin/bash
#SBATCH -A bif154
#SBATCH -J make_matrix
#SBATCH -o make_all_matrices.log
#SBATCH -N 20
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

#Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
###################################### Make Commands to Make the Matrices ######################################
################################################################################################################
#Change directory into ld matrices
cd $matrix_dir

#Send loci to file
ls $gwas_dir | grep "_causal_per_locus_" > $locus_file
#Loop over the directories and make the ld matrices
cmds=()
while read directory; do
	plink_file=$gwas_dir/$directory/$directory"_out_filt_sim"
	chr=$(echo "$directory" | sed 's/.*chr\([0-9]\+\).*/\1/')
	start=$(echo "$directory" | grep -oE '[0-9]+' | tail -2 | head -1)
	end=$(echo "$directory" | grep -oE '[0-9]+' | tail -1)
	if [ ! -f "$matrix_dir/$directory".unphased.vcor1" ]; then
		cmds+=("plink2 --bfile $plink_file --chr $chr --from-bp $start --to-bp $end --r-unphased square --out $matrix_dir/$directory")
	fi
done < $locus_file
rm $locus_file

#Launch commands
echo "JOB SIZE: ${#cmds[@]}"
for (( i=1; i<=${#cmds[@]}; i++ )); do
	srun -n1 -N1 --ntasks-per-node=1 --job-name=$i ${cmds[$i]} &
	if [ $(($i % $num_nodes)) == 0 ]; then 
	    wait
        fi
done
wait