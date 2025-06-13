#!/bin/bash
#SBATCH -A bif154
#SBATCH -J map_loci
#SBATCH -o map_all_loci.log
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
ld_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/ld_matrices"
out_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/finemapped_results"

#Set general script and file locations
locus_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/temp_files/make_matrix.loci.txt"
z_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/3_Fine-Mapping_Robustness_Testing/prep_z-score_files.py"
rsparsepro_script="/lustre/orion/bif154/proj-shared/mconery/tools/RSparsePro_LD/src/rsparsepro_ld.py"

#Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
############################################ Make Z-Files for All Loci #########################################
################################################################################################################
#Change directory intooutput directory
cd $out_dir

#Loop over the directories and make the z-score file commands
cmds=()
while read directory; do
	locus_vars=$ld_dir/$directory".unphased.vcor1.vars"
	gwas_file=$gwas_dir/$directory/$directory"_out_snps_sim_gwas.tsv"
	z_file=$gwas_dir/$directory/$directory".z-scores.tsv"
	if [ ! -f $z_file ]; then
		cmds+=("python $z_script -v $locus_vars -g $gwas_file -o $z_file --varcol=rsid --beta=beta --se=se")
	fi
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

################################################################################################################
################################################## Map All Loci ################################################
################################################################################################################
#Loop over the directories and make the mapping commands
cmds=()
shopt -s nullglob
while read directory; do
	z_file=$gwas_dir/$directory/$directory".z-scores.tsv"
	ld_files=($ld_dir/$directory.*unphased.vcor1)
	for ld_file in ${ld_files[@]}; do 
		temp=$(basename $ld_file)
		out_file=$out_dir/${temp%.unphased.vcor1}
		if [ ! -f $out_file.rsparsepro.txt ]; then
			cmds+=("python $rsparsepro_script --z $z_file --ld $ld_file --save $out_file")
		fi
	done
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
