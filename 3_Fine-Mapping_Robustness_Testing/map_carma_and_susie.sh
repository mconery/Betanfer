#!/bin/bash
#SBATCH -A bif154
#SBATCH -J map_carma_and_susie
#SBATCH -o map_carma_and_susie.log
#SBATCH -N 2000
#SBATCH -p batch
#SBATCH -t 1:00:00

#Activate conda directory
source ~/.bashrc

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set directories
gwas_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_gwas"
ld_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/ld_matrices"
out_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/finemapped_results"

#Set general script and file locations
locus_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/temp_files/make_matrix.loci.txt"
susie_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/3_Fine-Mapping_Robustness_Testing/finemap_SuSiE.R"
carma_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/3_Fine-Mapping_Robustness_Testing/finemap_CARMA.R"

#Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
############################################### SuSiE Fine-Mapping #############################################
################################################################################################################
#Change directory into output directory
cd $out_dir

#Activate the SuSiE Conda environment
conda activate susie

#Loop over the directories and make the mapping commands
cmds=()
shopt -s nullglob
while read directory; do
	summary_stats_file=$gwas_dir/$directory/$directory"_out_snps_sim_gwas.tsv"
	pheno_file=$gwas_dir/$directory/$directory"_out_snps_sim_phenos.tsv"
	ld_files=($ld_dir/$directory.*unphased.vcor1)
	for ld_file in ${ld_files[@]}; do 
		temp=$(basename $ld_file)
		out_file=$out_dir/${temp%.unphased.vcor1}".susie.txt"
		if [[ $ld_file == *"_traits"* ]]; then
			in_sample_flag=FALSE
		else
			in_sample_flag=TRUE
		fi
		if [ ! -f $out_file ]; then
			cmds+=("Rscript $susie_script $summary_stats_file $pheno_file $ld_file $out_file $in_sample_flag")
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

################################################################################################################
############################################### CARMA Fine-Mapping #############################################
################################################################################################################
#Activate the CARMA Conda environment
conda deactivate
conda activate carma

#Loop over the directories and make the mapping commands
cmds=()
shopt -s nullglob
while read directory; do
	summary_stats_file=$gwas_dir/$directory/$directory"_out_snps_sim_gwas.tsv"
	ld_files=($ld_dir/$directory.*unphased.vcor1)
	for ld_file in ${ld_files[@]}; do 
		temp=$(basename $ld_file)
		out_file=$out_dir/${temp%.unphased.vcor1}".carma.txt"
		if [ ! -f $out_file ]; then
			cmds+=("Rscript $carma_script $summary_stats_file $ld_file $out_file")
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
