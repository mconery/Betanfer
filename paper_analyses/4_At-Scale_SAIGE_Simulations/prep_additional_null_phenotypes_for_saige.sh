#!/bin/bash
#SBATCH -A bif154
#SBATCH -J prep_geno
#SBATCH -o prep_additional_null_phenotypes_for_saige.log
#SBATCH -N 1
#SBATCH -p batch
#SBATCH -q debug
#SBATCH -t 1:00:00

#Load plink tools and conda environment
module load miniforge3
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set directories
inp_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_genotypes"
geno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"

#Set script locations
null_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/simulate_null_phenotypes.py"
simulate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/simulate_phenotype.py"
decode_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/3_Fine-Mapping_Robustness_Testing/execute_command.py"

#Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
############################################ Simulate Null Phenotypes ##########################################
################################################################################################################
#Create empty command array
cmds=()
#Loop over the populations
for pop in AFR EUR; do
	cd $pheno_dir/$pop
	cmds+=("python $null_script -i $geno_dir/$pop/merged.plink.maf.ld.pcs.eigenvec -o $pheno_dir/$pop/merged2.maf.ld.continuous.txt -n 1000 -t continuous -r 510 -s 1001")
	cmds+=("python $null_script -i $geno_dir/$pop/merged.plink.maf.ld.pcs.eigenvec -o $pheno_dir/$pop/merged2.maf.ld.binary.txt -n 1000 -t binary -r 510 -s 1001")
done

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
################################################# Split Phenotypes #############################################
################################################################################################################
#Make commands to split the files
cmds=()
#Loop over the populations
for pop in AFR EUR; do for trait_type in binary continuous; do cd $pheno_dir/$pop; for i in {1001..2000}; do
		pheno_col=$(( $i + 12 - 1000 ))
		if [ ! -f "$pheno_dir"/"$pop"/Pheno"$i"."$trait_type".txt ]; then 
			cmds+=("python $decode_script -c cut___-f___1-12,"$pheno_col"___"$pheno_dir"/"$pop"/merged2.maf.ld."$trait_type".txt___>___"$pheno_dir"/"$pop"/Pheno"$i"."$trait_type".txt")
		fi
done; done; done

#Calculate number of possible jobs
job_size=${#cmds[@]}
num_parallel_jobs=$((($SLURM_CPUS_ON_NODE / 4) * num_nodes))
echo "JOB SIZE: $job_size"
#Launch commands
for (( i=0; i<${#cmds[@]}; i++ )); do
	j=$(( $i + 1 ))
	srun -n1 -c1 -N1 --exclusive --job-name=$j ${cmds[$i]} &
	if [ $(($j % $num_parallel_jobs)) == 0 ] || [ $j == $job_size ]; then 
	    wait
        fi
done
wait

################################################################################################################
################################## Replace the Phenotype Name in the Split File ################################
################################################################################################################
#Make commands to split the files
cmds=()
#Loop over the populations
for pop in AFR EUR; do for trait_type in binary continuous; do cd $pheno_dir/$pop; for i in {1001..2000}; do
			cmds+=("python $decode_script -c sed___-i___'s/Pheno$i/PHENO/g'___"$pheno_dir"/"$pop"/Pheno"$i"."$trait_type".txt")
done; done; done

#Calculate number of possible jobs
job_size=${#cmds[@]}
num_parallel_jobs=$((($SLURM_CPUS_ON_NODE / 4) * num_nodes))
echo "JOB SIZE: $job_size"
#Launch commands
for (( i=0; i<${#cmds[@]}; i++ )); do
	j=$(( $i + 1 ))
	srun -n1 -c1 -N1 --exclusive --job-name=$j ${cmds[$i]} &
	if [ $(($j % $num_parallel_jobs)) == 0 ] || [ $j == $job_size ]; then 
	    wait
        fi
done
wait
