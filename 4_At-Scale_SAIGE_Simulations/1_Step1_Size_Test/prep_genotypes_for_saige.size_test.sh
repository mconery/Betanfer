#!/bin/bash
#SBATCH -A bif154
#SBATCH -J prep_geno
#SBATCH -o prep_genotypes_for_saige.size_test.log
#SBATCH -N 20
#SBATCH -p batch
#SBATCH -t 2:00:00

#Load plink tools and conda environment
source ~/.bashrc
conda activate simulate_pheno

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set directories
inp_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_genotypes"
geno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"

#Set script locations
simulate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/simulate_null_phenotypes.py"
decode_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/3_Fine-Mapping_Robustness_Testing/execute_command.py"

#Get number of nodes
num_nodes="$SLURM_NNODES"

#Create and array of pruning degrees
prune_mults=(2 3 4)

################################################################################################################
############################## Filter for MAF and Make LD Pruning Filters ######################################
################################################################################################################
##Create empty command array
cmds=()
#Loop over two populations and different ld multiples for pruning
for pop in AFR EUR; do for prune_mult in ${prune_mults[@]}; do
	#Change into population output directory
	mkdir -p $geno_dir/$pop
	cd $geno_dir/$pop
	#Get keep file location
	keep_file=$inp_dir/$pop/$pop.400K_random_sample.tsv
	#Loop the chromosomes
	for chromo in {1..22}; do
		#Set plink file path
		plink_inp=$inp_dir/$pop/chr"$chromo".plink.tag
		if [ $pop == "EUR" ]; then 
			prune_degree=$(echo "0.1 * $prune_mult" | bc)
			plink_out=$geno_dir/$pop/chr"$chromo".plink.prune_mult_"$prune_mult".maf
			cmds+=("plink2 --bfile $plink_inp --keep $keep_file --set-all-var-ids "@:#:\$r:\$a" --rm-dup --maf 0.01 --indep-pairwise 1000 80 $prune_degree --make-bed --out $plink_out")
		else
			prune_degree=$(echo "0.05 * $prune_mult" | bc)
			plink_out=$geno_dir/$pop/chr"$chromo".plink.prune_mult_"$prune_mult".maf
			cmds+=("plink2 --bfile $plink_inp --keep $keep_file --set-all-var-ids "@:#:\$r:\$a" --rm-dup --maf 0.01 --indep-pairwise 1000 80 $prune_degree --make-bed --out $plink_out")
		fi
	done
done; done

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
################################################# Prune the Variants ###########################################
################################################################################################################
#Create empty command array
cmds=()
#Loop over two populations and different ld multiples for pruning
for pop in AFR EUR; do for prune_mult in ${prune_mults[@]}; do
	#Change into population output directory
	mkdir -p $geno_dir/$pop
	cd $geno_dir/$pop
	#Get keep file location
	keep_file=$inp_dir/$pop/$pop.400K_random_sample.tsv
	#Loop the chromosomes
	for chromo in {1..22}; do
		#Set plink file path
		plink_inp=$inp_dir/$pop/chr"$chromo".plink.tag
		plink_out=$geno_dir/$pop/chr"$chromo".plink.prune_mult_"$prune_mult".maf
		#Launch plink to prune variants
		cmds+=("plink2 --bfile $plink_out --extract $plink_out.prune.in --make-bed --out $plink_out.ld")
	done
done; done

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
###################################### Make Merged Files for LD Pruning Degrees ################################
################################################################################################################
#Create empty command array
cmds=()
#Loop over the populations and create commands for LD-pruned files
for pop in AFR EUR; do for prune_mult in ${prune_mults[@]}; do
	cd $geno_dir/$pop
	#Create merge list and merge files
	ls $geno_dir/$pop/chr*.plink.prune_mult_"$prune_mult".maf.ld.bim | sed 's/.bim//g' > $geno_dir/$pop/$pop.prune_mult_"$prune_mult".mergelist.txt
	cmds+=("plink2 --allow-no-sex --make-bed --pmerge-list $geno_dir/$pop/$pop.prune_mult_"$prune_mult".mergelist.txt bfile --out $geno_dir/$pop/merged.plink.prune_mult_"$prune_mult".maf.ld")
done; done

#Make commands for merged files with no LD pruning
#Loop over the populations
for pop in AFR EUR; do
	cd $geno_dir/$pop
	#Create merge list and merge files
	ls $geno_dir/$pop/chr*.plink.maf.bim | sed 's/.bim//g' > $geno_dir/$pop/$pop.no_ld_prune.mergelist.txt
	cmds+=("plink2 --allow-no-sex --make-bed --pmerge-list $geno_dir/$pop/$pop.no_ld_prune.mergelist.txt bfile --out $geno_dir/$pop/merged.plink.no_ld_prune.maf")
done;

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
