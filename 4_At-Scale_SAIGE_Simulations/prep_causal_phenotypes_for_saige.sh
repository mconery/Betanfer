#!/bin/bash
#SBATCH -A bif154
#SBATCH -J prep_geno
#SBATCH -o prep_causal_phenotypes_for_saige.log
#SBATCH -N 20
#SBATCH -p batch
#SBATCH -t 1:00:00
#SBATCH -C nvme

#Load plink tools and conda environment
source ~/.bashrc
conda activate simulate_pheno

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set directories
geno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"

#Set script locations
simulate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/simulate_phenotype.py"
selection_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/select_causal_variants.py"

#Set key variables for causal variant selection
set_numbers=(1 2 3)
m_beta=0.05
max_distance=500000
min_spacing=100000

#Set key variables for phenotype simulation
trait_types=(binary continuous)
prevalences=(0.1 0.25 0.5)

#Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
######################### Get List of All Variants that Pass MAF>5% in both Pops ###############################
################################################################################################################
##Create empty command array
#cmds=()
##Loop over two populations
#for pop in AFR EUR; do for chromo in {1..22}; do
#		plink_inp=$geno_dir/$pop/chr"$chromo".plink.maf
#		plink_out=$geno_dir/$pop/chr"$chromo".plink.maf_0.05
#		cmds+=("plink2 --bfile $plink_inp --maf 0.05 --make-just-bim --out $plink_out")
#done; done
#
##Launch commands
#job_size=${#cmds[@]}
#echo "JOB SIZE: $job_size"
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#        fi
#done
#wait

################################################################################################################
################################ Find the Common Variants Shared by Both Populations ###########################
################################################################################################################
##Loop over two populations
#for chromo in {1..22}; do
#	comm -12 <(sort $geno_dir/AFR/chr"$chromo".plink.maf_0.05.bim) <(sort $geno_dir/EUR/chr"$chromo".plink.maf_0.05.bim) > $geno_dir/shared_variants/chr"$chromo".plink.maf_0.05.shared.bim
#done
#
##Launch commands
#job_size=${#cmds[@]}
#echo "JOB SIZE: $job_size"
#for (( i=0; i<${#cmds[@]}; i++ )); do
#	j=$(( $i + 1 ))
#	srun -n1 -N1 --ntasks-per-node=1 --job-name=$j ${cmds[$i]} &
#	if [ $(($j % $num_nodes)) == 0 ] || [ $j == $job_size ]; then 
#	    wait
#        fi
#done
#wait

################################################################################################################
############################################ Select Causal Variants ############################################
################################################################################################################
##Loop over the number of variants
#for num_var in ${set_numbers[@]}; do
#	python $selection_script --input-dir $geno_dir/shared_variants --output-dir $pheno_dir/causal_variants --effect-size $m_beta --max-distance $max_distance --min-spacing $min_spacing --variants-per-set $num_var
#done
#Make a master list of phenotypes
#ls $pheno_dir/causal_variants/Pheno*.*_causal.*.txt > $pheno_dir/causal_variants/master_casual_list.txt


################################################################################################################
############################################## Simulate Phenotypes #############################################
################################################################################################################
#Loop over the populations and cast the pcs to the burst buffer
for pop in AFR EUR; do
	sbcast -pf $geno_dir/$pop/merged.plink.maf.ld.pcs.eigenvec /mnt/bb/${USER}/$pop.eigenvec
done

#Loop over the phenotype traits and create commands
cmds=()
while read file; do 
	#Get number of causal variants and heritability for simulation
	temp=${file%_causal*}
	num_var=${temp#$pheno_dir/causal_variants/Pheno*.}
	temp2=${temp%.*}
	pheno=${temp2#$pheno_dir/causal_variants/}
	gen_var=$(echo $num_var*0.0002 | bc)
	#Get chromosome
	chromo=$(cut -f 1 -d : $file | head -n 1)
	#Loop over the trait types and populations to create commands
	for trait_type in ${trait_types[@]}; do for pop in AFR EUR; do
		if [ $trait_type == "binary" ]; then
			for proportion in ${prevalences[@]}; do
				out_file=$pheno_dir/$pop/$pheno."$num_var"_causal."$trait_type"_"$proportion".txt
				if [ ! -f $out_file ]; then 
					cmds+=("python $simulate_script --plink $geno_dir/$pop/chr"$chromo".plink.maf --effects $file --trait-type $trait_type --heritability $gen_var --prevalence $proportion --pc-file /mnt/bb/${USER}/$pop.eigenvec --output $out_file")
				fi
			done
		else
			out_file=$pheno_dir/$pop/$pheno."$num_var"_causal."$trait_type".txt
			if [ ! -f $out_file ]; then 
				cmds+=("python $simulate_script --plink $geno_dir/$pop/chr"$chromo".plink.maf --effects $file --trait-type $trait_type --heritability $gen_var --pc-file /mnt/bb/${USER}/$pop.eigenvec --output $out_file")
			fi
		fi 
	done; done
done < $pheno_dir/causal_variants/master_casual_list.txt

#Get total number of cpus
total_cpus=$(($num_nodes * 56))
#Get maximum number of parallel jobs
max_parallel=$(($total_cpus / 7))

#Launch commands
job_size=${#cmds[@]}
echo "JOB SIZE: $job_size"
for (( i=0; i<${#cmds[@]}; i++ )); do
	j=$(( $i + 1 ))
	srun -N1 -n1 -c7 --job-name=$j ${cmds[$i]} &
	if [ $(($j % $max_parallel)) == 0 ] || [ $j == $job_size ]; then 
	    wait
        fi
done
wait