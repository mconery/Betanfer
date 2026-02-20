#!/bin/bash
#SBATCH -A bif154
#SBATCH -J prep_geno
#SBATCH -o prep_causal_phenotypes_for_saige.log
#SBATCH -N 20
#SBATCH -p batch
#SBATCH -t 2:00:00
#SBATCH -C nvme

#Load plink tools and conda environment
module load miniforge3
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set directories
geno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"

#Set script locations
simulate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/simulate_phenotype.py"
selection_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/select_causal_variants_ld.py"

#Set key variables for causal variant selection
set_numbers=(1 2 3)
heritability=0.00005
max_distance=500000
max_ld=0.05
replicates=1000

#Set Files
selection_bed="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Betanfer/paper_analyses/1_Locus_Defining/grch38_chr_bounds_mhc_only.bed"

#Set key variables for phenotype simulation
trait_types=(binary continuous)
prevalence=0.25

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
#	python $selection_script --input-dir $geno_dir/shared_variants --output-dir $pheno_dir/causal_variants --population-dirs $geno_dir/AFR $geno_dir/EUR --plink-naming-convention "chr*.plink.maf" --max-distance $max_distance --max-ld $max_ld --variants-per-set $num_var --num-sets $replicates --bed-file $selection_bed
#done
##Make a master list of phenotypes
#ls $pheno_dir/causal_variants/Pheno*.*_causal.*.txt > $pheno_dir/causal_variants/master_causal_list.txt


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
	#Get number of causal variants and effect size for simulation
	temp=${file%_causal*}
	num_var=${temp#$pheno_dir/causal_variants/Pheno*.}
	temp2=${temp%.*}
	pheno=${temp2#$pheno_dir/causal_variants/}
	#Get chromosome
	chromo=$(cut -f 1 -d : $file | head -n 1)
	#Loop over the trait types and populations to create commands
	for trait_type in ${trait_types[@]}; do for pop in AFR EUR; do
		if [ $trait_type == "binary" ]; then
			out_file=$pheno_dir/$pop/$pheno."$num_var"_causal.heritability_"$heritability"."$trait_type"_"$prevalence".txt
			if [ ! -f $out_file ]; then 
				cmds+=("python $simulate_script --plink $geno_dir/$pop/chr"$chromo".plink.maf --variants $file --trait-type $trait_type --per-snp-heritability $heritability --prevalence $prevalence --pc-file /mnt/bb/${USER}/$pop.eigenvec --output $out_file")
			fi
		else
			out_file=$pheno_dir/$pop/$pheno."$num_var"_causal.heritability_"$heritability"."$trait_type".txt
			if [ ! -f $out_file ]; then 
				cmds+=("python $simulate_script --plink $geno_dir/$pop/chr"$chromo".plink.maf --variants $file --trait-type $trait_type --per-snp-heritability $heritability --pc-file /mnt/bb/${USER}/$pop.eigenvec --output $out_file")
			fi
		fi 
	done; done
done < $pheno_dir/causal_variants/master_causal_list.txt

#Reset environment
conda deactivate
module unload miniforge3
conda deactivate

#Write commands to file
>$pheno_dir/simulate_commands.txt
for (( i=0; i<${#cmds[@]}; i++ )); do echo ${cmds[$i]} >> $pheno_dir/simulate_commands.txt; done

#Activate flux modules
module load python/3.10-miniforge3 r/4.4.1 amd/6.3.1 rocm/6.3.1 hwloc/2.11.1-gpu flux
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno
#Simulate phenotypes
mapfile -t cmds < $pheno_dir/simulate_commands.txt
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start "mapfile -t cmds < $pheno_dir/simulate_commands.txt; flux resource list; for ((i=0; i<${#cmds[@]}; i++)); do flux submit -N1 -n1 -c7 --output=flux_output/simulate_\$i.txt \${cmds[\${i}]}; done; flux queue drain;"
#Unload conda environment
conda deactivate