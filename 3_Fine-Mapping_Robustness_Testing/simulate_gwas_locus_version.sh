#!/bin/bash
#SBATCH -A bif154
#SBATCH -J runbiogwas 
#SBATCH -o simulate_gwas_locus_version.log
#SBATCH -N 400
#SBATCH -p batch
#SBATCH -t 1:00:00

#Activate conda directory
source ~/.bashrc
conda activate biogwas

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set directories
raw_geno_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_genotypes"
AFR_geno_dir="$raw_geno_dir/AFR"
EUR_geno_dir="$raw_geno_dir/EUR"
gwas_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_gwas"
image_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/biogwas_images"
gencode_dir="/lustre/orion/bif154/proj-shared/mconery/resources/gencode"
gencode_chr_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/gencode_chr_gtfs"

#Set script locations
causal_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/3_Fine-Mapping_Robustness_Testing/Select_Causal_Variants.py"
biogwas_script="/lustre/orion/bif154/proj-shared/mconery/tools/bioGWAS/biogwas.py"
gwas_locus_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/3_Fine-Mapping_Robustness_Testing/simulate_gwas_single_locus.sh"

#Set general file locations
gencode_gtf=$gencode_dir/"gencode.v19.chr_patch_hapl_scaff.annotation.gtf"
dependency_file="/lustre/orion/bif154/proj-shared/mconery/tools/bioGWAS/dependencies.yaml"
AFR_accept=$gwas_dir/AFR.grch37.acceptable_locus_ranges.bed #This file was made manually
EUR_accept=$gwas_dir/EUR.grch37.acceptable_locus_ranges.bed #This file was made manually

#Set parameters for defining causal variants and their effect sizes
num_loci=1000 #Number of loci to map
max_size=500000 #Maximum distance between causal variants in a locus
min_distance=1000000 #Minimum distance between causal variants falling into separate loci
cores_per_node="$SLURM_CPUS_ON_NODE"
num_nodes="$SLURM_NNODES"

################################################################################################################
############## Make plink files filtered to MAF > 0.01 and for a random sample of 400K individuals #############
################################################################################################################
##Create random samples of individuals
#cut -f 1-2 $AFR_geno_dir/chr1.plink.tag.fam | sort -R | head -n 400000 > $AFR_geno_dir/AFR.400K_random_sample.tsv
#cut -f 1-2 $EUR_geno_dir/chr1.plink.tag.fam | sort -R | head -n 400000 > $EUR_geno_dir/EUR.400K_random_sample.tsv
##Filter for maf and random sample
#for i in {1..22}; do
#	plink2 --bfile $AFR_geno_dir/chr"$i".plink.tag --keep $AFR_geno_dir/AFR.400K_random_sample.tsv --maf 0.01 --set-all-var-ids "@:#:\$r:\$a" --make-bed --out $AFR_geno_dir/AFR.400K_random_sample.chr"$i".plink.tag.maf 
#	plink2 --bfile $EUR_geno_dir/chr"$i".plink.tag --keep $EUR_geno_dir/EUR.400K_random_sample.tsv --maf 0.01 --set-all-var-ids "@:#:\$r:\$a" --make-bed --out $EUR_geno_dir/EUR.400K_random_sample.chr"$i".plink.tag.maf 
#done
##Create links for the files because they can't have periods in them
#for i in {1..22}; do for suffix in bed bim fam; do ln $AFR_geno_dir/AFR.400K_random_sample.chr$i.plink.tag.maf.$suffix $AFR_geno_dir/AFR_chr"$i"_biogwas.$suffix; ln $EUR_geno_dir/EUR.400K_random_sample.chr$i.plink.tag.maf.$suffix $EUR_geno_dir/EUR_chr"$i"_biogwas.$suffix; done; done


################################################################################################################
############ Prune Variants by LD to Create Variants in Approximate LD and Restrict to MAF > 0.05 ##############
################################################################################################################
##Prune the variants
#for i in {1..22}; do 
#	plink2 --bfile $AFR_geno_dir/AFR.400K_random_sample.chr"$i".plink.tag.maf --maf 0.05 --indep-pairwise 1000 5 0.01 --out $AFR_geno_dir/AFR.400K_random_sample.chr"$i".maf_0.05.1000_step_5.r2_lt_0.01
#	plink2 --bfile $EUR_geno_dir/EUR.400K_random_sample.chr"$i".plink.tag.maf --maf 0.05 --indep-pairwise 1000 5 0.01 --out $EUR_geno_dir/EUR.400K_random_sample.chr"$i".maf_0.05.1000_step_5.r2_lt_0.01
#done
##Combine lists of pruned variants into one file
#cat $AFR_geno_dir/AFR.400K_random_sample.chr*.maf_0.05.1000_step_5.r2_lt_0.01.prune.in > $AFR_geno_dir/AFR.400K_random_sample.maf_0.05.1000_step_5.r2_lt_0.01.prune.in
#cat $EUR_geno_dir/EUR.400K_random_sample.chr*.maf_0.05.1000_step_5.r2_lt_0.01.prune.in > $EUR_geno_dir/EUR.400K_random_sample.maf_0.05.1000_step_5.r2_lt_0.01.prune.in


################################################################################################################
################################ Define Causal Variants for Simulated Traits ###################################
################################################################################################################
##Set input file paths for step
#AFR_prune_in=$AFR_geno_dir/AFR.400K_random_sample.maf_0.05.1000_step_5.r2_lt_0.01.prune.in
#EUR_prune_in=$EUR_geno_dir/EUR.400K_random_sample.maf_0.05.1000_step_5.r2_lt_0.01.prune.in
##Loop over the number of causal variants per loci and select loci
#for i in 1 2 3; do
#	AFR_out=$gwas_dir/AFR."$i"_causal_per_locus.variants.txt
#	EUR_out=$gwas_dir/EUR."$i"_causal_per_locus.variants.txt
#	python $causal_script --prune-in $AFR_prune_in --filter-file $AFR_accept --num-loci $num_loci --max-size $max_size --min_distance $min_distance --var_per_locus $i --output $AFR_out
#	python $causal_script --prune-in $EUR_prune_in --filter-file $EUR_accept --num-loci $num_loci --max-size $max_size --min_distance $min_distance --var_per_locus $i --output $EUR_out
#done

################################################################################################################
####################### Split the Causal Variants up into Separate Files for Each Locus ########################
################################################################################################################
##Loop over populations and # causal variants to split up the causal variants into separate files
#for pop in AFR EUR; do for num_var in {1..3}; do 
#	temp=$(wc -l $gwas_dir/$pop."$num_var"_causal_per_locus.loci.bed | cut -d " " -f 1); for ((i=1; i<=temp; i++)); do 
#	#Create prefix for the locus
#	read -r -a array <<< "$(awk -v a=$i 'NR==a' $gwas_dir/$pop."$num_var"_causal_per_locus.loci.bed)"
#	file_prefix=$pop"_"$num_var"_causal_per_locus_chr"${array[0]}"_"${array[1]}"_"${array[2]}
#	#Get row numbers to use for pulling the locus's causal variantss
#	top_line_num=$(( ($i - 1) * $num_var + 1 ))
#	bottom_line_num=$(( ($i - 1) * $num_var + 1 + num_var ))
#	#Make the output directory for the mapping results
#	out_dir=$gwas_dir/$file_prefix
#	mkdir -p $out_dir
#	#Output the snps to the separated file
#	out_file=$out_dir/$file_prefix.variants.txt
#	srun -n1 -c1 -N 1 --exclusive awk -v t=$top_line_num -v b=$bottom_line_num 'NR>=t && NR<b' $gwas_dir/$pop."$num_var"_causal_per_locus.variants.txt > $out_file &
#done; wait; done; done


################################################################################################################
###################################### Split up the GTF File by Chromsome ######################################
################################################################################################################
##Loop over chromosomes and split file
#for i in {1..22}; do
#	srun -n1 -N1 -c1 --exclusive --ntasks-per-node=1 awk -v chr="chr${i}" 'BEGIN{OFS="\t"} /^##/ {print} $1==chr {print}' $gencode_gtf > $gencode_chr_dir/gencode_comprehensive_chr${i}.gtf &
#done
#wait

################################################################################################################
############################ Make Locus-Specific Plink Files and Simulate the GWAS #############################
################################################################################################################
#Loop over causal variant, populations, and then the loci
for pop in AFR EUR; do for num_var in {1..3}; do locus_file=$gwas_dir/$pop."$num_var"_causal_per_locus.loci.bed; cat $locus_file | while read line; do
	
	#Set locus file prefix
	read -r -a array <<< "$line"
	locus_prefix=$pop"_"$num_var"_causal_per_locus_chr"${array[0]}"_"${array[1]}"_"${array[2]}
	
	#Create locus-specific files
	out_dir="$gwas_dir/$locus_prefix"
	gencode_gtf=$gencode_chr_dir/"gencode_comprehensive_chr"$(echo $locus_prefix | grep -oP 'chr\K[0-9]+')".gtf"
	
	#Run the locus code
	srun --exclusive -n1 --nodes=1 --ntasks-per-node=1 $gwas_locus_script $locus_prefix $raw_geno_dir $out_dir $image_dir $biogwas_script $dependency_file $gencode_gtf & 
done; wait; done; done

