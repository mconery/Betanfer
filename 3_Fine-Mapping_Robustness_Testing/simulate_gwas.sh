#!/bin/bash
#SBATCH -A bif154
#SBATCH -J runbiogwas 
#SBATCH -o simulate_gwas.log
#SBATCH -N 1
#SBATCH -p extended
#SBATCH -t 24:00:00

#Activate conda directory
source ~/.bashrc
conda activate biogwas

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set directories
AFR_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_genotypes/AFR"
EUR_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_genotypes/EUR"
gwas_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_gwas"
image_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/biogwas_images"

#Set script and general file locations
causal_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/3_Fine-Mapping_Robustness_Testing/Select_Causal_Variants.py"
gencode_gtf="/lustre/orion/bif154/proj-shared/mconery/resources/gencode/gencode.v19.chr_patch_hapl_scaff.annotation.gtf"
biogwas_script="/lustre/orion/bif154/proj-shared/mconery/tools/bioGWAS/biogwas.py"
dependency_file="/lustre/orion/bif154/proj-shared/mconery/tools/bioGWAS/dependencies.yaml"

#Set parameters for defining causal variants
num_loci=1000 #Number of loci to map
max_size=250000 #Maximum distance between causal variants in a locus
min_distance=1000000 #Minimum distance between causal variants falling into separate loci

################################################################################################################
################ Make plink files filtered by MAF and for a random sample of 400K individuals ##################
################################################################################################################
##Create random samples of individuals
#cut -f 1-2 $AFR_dir/chr1.plink.tag.fam | sort -R | head -n 400000 > $AFR_dir/AFR.400K_random_sample.tsv
#cut -f 1-2 $EUR_dir/chr1.plink.tag.fam | sort -R | head -n 400000 > $EUR_dir/EUR.400K_random_sample.tsv
##Filter for maf and random sample
#for i in {1..22}; do
#	plink2 --bfile $AFR_dir/chr"$i".plink.tag --keep $AFR_dir/AFR.400K_random_sample.tsv --maf 0.01 --set-all-var-ids "@:#:\$r:\$a" --make-bed --out $AFR_dir/AFR.400K_random_sample.chr"$i".plink.tag.maf 
#	plink2 --bfile $EUR_dir/chr"$i".plink.tag --keep $EUR_dir/EUR.400K_random_sample.tsv --maf 0.01 --set-all-var-ids "@:#:\$r:\$a" --make-bed --out $EUR_dir/EUR.400K_random_sample.chr"$i".plink.tag.maf 
#done
##Create links for the files because they can't have periods in them
#for i in {1..22}; do for suffix in bed bim fam; do ln $AFR_dir/AFR.400K_random_sample.chr$i.plink.tag.maf.$suffix $AFR_dir/AFR_chr"$i"_biogwas.$suffix; ln $EUR_dir/EUR.400K_random_sample.chr$i.plink.tag.maf.$suffix $EUR_dir/EUR_chr"$i"_biogwas.$suffix; done; done


################################################################################################################
######################### Prune Variants by LD to Create Variants in Approximate LD ############################
################################################################################################################
##Prune the variants
#for i in {1..22}; do 
#	plink2 --bfile $AFR_dir/AFR.400K_random_sample.chr"$i".plink.tag.maf --indep-pairwise 1000 5 0.01 --out $AFR_dir/AFR.400K_random_sample.chr"$i".1000_step_5.r2_lt_0.01
#	plink2 --bfile $EUR_dir/EUR.400K_random_sample.chr"$i".plink.tag.maf --indep-pairwise 1000 5 0.01 --out $EUR_dir/EUR.400K_random_sample.chr"$i".1000_step_5.r2_lt_0.01
#done
##Combine lists of pruned variants into one file
#cat $AFR_dir/AFR.400K_random_sample.chr*.1000_step_5.r2_lt_0.01.prune.in > $AFR_dir/AFR.400K_random_sample.1000_step_5.r2_lt_0.01.prune.in
#cat $EUR_dir/EUR.400K_random_sample.chr*.1000_step_5.r2_lt_0.01.prune.in > $EUR_dir/EUR.400K_random_sample.1000_step_5.r2_lt_0.01.prune.in


################################################################################################################
################################ Define Causal Variants for Simulated Traits ###################################
################################################################################################################
##Set input file paths for step
#AFR_prune_in=$AFR_dir/AFR.400K_random_sample.1000_step_5.r2_lt_0.01.prune.in
#EUR_prune_in=$EUR_dir/EUR.400K_random_sample.1000_step_5.r2_lt_0.01.prune.in
#AFR_accept=$gwas_dir/AFR.grch37.acceptable_locus_ranges.bed
#EUR_accept=$gwas_dir/EUR.grch37.acceptable_locus_ranges.bed
##Loop over the number of causal variants per loci and select loci
#for i in 1 2 3; do
#	AFR_out=$gwas_dir/AFR."$i"_causal_per_locus.variants.txt
#	EUR_out=$gwas_dir/EUR."$i"_causal_per_locus.variants.txt
#	python $causal_script --prune-in $AFR_prune_in --filter-file $AFR_accept --num-loci 1000 --max-size 250000 --min_distance 1000000 --var_per_locus $i --output $AFR_out
#	python $causal_script --prune-in $EUR_prune_in --filter-file $EUR_accept --num-loci 1000 --max-size 250000 --min_distance 1000000 --var_per_locus $i --output $EUR_out
#done


################################################################################################################
###################################### Simulate the GWAS using bioGWAS #########################################
################################################################################################################
##Create the genotype files
#for pop in AFR EUR; do
#	pop_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/rsparsepro_robustness_testing/simulated_genotypes/$pop"
#	ls $pop_dir/"$pop"_chr*_biogwas.bed | sed 's/.bed/,/g' | awk -F "_" '{print $0$5}' | sort -V > $gwas_dir/$pop.genotypes.list
#done
#Run biogwas
for pop in EUR AFR; do for i in {1..3}; do
	python $biogwas_script --data_dir $gwas_dir --img_dir $image_dir -bfile --input_list $gwas_dir/$pop.genotypes.list --anno_file $gencode_gtf --use_causal_snps --causal_snps $gwas_dir/$pop."$i"_causal_per_locus.variants.txt -maf 0.01 -maf_min 0.01 --skip_simulation -d $dependency_file --pattern "$pop"_"$i"_causal
done; done
