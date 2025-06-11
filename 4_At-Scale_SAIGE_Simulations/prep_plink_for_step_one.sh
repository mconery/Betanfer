#!/bin/bash
#SBATCH -A bif154
#SBATCH -J prep_plink
#SBATCH -o prep_plink_for_step_one.log
#SBATCH -N 1
#SBATCH -p extended
#SBATCH -t 20:00:00

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

#Set simulation script location
simulate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/simulate_null_phenotypes.py"

################################################################################################################
################################################# Make Pruned Files ############################################
################################################################################################################
#Loop over two populations
for pop in AFR EUR; do 
	#Change into population output directory
	mkdir -p $geno_dir/$pop
	cd $geno_dir/$pop
	#Get keep file location
	keep_file=$inp_dir/$pop/$pop.400K_random_sample.tsv
	#Loop the chromosomes
	for chromo in {1..22}; do
		#Set plink file path
		plink_inp=$inp_dir/$pop/chr"$chromo".plink.tag
		plink_out=$geno_dir/$pop/chr"$chromo".plink.maf
		#Launch plink to prune variants
		plink2 --bfile $plink_inp --keep $keep_file --set-all-var-ids '@:#:$r:$a' --maf 0.05 --indep-pairwise 1000 80 0.1 --make-bed --out $plink_out
		plink2 --bfile $plink_out --extract $plink_out.prune.in --make-bed --out $plink_out.ld
	done
done

################################################################################################################
################################################# Make Merged Files ############################################
################################################################################################################
#Loop over the populations
for pop in AFR EUR; do
	cd $geno_dir/$pop
	#Create merge list and merge files
	ls $geno_dir/$pop/chr*.plink.maf.ld.bim | sed 's/.bim//g' > $pop.mergelist.txt
	plink2  --allow-no-sex --make-bed --pmerge-list $pop.mergelist.txt bfile --out $geno_dir/$pop/merged.plink.maf.ld
done

################################################################################################################
########################################### Calculate PCs and Traits ###########################################
################################################################################################################
#Loop over the populations
for pop in AFR EUR; do
	mkdir -p $pheno_dir/$pop
	cd $pheno_dir/$pop
	#Calculate PCs
	plink2  --bfile $geno_dir/$pop/merged.plink.maf.ld --pca approx 10 --out $geno_dir/$pop/merged.plink.maf.ld.pcs
	python $simulate_script -i $geno_dir/$pop/merged.plink.maf.ld.pcs.eigenvec -o $pheno_dir/$pop/merged.maf.ld.continuous.txt -n 1000 -t continuous
	python $simulate_script -i $geno_dir/$pop/merged.plink.maf.ld.pcs.eigenvec -o $pheno_dir/$pop/merged.maf.ld.continuous.txt -n 1000 -t binary
done

