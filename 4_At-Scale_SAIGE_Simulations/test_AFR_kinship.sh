#!/bin/bash
#SBATCH -A bif154
#SBATCH -J test_AFR_kin
#SBATCH -o test_AFR_kin.log
#SBATCH -N 1
#SBATCH -p extended
#SBATCH -t 4:00:00
#SBATCH -C nvme

#Run the plink command
plink2 --bfile \
	/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/AFR/merged.plink.maf.ld \
	--make-rel \
	--out /lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/AFR/kinship_results