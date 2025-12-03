#!/bin/bash

#SBATCH -A bif154
#SBATCH -J prep_geno
#SBATCH -o prep_causal_phenotypes_homo_effects.log
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
simulate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/5_Effect-Size_Homogeneity/simulate_phenotype_common_effect_sizes.py"

#Set key variables for phenotype simulation
trait_types=(binary continuous)
prevalence=0.25
heritability=0.00005

#Get number of nodes
num_nodes="$SLURM_NNODES"

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

  #Get number of causal variants for simulation
  temp=${file%_causal*}
  num_var=${temp#$pheno_dir/causal_variants/Pheno*.}
  temp2=${temp%.*}
  pheno=${temp2#$pheno_dir/causal_variants/}

  #Get chromosome
  chromo=$(cut -f 1 -d : $file | head -n 1)

  #Loop over the trait types to create commands
  for trait_type in ${trait_types[@]}; do

    # Scenario 1: AFR phenotypes using EUR as reference (effect sizes from EUR)
    if [ $trait_type == "binary" ]; then
      out_file=$pheno_dir/AFR/$pheno."$num_var"_causal.heritability_"$heritability"."$trait_type"_"$prevalence".matched_to_EUR.txt
      if [ ! -f $out_file ]; then
        cmds+=("python $simulate_script --plink-ref $geno_dir/EUR/chr"$chromo".plink.maf --plink-target $geno_dir/AFR/chr"$chromo".plink.maf --variants $file --trait-type $trait_type --per-snp-heritability $heritability --prevalence $prevalence --pc-file /mnt/bb/${USER}/AFR.eigenvec --output $out_file")
      fi
    else
      out_file=$pheno_dir/AFR/$pheno."$num_var"_causal.heritability_"$heritability"."$trait_type".matched_to_EUR.txt
      if [ ! -f $out_file ]; then
        cmds+=("python $simulate_script --plink-ref $geno_dir/EUR/chr"$chromo".plink.maf --plink-target $geno_dir/AFR/chr"$chromo".plink.maf --variants $file --trait-type $trait_type --per-snp-heritability $heritability --pc-file /mnt/bb/${USER}/AFR.eigenvec --output $out_file")
      fi
    fi

    # Scenario 2: EUR phenotypes using AFR as reference (effect sizes from AFR)
    if [ $trait_type == "binary" ]; then
      out_file=$pheno_dir/EUR/$pheno."$num_var"_causal.heritability_"$heritability"."$trait_type"_"$prevalence".matched_to_AFR.txt
      if [ ! -f $out_file ]; then
        cmds+=("python $simulate_script --plink-ref $geno_dir/AFR/chr"$chromo".plink.maf --plink-target $geno_dir/EUR/chr"$chromo".plink.maf --variants $file --trait-type $trait_type --per-snp-heritability $heritability --prevalence $prevalence --pc-file /mnt/bb/${USER}/EUR.eigenvec --output $out_file")
      fi
    else
      out_file=$pheno_dir/EUR/$pheno."$num_var"_causal.heritability_"$heritability"."$trait_type".matched_to_AFR.txt
      if [ ! -f $out_file ]; then
        cmds+=("python $simulate_script --plink-ref $geno_dir/AFR/chr"$chromo".plink.maf --plink-target $geno_dir/EUR/chr"$chromo".plink.maf --variants $file --trait-type $trait_type --per-snp-heritability $heritability --pc-file /mnt/bb/${USER}/EUR.eigenvec --output $out_file")
      fi
    fi

  done

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
