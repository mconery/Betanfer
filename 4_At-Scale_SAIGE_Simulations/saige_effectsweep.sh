#!/bin/bash
#SBATCH -A bif154
#SBATCH -J saige_effectsweep
#SBATCH -o saige_effectsweep.%j.log
#SBATCH -N 30                 
#SBATCH -q debug
#SBATCH -t 1:10:00
#SBATCH -C nvme

#############################  GLOBAL SETTINGS  ################################
# (edit only paths or SBATCH directives above this block if you move locations)
###############################################################################
# Directories
BASE="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/phenotype_tuning"
geno_dir="$BASE/simulated_genotypes"
pheno_dir="$BASE/simulated_phenotypes"
step1_out="$BASE/step1_results"
step2_out="$BASE/step2_results"

# Scripts
simulate_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/simulate_phenotype.py"
selection_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/select_causal_variants.py"
make_step1_cmd="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/make_step1_commands.py"
make_step2_cmd="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/make_step2_commands.py"
parsl_exec="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/parsl_run_1_command_per_node.py"
causal_extractor="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/causal_variant_results_extractor.py"

# Arrays
effect_sizes=(0.5 1 5 10 100)
set_numbers=(1 2 3)            # number of causal variants
trait_types=(binary continuous)
prevalence=0.25                # single value for all binary traits
replicates=5         # five phenotypes per combination
populations=(AFR EUR)
max_distance=500000
min_spacing=100000
gen_var=0.4

# Resources
num_nodes="$SLURM_NNODES"      # for parallel sims
cpus_per_node=56               # OLCF Crusher
gpus=3			       # GPUs for step 1

################################################################################
######################## 1. SELECT CAUSAL VARIANTS & SIMULATE ##################
################################################################################
echo "=== Selecting variants ==="

module load miniforge3
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno

mkdir -p "$pheno_dir/causal_variants"
mkdir -p "$pheno_dir/AFR" "$pheno_dir/EUR"

#Loop over the number of variants
for num_var in 1; do for m_beta in ${effect_sizes[@]}; do
	python $selection_script --input-dir $geno_dir/shared_variants --output-dir $pheno_dir/causal_variants --effect-size $m_beta --max-distance $max_distance --min-spacing $min_spacing --variants-per-set $num_var --num-sets 5
done; done
#Make a master list of phenotypes
ls $pheno_dir/causal_variants/Pheno*.*_causal.*.txt > $pheno_dir/causal_variants/master_casual_list.txt

echo "=== Generating phenotypes ==="

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
	temp3=${file%_space*}
	temp4=${temp3#$pheno_dir/causal_variants/Pheno*.effect_size_}
	effect_size=${temp4%.*}
	#Get chromosome
	chromo=$(cut -f 1 -d : $file | head -n 1)
	#Loop over the trait types and populations to create commands
	for trait_type in ${trait_types[@]}; do for pop in AFR EUR; do
		if [ $trait_type == "binary" ]; then
			out_file=$pheno_dir/$pop/$pheno."$num_var"_causal.effect_size_"$effect_size"."$trait_type"_"$prevalence".txt
			if [ ! -f $out_file ]; then 
				cmds+=("python $simulate_script --plink $geno_dir/$pop/chr"$chromo".plink.maf --effects $file --trait-type $trait_type --heritability $gen_var --prevalence $prevalence --pc-file /mnt/bb/${USER}/$pop.eigenvec --output $out_file")
			fi
		else
			out_file=$pheno_dir/$pop/$pheno."$num_var"_causal.effect_size_"$effect_size"."$trait_type".txt
			if [ ! -f $out_file ]; then 
				cmds+=("python $simulate_script --plink $geno_dir/$pop/chr"$chromo".plink.maf --effects $file --trait-type $trait_type --heritability $gen_var --pc-file /mnt/bb/${USER}/$pop.eigenvec --output $out_file")
			fi
		fi 
	done; done
done < $pheno_dir/causal_variants/master_casual_list.txt

#Reset environment
conda deactivate
module unload miniforge3

#Write commands to file
>$pheno_dir/simulate_commands.txt
for (( i=0; i<${#cmds[@]}; i++ )); do echo ${cmds[$i]} >> $pheno_dir/simulate_commands.txt; done

#Activate flux modules
module load python/3.10-miniforge3 r/4.4.1 amd/6.3.1 rocm/6.3.1 hwloc/2.11.1-gpu flux
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno
#Simulate phenotypes
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start "mapfile -t cmds < $pheno_dir/simulate_commands.txt; flux resource list; for ((i=0; i<${#cmds[@]}; i++)); do flux submit -N1 -n1 -c7 --output=flux_output/simulate_\$i.txt \${cmds[\${i}]}; done; flux queue drain;"
#Unload conda environment
conda deactivate


################################################################################
######################## 2. SAIGE STEP 1 (GRM / NULL) ##########################
################################################################################
echo "=== SAIGE Step 1 ==="

# cast merged plink once per pop to NVMe
for pop in "${populations[@]}"; do
  pre="$geno_dir/$pop/merged.plink.maf.ld"
  post="/mnt/bb/${USER}/${pop}"
  sbcast -pf "${pre}.bed" "${post}.bed"
  sbcast -pf "${pre}.bim" "${post}.bim"
  sbcast -pf "${pre}.fam" "${post}.fam"
done

# regenerate commands so new traits are picked up
python "$make_step1_cmd" --plink-dir $geno_dir --pheno-dir "$pheno_dir" --out-dir $step1_out --out "$step1_out/step1_commands.txt" --causal-flag

#Load conda environment
##### START OF SBCAST AND CONDA-UNPACK #####
#Set variables
ENV_NAME=RSAIGE_1.3.3_amd6.3.1_gpu
TAR_FILE=${ENV_NAME}.tar.gz
TAR_DIR=/lustre/orion/bif154/proj-shared/arodriguez/tools/conda_envs
NNODES=${SLURM_NNODES}
# Move a copy of the env to the NVMe on each node
echo "copying ${ENV_NAME} to each node in the job"
sbcast -pf ${TAR_DIR}/${TAR_FILE} /mnt/bb/${USER}/${TAR_FILE}
if [ ! "$?" == "0" ]; then
    # CHECK EXIT CODE. When SBCAST fails, it may leave partial files on the compute nodes, and if you continue to launch srun,
    # your application may pick up partially complete shared library files, which would give you confusing errors.
    echo "SBCAST failed!"
    exit 1
fi

# Untar the environment file (only need 1 task per node to do this)
srun -N ${NNODES} --ntasks-per-node 1 mkdir -p /mnt/bb/${USER}/${ENV_NAME}
echo "untaring ${ENV_NAME}"
srun -N ${NNODES} --ntasks-per-node 1 tar -xzf /mnt/bb/${USER}/${TAR_FILE} -C /mnt/bb/${USER}/${ENV_NAME}

# Unpack the env
source activate /mnt/bb/${USER}/${ENV_NAME}
srun -N ${NNODES} --ntasks-per-node 1 conda-unpack
##### END OF SBCAST AND CONDA-UNPACK #####

#Run Step 1
mapfile -t cmds < $step1_out/step1_commands.txt
srun -N ${num_nodes} -n ${num_nodes} -c 56 --gpus-per-node=8 flux start "mapfile -t cmds < $step1_out/step1_commands.txt; flux resource list; for ((i=0; i<${#cmds[@]}; i++)); do flux submit -N1 -n $gpus --gpus-per-task=1 --output=flux_output/step1_\$i.txt \${cmds[\${i}]}; done; flux queue drain;"

#Unload conflicting modules for step 2
module unload amd/6.3.1 rocm/6.3.1 hwloc/2.11.1-gpu flux

################################################################################
######################## 3. SAIGE STEP 2 (GWAS scan) ###########################
################################################################################
echo "=== SAIGE Step 2 ==="

#Activate needed modules
module load python/3.10-miniforge3
module load r/4.4.1
module load amd/6.4.0

##### START OF SBCAST AND CONDA-UNPACK #####
#Set variables
ENV_NAME=RSAIGE_1.3.3_amd_gpu
TAR_FILE=${ENV_NAME}.tar.gz
TAR_DIR=/lustre/orion/bif154/proj-shared/arodriguez/tools/conda_envs
NNODES=${SLURM_NNODES}
# Move a copy of the env to the NVMe on each node
echo "copying ${ENV_NAME} to each node in the job"
sbcast -pf ${TAR_DIR}/${TAR_FILE} /mnt/bb/${USER}/${TAR_FILE}
if [ ! "$?" == "0" ]; then
    # CHECK EXIT CODE. When SBCAST fails, it may leave partial files on the compute nodes, and if you continue to launch srun,
    # your application may pick up partially complete shared library files, which would give you confusing errors.
    echo "SBCAST failed!"
    exit 1
fi

# Untar the environment file (only need 1 task per node to do this)
srun -N ${NNODES} --ntasks-per-node 1 mkdir -p /mnt/bb/${USER}/${ENV_NAME}
echo "untaring ${ENV_NAME}"
srun -N ${NNODES} --ntasks-per-node 1 tar -xzf /mnt/bb/${USER}/${TAR_FILE} -C /mnt/bb/${USER}/${ENV_NAME}

# Unpack the env
source activate /mnt/bb/${USER}/${ENV_NAME}
srun -N ${NNODES} --ntasks-per-node 1 conda-unpack
##### END OF SBCAST AND CONDA-UNPACK #####

# create updated step2 commands
python "$make_step2_cmd" --plink-dir $geno_dir --pheno-dir "$pheno_dir" --out-dir $step2_out --step1-dir "$step1_out" --output-path "$step2_out/step2_commands.txt" --number-traits-per-job 1

# massive parallel sweep – one trait per node
python "$parsl_exec" "$step2_out/step2_commands.txt"

#Deactivate conda environment
conda deactivate

################################################################################
######################### 4. Extract Causal Variants ###########################
################################################################################

#Activate simulation environment
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno
#Run extractor script
python $causal_extractor --output-dir $step2_out --causal-var-dir $pheno_dir/causal_variants --out $BASE/combined_causal_results.csv

echo "SUCCESS: Pipeline complete."
