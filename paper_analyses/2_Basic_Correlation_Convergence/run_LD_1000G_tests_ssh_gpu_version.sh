#!/bin/bash -l
#PBS -A GeomicVar
#PBS -l walltime=45:00
#PBS -l filesystems=flare
#PBS -l select=2:ncpus=208:mpiprocs=1
#PBS -l place=scatter:excl
#PBS -q prod
#PBS -M mconery@anl.gov
#PBS -N onekg_test 
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/multinode_1000G_simulations.LD.gpu.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/multinode_1000G_simulations.LD.gpu.err

#Set up conda environment
source ~/.bashrc
conda activate finemap

#Set key directories and locations
onekg_script=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/2_Basic_Correlation_Convergence/OnekG_CausalLD.py
allocator_script=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Betanfer/paper_analyses/2_Basic_Correlation_Convergence/OnekG_Simulation_looping_Allocator.py
plink_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/plink_files
out_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/causal_simulation_results
keep_dir=/flare/GeomicVar/mconery/resources/1000G

#Set fixed variables
n_trait=1000
trait_type=continuous
chromosome=2

#Change directory into directory where script logs are stored
cd $(dirname $onekg_script)
cd ../../script_logs

#Create lists of commands
>multinode_1000G_simulations.gpu.commands.txt
for pop in EUR AFR ALL; do for n_trait in 50 100 200 300 400 500 1000; do
	temp="python $onekg_script --plink $plink_dir/chr"$chromosome"_10k_sample.bed --pairs=10000 --traits=$n_trait --outdir=$out_dir --keep=$keep_dir/1000G.$pop.txt --trait_type=$trait_type --low_ld_threshold=0.01 --high_ld_threshold=0.1"
	for effect_size in 0.1 0.5 1 2 10 100 1000; do
		temp_e="$temp"" --effect_size=$effect_size"
		echo $temp_e >> multinode_1000G_simulations.gpu.commands.txt
	done 
done; done

# Get the number of nodes
num_nodes=$(cat ${PBS_NODEFILE} | wc -l)

# Loop over nodes and launch the job with each config file
i=0
# get the host ip address (use the same for config files for all nodes)
for node in $(cat ${PBS_NODEFILE}); do
                # launch jobs for current node
                ssh $node "source ~/.bashrc; conda activate finemap; python $allocator_script --command_file $(pwd)/multinode_1000G_simulations.gpu.commands.txt --total_nodes=$num_nodes --node_num=$i" & 
                i=$((i+1))
done

wait