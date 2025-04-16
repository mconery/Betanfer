#!/bin/bash -l
#PBS -A GeomicVar
#PBS -l walltime=3:00:00
#PBS -l filesystems=flare
#PBS -l select=32
#PBS -q prod
#PBS -M mconery@anl.gov
#PBS -N onekg_test 
#PBS -o /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/multinode_1000G_simulations.log
#PBS -e /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs/multinode_1000G_simulations.err

#Set up conda environment
source ~/.bashrc
conda activate finemap_old

#Set key directories and locations
onekg_script=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/2_Correlation_Inference/1000G_Simulation_Test.py
plink_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/plink_files
out_dir=/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/simulation_results
keep_dir=/flare/GeomicVar/mconery/resources/1000G

#Change directory into directory where script logs are stored
cd $(dirname $onekg_script)
cd ../../script_logs

#Create lists of commands
>multinode_1000G_simulations.commands.txt
for chromosome in 2 3 4; do
	for pop in EUR AFR ALL; do
		for trait_type in continuous binary; do
			for n_trait in 50 100 200 300 400 500 1000; do 
				for n_causal in 0 1 2 3 4 5; do 
					temp="python $onekg_script --plink $plink_dir/chr"$chromosome"_10k_sample.bed --pairs=10000 --traits=$n_trait --outdir=$out_dir --keep=$keep_dir/1000G.$pop.txt --trait_type=$trait_type --n_causal=$n_causal --min_distance 30000000"
					if [ "$n_causal" != 0 ]; then
						for effect_size in 0.1 0.5 1 2 5; do
							temp_e="$temp"" --effect_size=$effect_size"
							if [ "$trait_type" == "binary" ]; then
								for case_prob in 0.05 0.1 0.2 0.3 0.4 0.5; do
									temp_eb="$temp_e"" --case_prob=$case_prob"
									echo $temp_eb >> multinode_1000G_simulations.commands.txt
							 	done
							else
								echo $temp_e >> multinode_1000G_simulations.commands.txt
							fi
						done 
					else
						if [ "$trait_type" == "binary" ]; then
							for case_prob in 0.1 0.2 0.3 0.4 0.5; do
								temp_b="$temp"" --case_prob=$case_prob"
								echo $temp_b >> multinode_1000G_simulations.commands.txt
							done
						else
							echo $temp >> multinode_1000G_simulations.commands.txt
						fi
					fi
				done
			done
		done
	done
done

#Get the number of nodes
NUM_NODES=$(cat ${PBS_NODEFILE} | wc -l)

# Split commands into equal parts per node
split -d -n l/${NUM_NODES} --additional-suffix=".cmds" multinode_1000G_simulations.commands.txt onekg_cmds_part_

# Execute commands using MPI
mpirun -np ${NUM_NODES} bash -c '
    RANK=${OMPI_COMM_WORLD_RANK:-${PMI_RANK}}  # Handle different MPI implementations
    FILE_NAME=$(printf "onekg_cmds_part_%02d.cmds" ${RANK})
    while IFS= read -r cmd; do
        echo "Rank ${RANK} executing: ${cmd}"  # For debugging
        eval "${cmd}"
    done < "${FILE_NAME}"
'