#!/bin/bash
#SBATCH -A bif154
#SBATCH -J make_matrix
#SBATCH -o make_raw_in-sample_matrices.log
#SBATCH -N 10
#SBATCH -p batch
#SBATCH -t 1:00:00

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set input directories and files
var_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants/master_causal_list.txt"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample"
plink_dirs="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/AFR,/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/EUR"
cmd_out_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/raw_in-sample/make_raw_in-sample_ld.cmds.txt"

#Set script locations
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/make_raw_in-sample_matrix_command_writer.py"
parsl_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/parsl_run_1_command_per_node.py"

#set populations
populations=AFR,EUR

#Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
############################################## Make Basic LD Matrices ##########################################
################################################################################################################

#Load modules and activate conda directory
module load flux
module load miniforge3
source activate /lustre/orion/bif154/proj-shared/mconery/tools/environments/simulate_pheno

#Make command file
python $command_script --cmd_out $cmd_out_file --matrix_dir $matrix_dir --plink_dirs $plink_dirs --populations $populations --variant_lists_file $var_file --variant-file-suffix ".variant_set"
#Us flux to run matrix creation
mapfile -t cmds < ${cmd_out_file}
srun -N ${num_nodes} -n ${num_nodes} -c 56 flux start \
    "flux resource list; mapfile -t cmds < ${cmd_out_file}; for ((i=0; i<${#cmds[@]}; i++)); do flux submit -N1 -n1 -c4 --output=flux_output/make_raw_matrix_\$i.txt \${cmds[\$i]}; done; flux queue drain;" 
echo 'SUCCESS: All raw in-sample matrices created.'
