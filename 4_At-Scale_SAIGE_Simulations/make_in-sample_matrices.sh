#!/bin/bash
#SBATCH -A bif154
#SBATCH -J make_matrix
#SBATCH -o make_in-sample_matrices.log
#SBATCH -N 10
#SBATCH -p batch
#SBATCH -t 1:00:00

#Activate conda directory
conda activate parsl_bash

################################################################################################################
########################## Define directories and other key file/script locations ##############################
################################################################################################################
#Set input directories and files
var_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants/master.locus_list.txt"
matrix_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
plink_dirs="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/AFR,/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes/EUR"
cmd_out_file="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices/make_ld_matrices.cmds.txt"

#Set script locations
command_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/make_matrix_command_writer.py"
parsl_script="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/parsl_run_1_command_per_node.py"

#set populations
populations=AFR,EUR

#Get number of nodes
num_nodes="$SLURM_NNODES"

################################################################################################################
############################################## Make Basic LD Matrices ##########################################
################################################################################################################

#Make command file
python $command_script --cmd_out $cmd_out_file --matrix_dir $matrix_dir --plink_dirs $plink_dirs --populations $populations --variant_lists_file $var_file
#Us parsl to run matrix creation
python $parsl_script $cmd_out_file
