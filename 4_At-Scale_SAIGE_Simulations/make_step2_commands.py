#!/usr/bin/env python3

#Set file of commands
output_path = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/step2_commands.txt"

# Set directory paths and variables
# Directories (adjust as needed)
plink_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
pheno_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes"
causal_var_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes/causal_variants"
inp_dir="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step1_results"
out_dir = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
saige_step2_script = "/lustre/orion/bif154/proj-shared/arodriguez/tools/SAIGE-GPU/src/SAIGE/extdata/step2_SPAtests.R"

populations = ["AFR", "EUR"]
trait_types = ["binary", "continuous"]
username = "mconery"
number_traits_per_job = 15
causal_var_suffix = ".effect_size_0.05.500000_space100000.txt"

#Set variables specific to SAIGE Step 2
minMAF=0
minMAC=0.5
LOCO="TRUE"
is_Firth_beta="TRUE"
allele_order="ref-first"
is_fastTest="TRUE"
pCutoffforFirth=0.05 

############################################################################################################################################################

#Library imports
import os
import glob
from concurrent.futures import ThreadPoolExecutor
import csv

num_cores = os.cpu_count()  # Automatically detect number of cores

############################################################################################################################################################
###################################################################### Useful Helper Functions #############################################################
############################################################################################################################################################

#Set function for extracting chromosomes for causal traits
def get_chromosome_from_file(file_path):
    """
    Reads the first line of a file and returns the chromosome part of the variant ID.
    Assumes the first column is 'chr:pos:ref:alt'.
    """
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        if not first_line:
            return None  # Handle empty files gracefully
        variant_id = first_line.split()[0]  # Adapt if your file is comma/tab separated
        chromosome = variant_id.split(':')[0]
        return [os.path.basename(file_path).replace(causal_var_suffix, ""), chromosome]

#Make a function that splits up a list into sections
def chunked_generator(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

############################################################################################################################################################
################################################## Make Critical Output Directories if They Don't Exist Yet ################################################
############################################################################################################################################################

for pop in populations:
    os.makedirs(os.path.join(out_dir, pop), exist_ok=True)
    os.makedirs(os.path.join(out_dir, pop + '_manifests'), exist_ok=True)

#Create commands list
cmds = []
#Loop over the populations and trait types
for pop in populations:
    #Cast the sample file to the burst buffer (Consider making this more robust in future)
    sample_file=f"{plink_dir}/{pop}/sample_list.txt"
    ###########################################################################################################################
    ######################################### Identify all the output files to create #########################################
    ###########################################################################################################################
    # Find all rda and varianceRatio files for the population 
    rda_files = glob.glob(f"{inp_dir}/{pop}/*.step1_output.rda")
    ratio_files = glob.glob(f"{inp_dir}/{pop}/*.step1_output.varianceRatio.txt")
    causal_var_files = glob.glob(f"{causal_var_dir}/Pheno*.txt")
    #Get base file names
    rda_files = [os.path.basename(x).replace(".step1_output.rda", "") for x in rda_files]
    ratio_files = [os.path.basename(x).replace(".step1_output.varianceRatio.txt", "") for x in ratio_files]
    causal_var_files = {".".join(os.path.basename(x).split(".")[:2]):x for x in causal_var_files}
    #Find intersection of sets
    shared_files = list(set(rda_files) & set(ratio_files))
    causal_dict = {x : os.path.join(causal_var_dir, ".".join(x.split(".")[:2]) + causal_var_suffix) for x in shared_files if ".".join(x.split(".")[:2]) in causal_var_files} #Make sure causal variants file exists too
    #Loop the files and identify the chromosomes for each output file
    with ThreadPoolExecutor(max_workers=num_cores) as executor:
        causal_chromosomes = list(executor.map(get_chromosome_from_file, causal_dict.values()))
    #Get list of all chromosomes and make a dictionary linking prefixes to their chromosomes
    chromosomes = list(set([x[1] for x in causal_chromosomes]))
    chromo_dict = {os.path.basename(".".join(os.path.basename(x[0]).split(".")[:2])) : x[1] for x in causal_chromosomes}
    #Extract the list of null shared files and those with known causal variants
    null_files = [x for x in shared_files if "causal" not in x]
    causal_files = [x for x in shared_files if "causal" in x]
    #Map the shared files to their chromosomes
    chromo_to_file_dict = {chromo : [x for x in causal_files if chromo_dict[".".join(x.split(".")[:2])] == chromo] + null_files for chromo in chromosomes}
    #Prune for output files that have already been made
    chromo_to_file_dict = {chromo : [prefix for prefix in chromo_to_file_dict[chromo] if not os.path.exists(os.path.join(out_dir,f"{prefix}.chr{chromo}.step2.txt"))] for chromo in chromosomes}
    ###########################################################################################################################
    ################################################## Prepare the File Commands ##############################################
    ###########################################################################################################################
    #Get bgen and variant list file locations
    bgen_files = glob.glob(f"{plink_dir}/{pop}/chr*.plink.maf.bgen")
    var_files = glob.glob(f"{plink_dir}/{pop}/chr*.plink.maf.var_list.txt")
    #Get chromosomes and set as bgen file names
    chromosomes_bgen = [os.path.basename(x).rstrip(".plink.maf.bgen").lstrip("chr") for x in bgen_files]
    chromosomes_var = [os.path.basename(x).rstrip(".plink.maf.var_list.txt").lstrip("chr") for x in var_files]
    chromosomes= list(set(chromosomes_bgen) & set(chromosomes_var))
    bgen_files = {chromosomes[x]:f"{plink_dir}/{pop}/chr{chromosomes[x]}.plink.maf.bgen" for x in range(len(chromosomes))}
    var_files = {chromosomes[x]:f"{plink_dir}/{pop}/chr{chromosomes[x]}.plink.maf.var_list.txt" for x in range(len(chromosomes))}
    #Loop over the chromosomes and combine the traits together into groups of the given size
    for chromo in chromosomes:
        #Get location of corresponding bgen file
        bgen_file = bgen_files[chromo]
        var_list = var_files[chromo]
        #Get the prefixes corresponding to the current chromosome and then chunk themp 
        chromo_file_prefixes = chromo_to_file_dict[chromo]
        chunk_counter = 1
        for sublist in chunked_generator(chromo_file_prefixes, number_traits_per_job):
            #Make manifest file
            manifest_file = os.path.join(out_dir, pop + '_manifests', 'chr' + chromo + '.chunk_' + str(chunk_counter) + '.manifest.txt')
            #Create list of files to go into manifest and then write it to the manifest
            manifest_list = [[f"{inp_dir}/{pop}/{prefix}.step1_output.rda", f"{inp_dir}/{pop}/{prefix}.step1_output.varianceRatio.txt", f"{out_dir}/{pop}/{prefix}.chr{chromo}.step2.txt"] for prefix in sublist]
            with open(manifest_file, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(manifest_list)
            #Create command
            chunk_command = f"Rscript {saige_step2_script} --bgenFile={bgen_file} --bgenFileIndex={bgen_file}.bgi \
                --manifestFile={manifest_file} --sampleFile={sample_file} --chrom={chromo} --AlleleOrder=ref-first \
                    --minMAF={minMAF} --minMAC={minMAC} --LOCO={LOCO} --is_fastTest={is_fastTest} --is_Firth_beta={is_Firth_beta} \
                        --idstoIncludeFile={var_list} \
                        --AlleleOrder={allele_order} --pCutoffforFirth={pCutoffforFirth} --nThreads={num_cores/2}"
            cmds.append(chunk_command)
            #Increment the counter
            chunk_count = chunk_counter + 1
                        
#Write to file 
with open(output_path, 'w') as f:
    for cmd in cmds:
        f.write(cmd + '\n')
