#!/usr/bin/env python3

import os
import glob
import csv
import argparse
from concurrent.futures import ThreadPoolExecutor

# ------------------ Argument Parsing -----------------------
parser = argparse.ArgumentParser(
    description="Generate SAIGE Step 2 commands."
)
parser.add_argument("--output-path", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/step2_commands.txt", help="Output file for step2 commands.")
parser.add_argument("--plink-dir", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes", help="Directory for genotype files.")
parser.add_argument("--pheno-dir", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes", help="Directory for phenotype files.")
parser.add_argument("--step1-dir", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step1_results", help="Step 1 results directory.")
parser.add_argument("--out-dir", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results", help="Output directory for step 2 results.")
parser.add_argument("--saige-step2-script", default="/lustre/orion/bif154/proj-shared/arodriguez/tools/SAIGE-GPU/src/SAIGE/extdata/step2_SPAtests.R", help="SAIGE Step 2 R script path.")
parser.add_argument("--timing-file", default="timing2.csv", help="Timing file name.")
parser.add_argument("--populations", nargs="+", default=["AFR", "EUR"], help="List of populations.")
parser.add_argument("--number-traits-per-job", type=int, default=15, help="Number of traits per job.")
parser.add_argument("--causal-var-suffix", default=".500000_space100000.txt", help="Suffix for causal variant files.")
parser.add_argument("--minMAF", type=float, default=0, help="Minimum MAF.")
parser.add_argument("--minMAC", type=float, default=0.5, help="Minimum MAC.")
parser.add_argument("--LOCO", default="TRUE", help="Leave One Chromosome Out setting.")
parser.add_argument("--is-Firth-beta", default="TRUE", help="Firth beta setting.")
parser.add_argument("--allele-order", default="ref-first", help="Allele order.")
parser.add_argument("--is-fastTest", default="FALSE", help="Use fast test.")
parser.add_argument("--pCutoffforFirth", type=float, default=0.05, help="p-value cutoff for Firth correction.")
parser.add_argument("--num-cores", type=int, default=112, help="Number of cores.")

args = parser.parse_args()

# ------------------ Adjustable Parameters -----------------------
output_path = args.output_path
plink_dir = args.plink_dir
pheno_dir = args.pheno_dir
causal_var_dir = os.path.join(pheno_dir, "causal_variants")
inp_dir = args.step1_dir
out_dir = args.out_dir
saige_step2_script = args.saige_step2_script
timing_file = args.timing_file
populations = args.populations
number_traits_per_job = args.number_traits_per_job
causal_var_suffix = args.causal_var_suffix
minMAF = args.minMAF
minMAC = args.minMAC
LOCO = args.LOCO
is_Firth_beta = args.is_Firth_beta
allele_order = args.allele_order
is_fastTest = args.is_fastTest
pCutoffforFirth = args.pCutoffforFirth
num_cores = args.num_cores

# ------------------ Helper Functions -----------------------

def chunked_generator(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def get_chromosome_from_file(file_path):
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        if not first_line:
            return None
        variant_id = first_line.split()[0]
        chromosome = variant_id.split(':')[0]
        return [os.path.basename(file_path).replace(causal_var_suffix, ""), chromosome]

def create_commands_for_chromosome(chromo, chromo_file_prefixes, bgen_file, var_list, pop, inp_dir, out_dir, saige_step2_script, plink_dir, num_cores, number_traits_per_job):
    cmds = []
    os.makedirs(os.path.join(out_dir, f"{pop}_manifests"), exist_ok=True)
    chunk_counter = 1
    sample_file = f"{plink_dir}/{pop}/sample_list.txt"
    for sublist in chunked_generator(chromo_file_prefixes, number_traits_per_job):
        #Make chunk directory if it doesn't exist yet
        os.makedirs(f"{out_dir}/{pop}/chunk_{chunk_counter}",  exist_ok=True)
        manifest_file = os.path.join(out_dir, f"{pop}_manifests", f"chr{chromo}.chunk_{chunk_counter}.manifest.txt")
        manifest_list = [
            [
                f"{inp_dir}/{pop}/{prefix}.step1_output.rda",
                f"{inp_dir}/{pop}/{prefix}.step1_output.varianceRatio.txt",
                f"{out_dir}/{pop}/chunk_{chunk_counter}/{prefix}.chr{chromo}.step2.txt"
            ] for prefix in sublist
        ]
        with open(manifest_file, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(manifest_list)
        chunk_command = (
            f"start_time=$(date +%s.%N); Rscript {saige_step2_script} --bgenFile={bgen_file} --bgenFileIndex={bgen_file}.bgi "
            f"--manifestFile={manifest_file} --sampleFile={sample_file} --chrom={chromo} --AlleleOrder=ref-first "
            f"--minMAF={minMAF} --minMAC={minMAC} --LOCO={LOCO} --is_fastTest={is_fastTest} --is_Firth_beta={is_Firth_beta} "
            f"--idstoIncludeFile={var_list} "
            f"--AlleleOrder={allele_order} --pCutoffforFirth={pCutoffforFirth} --nThreads={int(num_cores/2)}; end_time=$(date +%s.%N); "
            f'duration=$(echo "$end_time - $start_time" | bc -l); echo "chunk_{chunk_counter},chr{chromo},$duration" >> {out_dir}/{pop}/timing2.csv; '
            f'mv {out_dir}/{pop}/chunk_{chunk_counter}/*.chr{chromo}.step2.txt {out_dir}/{pop}'
        )
        #Add timing command
        cmds.append(chunk_command)
        chunk_counter += 1
    return cmds

# ------------------ Output Directories -----------------------

for pop in populations:
    os.makedirs(os.path.join(out_dir, pop), exist_ok=True)
    os.makedirs(os.path.join(out_dir, pop + '_manifests'), exist_ok=True)

cmds = []

# ------------------ Main Pipeline -----------------------

for pop in populations:
    # Collect file and chromosome info
    rda_files = glob.glob(f"{inp_dir}/{pop}/*.step1_output.rda")
    ratio_files = glob.glob(f"{inp_dir}/{pop}/*.step1_output.varianceRatio.txt")
    causal_var_files = glob.glob(f"{causal_var_dir}/Pheno*.txt")
    rda_files = [os.path.basename(x).replace(".step1_output.rda", "") for x in rda_files]
    ratio_files = [os.path.basename(x).replace(".step1_output.varianceRatio.txt", "") for x in ratio_files]
    causal_var_files = {".".join(os.path.basename(x).split(".")[:3]):x for x in causal_var_files}
    shared_files = list(set(rda_files) & set(ratio_files))
    causal_dict = {
        x: os.path.join(causal_var_dir, ".".join(x.split(".")[:3]) + causal_var_suffix)
        for x in shared_files if ".".join(x.split(".")[:3]) in causal_var_files
    }
    with ThreadPoolExecutor(max_workers=num_cores) as executor:
        causal_chromosomes = list(executor.map(get_chromosome_from_file, causal_dict.values()))
        causal_chromosomes = [x for x in causal_chromosomes if x is not None]
    chromosomes = sorted(list(set([x[1] for x in causal_chromosomes])))
    chromo_dict = {os.path.basename(".".join(x[0].split(".")[:2])): x[1] for x in causal_chromosomes}
    null_files = [x for x in shared_files if "causal" not in x]
    causal_files = [x for x in shared_files if "causal" in x]
    if null_files != []:
        chromosomes = [str(x) for x in range(1,23,1)]
    chromo_to_file_dict = {
        chromo: [x for x in causal_files if chromo_dict.get(".".join(x.split(".")[:2])) == chromo] + null_files
        for chromo in chromosomes
    }
    chromo_to_file_dict = {
        chromo: [
            prefix for prefix in chromo_to_file_dict[chromo]
            if not os.path.exists(os.path.join(out_dir, pop, f"{prefix}.chr{chromo}.step2.txt"))
        ]
        for chromo in chromosomes
    }
    bgen_files_raw = glob.glob(f"{plink_dir}/{pop}/chr*.plink.maf.bgen")
    var_files_raw = glob.glob(f"{plink_dir}/{pop}/chr*.plink.maf.var_list.txt")
    chromosomes_bgen = [os.path.basename(x).replace(".plink.maf.bgen", "").replace("chr", "") for x in bgen_files_raw]
    chromosomes_var = [os.path.basename(x).replace(".plink.maf.var_list.txt", "").replace("chr", "") for x in var_files_raw]
    chr_intersection = sorted(list(set(chromosomes_bgen) & set(chromosomes_var) & set(chromosomes)))
    bgen_files = {chrom: f"{plink_dir}/{pop}/chr{chrom}.plink.maf.bgen" for chrom in chr_intersection}
    var_files = {chrom: f"{plink_dir}/{pop}/chr{chrom}.plink.maf.var_list.txt" for chrom in chr_intersection}

    # ------------ Parallel Command Generation ------------
    if len(chr_intersection) == 0:
        continue #Case where nothing needs mapping for the population
    with ThreadPoolExecutor(max_workers=min(num_cores, len(chr_intersection))) as executor:
        future_cmds = []
        for chromo in chr_intersection:
            chromo_file_prefixes = chromo_to_file_dict.get(chromo, [])
            bgen_file = bgen_files[chromo]
            var_list  = var_files[chromo]
            if chromo_file_prefixes:
                future_cmds.append(
                    executor.submit(
                        create_commands_for_chromosome,
                        chromo, chromo_file_prefixes, bgen_file, var_list, pop, inp_dir, out_dir,
                        saige_step2_script, plink_dir, num_cores, number_traits_per_job
                    )
                )
        for future in future_cmds:
            cmds.extend(future.result())

# ------------------ Write Commands to Output -----------------------
with open(output_path, 'w') as f:
    for cmd in cmds:
        f.write(cmd + '\n')
