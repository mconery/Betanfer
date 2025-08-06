#!/usr/bin/env python3

import os
import glob
import argparse

# --- Argument parsing ---
parser = argparse.ArgumentParser(
    description="Generate SAIGE Step 1 commands"
)
parser.add_argument(
    "--plink-dir", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes",
    help="Directory containing plink genotype files"
)
parser.add_argument(
    "--pheno-dir", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_phenotypes",
    help="Directory containing phenotype files"
)
parser.add_argument(
    "--out", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step1_results/step1_commands.txt",
    help="Path for output commands file"
)
parser.add_argument(
    "--out-dir", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step1_results",
    help="Directory for output results"
)
parser.add_argument(
    "--saige-path", default="/lustre/orion/bif154/proj-shared/arodriguez/tools/SAIGE_AMD_6.3.1/SAIGE-GPU/src/SAIGE",
    help="Full path to SAIGE"
)
parser.add_argument(
    "--saige-step1-script", default="/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/run_step1_saige_single_trait.sh",
    help="SAIGE step 1 runner script"
)
parser.add_argument(
    "--populations", nargs="+", default=["AFR", "EUR"], help="List of populations"
)
parser.add_argument(
    "--trait-types", nargs="+", default=["binary", "continuous"], help="List of trait types"
)
parser.add_argument(
    "--causal-flag", action="store_true", help="Include causal phenotypes"
)
args = parser.parse_args()

# Variables from args
plink_dir = args.plink_dir
pheno_dir = args.pheno_dir
out_dir = args.out_dir
path_to_saige = args.saige_path
saige_step1_script = args.saige_step1_script
output_path = args.out
populations = args.populations
trait_types = args.trait_types
causal_flag = args.causal_flag

cmds = []

for pop in populations:
    plink_file_post = f"/mnt/bb/mconery/{pop}"
    pop_out_dir = os.path.join(out_dir, pop)
    for trait_type in trait_types:
        pattern = os.path.join(pheno_dir, pop, f"Pheno*.{trait_type}*txt")
        phenotype_files = glob.glob(pattern)
        if not causal_flag:
            phenotype_files = [x for x in phenotype_files if 'causal' not in x]
        for file in phenotype_files:
            temp = os.path.basename(file)
            prefix = temp[:-4] if temp.endswith('.txt') else temp
            step1_output_rda = os.path.join(pop_out_dir, f"{prefix}.step1_output.rda")
            step1_output_variance = os.path.join(pop_out_dir, f"{prefix}.step1_output.varianceRatio.txt")
            if not (os.path.exists(step1_output_rda) and os.path.exists(step1_output_variance)):
                cmd = f"{saige_step1_script} {plink_file_post} {file} {prefix} {pop_out_dir} {path_to_saige} {trait_type}"
                cmds.append(cmd)

with open(output_path, 'w') as f:
    for cmd in cmds:
        f.write(cmd + '\n')
