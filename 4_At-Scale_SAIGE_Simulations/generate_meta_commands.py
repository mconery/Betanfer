#!/usr/bin/env python3

"""

Flexible script to identify unique traits for meta-analysis across multiple populations
and generate GWAMA commands accordingly.
Author: Computational Genomics Pipeline
Date: Updated for multi-population support

"""

import os
import glob
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Generate GWAMA meta-analysis commands for multiple populations'
    )
    parser.add_argument(
        '--pop_names',
        required=True,
        help='Comma-separated list of population names (e.g. AFR,EUR)'
    )
    parser.add_argument(
        '--pop_dirs',
        required=True,
        help='Comma-separated list of directories containing population GWAS results; must correspond in order to pop_names'
    )
    parser.add_argument(
        '--output_dir',
        required=True,
        help='Directory for meta-analysis output files'
    )
    parser.add_argument(
        '--commands_file',
        default='meta_analysis_commands.txt',
        help='File to write GWAMA commands (default: meta_analysis_commands.txt)'
    )
    parser.add_argument(
        '--trait_types',
        nargs='+',
        default=['binary', 'continuous'],
        help='Trait types to process (default: binary continuous)'
    )
    parser.add_argument(
        '--analysis_type',
        choices=['random', 'fixed'],
        default='random',
        help='Type of meta-analysis to perform: random (default) or fixed effects'
    )
    return parser.parse_args()

def find_gwas_files(directory):
    """
    Find all GWAS result files matching the pattern:
    ...step2.txt
    Returns dictionary {(phenotype, trait_type, chromosome): filepath}
    """
    files_dict = {}
    pattern = os.path.join(directory, "*.step2.txt")
    for filepath in glob.glob(pattern):
        filename = os.path.basename(filepath)
        parts = filename.split('.')
        if len(parts) >= 4 and parts[-2] == 'step2' and parts[-1] == 'txt':
            chromosome = parts[-3]
            trait_type = parts[-4]
            phenotype = '.'.join(parts[:-4])
            files_dict[(phenotype, trait_type, chromosome)] = filepath
    return files_dict

def find_common_traits_across_pops(pop_files_dict):
    """
    Given a dict mapping pop_name -> {(phenotype, trait_type, chromosome): filepath},
    find trait keys common to all populations.
    Returns:
    dict {(phenotype, trait_type, chromosome): {pop_name: filepath, ...}}
    """
    # Find set of keys present in all populations
    all_keys_sets = [set(d.keys()) for d in pop_files_dict.values()]
    common_keys = set.intersection(*all_keys_sets)
    matched = {}
    for key in common_keys:
        matched[key] = {pop: pop_files_dict[pop][key] for pop in pop_files_dict.keys()}
    return matched

def check_file_exists(filepath):
    return os.path.exists(filepath) and os.path.getsize(filepath) > 0

def prepare_meta_task(args):
    # Unpack arguments
    phenotype, trait_type, chromosome, pop_files, output_dir, analysis_type = args
    pheno_out_dir = os.path.join(output_dir, f"{phenotype}.{trait_type}")
    output_prefix = os.path.join(pheno_out_dir, f"{phenotype}.{trait_type}.{chromosome}.{analysis_type}_meta")
    input_list_file = f"{output_prefix}.input_files.txt"
    expected_output = f"{output_prefix}.out"
    # If meta-analysis already run, skip
    if os.path.exists(expected_output) and os.path.getsize(expected_output) > 0:
        return None
    os.makedirs(pheno_out_dir, exist_ok=True)
    with open(input_list_file, 'w') as f:
        for fpath in pop_files.values():
            f.write(f"{fpath}\n")
    # Build GWAMA command based on analysis type
    if analysis_type == 'random':
        gwama_cmd = (
            f"GWAMA --filelist {input_list_file} --output {output_prefix} --random "
            "--quantitative --indel_alleles --name_marker MarkerID --name_ea Allele2 --name_nea Allele1 "
            "--name_n N --name_eaf AF_Allele2 --name_beta BETA --name_se SE"
        )
    else:  # fixed effects
        gwama_cmd = (
            f"GWAMA --filelist {input_list_file} --output {output_prefix} "
            "--quantitative --indel_alleles --name_marker MarkerID --name_ea Allele2 --name_nea Allele1 "
            "--name_n N --name_eaf AF_Allele2 --name_beta BETA --name_se SE"
        )
    return gwama_cmd

def generate_gwama_commands(matched_traits, output_dir, analysis_type, max_workers=32):
    """
    Prepare all directories and input files in parallel,
    return pure GWAMA command lines only.
    """
    # Prepare argument tuples for parallelization
    args_list = [
        (phenotype, trait_type, chromosome, pop_files, output_dir, analysis_type)
        for (phenotype, trait_type, chromosome), pop_files in matched_traits.items()
    ]
    commands = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for res in executor.map(prepare_meta_task, args_list):
            if res is not None:
                commands.append(res)
    return commands

def write_commands(commands, commands_file):
    with open(commands_file, 'w') as f:
        for cmd in commands:
            f.write(cmd + '\n')
    print(f"Generated {len(commands)} commands written to {commands_file}")

def main():
    args = parse_arguments()
    pop_names_list = [p.strip() for p in args.pop_names.split(',')]
    pop_dirs_list = [d.strip() for d in args.pop_dirs.split(',')]
    if len(pop_names_list) != len(pop_dirs_list):
        print("ERROR: The number of populations and directories must match!")
        exit(1)
    # Build population -> files_dict
    pop_files_dict = {}
    for pop_name, pop_dir in zip(pop_names_list, pop_dirs_list):
        pop_files_dict[pop_name] = find_gwas_files(pop_dir)
    matched_traits = find_common_traits_across_pops(pop_files_dict)
    print(f"Found {len(matched_traits)} matching trait-chromosome combinations across populations: {', '.join(pop_names_list)}")
    print(f"Analysis type: {args.analysis_type} effects meta-analysis")
    commands = generate_gwama_commands(matched_traits, args.output_dir, args.analysis_type, max_workers=os.cpu_count())
    write_commands(commands, args.commands_file)
    # Summary
    traits_summary = defaultdict(set)
    for (phenotype, trait_type, chromosome) in matched_traits.keys():
        traits_summary[(phenotype, trait_type)].add(chromosome)
    print("\nSummary:")
    print(f"Unique phenotype-trait combinations: {len(traits_summary)}")
    for (phenotype, trait_type), chromosomes in traits_summary.items():
        print(f" {phenotype} ({trait_type}): {len(chromosomes)} chromosomes")

if __name__ == "__main__":
    main()