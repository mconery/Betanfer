#!/usr/bin/env python3

"""
Flexible script to identify unique traits for meta-analysis across multiple populations
and generate GWAMA commands accordingly.

Author: Computational Genomics Pipeline
Date: Updated for multi-population support with matched effects, drop_causal, and dual analysis types
"""

import os
import glob
import argparse
import re
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
        choices=['random', 'fixed', 'both'],
        default='both',
        help='Type of meta-analysis to perform: random, fixed, or both (default: both)'
    )

    return parser.parse_args()


def find_gwas_files(directory):
    """
    Find all GWAS result files matching the pattern:
    ...step2.txt
    
    Categorizes files into three types:
    1. Default: standard files (no matched_to or drop_causal)
    2. Matched effects: files with 'matched_to_POP' in the name
    3. Drop causal: files with 'drop_causal' in the name
    
    Returns dictionary {(phenotype, trait_type, chromosome, file_type, modifier): filepath}
    where file_type is 'default', 'matched_to', or 'drop_causal'
    and modifier is None, target_pop, or source_pop respectively
    """
    files_dict = {}
    pattern = os.path.join(directory, "*.step2.txt")
    
    for filepath in glob.glob(pattern):
        filename = os.path.basename(filepath)
        parts = filename.split('.')
        
        if len(parts) >= 4 and parts[-2] == 'step2' and parts[-1] == 'txt':
            chromosome = parts[-3]
            trait_type = parts[-4]
            phenotype_base = '.'.join(parts[:-4])
            
            # Determine file type and extract modifier
            if 'matched_to_' in phenotype_base:
                # Extract the population being matched to
                match_obj = re.search(r'matched_to_([A-Z]+)$', phenotype_base)
                if match_obj:
                    target_pop = match_obj.group(1)
                    phenotype = phenotype_base.replace(f'.matched_to_{target_pop}', '')
                    key = (phenotype, trait_type, chromosome, 'matched_to', target_pop)
                    files_dict[key] = filepath
                    
            elif 'drop_causal' in phenotype_base:
                # This file has causal variants dropped
                phenotype = phenotype_base.replace('.drop_causal', '')
                key = (phenotype, trait_type, chromosome, 'drop_causal', None)
                files_dict[key] = filepath
                
            else:
                # Default file
                key = (phenotype_base, trait_type, chromosome, 'default', None)
                files_dict[key] = filepath
    
    return files_dict


def find_common_traits_across_pops(pop_files_dict, pop_names_list):
    """
    Given a dict mapping pop_name -> {(phenotype, trait_type, chromosome, file_type, modifier): filepath},
    identify meta-analysis pairs:
    
    1. Standard pairs: (phenotype, trait_type, chromosome, 'default') from all populations
    2. Matched pairs: matched_to_X file from population A paired with default file from population X
    3. Drop causal pairs: drop_causal file from population A paired with default file from population B
    
    Returns list of tuples:
    (phenotype, trait_type, chromosome, analysis_variant, {pop_name: filepath, ...})
    where analysis_variant describes the type of meta-analysis
    """
    meta_tasks = []
    
    # First, extract all unique (phenotype, trait_type, chromosome) combinations from default files
    default_files = {}  # {(pheno, trait, chr): {pop: filepath}}
    matched_files = {}  # {(pheno, trait, chr, target_pop): {pop: filepath}}
    drop_causal_files = {}  # {(pheno, trait, chr): {pop: filepath}}
    
    for pop_name in pop_names_list:
        if pop_name not in pop_files_dict:
            continue
            
        pop_dict = pop_files_dict[pop_name]
        
        for (phenotype, trait_type, chromosome, file_type, modifier), filepath in pop_dict.items():
            if file_type == 'default':
                key = (phenotype, trait_type, chromosome)
                if key not in default_files:
                    default_files[key] = {}
                default_files[key][pop_name] = filepath
                
            elif file_type == 'matched_to':
                key = (phenotype, trait_type, chromosome, modifier)  # modifier is target_pop
                if key not in matched_files:
                    matched_files[key] = {}
                matched_files[key][pop_name] = filepath
                
            elif file_type == 'drop_causal':
                key = (phenotype, trait_type, chromosome)
                if key not in drop_causal_files:
                    drop_causal_files[key] = {}
                drop_causal_files[key][pop_name] = filepath
    
    # 1. Standard meta-analyses: all populations have default files
    for (phenotype, trait_type, chromosome), pop_files in default_files.items():
        if len(pop_files) == len(pop_names_list):
            meta_tasks.append((phenotype, trait_type, chromosome, 'standard', pop_files))
    
    # 2. Matched effects meta-analyses
    # For matched_to_X file from pop A, pair with default file from pop X
    for (phenotype, trait_type, chromosome, target_pop), pop_files in matched_files.items():
        # pop_files has the population(s) that have matched_to_target_pop files
        # We need to find default files from target_pop
        if target_pop in pop_files_dict:
            target_pop_dict = pop_files_dict[target_pop]
            default_key = (phenotype, trait_type, chromosome, 'default', None)
            if default_key in target_pop_dict:
                target_default_file = target_pop_dict[default_key]
                
                # For each population with matched_to_target_pop file
                for source_pop, matched_file in pop_files.items():
                    if source_pop != target_pop:
                        paired_files = {
                            source_pop: matched_file,
                            target_pop: target_default_file
                        }
                        variant_name = f"matched_to_{target_pop}"
                        meta_tasks.append(
                            (phenotype, trait_type, chromosome, variant_name, paired_files)
                        )
    
    # 3. Drop causal meta-analyses
    # For drop_causal file from pop A, pair with default file from pop B
    for (phenotype, trait_type, chromosome), pop_files in drop_causal_files.items():
        # pop_files has population(s) with drop_causal files
        # Pair each with default files from other populations
        for source_pop, drop_file in pop_files.items():
            # Find all populations that have default files for this trait
            default_key = (phenotype, trait_type, chromosome)
            if default_key in default_files:
                for target_pop, default_file in default_files[default_key].items():
                    if source_pop != target_pop:
                        paired_files = {
                            source_pop: drop_file,
                            target_pop: default_file
                        }
                        variant_name = f"drop_causal_{source_pop}"
                        meta_tasks.append(
                            (phenotype, trait_type, chromosome, variant_name, paired_files)
                        )
    
    return meta_tasks


def check_file_exists(filepath):
    return os.path.exists(filepath) and os.path.getsize(filepath) > 0


def prepare_meta_task(args):
    """
    Prepare meta-analysis task: create output directory and input file list.
    Returns GWAMA command or None if output already exists.
    """
    phenotype, trait_type, chromosome, analysis_variant, pop_files, output_dir, analysis_type = args
    
    pheno_out_dir = os.path.join(output_dir, f"{phenotype}.{trait_type}")
    
    # Build output filename incorporating the analysis variant
    if analysis_variant == 'standard':
        output_prefix = os.path.join(
            pheno_out_dir,
            f"{phenotype}.{trait_type}.{chromosome}.{analysis_type}_meta"
        )
    else:
        # For matched or drop_causal variants, include them in the filename
        output_prefix = os.path.join(
            pheno_out_dir,
            f"{phenotype}.{trait_type}.{chromosome}.{analysis_variant}.{analysis_type}_meta"
        )
    
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


def generate_gwama_commands(meta_tasks, output_dir, analysis_types, max_workers=32):
    """
    Prepare all directories and input files in parallel,
    return pure GWAMA command lines only.
    
    analysis_types: list of 'random' and/or 'fixed'
    """
    # Expand meta_tasks to include all combinations of analysis types
    expanded_args_list = []
    
    for phenotype, trait_type, chromosome, analysis_variant, pop_files in meta_tasks:
        for analysis_type in analysis_types:
            args = (phenotype, trait_type, chromosome, analysis_variant, pop_files, output_dir, analysis_type)
            expanded_args_list.append(args)
    
    commands = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for res in executor.map(prepare_meta_task, expanded_args_list):
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
    
    # Find meta-analysis tasks
    meta_tasks = find_common_traits_across_pops(pop_files_dict, pop_names_list)
    
    print(f"Found {len(meta_tasks)} meta-analysis tasks across populations: {', '.join(pop_names_list)}")
    
    # Determine analysis types to run
    if args.analysis_type == 'both':
        analysis_types = ['random', 'fixed']
        analysis_desc = "random and fixed effects"
    else:
        analysis_types = [args.analysis_type]
        analysis_desc = f"{args.analysis_type} effects"
    
    print(f"Analysis type(s): {analysis_desc} meta-analysis")
    
    commands = generate_gwama_commands(meta_tasks, args.output_dir, analysis_types, max_workers=os.cpu_count())
    
    write_commands(commands, args.commands_file)
    
    # Summary
    tasks_summary = defaultdict(lambda: defaultdict(set))
    
    for (phenotype, trait_type, chromosome, analysis_variant, pop_files) in meta_tasks:
        tasks_summary[(phenotype, trait_type)][analysis_variant].add(chromosome)
    
    print("\nSummary:")
    print(f"Unique phenotype-trait combinations: {len(tasks_summary)}")
    
    for (phenotype, trait_type), variants in tasks_summary.items():
        print(f" {phenotype} ({trait_type}):")
        for variant, chromosomes in variants.items():
            print(f"   {variant}: {len(chromosomes)} chromosomes")


if __name__ == "__main__":
    main()
