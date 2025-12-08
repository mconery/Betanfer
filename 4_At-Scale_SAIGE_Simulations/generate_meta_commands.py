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
    Find all GWAS result files matching the SAIGE step2 output patterns.
    
    File format categories:
    1. Basic:        Pheno#.continuous/binary.chr#.step2.txt
    2. Causal:       Pheno#.#_causal.heritability_#.continuous/binary.chr#.step2.txt
    3. Drop causal:  Pheno#.#_causal.heritability_#.continuous/binary.chr#.drop_causal.step2.txt
    4. Matched:      Pheno#.#_causal.heritability_#.continuous/binary.matched_to_POP.chr#.step2.txt
    
    Returns dictionary {(phenotype, trait_type, chromosome, file_type, modifier): filepath}
    where:
    - file_type: 'default', 'causal', 'drop_causal', or 'matched_to'
    - modifier: None for default/causal/drop_causal, target_pop for matched_to
    """
    files_dict = {}
    pattern = os.path.join(directory, "*.step2.txt")
    
    for filepath in glob.glob(pattern):
        filename = os.path.basename(filepath)
        
        # All valid files must end with .step2.txt
        if not filename.endswith('.step2.txt'):
            continue
        
        # Remove the .step2.txt extension
        base_name = filename[:-10]  # Remove '.step2.txt'
        
        # Determine file type and extract components
        file_type = None
        modifier = None
        chromosome = None
        trait_type = None
        phenotype = None
        
        # Pattern 3: Drop causal (must check before pattern 4 since 4 also has matched_to)
        if '.drop_causal' in base_name:
            # Format: Pheno#.#_causal.heritability_#.continuous/binary.chr#.drop_causal
            # Split on .drop_causal to separate before and after
            minus_drop = base_name.rstrip('.drop_causal')
            
            # Last part before drop should be chromosome number
            parts = minus_drop.split(".")
            if parts[-1].startswith("chr"):
                chromosome = parts[-1].lstrip("chr")  # Extract chromosome number
                file_type = 'drop_causal'
                phenotype = ".".join(parts[:-1])
                for i in range(len(parts)):
                    if parts[i] == "continuous":
                        trait_type = "continuous"
                        break
                    elif parts[i] == "binary":
                        trait_type = "binary"
                        break
                    elif parts[i] == "binary_0":
                        trait_type = ".".join(parts[i:i+2])           
                        break
        
        # Pattern 4: Matched to population
        elif '.matched_to_' in base_name:
            # Format: Pheno#.#_causal.heritability_#.continuous/binary.matched_to_POP.chr#
            # Split on .matched_to_ to extract population and chromosome
            parts = base_name.split('.matched_to_')
            if len(parts) == 2:
                before_matched = parts[0]
                after_matched = parts[1]
                
                # after_matched is "POP.chr#"
                matched_parts = after_matched.split('.')
                if len(matched_parts) >= 2:
                    modifier = matched_parts[0]  # Population name
                    
                    # Find chromosome part
                    for part in matched_parts[1:]:
                        if part.startswith('chr'):
                            chromosome = part[3:]
                            break
                    
                    if chromosome:
                        file_type = 'matched_to'
                        phenotype = before_matched
        
        # Pattern 2: Causal (with heritability)
        elif '_causal.heritability_' in base_name:
            # Format: Pheno#.#_causal.heritability_#.continuous/binary.chr#
            # Split on .chr to separate phenotype info from chromosome
            if '.chr' in base_name:
                before_chr = base_name.split('.chr')[0]
                chr_part = base_name.split('.chr')[1]
                
                # chr_part should be just the chromosome number
                if chr_part and chr_part[0].isdigit():
                    chromosome = chr_part
                    file_type = 'causal'
                    phenotype = before_chr
                    
                    # Extract trait type (last part before .chr)
                    parts = before_chr.split('.')
                    if parts:
                        potential_trait = parts[-1]
                        if potential_trait in ['continuous', 'binary']:
                            trait_type = potential_trait
                        elif "binary_0" == parts[-2]:
                            trait_type = ".".join(parts[-2:])
        
        # Pattern 1: Basic (no causal, no heritability)
        else:
            # Format: Pheno#.continuous/binary.chr#
            if '.chr' in base_name:
                before_chr = base_name.split('.chr')[0]
                chr_part = base_name.split('.chr')[1]
                
                # chr_part should be just the chromosome number
                if chr_part and chr_part[0].isdigit():
                    chromosome = chr_part
                    file_type = 'default'
                    phenotype = before_chr
                    
                    # Extract trait type (last part before .chr)
                    parts = before_chr.split('.')
                    if parts:
                        potential_trait = parts[-1]
                        if potential_trait in ['continuous', 'binary']:
                            trait_type = potential_trait
                        elif "binary_0" == parts[-2]:
                            trait_type = ".".join(parts[-2:])
        
        # If we successfully parsed the file, add it to the dictionary
        if file_type and phenotype and chromosome:
            # If trait_type wasn't extracted, try to find it in phenotype
            if not trait_type:
                parts = phenotype.split('.')
                for i in range(len(parts)):
                    if parts[i] == "continuous":
                        trait_type = "continuous"
                        break
                    elif parts[i] == "binary":
                        trait_type = "binary"
                        break
                    elif parts[i] == "binary_0":
                        trait_type = ".".join(parts[i:i+2])           
                        break
            
            if trait_type:
                key = (phenotype, trait_type, chromosome, file_type, modifier)
                files_dict[key] = filepath
    
    return files_dict


def find_common_traits_across_pops(pop_files_dict, pop_names_list):
    """
    Given a dict mapping pop_name -> {(phenotype, trait_type, chromosome, file_type, modifier): filepath},
    identify meta-analysis pairs:
    
    1. Standard pairs: (phenotype, trait_type, chromosome, 'default'/'causal') from all populations
    2. Matched pairs: matched_to_X file from population A paired with default/causal file from population X
    3. Drop causal pairs: drop_causal file from population A paired with default/causal file from population B
    
    Returns list of tuples:
    (phenotype, trait_type, chromosome, analysis_variant, {pop_name: filepath, ...})
    where analysis_variant describes the type of meta-analysis
    """
    meta_tasks = []
    
    # First, extract all unique (phenotype, trait_type, chromosome) combinations from default/causal files
    standard_files = {}  # {(pheno, trait, chr, file_type): {pop: filepath}} for 'default' or 'causal'
    matched_files = {}  # {(pheno, trait, chr, target_pop): {pop: filepath}}
    drop_causal_files = {}  # {(pheno, trait, chr): {pop: filepath}}
    
    for pop_name in pop_names_list:
        if pop_name not in pop_files_dict:
            continue
            
        pop_dict = pop_files_dict[pop_name]
        
        for (phenotype, trait_type, chromosome, file_type, modifier), filepath in pop_dict.items():
            if file_type in ['default', 'causal']:
                key = (phenotype, trait_type, chromosome, file_type)
                if key not in standard_files:
                    standard_files[key] = {}
                standard_files[key][pop_name] = filepath
                
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
    
    # 1. Standard meta-analyses: all populations have the same file type for a given trait
    for (phenotype, trait_type, chromosome, file_type), pop_files in standard_files.items():
        if len(pop_files) == len(pop_names_list):
            variant_name = 'standard' if file_type == 'default' else f'standard_{file_type}'
            meta_tasks.append((phenotype, trait_type, chromosome, variant_name, pop_files))
    
    # 2. Matched effects meta-analyses
    # For matched_to_X file from pop A, pair with standard (default/causal) file from pop X
    for (phenotype, trait_type, chromosome, target_pop), pop_files in matched_files.items():
        # pop_files has the population(s) that have matched_to_target_pop files
        # We need to find standard files from target_pop
        if target_pop in pop_files_dict:
            target_pop_dict = pop_files_dict[target_pop]
            
            # Look for either default or causal file with same phenotype/trait/chromosome in target pop
            standard_file = None
            for (tp, tt, tc, ft, mod), fpath in target_pop_dict.items():
                if tp == phenotype and tt == trait_type and tc == chromosome and ft in ['default', 'causal']:
                    standard_file = fpath
                    break
            
            if standard_file:
                # For each population with matched_to_target_pop file
                for source_pop, matched_file in pop_files.items():
                    if source_pop != target_pop:
                        paired_files = {
                            source_pop: matched_file,
                            target_pop: standard_file
                        }
                        variant_name = f"matched_to_{target_pop}"
                        meta_tasks.append(
                            (phenotype, trait_type, chromosome, variant_name, paired_files)
                        )
    
    # 3. Drop causal meta-analyses
    # For drop_causal file from pop A, pair with standard file from pop B
    for (phenotype, trait_type, chromosome), pop_files in drop_causal_files.items():
        # pop_files has population(s) with drop_causal files
        # Pair each with standard files from other populations
        for source_pop, drop_file in pop_files.items():
            # Find all populations that have standard (default/causal) files for this trait
            for (tp, tt, tc, ft), std_pop_files in standard_files.items():
                if tp == phenotype and tt == trait_type and tc == chromosome and ft in ['default', 'causal']:
                    for target_pop, standard_file in std_pop_files.items():
                        if source_pop != target_pop:
                            paired_files = {
                                source_pop: drop_file,
                                target_pop: standard_file
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
    
    pheno_out_dir = os.path.join(output_dir, f"{phenotype}")
    
    # Build output filename incorporating the analysis variant
    if analysis_variant == 'standard' or analysis_variant.startswith('standard_'):
        output_prefix = os.path.join(
            pheno_out_dir,
            f"{phenotype}.chr{chromosome}.{analysis_type}_meta"
        )
    else:
        # For matched or drop_causal variants, include them in the filename
        output_prefix = os.path.join(
            pheno_out_dir,
            f"{phenotype}.chr{chromosome}.{analysis_variant}.{analysis_type}_meta"
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
