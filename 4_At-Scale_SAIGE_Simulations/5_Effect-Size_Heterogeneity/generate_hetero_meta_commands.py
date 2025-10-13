#!/usr/bin/env python3

import os
import sys
import argparse
import re
import glob
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

def parse_filename(filename, default_case_prop=0.25):
    """
    Parse GWAS filename to extract components.
    
    Filename formats:
    - With causal variants: "Pheno#.#_causal.heritability_#.continuous/binary_#.chr#.step2.txt"
    - Without causal variants: "Pheno#.continuous/binary_#.chr#.step2.txt"
    
    Returns dict with keys: phenotype, n_causal, heritability, trait_type, case_prop, chromosome
    """
    basename = os.path.basename(filename)
    
    # Remove .step2.txt suffix
    if not basename.endswith('.step2.txt'):
        return None
    basename = basename[:-10]  # Remove '.step2.txt'
    
    parsed = {}
    
    # Match the phenotype at the beginning
    pheno_match = re.match(r'^(Pheno\d+)', basename)
    if not pheno_match:
        return None
    parsed['phenotype'] = pheno_match.group(1)
    
    # Look for causal variants pattern: #_causal or #.#_causal
    causal_match = re.search(r'(\d+)_causal', basename)
    if causal_match:
        parsed['n_causal'] = int(causal_match.group(1))
    else:
        parsed['n_causal'] = None
    
    # Look for heritability pattern: heritability_NUMBER (where NUMBER can have decimals)
    herit_match = re.search(r'heritability_(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)', basename)
    if herit_match:
        parsed['heritability'] = float(herit_match.group(1))
    else:
        parsed['heritability'] = 0
    
    # Look for trait type
    if 'continuous' in basename:
        parsed['trait_type'] = 'continuous'
        parsed['case_prop'] = None
    elif 'binary' in basename:
        parsed['trait_type'] = 'binary'
        # Look for case proportion after binary_
        case_match = re.search(r'binary_(\d+(?:\.\d+)?)', basename)
        if case_match:
            parsed['case_prop'] = float(case_match.group(1))
        else:
            parsed['case_prop'] = default_case_prop
    else:
        parsed['trait_type'] = None
        parsed['case_prop'] = None
    
    # Look for chromosome: chr#
    chr_match = re.search(r'chr(\d+)', basename)
    if chr_match:
        parsed['chromosome'] = int(chr_match.group(1))
    else:
        parsed['chromosome'] = None
    
    return parsed

def find_gwas_files(pop_dirs):
    """
    Find all GWAS files in the population directories and parse their filenames.
    
    Returns dict: {pop_name: [parsed_file_info, ...]}
    """
    files_by_pop = {}
    
    for pop_name, pop_dir in pop_dirs.items():
        files_by_pop[pop_name] = []
        
        # Find all .step2.txt files
        pattern = os.path.join(pop_dir, "*.step2.txt")
        gwas_files = glob.glob(pattern)
        
        for filepath in gwas_files:
            parsed = parse_filename(filepath)
            if parsed:
                parsed['filepath'] = filepath
                parsed['pop_name'] = pop_name
                files_by_pop[pop_name].append(parsed)
    
    return files_by_pop

def match_files_for_meta_analysis(files_by_pop, target_heritabilities):
    """
    Match files across populations for meta-analysis based on target heritabilities.
    
    Args:
        files_by_pop: Dict of {pop_name: [parsed_file_info, ...]}
        target_heritabilities: List of target heritability values for each population
    
    Returns list of matched file groups for meta-analysis
    """
    pop_names = list(files_by_pop.keys())
    
    if len(pop_names) != len(target_heritabilities):
        raise ValueError(f"Number of populations ({len(pop_names)}) must match number of heritabilities ({len(target_heritabilities)})")
    
    #Get unique numbers of causal variants (used to handle no causal variant case)
    unique_causal_counts = []
    for pop in files_by_pop.keys():
        temp = files_by_pop[pop]
        unique_causal_counts.extend(list(set([x['n_causal'] for x in temp if x['n_causal'] != None])))
    unique_causal_counts = list(set(unique_causal_counts))
    
    # Group files by phenotype, trait type, and chromosome for each population
    grouped_files = {}
    for pop_name in pop_names:
        grouped_files[pop_name] = defaultdict(list)
        
        for file_info in files_by_pop[pop_name]:
            if file_info['n_causal'] != None:
                # Create key for grouping: (phenotype, trait_type, case_prop, chromosome)
                key = (
                    file_info['phenotype'],
                    file_info['n_causal'],
                    file_info['trait_type'],
                    file_info.get('case_prop'),
                    file_info['chromosome']
                )
                grouped_files[pop_name][key].append(file_info)
            else:
                for causal_count in unique_causal_counts:
                    # Create key for grouping: (phenotype, trait_type, case_prop, chromosome)
                    key = (
                        file_info['phenotype'],
                        causal_count,
                        file_info['trait_type'],
                        file_info.get('case_prop'),
                        file_info['chromosome']
                    )
                    grouped_files[pop_name][key].append(file_info)
    
    # Find matches across populations
    matches = []
    
    # Get all possible keys from first population
    first_pop = [pop_names[x] for x in range(len(target_heritabilities)) if target_heritabilities[x] != 0 ][0]
    for key in grouped_files[first_pop]:
        # Check if this key exists in all populations
        if all(key in grouped_files[pop] for pop in pop_names):
            # Find matching files with target heritabilities
            matched_files = []
            
            for i, (pop_name, target_herit) in enumerate(zip(pop_names, target_heritabilities)):
                # Find file with matching heritability in this population
                matching_file = None
                
                for file_info in grouped_files[pop_name][key]:
                    file_herit = file_info.get('heritability')
                    if file_herit is None:
                        file_herit = 0.0  # Treat None as 0 for traits with no causal variants
                    
                    if abs(file_herit - target_herit) < 1e-10:  # Account for floating point precision
                        matching_file = file_info
                        break
                
                if matching_file:
                    matched_files.append(matching_file)
                else:
                    # No matching file found for this population/heritability combination
                    break
            
            # Only add to matches if we found files for all populations
            if len(matched_files) == len(pop_names):
                matches.append(matched_files)
        #Print keys missing in one or more populations
        else:
            print(f"Missing SAIGE Step2 results for chr{key[4]} of {key[2]} {key[0]}.{key[1]}_causal")
    
    return matches

def write_single_gwama_command(args):
    #Extract the important variables
    matched_files, output_dir, analysis_type = args
    # Create output filename based on first file's properties
    first_file = matched_files[0]
    
    # Build output name
    output_parts = [first_file['phenotype']]
    
    # Add causal info if present (use first file's info as reference)
    if first_file.get('n_causal') is not None:
        # Format the number to remove unnecessary decimals
        n_causal = first_file['n_causal']
        if n_causal == int(n_causal):
            output_parts.append(f"{int(n_causal)}_causal")
        else:
            output_parts.append(f"{n_causal}_causal")
    
    # Add heritabilities from all populations
    herit_parts = []
    for file_info in matched_files:
        herit = file_info.get('heritability')
        if herit is None:
            herit = 0.0
        # Format heritability to avoid scientific notation issues
        if herit == 0:
            herit_parts.append("0")
        elif herit == int(herit):
            herit_parts.append(str(int(herit)))
        else:
            herit_parts.append(str(herit))
    output_parts.append(f"heritability_{'_'.join(herit_parts)}")
    
    # Add trait type and case proportion
    if first_file.get('case_prop') is not None:
        case_prop = first_file['case_prop']
        if case_prop == int(case_prop):
            output_parts.append(f"{first_file['trait_type']}_{int(case_prop)}")
        else:
            output_parts.append(f"{first_file['trait_type']}_{case_prop}")
    else:
        output_parts.append(first_file['trait_type'])
    
    #Make phenotype output directory variable
    pheno_out_dir = os.path.join(output_dir, '.'.join(output_parts))
    
    # Add chromosome
    output_parts.append(f"chr{first_file['chromosome']}")
    
    output_filename = '.'.join(output_parts) + f'.{analysis_type}_meta'
    output_path = os.path.join(pheno_out_dir, output_filename)
    expected_output = f"{output_path}.out"
    input_list_file = f"{output_path}.input_files.txt"
    # If meta-analysis already run, skip
    if os.path.exists(expected_output) and os.path.getsize(expected_output) > 0:
        return None
    
    # Create input file list for GWAMA
    input_files = []
    for file_info in matched_files:
        input_files.append(file_info['filepath'])
    #Make phenotype output directory and gwama input file list if necessary
    os.makedirs(pheno_out_dir, exist_ok=True)
    with open(input_list_file, 'w') as f:
        for fpath in input_files:
            f.write(f"{fpath}\n")
    
    # Generate GWAMA command
    cmd_parts = [
        "GWAMA",
        f"--filelist {input_list_file}",
        f"--output {output_path}",
        "--quantitative", "--indel_alleles", "--name_marker MarkerID", "--name_ea Allele2", "--name_nea Allele1",
        "--name_n N", "--name_eaf AF_Allele2", "--name_beta BETA", "--name_se SE"
    ]
    if analysis_type == "random":
        cmd_parts.append("--random")
    #Return the command
    return ' '.join(cmd_parts)

def generate_gwama_commands(matches, output_dir, analysis_type, max_workers=os.cpu_count()):
    """
    Generate GWAMA meta-analysis commands for matched file groups.
    
    Returns list of command strings
    """
    commands = []
    # Prepare argument tuples for parallelization
    args_list = [(matched_files, output_dir, analysis_type) for matched_files in matches]
    
    #Execute command generation in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for res in executor.map(write_single_gwama_command, args_list):
            if res is not None:
                commands.append(res)
    
    return commands

def main():
    parser = argparse.ArgumentParser(description='Generate GWAMA meta-analysis commands for GWAS files with heritability filtering')
    
    parser.add_argument('--pop_names', required=True, help='Comma-separated population names')
    parser.add_argument('--pop_dirs', required=True, help='Comma-separated population directories')
    parser.add_argument('--heritabilities', required=True, help='Comma-separated per-variant heritabilities for each population')
    parser.add_argument('--output_dir', required=True, help='Output directory for meta-analysis results')
    parser.add_argument('--commands_file', required=True, help='File to write GWAMA commands')
    parser.add_argument('--analysis_type', default='random', choices=['fixed', 'random'], help='Meta-analysis type')
    
    args = parser.parse_args()
    
    # Parse inputs
    pop_names = [name.strip() for name in args.pop_names.split(',')]
    pop_dirs_list = [dir.strip() for dir in args.pop_dirs.split(',')]
    heritabilities = [float(h.strip()) for h in args.heritabilities.split(',')]
    
    if len(pop_names) != len(pop_dirs_list):
        print(f"Error: Number of population names ({len(pop_names)}) must match number of directories ({len(pop_dirs_list)})")
        sys.exit(1)
    
    if len(pop_names) != len(heritabilities):
        print(f"Error: Number of populations ({len(pop_names)}) must match number of heritabilities ({len(heritabilities)})")
        sys.exit(1)
    
    # Create pop_dirs dictionary
    pop_dirs = dict(zip(pop_names, pop_dirs_list))
    
    print(f"Populations: {pop_names}")
    print(f"Target heritabilities: {heritabilities}")
    
    # Find and parse GWAS files
    print("Finding and parsing GWAS files...")
    files_by_pop = find_gwas_files(pop_dirs)
    
    total_files = sum(len(files) for files in files_by_pop.values())
    print(f"Found {total_files} total GWAS files across {len(pop_names)} populations")
    
    for pop_name, files in files_by_pop.items():
        print(f"  {pop_name}: {len(files)} files")
    
    # Match files for meta-analysis
    print("Matching files for meta-analysis...")
    matches = match_files_for_meta_analysis(files_by_pop, heritabilities)
    
    print(f"Found {len(matches)} matching file groups for meta-analysis")
    
    if len(matches) == 0:
        print("Warning: No matching files found for the specified heritabilities")
        print("Available heritabilities by population:")
        for pop_name, files in files_by_pop.items():
            herit_values = set()
            for file_info in files:
                herit = file_info.get('heritability')
                if herit is None:
                    herit = 0.0
                herit_values.add(herit)
            print(f"  {pop_name}: {sorted(herit_values)}")
        sys.exit(1)
    
    # Generate GWAMA commands
    print("Generating GWAMA commands...")
    commands = generate_gwama_commands(matches, args.output_dir, args.analysis_type)
    
    # Write commands to file
    os.makedirs(args.output_dir, exist_ok=True)
    with open(args.commands_file, 'w') as f:
        for cmd in commands:
            f.write(cmd + '\n')
    
    print(f"Generated {len(commands)} commands and saved to {args.commands_file}")
    
    # Show sample matches for verification
    print("\nSample file matches:")
    for i, matched_files in enumerate(matches[:3]):  # Show first 3 matches
        print(f"  Match {i+1}:")
        for file_info in matched_files:
            herit = file_info.get('heritability')
            if herit is None:
                herit = 0.0
            print(f"    {file_info['pop_name']}: {os.path.basename(file_info['filepath'])} (h²={herit})")
        if i < len(matches) - 1:
            print()

if __name__ == '__main__':
    main()