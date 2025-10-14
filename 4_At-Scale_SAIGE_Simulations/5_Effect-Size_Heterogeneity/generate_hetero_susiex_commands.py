#!/usr/bin/env python3

"""
Script to generate SuSiEx multi-population fine-mapping commands for all loci.

This script creates a command file where each line runs SuSiEx fine-mapping
for a different locus across all available matrix types and populations.

Optimized for parallel processing on HPC systems with file existence checks.

UPDATED:
- GWAS files are located directly in population directories
- LD matrices are located in matrix_type/population/ structure (no weighted/Meta subdirs)
- Support for heterogeneous effect sizes across populations
- Support for null traits (no causal variants)
- Per-population heritability filtering
"""

import os
import glob
import argparse
import re
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial

def detect_matrix_subdirs(ld_base_dir, matrix_filter=None):
    """
    Automatically detect matrix subdirectories that contain population subdirectories.
    
    Parameters:
    - ld_base_dir: Base directory containing LD matrices
    - matrix_filter: Optional list of matrix subdirectory names to include
    
    Returns:
    - list: List of matrix type directory names
    """
    matrix_types = []
    
    # Look for all subdirectories in ld_base_dir
    if not os.path.exists(ld_base_dir):
        print(f"Warning: LD base directory does not exist: {ld_base_dir}")
        return matrix_types
    
    for item in os.listdir(ld_base_dir):
        matrix_dir = os.path.join(ld_base_dir, item)
        if not os.path.isdir(matrix_dir):
            continue
        
        # Apply matrix filter if specified
        if matrix_filter and item not in matrix_filter:
            print(f"Skipping matrix type (not in filter): {item}")
            continue
        
        # Check if this directory contains population subdirectories
        has_population_subdirs = False
        common_populations = ['AFR', 'EUR', 'EAS', 'SAS', 'AMR']
        
        for pop in common_populations:
            pop_dir = os.path.join(matrix_dir, pop)
            if os.path.isdir(pop_dir):
                has_population_subdirs = True
                break
        
        if has_population_subdirs:
            matrix_types.append(item)
            print(f"Found matrix type: {item}")
        else:
            print(f"Warning: No population subdirectories found in {matrix_dir}")
    
    return matrix_types


def find_matrix_files_by_population(ld_base_dir, matrix_type, pheno_pattern, populations):
    """
    Search for matrix files for all populations that match the given phenotype pattern.
    
    Parameters:
    - ld_base_dir: Base directory containing LD matrices
    - matrix_type: Matrix type directory name (e.g., "raw_in-sample")
    - pheno_pattern: Phenotype pattern (e.g., "Pheno1.2_causal" or "Pheno1" for null)
    - populations: List of population names (e.g., ["AFR", "EUR"])
    
    Returns:
    - dict: population -> full path to matching matrix file, or None if not found
    """
    matrix_files = {}
    
    for pop in populations:
        # Check population-specific directory directly under matrix_type
        search_dir = os.path.join(ld_base_dir, matrix_type, pop)
        
        if not os.path.exists(search_dir):
            matrix_files[pop] = None
            continue
        
        # Search for files that start with the pheno_pattern and end with .vcor1
        pattern = os.path.join(search_dir, f"{pheno_pattern}.*.vcor1")
        matching_files = glob.glob(pattern)
        
        if len(matching_files) == 1:
            matrix_files[pop] = matching_files[0]
        elif len(matching_files) > 1:
            # For parallel processing, we'll just use first match
            matrix_files[pop] = matching_files[0]
        else:
            # No matches found
            matrix_files[pop] = None
    
    return matrix_files


def parse_gwas_filename(filename):
    """
    Parse GWAS filename to extract components, handling both null and causal traits.
    
    Expected formats:
    - With causal variants: Pheno#.#_causal.heritability_#.continuous/binary_#.chr#.step2.txt
    - Without causal variants: Pheno#.continuous/binary_#.chr#.step2.txt
    
    Parameters:
    - filename: The GWAS filename
    
    Returns:
    - dict: Dictionary with parsed components or None if parsing fails
        Components: pheno_num, n_causal, heritability, trait_type, chromosome, pheno_pattern, has_causal
    """
    # Remove the .step2.txt extension
    base_name = filename.replace('.step2.txt', '')
    
    # Split by dots
    parts = base_name.split('.')
    
    if len(parts) < 3:  # Need at least: Pheno#, trait_type, chr#
        return None
    
    try:
        # Extract phenotype number (e.g., "Pheno1" -> 1)
        pheno_match = re.match(r'^Pheno(\d+)$', parts[0])
        if not pheno_match:
            return None
        pheno_num = int(pheno_match.group(1))
        
        # Check if this is a causal trait or null trait
        has_causal = False
        n_causal = 0
        heritability = 0.0
        
        if len(parts) > 1 and '_causal' in parts[1]:
            # This is a causal trait: Pheno#.#_causal.heritability_#.trait_type.chr#
            has_causal = True
            
            # Extract number of causal variants
            causal_match = re.match(r'^(\d+)_causal$', parts[1])
            if not causal_match:
                return None
            n_causal = int(causal_match.group(1))
            
            # Find heritability - it starts with "heritability_"
            if len(parts) < 5:  # Need: Pheno#, #_causal, heritability_#, trait_type, chr#
                return None
            
            # The heritability can span multiple parts if it has a decimal point
            # e.g., ["heritability_0", "00005"] or ["heritability_0", "0003"]
            herit_idx = 2
            if not parts[herit_idx].startswith('heritability_'):
                return None
            
            # Reconstruct heritability value, handling decimal points
            herit_parts = [parts[herit_idx].replace('heritability_', '')]
            herit_idx += 1
            
            # Continue collecting parts until we hit a trait type or 'chr'
            while herit_idx < len(parts):
                part = parts[herit_idx]
                if part in ['continuous', 'binary_0'] or part.startswith('chr'):
                    break
                herit_parts.append(part)
                herit_idx += 1
            
            # Reconstruct heritability as float
            herit_str = '.'.join(herit_parts)
            heritability = float(herit_str)
            
            # Now find trait type starting from herit_idx
            trait_type_idx = herit_idx
        else:
            # This is a null trait: Pheno#.trait_type.chr#
            has_causal = False
            n_causal = 0
            heritability = 0.0
            trait_type_idx = 1
        
        # Parse trait type (continuous or binary_#)
        if trait_type_idx >= len(parts):
            return None
        
        if parts[trait_type_idx] == 'continuous':
            trait_type = 'continuous'
            chromosome_idx = trait_type_idx + 1
        elif parts[trait_type_idx].startswith('binary_'):
            # Handle binary_# where # can have decimal (e.g., binary_0.25)
            # This can be split as ["binary_0", "25"] or just ["binary_0"]
            trait_parts = [parts[trait_type_idx]]
            next_idx = trait_type_idx + 1
            
            # Check if next part is part of the binary proportion
            if next_idx < len(parts) and not parts[next_idx].startswith('chr'):
                # This is the decimal continuation
                trait_parts.append(parts[next_idx])
                chromosome_idx = next_idx + 1
            else:
                chromosome_idx = next_idx
            
            trait_type = '.'.join(trait_parts)
        else:
            return None
        
        # Get chromosome
        if chromosome_idx >= len(parts):
            return None
        
        chromosome = parts[chromosome_idx]
        
        # Construct the pattern identifier for matching
        # For causal traits: Pheno#.#_causal
        # For null traits: Pheno#
        if has_causal:
            pheno_pattern = f"Pheno{pheno_num}.{n_causal}_causal"
        else:
            pheno_pattern = f"Pheno{pheno_num}"
        
        return {
            'pheno_num': pheno_num,
            'n_causal': n_causal,
            'heritability': heritability,
            'trait_type': trait_type,
            'chromosome': chromosome,
            'pheno_pattern': pheno_pattern,
            'has_causal': has_causal
        }
    
    except (IndexError, ValueError) as e:
        return None


def check_output_completeness(output_path, min_file_size=100):
    """
    Check if output file exists and appears complete.
    
    Parameters:
    - output_path: Path to output file
    - min_file_size: Minimum file size in bytes to consider complete
    
    Returns:
    - bool: True if file exists and appears complete
    """
    if not os.path.exists(output_path):
        return False
    
    try:
        # Check file size
        file_size = os.path.getsize(output_path)
        if file_size < min_file_size:
            return False
        
        # For SuSiEx output files, check if they have a header and at least one data row
        with open(output_path, 'r') as f:
            lines = f.readlines()
        
        # Should have at least 2 lines (header + 1 data row)
        if len(lines) < 2:
            return False
        
        # Check if header contains expected columns for SuSiEx output
        header = lines[0].strip().split('\t')
        expected_cols = {'VARIANT_ID', 'PIP', 'converged', 'cs'}
        if not expected_cols.issubset(set(header)):
            return False
        
        return True
    
    except (IOError, OSError, UnicodeDecodeError):
        # If we can't read the file, consider it incomplete
        return False


def get_population_sample_sizes(populations, default_n=40000):
    """
    Get sample sizes for each population. For now, using default values.
    In the future, this could read from a config file or be passed as arguments.
    
    Parameters:
    - populations: List of population names
    - default_n: Default sample size to use
    
    Returns:
    - dict: population -> sample size
    """
    # Default sample sizes - could be made configurable
    default_sizes = {
        'AFR': 40000,
        'EUR': 40000,
        'EAS': 40000,
        'SAS': 40000,
        'AMR': 40000
    }
    
    sample_sizes = {}
    for pop in populations:
        sample_sizes[pop] = default_sizes.get(pop, default_n)
    
    return sample_sizes


def find_matching_gwas_files(gwas_base_dir, populations, pheno_num, n_causal, heritabilities, trait_type, chromosome):
    """
    Find matching GWAS files for each population with the specified heritability values.
    
    Parameters:
    - gwas_base_dir: Base directory containing GWAS files
    - populations: List of population names
    - pheno_num: Phenotype number
    - n_causal: Number of causal variants (0 for null traits)
    - heritabilities: List of heritability values (one per population)
    - trait_type: Trait type (continuous or binary_#)
    - chromosome: Chromosome identifier
    
    Returns:
    - dict: population -> gwas file path, or None if not found
    """
    gwas_files = {}
    
    for i, pop in enumerate(populations):
        pop_gwas_dir = os.path.join(gwas_base_dir, pop)
        
        if not os.path.exists(pop_gwas_dir):
            gwas_files[pop] = None
            continue
        
        target_heritability = heritabilities[i]
        
        # Construct expected filename based on heritability and n_causal
        if target_heritability == 0.0 or n_causal == 0:
            # Null trait - no causal variants
            expected_filename = f"Pheno{pheno_num}.{trait_type}.{chromosome}.step2.txt"
        else:
            # Causal trait
            expected_filename = f"Pheno{pheno_num}.{n_causal}_causal.heritability_{target_heritability}.{trait_type}.{chromosome}.step2.txt"
        
        expected_path = os.path.join(pop_gwas_dir, expected_filename)
        
        if os.path.exists(expected_path):
            gwas_files[pop] = expected_path
        else:
            # Try to find matching file by parsing all files
            all_files = glob.glob(os.path.join(pop_gwas_dir, f"Pheno{pheno_num}*.{chromosome}.step2.txt"))
            
            found = False
            for candidate_file in all_files:
                parsed = parse_gwas_filename(os.path.basename(candidate_file))
                
                if parsed and \
                   parsed['pheno_num'] == pheno_num and \
                   parsed['n_causal'] == n_causal and \
                   abs(parsed['heritability'] - target_heritability) < 1e-10 and \
                   parsed['trait_type'] == trait_type and \
                   parsed['chromosome'] == chromosome:
                    gwas_files[pop] = candidate_file
                    found = True
                    break
            
            if not found:
                gwas_files[pop] = None
    
    return gwas_files


def process_trait_combination(args_tuple):
    """
    Process a specific trait combination across populations with specified heritabilities.
    
    This function is designed to be called in parallel.
    
    Parameters:
    - args_tuple: Tuple containing all necessary arguments
    
    Returns:
    - tuple: (commands, warnings, skipped_info)
    """
    (pheno_num, n_causal, trait_type, chromosome, heritabilities, 
     gwas_base_dir, ld_base_dir, geno_base_dir, matrix_types, 
     output_base_dir, wrapper_script, populations, plink_conv, skip_existing) = args_tuple
    
    commands = []
    warnings = []
    skipped_info = []
    
    # Construct pheno_pattern for LD matrix matching
    if n_causal > 0:
        pheno_pattern = f"Pheno{pheno_num}.{n_causal}_causal"
    else:
        pheno_pattern = f"Pheno{pheno_num}"
    
    # Find matching GWAS files for each population
    gwas_files = find_matching_gwas_files(
        gwas_base_dir, populations, pheno_num, n_causal, 
        heritabilities, trait_type, chromosome
    )
    
    # Check if we have GWAS files for all populations
    missing_pops = [pop for pop, path in gwas_files.items() if path is None]
    if missing_pops:
        herit_str = ','.join([str(h) for h in heritabilities])
        warnings.append(
            f"Missing GWAS files for populations {missing_pops} "
            f"(Pheno{pheno_num}, {n_causal} causal, h²=[{herit_str}], {trait_type}, {chromosome})"
        )
        return commands, warnings, skipped_info
    
    # Get plink paths
    plink_paths = [os.path.join(geno_base_dir, pop, plink_conv.replace("CHR#", chromosome.replace('chr', ''))) 
                   for pop in populations]
    
    # Verify plink files exist
    for extension in [".bed", ".bim", ".fam"]:
        for plink_path in plink_paths:
            if not os.path.exists(str(plink_path) + extension):
                warnings.append(f"Missing Plink Ref File: {str(plink_path)}{extension}")
                return commands, warnings, skipped_info
    
    # Get sample sizes for populations
    sample_sizes = get_population_sample_sizes(populations)
    
    # For each detected matrix type
    for matrix_type in matrix_types:
        # Find matrix files for all populations
        matrix_files = find_matrix_files_by_population(ld_base_dir, matrix_type, pheno_pattern, populations)
        
        # Check if we have matrix files for all populations
        missing_matrices = [pop for pop, path in matrix_files.items() if path is None]
        if missing_matrices:
            warnings.append(
                f"Missing matrix files for populations {missing_matrices} "
                f"in {pheno_pattern}, {matrix_type}"
            )
            continue
        
        # Construct output directory and file
        output_dir = os.path.join(output_base_dir, matrix_type)
        
        # Output filename format: Pheno#.#_causal.h2_X_Y.trait_type.matrix_type.txt
        # where X_Y represents the heritability values for each population
        herit_str = '_'.join([str(h).replace('.', 'p') for h in heritabilities])
        if n_causal > 0:
            output_filename = f"Pheno{pheno_num}.{n_causal}_causal.h2_{herit_str}.{trait_type}.{matrix_type}.txt"
        else:
            output_filename = f"Pheno{pheno_num}.h2_{herit_str}.{trait_type}.{matrix_type}.txt"
        
        output_path = os.path.join(output_dir, output_filename)
        
        # Check if output file already exists and is complete
        if skip_existing and check_output_completeness(output_path):
            skipped_info.append({
                'pheno_num': pheno_num,
                'n_causal': n_causal,
                'heritabilities': heritabilities,
                'matrix_type': matrix_type,
                'trait_type': trait_type,
                'chromosome': chromosome,
                'output_path': output_path,
                'populations': populations
            })
            continue
        
        # Create comma-delimited lists for SuSiEx
        gwas_files_list = ",".join([gwas_files[pop] for pop in populations])
        ld_files_list = ",".join([matrix_files[pop] for pop in populations])
        n_gwas_list = ",".join([str(sample_sizes[pop]) for pop in populations])
        
        # Create the SuSiEx command
        susiex_cmd = (f"python {wrapper_script} "
                     f"--gwas_files {gwas_files_list} "
                     f"--ld_files {ld_files_list} "
                     f"--n_gwas {n_gwas_list} "
                     f"--ref_panels {','.join(plink_paths)} "
                     f"--output_file {output_path}")
        
        # Store command with metadata
        commands.append((susiex_cmd, output_dir, matrix_type, pheno_pattern, trait_type,
                        chromosome, heritabilities, populations, output_filename))
    
    return commands, warnings, skipped_info


def parse_matrix_filter(matrix_string):
    """
    Parse comma-delimited matrix filter string.
    
    Parameters:
    - matrix_string: Comma-delimited string of matrix names
    
    Returns:
    - list: List of matrix names, or None if empty
    """
    if not matrix_string:
        return None
    
    # Split by comma and strip whitespace
    matrix_list = [name.strip() for name in matrix_string.split(',')]
    
    # Remove empty strings
    matrix_list = [name for name in matrix_list if name]
    
    return matrix_list if matrix_list else None


def parse_populations(populations_string):
    """
    Parse comma-delimited populations string.
    
    Parameters:
    - populations_string: Comma-delimited string of population names
    
    Returns:
    - list: List of population names
    """
    if not populations_string:
        return ["AFR", "EUR"]  # Default populations
    
    # Split by comma and strip whitespace
    pop_list = [name.strip().upper() for name in populations_string.split(',')]
    
    # Remove empty strings
    pop_list = [name for name in pop_list if name]
    
    return pop_list if pop_list else ["AFR", "EUR"]


def parse_heritabilities(heritabilities_string):
    """
    Parse comma-delimited heritabilities string.
    
    Parameters:
    - heritabilities_string: Comma-delimited string of heritability values
    
    Returns:
    - list: List of heritability values as floats
    """
    if not heritabilities_string:
        return None
    
    # Split by comma and convert to float
    try:
        herit_list = [float(h.strip()) for h in heritabilities_string.split(',')]
        return herit_list
    except ValueError as e:
        print(f"Error parsing heritabilities: {e}")
        return None


def collect_unique_trait_combinations(gwas_base_dir, populations, heritabilities):
    """
    Collect all unique trait combinations from GWAS files across populations.
    
    A trait combination is defined by:
    - Phenotype number
    - Number of causal variants
    - Trait type
    - Chromosome
    
    Parameters:
    - gwas_base_dir: Base directory containing GWAS files
    - populations: List of population names
    - heritabilities: List of heritability values (one per population)
    
    Returns:
    - list: List of unique trait combination tuples
    """
    trait_combinations = set()
    
    # Scan all populations to find all unique trait combinations
    for pop in populations:
        pop_gwas_dir = os.path.join(gwas_base_dir, pop)
        
        if not os.path.exists(pop_gwas_dir):
            print(f"Warning: Population directory not found: {pop_gwas_dir}")
            continue
        
        # Find all GWAS files
        all_gwas_files = glob.glob(os.path.join(pop_gwas_dir, "*.step2.txt"))
        
        for gwas_file in all_gwas_files:
            filename = os.path.basename(gwas_file)
            parsed = parse_gwas_filename(filename)
            
            if parsed:
                # Create a unique identifier for this trait combination
                trait_key = (
                    parsed['pheno_num'],
                    parsed['n_causal'],
                    parsed['trait_type'],
                    parsed['chromosome']
                )
                trait_combinations.add(trait_key)
    
    print(f"\nFound {len(trait_combinations)} unique trait combinations across all populations")
    
    return sorted(list(trait_combinations))


def create_susiex_commands(gwas_base_dir, ld_base_dir, geno_base_dir, output_base_dir, 
                          wrapper_script, output_commands_file, populations, plink_conv, 
                          heritabilities, matrix_filter=None, n_processes=None, skip_existing=False):
    """
    Generate SuSiEx multi-population fine-mapping commands for all trait combinations 
    using parallel processing.
    
    Parameters:
    - gwas_base_dir: Base directory containing GWAS summary statistics
    - ld_base_dir: Base directory containing LD matrices
    - geno_base_dir: Base directory containing genotype files
    - output_base_dir: Base directory for fine-mapping results
    - wrapper_script: Path to finemap_SuSiEx.py script
    - output_commands_file: File to write commands to
    - populations: List of population names
    - plink_conv: PLINK file naming convention
    - heritabilities: List of heritability values (one per population)
    - matrix_filter: Optional list of matrix subdirectory names to include
    - n_processes: Number of processes to use (None = use all available cores)
    - skip_existing: Whether to skip generating commands for existing complete output files
    """
    # Print heritability information
    print(f"\nPer-population heritabilities: {heritabilities}")
    for i, pop in enumerate(populations):
        print(f"  {pop}: {heritabilities[i]}")
    
    # Automatically detect matrix types
    print(f"\nDetecting matrix subdirectories in: {ld_base_dir}")
    if matrix_filter:
        print(f"Filtering for matrix types: {matrix_filter}")
    
    matrix_types = detect_matrix_subdirs(ld_base_dir, matrix_filter)
    
    if not matrix_types:
        if matrix_filter:
            print(f"Error: No valid matrix configurations found matching filter: {matrix_filter}")
        else:
            print("Error: No valid matrix configurations found!")
        return []
    
    print(f"Matrix types to process: {matrix_types}")
    print(f"Populations to process: {populations}")
    
    if skip_existing:
        print(f"File existence checking: ENABLED - will skip existing complete output files")
    else:
        print(f"File existence checking: DISABLED - will generate all commands")
    
    # Create output directories now
    print(f"\nCreating output directories...")
    for matrix_type in matrix_types:
        output_dir = os.path.join(output_base_dir, matrix_type)
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, "timing"), exist_ok=True)
        print(f"  Created: {output_dir}")
    
    # Collect all unique trait combinations
    print(f"\nCollecting trait combinations from: {gwas_base_dir}")
    trait_combinations = collect_unique_trait_combinations(gwas_base_dir, populations, heritabilities)
    
    print(f"Found {len(trait_combinations)} unique trait combinations to process")
    
    # Determine number of processes
    if n_processes is None:
        n_processes = cpu_count()
    print(f"Using {n_processes} processes for parallel processing")
    
    # Prepare arguments for parallel processing
    process_args = []
    for pheno_num, n_causal, trait_type, chromosome in trait_combinations:
        args_tuple = (
            pheno_num, n_causal, trait_type, chromosome, heritabilities,
            gwas_base_dir, ld_base_dir, geno_base_dir, matrix_types,
            output_base_dir, wrapper_script, populations, plink_conv, skip_existing
        )
        process_args.append(args_tuple)
    
    # Process trait combinations in parallel
    print(f"\nProcessing {len(process_args)} trait combinations in parallel...")
    all_commands = []
    all_warnings = []
    all_skipped = []
    
    with Pool(processes=n_processes) as pool:
        results = pool.map(process_trait_combination, process_args)
    
    # Collect results
    processed_traits = set()
    for commands, warnings, skipped_info in results:
        all_commands.extend(commands)
        all_warnings.extend(warnings)
        all_skipped.extend(skipped_info)
        
        # Track processed traits for summary
        for cmd_info in commands:
            if len(cmd_info) >= 4:
                processed_traits.add(cmd_info[3])  # pheno_pattern is at index 3
    
    # Print warnings
    if all_warnings:
        print(f"\nWarnings encountered ({len(all_warnings)}):")
        for warning in all_warnings[:10]:  # Show first 10 warnings
            print(f"  {warning}")
        if len(all_warnings) > 10:
            print(f"  ... and {len(all_warnings) - 10} more warnings")
    
    # Print skipped files summary
    if all_skipped:
        print(f"\nExisting files skipped ({len(all_skipped)}):")
        skipped_by_matrix = {}
        for skipped in all_skipped:
            matrix_type = skipped['matrix_type']
            if matrix_type not in skipped_by_matrix:
                skipped_by_matrix[matrix_type] = 0
            skipped_by_matrix[matrix_type] += 1
        
        for matrix_type, count in skipped_by_matrix.items():
            print(f"  {matrix_type}: {count} files")
    
    # Print processing summary
    print(f"\nProcessing Summary:")
    print(f"  Total commands generated: {len(all_commands)}")
    if skip_existing:
        print(f"  Existing files skipped: {len(all_skipped)}")
        print(f"  New commands needed: {len(all_commands)}")
    
    # Show examples of commands
    if all_commands:
        print(f"\nExample commands (first 3):")
        for cmd_info in all_commands[:3]:
            if len(cmd_info) >= 9:
                _, _, matrix_type, pheno_pattern, trait_type, chromosome, heritabilities_vals, pops, output_filename = cmd_info
                herit_str = ','.join([str(h) for h in heritabilities_vals])
                print(f"  {pheno_pattern} ({trait_type}, {chromosome}, h²=[{herit_str}]) -> {matrix_type}")
                print(f"    Populations: {','.join(pops)}")
                print(f"    Output: {output_filename}")
        
        if len(all_commands) > 3:
            print(f"  ... and {len(all_commands) - 3} more commands")
    
    # Extract just the commands for writing to file
    commands_only = [cmd_info[0] for cmd_info in all_commands]
    
    # Write commands to file
    if output_commands_file:
        print(f"\nWriting {len(commands_only)} commands to: {output_commands_file}")
        with open(output_commands_file, 'w') as f:
            for cmd in commands_only:
                f.write(cmd + '\n')
    
    return commands_only


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description='Generate SuSiEx multi-population fine-mapping commands with heterogeneous effects (optimized for HPC)'
    )
    
    # Define default paths
    default_wrapper = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/finemap_SuSiEx.py"
    default_gwas = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
    default_geno = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/simulated_genotypes"
    default_ld = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
    default_output = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results/susiex_results"
    default_commands = f"{default_output}/susiex_commands.txt"
    
    parser.add_argument('--wrapper-script', '-w', default=default_wrapper,
                       help=f'Path to finemap_SuSiEx.py script (default: {default_wrapper})')
    
    parser.add_argument('--gwas-dir', '-g', default=default_gwas,
                       help=f'Base directory containing GWAS summary statistics with population subdirectories (default: {default_gwas})')
    
    parser.add_argument('--ld-dir', '-l', default=default_ld,
                       help=f'Base directory containing LD matrices (default: {default_ld})')
    
    parser.add_argument('--geno-dir', '-n', default=default_geno,
                       help=f'Base directory of Genotype (Plink) Files (default: {default_geno})')
    
    parser.add_argument('--output-dir', '-o', default=default_output,
                       help=f'Base directory for fine-mapping results (default: {default_output})')
    
    parser.add_argument('--commands-file', '-c', default=default_commands,
                       help=f'Output file for commands (default: {default_commands})')
    
    parser.add_argument('--populations', '-pop', type=str, default="AFR,EUR",
                       help='Comma-delimited list of population names (default: "AFR,EUR")')
    
    parser.add_argument('--plink_conv', '-pc', type=str, default="CHR#.plink.maf",
                       help='Frequency ref file naming convention (default: "CHR#.plink.maf")')
    
    parser.add_argument('--heritabilities', '-h2', type=str, required=True,
                       help='REQUIRED: Comma-delimited list of per-variant heritability values, one per population (e.g., "0.00005,0.0003" or "0,0.0003" for null in first pop)')
    
    parser.add_argument('--matrix-types', '-m', type=str, default=None,
                       help='Comma-delimited list of matrix subdirectory names to process (e.g., "raw_in-sample,PC-corrected"). If not specified, all detected matrix types will be processed.')
    
    parser.add_argument('--processes', '-p', type=int, default=None,
                       help='Number of processes to use for parallel processing (default: use all available cores)')
    
    parser.add_argument('--skip-existing', '-s', action='store_true',
                       help='Skip generating commands for output files that already exist and appear complete')
    
    parser.add_argument('--dry-run', '-d', action='store_true',
                       help='Print first 10 commands without writing to file')
    
    parser.add_argument('--list-matrices', action='store_true',
                       help='List available matrix subdirectories and exit')
    
    args = parser.parse_args()
    
    # Parse populations
    populations = parse_populations(args.populations)
    print(f"Populations: {populations}")
    
    # Parse heritabilities
    heritabilities = parse_heritabilities(args.heritabilities)
    
    if heritabilities is None:
        print("Error: Invalid heritabilities argument")
        return 1
    
    if len(heritabilities) != len(populations):
        print(f"Error: Number of heritabilities ({len(heritabilities)}) must match number of populations ({len(populations)})")
        return 1
    
    # Handle list-matrices option
    if args.list_matrices:
        print(f"Scanning for matrix subdirectories in: {args.ld_dir}")
        
        if not os.path.exists(args.ld_dir):
            print(f"Error: LD directory does not exist: {args.ld_dir}")
            return 1
        
        matrix_types = detect_matrix_subdirs(args.ld_dir)
        
        if matrix_types:
            print("\nAvailable matrix types:")
            for matrix_type in matrix_types:
                print(f"  {matrix_type}")
                
                # Show sample matrix files for this type from first population
                first_pop = populations[0]
                search_dir = os.path.join(args.ld_dir, matrix_type, first_pop)
                
                if os.path.exists(search_dir):
                    sample_files = glob.glob(os.path.join(search_dir, "*.vcor1"))[:3]
                    if sample_files:
                        print(f"    Sample files ({first_pop}):")
                        for f in sample_files:
                            print(f"      {os.path.basename(f)}")
                        if len(glob.glob(os.path.join(search_dir, "*.vcor1"))) > 3:
                            print("      ... and more")
            
            print("\nExample usage:")
            print("  --matrix-types \"raw_in-sample,PC-corrected\"")
            print(f"  --matrix-types \"{','.join(matrix_types)}\"")
        else:
            print("No valid matrix configurations found.")
        
        return 0
    
    # Check if directories exist
    for dir_path, name in [(args.gwas_dir, 'GWAS base'), (args.ld_dir, 'LD matrices')]:
        if not os.path.exists(dir_path):
            print(f"Error: {name} directory does not exist: {dir_path}")
            return 1
    
    # Check population directories exist
    for pop in populations:
        pop_dir = os.path.join(args.gwas_dir, pop)
        if not os.path.exists(pop_dir):
            print(f"Error: Population directory not found: {pop_dir}")
            return 1
    
    if not os.path.exists(args.wrapper_script):
        print(f"Error: Wrapper script does not exist: {args.wrapper_script}")
        return 1
    
    # Parse matrix filter
    matrix_filter = parse_matrix_filter(args.matrix_types)
    if matrix_filter:
        print(f"Matrix filter applied: {matrix_filter}")
    
    # Validate number of processes
    if args.processes is not None and args.processes < 1:
        print(f"Error: Number of processes must be >= 1, got {args.processes}")
        return 1
    
    # Generate commands
    commands = create_susiex_commands(
        gwas_base_dir=args.gwas_dir,
        ld_base_dir=args.ld_dir,
        geno_base_dir=args.geno_dir,
        output_base_dir=args.output_dir,
        wrapper_script=args.wrapper_script,
        output_commands_file=args.commands_file if not args.dry_run else None,
        populations=populations,
        plink_conv=args.plink_conv,
        heritabilities=heritabilities,
        matrix_filter=matrix_filter,
        n_processes=args.processes,
        skip_existing=args.skip_existing
    )
    
    if args.dry_run:
        print("\nDry run - first 10 commands:")
        for i, cmd in enumerate(commands[:10]):
            print(f"{i+1:3d}: {cmd}")
        if len(commands) > 10:
            print(f"... and {len(commands) - 10} more commands")
    else:
        print(f"\nCommands written to: {args.commands_file}")
        print(f"Total commands: {len(commands)}")
        print("\nTo run all commands, you can use:")
        print(f"  bash {args.commands_file}")
        print("\nOr submit as array job to SLURM (see script comments for example)")
    
    return 0


if __name__ == "__main__":
    exit(main())
