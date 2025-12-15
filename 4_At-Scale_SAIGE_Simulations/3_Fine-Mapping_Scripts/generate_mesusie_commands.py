#!/usr/bin/env python3

"""
Script to generate MESuSiE multi-population fine-mapping commands for all loci.

This script creates a command file where each line runs MESuSiE fine-mapping
for a different locus across all available matrix types and populations.

Adapted from generate_susiex_commands.py to work with finemap_MESuSiE.R wrapper.

Key differences from SuSiEx version:
- No PLINK reference files needed (--ref_panels, --geno-dir, --plink_conv removed)
- Uses --ancestry_labels parameter instead of inferring from paths
- Calls R script instead of Python script
- Different output file structure
"""

import os
import glob
import argparse
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

def find_matrix_files_by_population(ld_base_dir, matrix_type, pheno_causal, populations):
    """
    Search for matrix files for all populations that match the given phenotype/causal pattern.
    
    Parameters:
    - ld_base_dir: Base directory containing LD matrices
    - matrix_type: Matrix type directory name (e.g., "raw_in-sample")
    - pheno_causal: Phenotype pattern (e.g., "Pheno1.2_causal")
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
        
        # Search for files that start with the pheno_causal pattern and end with .vcor1
        pattern = os.path.join(search_dir, f"{pheno_causal}.*.vcor1")
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

def normalize_heritability(heritability_str):
    """
    Normalize heritability string for comparison.
    Converts decimal representations to consistent format.
    
    Parameters:
    - heritability_str: Heritability string from filename (e.g., "heritability_0.0003")
    
    Returns:
    - float: Normalized heritability value
    """
    try:
        # Remove "heritability_" prefix if present
        if heritability_str.startswith("heritability_"):
            heritability_str = heritability_str.replace("heritability_", "")
        
        # Convert to float for consistent comparison
        return float(heritability_str)
    except (ValueError, TypeError):
        return None

def parse_gwas_filename(meta_filename, target_heritability=None):
    """
    Parse GWAS filename to extract components, handling decimal heritability values.
    Expected format: Pheno#.#_causal.heritability_#.#.continuous/binary_0.25.chr#.step2.txt
    
    Parameters:
    - meta_filename: The GWAS filename
    - target_heritability: Target heritability value to filter for (None = no filter)
    
    Returns:
    - dict: Dictionary with parsed components or None if parsing fails
    """
    # Remove the .step2.txt extension
    base_name = meta_filename.replace('.step2.txt', '')
    
    # Split by dots
    parts = base_name.split('.')
    if len(parts) < 5:  # Need at least: Pheno#, #_causal, heritability_#, #, trait_type, chr#
        return None
    
    try:
        # First two parts: Pheno# and #_causal
        pheno_causal = parts[0] + "." + parts[1]  # e.g., "Pheno1.2_causal"
        
        # Find heritability parts - they start with "heritability_" and continue until we hit trait type
        heritability_parts = []
        trait_type_index = None
        
        for i in range(2, len(parts)):
            part = parts[i]
            # Check if this is a trait type (continuous or binary_0.25)
            if part in ['continuous', 'binary_0'] or (part == '25' and i > 2 and parts[i-1] == 'binary_0'):
                trait_type_index = i
                break
            else:
                heritability_parts.append(part)
        
        if trait_type_index is None:
            return None
        
        # Reconstruct heritability (e.g., "heritability_0.01" from ["heritability_0", "01"])
        heritability = ".".join(heritability_parts)
        
        # Filter by heritability if specified
        if target_heritability is not None:
            parsed_heritability = normalize_heritability(heritability)
            if parsed_heritability is None:
                return None
            
            # Compare heritability values with small tolerance for floating point comparison
            if abs(parsed_heritability - target_heritability) > 1e-10:
                return None
        
        # Handle trait type - could be "continuous" or "binary_0.25"
        if parts[trait_type_index] == 'continuous':
            trait_type = 'continuous'
            chromosome_index = trait_type_index + 1
        elif parts[trait_type_index] == 'binary_0' and trait_type_index + 1 < len(parts) and parts[trait_type_index + 1] == '25':
            trait_type = 'binary_0.25'
            chromosome_index = trait_type_index + 2
        else:
            return None
        
        # Get chromosome
        if chromosome_index >= len(parts):
            return None
        chromosome = parts[chromosome_index]
        
        return {
            'pheno_causal': pheno_causal,
            'heritability': heritability,
            'trait_type': trait_type,
            'chromosome': chromosome
        }
    
    except (IndexError, ValueError):
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
        
        # For MESuSiE output files, check if they have a header and at least one data row
        with open(output_path, 'r') as f:
            lines = f.readlines()
        
        # Should have at least 2 lines (header + 1 data row)
        if len(lines) < 2:
            return False
        
        # Check if header contains expected columns for MESuSiE output
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

def process_gwas_file_pattern(args_tuple):
    """
    Process a GWAS file pattern to find matching files across populations and generate MESuSiE commands.
    This function is designed to be called in parallel.
    
    Parameters:
    - args_tuple: Tuple containing (gwas_file_pattern, gwas_base_dir, ld_base_dir, matrix_types,
                  output_base_dir, wrapper_script, populations, skip_existing)
    
    Returns:
    - tuple: (commands, warnings, skipped_info)
    """
    gwas_file_pattern, gwas_base_dir, ld_base_dir, matrix_types, output_base_dir, wrapper_script, populations, skip_existing, target_heritability = args_tuple
    
    commands = []
    warnings = []
    skipped_info = []
    
    # Parse the filename pattern from the sample file
    sample_filename = os.path.basename(gwas_file_pattern)
    parsed = parse_gwas_filename(sample_filename)
    if parsed is None:
        warnings.append(f"Could not parse filename: {sample_filename}")
        return commands, warnings, skipped_info
    
    pheno_causal = parsed['pheno_causal']
    heritability = parsed['heritability']
    trait_type = parsed['trait_type']
    chromosome = parsed['chromosome']
    
    # Find matching GWAS files in each population directory
    gwas_files = {}
    for pop in populations:
        pop_gwas_dir = os.path.join(gwas_base_dir, pop)
        if not os.path.exists(pop_gwas_dir):
            warnings.append(f"Population directory not found: {pop_gwas_dir}")
            continue
        
        # Look for files that match this exact pattern
        expected_filename = f"{pheno_causal}.{heritability}.{trait_type}.{chromosome}.step2.txt"
        expected_path = os.path.join(pop_gwas_dir, expected_filename)
        
        if os.path.exists(expected_path):
            gwas_files[pop] = expected_path
        else:
            # Try broader pattern matching
            broad_pattern = f"{pheno_causal}.*.{trait_type}.{chromosome}.step2.txt"
            matching_files = glob.glob(os.path.join(pop_gwas_dir, broad_pattern))
            
            if not matching_files:
                # Even broader - just look for files with the phenotype and check details
                all_causal_files = glob.glob(os.path.join(pop_gwas_dir, f"{pheno_causal}*.step2.txt"))
                matching_files = []
                for f in all_causal_files:
                    fname = os.path.basename(f)
                    parsed_f = parse_gwas_filename(fname)
                    if (parsed_f and 
                        parsed_f['pheno_causal'] == pheno_causal and 
                        parsed_f['chromosome'] == chromosome and 
                        parsed_f['trait_type'] == trait_type):
                        matching_files.append(f)
            
            if len(matching_files) == 1:
                gwas_files[pop] = matching_files[0]
            elif len(matching_files) > 1:
                gwas_files[pop] = matching_files[0]
                warnings.append(f"Multiple GWAS files matched for {pop}: {[os.path.basename(f) for f in matching_files]}, using {os.path.basename(matching_files[0])}")
            else:
                warnings.append(f"No matching GWAS file found for {pop} with pattern {expected_filename}")
                continue
    
    # If we don't have GWAS files for all populations, skip
    if len(gwas_files) != len(populations):
        warnings.append(f"Missing GWAS files for some populations - found {len(gwas_files)}/{len(populations)} populations for {pheno_causal}")
        return commands, warnings, skipped_info
    
    # Get sample sizes for populations
    sample_sizes = get_population_sample_sizes(populations)
    
    # For each detected matrix type
    for matrix_type in matrix_types:
        # Find matrix files for all populations
        matrix_files = find_matrix_files_by_population(ld_base_dir, matrix_type, pheno_causal, populations)
        
        # Check if we have matrix files for all populations
        missing_matrices = [pop for pop, path in matrix_files.items() if path is None]
        if missing_matrices:
            warnings.append(f"Missing matrix files for populations {missing_matrices} in {pheno_causal}, {matrix_type}")
            continue
        
        # Construct output directory and file
        output_dir = os.path.join(output_base_dir, matrix_type)
        
        # Output filename format: Pheno#.#_causal.continuous/binary_0.25.matrix_type.txt
        output_filename = f"{pheno_causal}.{trait_type}.{matrix_type}.txt"
        output_path = os.path.join(output_dir, output_filename)
        
        # Check if output file already exists and is complete
        if skip_existing and check_output_completeness(output_path):
            skipped_info.append({
                'pheno_causal': pheno_causal,
                'matrix_type': matrix_type,
                'trait_type': trait_type,
                'chromosome': chromosome,
                'heritability': heritability,
                'output_path': output_path,
                'populations': populations
            })
            continue
        
        # Create comma-delimited lists for MESuSiE
        gwas_files_list = ",".join([gwas_files[pop] for pop in populations])
        ld_files_list = ",".join([matrix_files[pop] for pop in populations])
        n_gwas_list = ",".join([str(sample_sizes[pop]) for pop in populations])
        ancestry_labels_list = ",".join(populations)
        
        # Create the MESuSiE command
        # Rscript finemap_MESuSiE.R --gwas_files file1,file2 --ld_files matrix1,matrix2 --n_gwas n1,n2 --ancestry_labels pop1,pop2 --output_file output.txt
        mesusie_cmd = (f"Rscript {wrapper_script} "
                      f"--gwas_files {gwas_files_list} "
                      f"--ld_files {ld_files_list} "
                      f"--n_gwas {n_gwas_list} "
                      f"--ancestry_labels {ancestry_labels_list} "
                      f"--output_file {output_path}")
        
        # Store command with metadata
        commands.append((mesusie_cmd, output_dir, matrix_type, pheno_causal, trait_type,
                        chromosome, heritability, populations, output_filename))
    
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

def create_mesusie_commands(gwas_base_dir, ld_base_dir, output_base_dir, wrapper_script, output_commands_file,
                           populations, matrix_filter=None, n_processes=None, skip_existing=False, target_heritability=None):
    """
    Generate MESuSiE multi-population fine-mapping commands for all loci using parallel processing.
    
    Parameters:ls 
    - gwas_base_dir: Base directory containing GWAS summary statistics (with population subdirs)
    - ld_base_dir: Base directory containing LD matrices
    - output_base_dir: Base directory for fine-mapping results
    - wrapper_script: Path to finemap_MESuSiE.R script
    - output_commands_file: File to write commands to
    - populations: List of population names
    - matrix_filter: Optional list of matrix subdirectory names to include
    - n_processes: Number of processes to use (None = use all available cores)
    - skip_existing: Whether to skip generating commands for existing complete output files
    - target_heritability: Target heritability value to filter for (None = no filter)
    """
    
    # Print heritability filter information
    if target_heritability is not None:
        print(f"Heritability filter: {target_heritability}")
        print(f"Only processing traits with heritability = {target_heritability}")
    else:
        print(f"Heritability filter: DISABLED - processing all heritability values")
    
    # Automatically detect matrix types
    print(f"Detecting matrix subdirectories in: {ld_base_dir}")
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
        os.makedirs(os.path.join(output_dir, "rds_files"), exist_ok=True)  # MESuSiE creates RDS files
        print(f" Created: {output_dir}")
    
    # Find all unique GWAS file patterns from the first population directory
    first_pop = populations[0]
    pop_gwas_dir = os.path.join(gwas_base_dir, first_pop)
    print(f"\nSearching for causal files in: {pop_gwas_dir}")
    
    if not os.path.exists(pop_gwas_dir):
        print(f"Error: Population GWAS directory does not exist: {pop_gwas_dir}")
        return []
    
    causal_files = glob.glob(os.path.join(pop_gwas_dir, "*_causal*.step2.txt"))
    # Filter out files with unwanted patterns
    causal_files = [f for f in causal_files 
                    if 'drop_causal' not in os.path.basename(f) 
                    and 'matched_to_' not in os.path.basename(f)]
    print(f"Found {len(causal_files)} standard causal files")
    
    # Group files by phenotype pattern to avoid duplicates
    unique_patterns = set()
    gwas_file_patterns = []
    
    for causal_file in causal_files:
        filename = os.path.basename(causal_file)
        parsed = parse_gwas_filename(filename, target_heritability=target_heritability)
        if parsed:
            pattern_key = f"{parsed['pheno_causal']}.{parsed['trait_type']}.{parsed['chromosome']}"
            if pattern_key not in unique_patterns:
                unique_patterns.add(pattern_key)
                gwas_file_patterns.append(causal_file)
    
    print(f"Found {len(gwas_file_patterns)} unique phenotype patterns")
    
    # Determine number of processes
    if n_processes is None:
        n_processes = cpu_count()
    print(f"Using {n_processes} processes for parallel processing")
    
    # Prepare arguments for parallel processing
    process_args = []
    for gwas_file_pattern in gwas_file_patterns:
        args_tuple = (gwas_file_pattern, gwas_base_dir, ld_base_dir, matrix_types,
                     output_base_dir, wrapper_script, populations, skip_existing, target_heritability)
        process_args.append(args_tuple)
    
    # Process file patterns in parallel
    print(f"\nProcessing {len(process_args)} file patterns in parallel...")
    
    all_commands = []
    all_warnings = []
    all_skipped = []
    
    with Pool(processes=n_processes) as pool:
        results = pool.map(process_gwas_file_pattern, process_args)
        
    # Collect results
    processed_phenotypes = set()
    for commands, warnings, skipped_info in results:
        all_commands.extend(commands)
        all_warnings.extend(warnings)
        all_skipped.extend(skipped_info)
        
        # Track processed phenotypes for summary
        for cmd_info in commands:
            if len(cmd_info) >= 4:  # Check if we have the expected tuple structure
                processed_phenotypes.add(cmd_info[3])  # pheno_causal is at index 3
                
    # Print warnings
    if all_warnings:
        print(f"\nWarnings encountered ({len(all_warnings)}):")
        for warning in all_warnings:
            print(f"  {warning}")
            
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
            
        # Show a few examples
        print(f"  Examples:")
        for i, skipped in enumerate(all_skipped[:3]):
            print(f"    {skipped['pheno_causal']}.{skipped['trait_type']} (h²={skipped['heritability']}) -> {skipped['matrix_type']} (already complete)")
        if len(all_skipped) > 3:
            print(f"    ... and {len(all_skipped) - 3} more")
            
    # Print processing summary
    print(f"\nProcessing Summary:")
    print(f"  Total commands generated: {len(all_commands)}")
    if skip_existing:
        print(f"  Existing files skipped: {len(all_skipped)}")
        print(f"  New commands needed: {len(all_commands)}")
        
    if target_heritability is not None:
        print(f"  Heritability filter applied: only traits with h² = {target_heritability}")
        
    # Show examples of commands
    for cmd_info in all_commands[:3]:  # Show first 3 as examples
        if len(cmd_info) >= 9:
            _, _, matrix_type, pheno_causal, trait_type, chromosome, heritability, pops, output_filename = cmd_info
            print(f"    {pheno_causal} ({trait_type}, {chromosome}, h²={heritability}) -> {matrix_type}")
            print(f"      Populations: {','.join(pops)}")
            print(f"      Output: {output_filename}")
    if len(all_commands) > 3:
        print(f"    ... and {len(all_commands) - 3} more commands")
        
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
    parser = argparse.ArgumentParser(description='Generate MESuSiE multi-population fine-mapping commands (optimized for HPC)')

    # Define default paths (adapted for MESuSiE)
    default_wrapper = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/finemap_MESuSiE.R"
    default_gwas = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results"
    default_ld = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
    default_output = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results/mesusie_results"
    default_commands = f"{default_output}/mesusie_commands.txt"

    parser.add_argument('--wrapper-script', '-w', default=default_wrapper,
                       help=f'Path to finemap_MESuSiE.R script (default: {default_wrapper})')
    parser.add_argument('--gwas-dir', '-g', default=default_gwas,
                       help=f'Base directory containing GWAS summary statistics with population subdirectories (default: {default_gwas})')
    parser.add_argument('--ld-dir', '-l', default=default_ld,
                       help=f'Base directory containing LD matrices (default: {default_ld})')
    parser.add_argument('--output-dir', '-o', default=default_output,
                       help=f'Base directory for fine-mapping results (default: {default_output})')
    parser.add_argument('--commands-file', '-c', default=default_commands,
                       help=f'Output file for commands (default: {default_commands})')
    parser.add_argument('--populations', '-pop', type=str, default="AFR,EUR",
                       help='Comma-delimited list of population names (default: "AFR,EUR")')
    parser.add_argument('--matrix-types', '-m', type=str, default=None,
                       help='Comma-delimited list of matrix subdirectory names to process (e.g., "raw_in-sample,PC-corrected"). If not specified, all detected matrix types will be processed.')
    parser.add_argument('--processes', '-p', type=int, default=None,
                       help='Number of processes to use for parallel processing (default: use all available cores)')
    parser.add_argument('--skip-existing', '-s', action='store_true',
                       help='Skip generating commands for output files that already exist and appear complete')
    parser.add_argument('--heritability', type=float, default=0.0003,
                        help='Target heritability value to filter traits (default: 0.0003)')
    parser.add_argument('--dry-run', '-d', action='store_true',
                       help='Print first 10 commands without writing to file')
    parser.add_argument('--list-matrices', action='store_true',
                       help='List available matrix subdirectories and exit')

    args = parser.parse_args()

    # Parse populations
    populations = parse_populations(args.populations)
    print(f"Populations: {populations}")
    
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
            print(f"Error: Population directory does not exist: {pop_dir}")
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
    commands = create_mesusie_commands(
        gwas_base_dir=args.gwas_dir,
        ld_base_dir=args.ld_dir,
        output_base_dir=args.output_dir,
        wrapper_script=args.wrapper_script,
        output_commands_file=args.commands_file if not args.dry_run else None,
        populations=populations,
        matrix_filter=matrix_filter,
        n_processes=args.processes,
        skip_existing=args.skip_existing,
        target_heritability=args.heritability
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
