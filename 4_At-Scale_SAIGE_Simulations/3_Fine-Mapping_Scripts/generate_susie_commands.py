#!/usr/bin/env python3
"""
Script to generate SuSiE fine-mapping commands for all loci.

This script creates a command file where each line runs SuSiE fine-mapping
for a different locus across all available matrix types. 

Optimized for parallel processing on HPC systems.
"""

import os
import glob
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial

def detect_matrix_subdirs(ld_base_dir, matrix_filter=None):
    """
    Automatically detect matrix subdirectories and their subfolder structure.

    Parameters:
    - ld_base_dir: Base directory containing LD matrices
    - matrix_filter: Optional list of matrix subdirectory names to include

    Returns:
    - dict: mapping of matrix_type -> subfolder_name
    """
    matrix_configs = {}

    # Look for all subdirectories in ld_base_dir
    if not os.path.exists(ld_base_dir):
        print(f"Warning: LD base directory does not exist: {ld_base_dir}")
        return matrix_configs

    for item in os.listdir(ld_base_dir):
        matrix_dir = os.path.join(ld_base_dir, item)
        if not os.path.isdir(matrix_dir):
            continue

        # Apply matrix filter if specified
        if matrix_filter and item not in matrix_filter:
            print(f"Skipping matrix type (not in filter): {item}")
            continue

        # Check for 'weighted' or 'Meta' subdirectory
        weighted_path = os.path.join(matrix_dir, "weighted")
        meta_path = os.path.join(matrix_dir, "Meta")

        if os.path.isdir(weighted_path):
            matrix_configs[item] = "weighted"
            print(f"Found matrix type: {item} -> weighted")
        elif os.path.isdir(meta_path):
            matrix_configs[item] = "Meta"  
            print(f"Found matrix type: {item} -> Meta")
        else:
            print(f"Warning: No 'weighted' or 'Meta' subdirectory found in {matrix_dir}")

    return matrix_configs

def find_matrix_file(ld_base_dir, matrix_type, subfolder, pheno_causal):
    """
    Search for a matrix file that matches the given phenotype/causal pattern.

    Parameters:
    - ld_base_dir: Base directory containing LD matrices
    - matrix_type: Matrix type directory name (e.g., "beta-inferred_2000")
    - subfolder: Subfolder name ("weighted" or "Meta")
    - pheno_causal: Phenotype pattern (e.g., "Pheno1.2_causal")

    Returns:
    - str: Full path to matching matrix file, or None if not found
    """
    search_dir = os.path.join(ld_base_dir, matrix_type, subfolder)

    if not os.path.exists(search_dir):
        return None

    # Search for files that start with the pheno_causal pattern and end with .vcor1
    pattern = os.path.join(search_dir, f"{pheno_causal}.*.vcor1")
    matching_files = glob.glob(pattern)

    if len(matching_files) == 1:
        return matching_files[0]
    elif len(matching_files) > 1:
        # For parallel processing, we'll handle warnings differently
        return matching_files[0]  # Just use first match
    else:
        # No matches found
        return None

def process_causal_dir(args_tuple):
    """
    Process a single causal directory to generate commands.
    This function is designed to be called in parallel.

    Parameters:
    - args_tuple: Tuple containing (causal_dir, ld_base_dir, matrix_configs, output_base_dir, wrapper_script)

    Returns:
    - list: List of commands generated for this directory
    """
    causal_dir, ld_base_dir, matrix_configs, output_base_dir, wrapper_script = args_tuple

    commands = []
    warnings = []

    if not os.path.isdir(causal_dir):
        return commands, warnings

    # Find the .meta.out file in this directory
    meta_files = glob.glob(os.path.join(causal_dir, "*.meta.out"))

    if len(meta_files) != 1:
        warning = f"Expected 1 .meta.out file in {causal_dir}, found {len(meta_files)}"
        if len(meta_files) > 1:
            warning += f" - Files: {[os.path.basename(f) for f in meta_files]}"
        warnings.append(warning)
        return commands, warnings

    meta_file = meta_files[0]
    meta_filename = os.path.basename(meta_file)

    # Parse the filename to extract components
    # Format: Pheno#.#_causal.heritability_#.continuous/binary_0.25.chr#.meta.out
    filename_parts = meta_filename.replace('.meta.out', '').split('.')

    if len(filename_parts) < 4:
        warnings.append(f"Unexpected filename format: {meta_filename} - Parts: {filename_parts}")
        return commands, warnings

    try:
        pheno_causal = filename_parts[0] + "." + filename_parts[1]  # e.g., "Pheno1.2_causal"
        heritability = filename_parts[2]  # e.g., "heritability_0.01"
        trait_type = filename_parts[3]    # e.g., "continuous" or "binary_0.25"
        chromosome = filename_parts[4]    # e.g., "chr1"

    except IndexError:
        warnings.append(f"Could not parse filename: {meta_filename}")
        return commands, warnings

    # For each detected matrix configuration
    for matrix_type, subfolder in matrix_configs.items():
        # Search for the actual matrix file that matches this phenotype
        matrix_path = find_matrix_file(ld_base_dir, matrix_type, subfolder, pheno_causal)

        if matrix_path is None:
            warnings.append(f"No matrix file found for {pheno_causal} in {matrix_type}/{subfolder}")
            continue

        # Construct output directory and file
        output_dir = os.path.join(output_base_dir, matrix_type)
        output_filename = f"{pheno_causal}.{subfolder}.txt"
        output_path = os.path.join(output_dir, output_filename)

        # Create the SuSiE command
        # Usage: Rscript finemap_SuSiE.R summary_stats_file ld_file out_file [optional_params]
        susie_cmd = f"Rscript {wrapper_script} {meta_file} {matrix_path} {output_path}"

        commands.append((susie_cmd, output_dir, matrix_type, pheno_causal, trait_type, chromosome, os.path.basename(matrix_path)))

    return commands, warnings

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

def create_susie_commands(gwas_dir, ld_base_dir, output_base_dir, wrapper_script, output_commands_file, matrix_filter=None, n_processes=None):
    """
    Generate SuSiE fine-mapping commands for all loci using parallel processing.

    Parameters:
    - gwas_dir: Directory containing GWAS summary statistics
    - ld_base_dir: Base directory containing LD matrices  
    - output_base_dir: Base directory for fine-mapping results
    - wrapper_script: Path to finemap_SuSiE.R script
    - output_commands_file: File to write commands to
    - matrix_filter: Optional list of matrix subdirectory names to include
    - n_processes: Number of processes to use (None = use all available cores)
    """

    # Automatically detect matrix configurations
    print(f"Detecting matrix subdirectories in: {ld_base_dir}")
    if matrix_filter:
        print(f"Filtering for matrix types: {matrix_filter}")

    matrix_configs = detect_matrix_subdirs(ld_base_dir, matrix_filter)

    if not matrix_configs:
        if matrix_filter:
            print(f"Error: No valid matrix configurations found matching filter: {matrix_filter}")
        else:
            print("Error: No valid matrix configurations found!")
        return []

    print(f"Matrix configurations to process: {matrix_configs}")

    # Create output directories now
    print("\nCreating output directories...")
    for matrix_type in matrix_configs.keys():
        output_dir = os.path.join(output_base_dir, matrix_type)
        os.makedirs(output_dir, exist_ok=True)
        print(f"  Created: {output_dir}")

    # Find all directories with "_causal" in the name
    print(f"\nSearching for causal directories in: {gwas_dir}")
    causal_dirs = glob.glob(os.path.join(gwas_dir, "*_causal*"))
    print(f"Found {len(causal_dirs)} causal directories")

    # Determine number of processes
    if n_processes is None:
        n_processes = cpu_count()
    print(f"Using {n_processes} processes for parallel processing")

    # Prepare arguments for parallel processing
    process_args = []
    for causal_dir in causal_dirs:
        args_tuple = (causal_dir, ld_base_dir, matrix_configs, output_base_dir, wrapper_script)
        process_args.append(args_tuple)

    # Process directories in parallel
    print(f"\nProcessing {len(process_args)} directories in parallel...")

    all_commands = []
    all_warnings = []

    with Pool(processes=n_processes) as pool:
        results = pool.map(process_causal_dir, process_args)

    # Collect results
    processed_phenotypes = set()
    for commands, warnings in results:
        all_commands.extend(commands)
        all_warnings.extend(warnings)

        # Track processed phenotypes for summary
        for cmd_info in commands:
            if len(cmd_info) >= 4:  # Check if we have the expected tuple structure
                processed_phenotypes.add(cmd_info[3])  # pheno_causal is at index 3

    # Print warnings
    if all_warnings:
        print(f"\nWarnings encountered ({len(all_warnings)}):")
        for warning in all_warnings:
            print(f"  {warning}")

    # Print processing summary
    print(f"\nProcessing Summary:")
    for cmd_info in all_commands[:5]:  # Show first 5 as examples
        if len(cmd_info) >= 7:
            _, _, matrix_type, pheno_causal, trait_type, chromosome, matrix_file = cmd_info
            print(f"  {pheno_causal} ({trait_type}, {chromosome}) -> {matrix_type}")
            print(f"    Matrix: {matrix_file}")

    if len(all_commands) > 5:
        print(f"    ... and {len(all_commands) - 5} more commands")

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

    parser = argparse.ArgumentParser(description='Generate SuSiE fine-mapping commands (optimized for HPC)')

    # Define default paths
    default_wrapper = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/finemap_SuSiE.R"
    default_gwas = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"
    default_ld = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
    default_output = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results/susie_results"
    default_commands = f"{default_output}/susie_commands.txt"

    parser.add_argument('--wrapper-script', '-w', default=default_wrapper,
                        help=f'Path to finemap_SuSiE.R script (default: {default_wrapper})')
    parser.add_argument('--gwas-dir', '-g', default=default_gwas,
                        help=f'Directory containing GWAS summary statistics (default: {default_gwas})')
    parser.add_argument('--ld-dir', '-l', default=default_ld,
                        help=f'Base directory containing LD matrices (default: {default_ld})')
    parser.add_argument('--output-dir', '-o', default=default_output,
                        help=f'Base directory for fine-mapping results (default: {default_output})')
    parser.add_argument('--commands-file', '-c', default=default_commands,
                        help=f'Output file for commands (default: {default_commands})')
    parser.add_argument('--matrix-types', '-m', type=str, default=None,
                        help='Comma-delimited list of matrix subdirectory names to process (e.g., "raw_in-sample,PC-corrected"). If not specified, all detected matrix types will be processed.')
    parser.add_argument('--processes', '-p', type=int, default=None,
                        help='Number of processes to use for parallel processing (default: use all available cores)')
    parser.add_argument('--dry-run', '-d', action='store_true',
                        help='Print first 10 commands without writing to file')
    parser.add_argument('--list-matrices', action='store_true',
                        help='List available matrix subdirectories and exit')

    args = parser.parse_args()

    # Handle list-matrices option
    if args.list_matrices:
        print(f"Scanning for matrix subdirectories in: {args.ld_dir}")
        if not os.path.exists(args.ld_dir):
            print(f"Error: LD directory does not exist: {args.ld_dir}")
            return 1

        matrix_configs = detect_matrix_subdirs(args.ld_dir)
        if matrix_configs:
            print("\nAvailable matrix types:")
            for matrix_type, subfolder in matrix_configs.items():
                print(f"  {matrix_type} ({subfolder})")

                # Show sample matrix files for this type
                search_dir = os.path.join(args.ld_dir, matrix_type, subfolder)
                if os.path.exists(search_dir):
                    sample_files = glob.glob(os.path.join(search_dir, "*.vcor1"))[:3]
                    if sample_files:
                        print(f"    Sample files:")
                        for f in sample_files:
                            print(f"      {os.path.basename(f)}")
                        if len(glob.glob(os.path.join(search_dir, "*.vcor1"))) > 3:
                            print(f"      ... and more")

            print("\nExample usage:")
            print(f"  --matrix-types \"raw_in-sample,PC-corrected\"")
            print(f"  --matrix-types \"{','.join(matrix_configs.keys())}\"")
        else:
            print("No valid matrix configurations found.")
        return 0

    # Check if directories exist
    for dir_path, name in [(args.gwas_dir, 'GWAS'), (args.ld_dir, 'LD matrices')]:
        if not os.path.exists(dir_path):
            print(f"Error: {name} directory does not exist: {dir_path}")
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
    commands = create_susie_commands(
        gwas_dir=args.gwas_dir,
        ld_base_dir=args.ld_dir,
        output_base_dir=args.output_dir,
        wrapper_script=args.wrapper_script,
        output_commands_file=args.commands_file if not args.dry_run else None,
        matrix_filter=matrix_filter,
        n_processes=args.processes
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
