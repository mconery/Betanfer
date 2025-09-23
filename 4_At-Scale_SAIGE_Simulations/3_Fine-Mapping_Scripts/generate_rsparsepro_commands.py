#!/usr/bin/env python3

"""
generate_rsparsepro_commands.py
Generate commands for running RSparsePro fine-mapping on multiple loci

This script creates a command file where each line contains a command to run
RSparsePro fine-mapping on a single locus using the wrapper script.

Optimized for parallel processing on HPC systems with file existence checks.
"""

import os
import glob
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial

def parse_gwas_filename(meta_filename):
    """
    Parse GWAS filename to extract components, handling decimal heritability values.
    Expected format: Pheno#.#_causal.heritability_#.#.continuous/binary_0.25.chr#.meta.out
    
    Parameters:
    - meta_filename: The GWAS filename
    
    Returns:
    - dict: Dictionary with parsed components or None if parsing fails
    """
    # Remove the .meta.out extension
    base_name = meta_filename.replace('.meta.out', '')
    
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

def find_ld_file(ld_dir, pheno_causal):
    """
    Search for an LD file that matches the given phenotype/causal pattern.
    
    Parameters:
    - ld_dir: Directory containing LD files
    - pheno_causal: Phenotype pattern (e.g., "Pheno1.2_causal")
    
    Returns:
    - str: Full path to matching LD file, or None if not found
    """
    if not os.path.exists(ld_dir):
        return None
    
    # Search for files that start with the pheno_causal pattern and end with .vcor1
    pattern = os.path.join(ld_dir, f"{pheno_causal}.*.vcor1")
    matching_files = glob.glob(pattern)
    
    if len(matching_files) == 1:
        # Check if corresponding .vars file exists
        vars_file = matching_files[0] + ".vars"
        if os.path.exists(vars_file):
            return matching_files[0]
        else:
            return None
    elif len(matching_files) > 1:
        # For parallel processing, use first match that has a .vars file
        for ld_file in matching_files:
            vars_file = ld_file + ".vars"
            if os.path.exists(vars_file):
                return ld_file
        return None
    else:
        # No matches found
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
        
        # For RSparsePro output files, check if they have expected content
        with open(output_path, 'r') as f:
            lines = f.readlines()
        
        # Should have at least 2 lines (header + 1 data row)
        if len(lines) < 2:
            return False
        
        # Check if header contains expected columns for RSparsePro
        header = lines[0].strip().split('\t')
        expected_cols = {'PIP', 'beta', 'z_estimated'}  # Common RSparsePro output columns
        
        # Check if at least one expected column is present
        if not any(col in header for col in expected_cols):
            return False
        
        return True
    
    except (IOError, OSError, UnicodeDecodeError):
        # If we can't read the file, consider it incomplete
        return False

def process_causal_dir(args_tuple):
    """
    Process a single causal directory to generate commands.
    This function is designed to be called in parallel.
    
    Parameters:
    - args_tuple: Tuple containing (causal_dir, ld_dir, output_dir, wrapper_script, varcol, beta_col, se_col, skip_existing)
    
    Returns:
    - tuple: (commands, warnings, skipped_info)
    """
    causal_dir, ld_dir, output_dir, wrapper_script, varcol, beta_col, se_col, skip_existing = args_tuple
    
    commands = []
    warnings = []
    skipped_info = []
    
    if not os.path.isdir(causal_dir):
        return commands, warnings, skipped_info
    
    # Find the .meta.out file in this directory
    meta_files = glob.glob(os.path.join(causal_dir, "*.meta.out"))
    
    if len(meta_files) != 1:
        warning = f"Expected 1 .meta.out file in {causal_dir}, found {len(meta_files)}"
        if len(meta_files) > 1:
            warning += f" - Files: {[os.path.basename(f) for f in meta_files]}"
        warnings.append(warning)
        return commands, warnings, skipped_info
    
    meta_file = meta_files[0]
    meta_filename = os.path.basename(meta_file)
    
    # Parse the filename
    parsed = parse_gwas_filename(meta_filename)
    if parsed is None:
        warnings.append(f"Could not parse filename: {meta_filename}")
        return commands, warnings, skipped_info
    
    pheno_causal = parsed['pheno_causal']
    heritability = parsed['heritability']
    trait_type = parsed['trait_type']
    chromosome = parsed['chromosome']
    
    # Search for the LD file that matches this phenotype
    ld_file = find_ld_file(ld_dir, pheno_causal)
    
    if ld_file is None:
        warnings.append(f"No LD file (.vcor1 + .vars) found for {pheno_causal} in {ld_dir}")
        return commands, warnings, skipped_info
    
    # Construct output file path
    output_filename = f"{pheno_causal}.{trait_type}.rsparsepro.txt"
    output_path = os.path.join(output_dir, output_filename)
    
    # Check if output file already exists and is complete
    if skip_existing and check_output_completeness(output_path):
        skipped_info.append({
            'pheno_causal': pheno_causal,
            'trait_type': trait_type,
            'chromosome': chromosome,
            'heritability': heritability,
            'output_path': output_path,
            'ld_file': os.path.basename(ld_file)
        })
        return commands, warnings, skipped_info
    
    # Create the RSparsePro command
    # Usage: ./finemap_RSparsePro.sh [gwas_file] [ld_file] [output_file] [varcol] [beta_col] [se_col]
    rsparsepro_cmd = f"bash {wrapper_script} {meta_file} {ld_file} {output_path} {varcol} {beta_col} {se_col}"
    
    commands.append((rsparsepro_cmd, output_dir, pheno_causal, trait_type, chromosome, heritability, os.path.basename(ld_file), output_filename))
    
    return commands, warnings, skipped_info

def create_rsparsepro_commands(gwas_dir, ld_dir, output_dir, wrapper_script, output_commands_file,
                              varcol='rs_number', beta_col='beta', se_col='se', 
                              n_processes=None, skip_existing=False):
    """
    Generate RSparsePro fine-mapping commands for all loci using parallel processing.
    
    Parameters:
    - gwas_dir: Directory containing GWAS summary statistics
    - ld_dir: Directory containing LD matrices
    - output_dir: Directory for fine-mapping results
    - wrapper_script: Path to finemap_RSparsePro.sh script
    - output_commands_file: File to write commands to
    - varcol: Column name for variant IDs
    - beta_col: Column name for beta coefficients
    - se_col: Column name for standard errors
    - n_processes: Number of processes to use (None = use all available cores)
    - skip_existing: Whether to skip generating commands for existing complete output files
    """
    
    print(f"GWAS directory: {gwas_dir}")
    print(f"LD directory: {ld_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Column settings: varcol={varcol}, beta_col={beta_col}, se_col={se_col}")
    
    if skip_existing:
        print(f"File existence checking: ENABLED - will skip existing complete output files")
    else:
        print(f"File existence checking: DISABLED - will generate all commands")
    
    # Create output directory
    print(f"\nCreating output directories...")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "timing"), exist_ok=True)
    print(f" Created: {output_dir}")
    
    # Find all directories with "_causal" in the name
    print(f"\nSearching for causal directories in: {gwas_dir}")
    causal_dirs = glob.glob(os.path.join(gwas_dir, "*_causal*"))
    print(f"Found {len(causal_dirs)} causal directories")
    
    if not causal_dirs:
        print("No causal directories found!")
        return []
    
    # Determine number of processes
    if n_processes is None:
        n_processes = cpu_count()
    print(f"Using {n_processes} processes for parallel processing")
    
    # Prepare arguments for parallel processing
    process_args = []
    for causal_dir in causal_dirs:
        args_tuple = (causal_dir, ld_dir, output_dir, wrapper_script, varcol, beta_col, se_col, skip_existing)
        process_args.append(args_tuple)
    
    # Process directories in parallel
    print(f"\nProcessing {len(process_args)} directories in parallel...")
    
    all_commands = []
    all_warnings = []
    all_skipped = []
    
    with Pool(processes=n_processes) as pool:
        results = pool.map(process_causal_dir, process_args)
    
    # Collect results
    processed_phenotypes = set()
    for commands, warnings, skipped_info in results:
        all_commands.extend(commands)
        all_warnings.extend(warnings)
        all_skipped.extend(skipped_info)
        
        # Track processed phenotypes for summary
        for cmd_info in commands:
            if len(cmd_info) >= 3:  # Check if we have the expected tuple structure
                processed_phenotypes.add(cmd_info[2])  # pheno_causal is at index 2
    
    # Print warnings
    if all_warnings:
        print(f"\nWarnings encountered ({len(all_warnings)}):")
        for warning in all_warnings:
            print(f" {warning}")
    
    # Print skipped files summary
    if all_skipped:
        print(f"\nExisting files skipped ({len(all_skipped)}):")
        print(f" Total: {len(all_skipped)} files")
        
        # Show a few examples
        print(f" Examples:")
        for i, skipped in enumerate(all_skipped[:3]):
            print(f" {skipped['pheno_causal']}.{skipped['trait_type']} (h²={skipped['heritability']}) -> (already complete)")
        
        if len(all_skipped) > 3:
            print(f" ... and {len(all_skipped) - 3} more")
    
    # Print processing summary
    print(f"\nProcessing Summary:")
    print(f" Total commands generated: {len(all_commands)}")
    if skip_existing:
        print(f" Existing files skipped: {len(all_skipped)}")
        print(f" New commands needed: {len(all_commands)}")
    
    # Show examples of commands
    for cmd_info in all_commands[:5]:  # Show first 5 as examples
        if len(cmd_info) >= 8:
            _, _, pheno_causal, trait_type, chromosome, heritability, ld_file, output_filename = cmd_info
            print(f" {pheno_causal} ({trait_type}, {chromosome}, h²={heritability})")
            print(f" LD: {ld_file}")
            print(f" Output: {output_filename}")
    
    if len(all_commands) > 5:
        print(f" ... and {len(all_commands) - 5} more commands")
    
    # Extract just the commands for writing to file
    commands_only = [cmd_info[0] for cmd_info in all_commands]
    
    # Write commands to file
    if output_commands_file:
        print(f"\nWriting {len(commands_only)} commands to: {output_commands_file}")
        with open(output_commands_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# RSparsePro fine-mapping commands\n")
            f.write(f"# Generated by {os.path.basename(__file__)}\n\n")
            for cmd in commands_only:
                f.write(cmd + '\n')
        
        # Make command file executable
        os.chmod(output_commands_file, 0o755)
    
    return commands_only

def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(description='Generate RSparsePro fine-mapping commands (optimized for HPC)')
    
    # Define default paths - update these to point to RSparsePro scripts
    default_wrapper = "./finemap_RSparsePro.sh"
    default_gwas = "/path/to/gwas/results"
    default_ld = "/path/to/ld/matrices"
    default_output = "/path/to/output/rsparsepro_results"
    default_commands = "rsparsepro_commands.sh"
    
    parser.add_argument('--wrapper-script', '-w', default=default_wrapper,
                       help=f'Path to finemap_RSparsePro.sh script (default: {default_wrapper})')
    
    parser.add_argument('--gwas-dir', '-g', required=True,
                       help='Directory containing GWAS summary statistics')
    
    parser.add_argument('--ld-dir', '-l', required=True,
                       help='Directory containing LD matrices (.vcor1 files)')
    
    parser.add_argument('--output-dir', '-o', required=True,
                       help='Directory for fine-mapping results')
    
    parser.add_argument('--commands-file', '-c', default=default_commands,
                       help=f'Output file for commands (default: {default_commands})')
    
    parser.add_argument('--varcol', default='rs_number',
                       help='Column name for variant IDs (default: rs_number)')
    
    parser.add_argument('--beta-col', default='beta',
                       help='Column name for beta coefficients (default: beta)')
    
    parser.add_argument('--se-col', default='se',
                       help='Column name for standard errors (default: se)')
    
    parser.add_argument('--processes', '-p', type=int, default=None,
                       help='Number of processes to use for parallel processing (default: use all available cores)')
    
    parser.add_argument('--skip-existing', '-s', action='store_true',
                       help='Skip generating commands for output files that already exist and appear complete')
    
    parser.add_argument('--force', '-f', action='store_true',
                       help='Force regeneration of all commands, even if output files exist')
    
    parser.add_argument('--dry-run', '-d', action='store_true',
                       help='Print first 10 commands without writing to file')
    
    args = parser.parse_args()
    
    # Validate mutually exclusive options
    if args.skip_existing and args.force:
        print("Error: --skip-existing and --force options are mutually exclusive")
        return 1
    
    # Check if directories exist
    for dir_path, name in [(args.gwas_dir, 'GWAS'), (args.ld_dir, 'LD matrices')]:
        if not os.path.exists(dir_path):
            print(f"Error: {name} directory does not exist: {dir_path}")
            return 1
    
    if not os.path.exists(args.wrapper_script):
        print(f"Error: Wrapper script does not exist: {args.wrapper_script}")
        return 1
    
    # Validate number of processes
    if args.processes is not None and args.processes < 1:
        print(f"Error: Number of processes must be >= 1, got {args.processes}")
        return 1
    
    # Set skip_existing based on force flag
    skip_existing = args.skip_existing and not args.force
    
    # Generate commands
    commands = create_rsparsepro_commands(
        gwas_dir=args.gwas_dir,
        ld_dir=args.ld_dir,
        output_dir=args.output_dir,
        wrapper_script=args.wrapper_script,
        output_commands_file=args.commands_file if not args.dry_run else None,
        varcol=args.varcol,
        beta_col=args.beta_col,
        se_col=args.se_col,
        n_processes=args.processes,
        skip_existing=skip_existing
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
        print(f" bash {args.commands_file}")
        print("\nOr submit as array job to SLURM")
    
    return 0

if __name__ == "__main__":
    exit(main())