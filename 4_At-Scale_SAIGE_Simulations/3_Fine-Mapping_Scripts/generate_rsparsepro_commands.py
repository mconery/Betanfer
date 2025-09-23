#!/usr/bin/env python3

"""
generate_rsparsepro_commands.py
Generate commands for running RSparsePro fine-mapping on multiple loci

This script creates a command file where each line contains a command to run
RSparsePro fine-mapping on a single locus using the wrapper script.

Usage:
    python generate_rsparsepro_commands.py --input_dir <input_directory> --output_dir <output_directory> --command_file <command_file> [options]

Required arguments:
    --input_dir: Directory containing input GWAS and LD files
    --output_dir: Directory where RSparsePro results will be saved
    --command_file: Output file containing commands for each locus

Optional arguments:
    --gwas_file: Single GWAS file to use for all loci (default: look for individual files)
    --ld_suffix: Suffix for LD files (default: .vcor1)
    --varcol: Column name for variant IDs (default: rs_number)
    --beta_col: Column name for beta coefficients (default: beta)
    --se_col: Column name for standard errors (default: se)
    --wrapper_script: Path to RSparsePro wrapper script (default: ./finemap_RSparsePro.sh)
    --skip_existing: Skip loci where output files already exist (default: False)
    --force: Force regeneration of all commands, even if output files exist (default: False)
    --help: Show this help message
"""

import os
import argparse
import glob
import sys
from pathlib import Path

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Generate commands for RSparsePro fine-mapping",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        '--input_dir', 
        type=str, 
        required=True,
        help='Directory containing input GWAS and LD files'
    )
    
    parser.add_argument(
        '--output_dir', 
        type=str, 
        required=True,
        help='Output directory for RSparsePro results'
    )
    
    parser.add_argument(
        '--command_file', 
        type=str, 
        required=True,
        help='Output file to write commands'
    )
    
    parser.add_argument(
        '--gwas_file', 
        type=str, 
        default=None,
        help='Single GWAS file to use for all loci (if not provided, looks for individual files)'
    )
    
    parser.add_argument(
        '--ld_suffix', 
        type=str, 
        default='.vcor1',
        help='Suffix for LD files (default: .vcor1)'
    )
    
    parser.add_argument(
        '--varcol', 
        type=str, 
        default='rs_number',
        help='Column name for variant IDs (default: rs_number)'
    )
    
    parser.add_argument(
        '--beta_col', 
        type=str, 
        default='beta',
        help='Column name for beta coefficients (default: beta)'
    )
    
    parser.add_argument(
        '--se_col', 
        type=str, 
        default='se',
        help='Column name for standard errors (default: se)'
    )
    
    parser.add_argument(
        '--wrapper_script', 
        type=str, 
        default='./finemap_RSparsePro.sh',
        help='Path to RSparsePro wrapper script (default: ./finemap_RSparsePro.sh)'
    )
    
    parser.add_argument(
        '--skip_existing', 
        action='store_true',
        help='Skip loci where output files already exist (default: False)'
    )
    
    parser.add_argument(
        '--force', 
        action='store_true',
        help='Force regeneration of all commands, even if output files exist (default: False)'
    )
    
    return parser.parse_args()

def find_locus_files(input_dir, ld_suffix, gwas_file=None):
    """
    Find matching GWAS and LD files for each locus
    
    Returns:
        list: List of tuples (locus_name, gwas_file, ld_file)
    """
    input_path = Path(input_dir)
    
    # Find all LD files
    ld_files = list(input_path.glob(f'*{ld_suffix}'))
    
    locus_files = []
    
    for ld_file in ld_files:
        # Extract locus name by removing suffix
        locus_name = ld_file.stem.replace(ld_suffix.replace('.vcor1', ''), '')
        
        # Check if corresponding .vars file exists
        vars_file = input_path / f"{ld_file.name}.vars"
        if not vars_file.exists():
            print(f"Warning: Variants file not found for locus {locus_name}: {vars_file}")
            continue
        
        # Determine GWAS file to use
        if gwas_file:
            # Use single GWAS file for all loci
            gwas_path = gwas_file
        else:
            # Look for locus-specific GWAS file
            # Try common GWAS file patterns
            gwas_patterns = [
                f"{locus_name}.gwas.txt",
                f"{locus_name}.gwas.tsv",
                f"{locus_name}.tsv",
                f"{locus_name}.txt",
                f"{locus_name}_gwas.txt",
                f"{locus_name}_gwas.tsv"
            ]
            
            gwas_path = None
            for pattern in gwas_patterns:
                potential_gwas = input_path / pattern
                if potential_gwas.exists():
                    gwas_path = str(potential_gwas)
                    break
            
            if not gwas_path:
                print(f"Warning: GWAS file not found for locus {locus_name}. Tried patterns: {gwas_patterns}")
                continue
        
        # Validate GWAS file exists
        if not os.path.exists(gwas_path):
            print(f"Warning: GWAS file does not exist: {gwas_path}")
            continue
            
        locus_files.append((locus_name, gwas_path, str(ld_file)))
    
    return locus_files

def check_output_exists(output_file):
    """
    Check if output file exists and appears to be a valid RSparsePro result
    
    Returns:
        bool: True if output exists and appears valid, False otherwise
    """
    if not os.path.exists(output_file):
        return False
    
    # Check if file is not empty
    if os.path.getsize(output_file) == 0:
        return False
    
    # Check if file contains expected RSparsePro output format
    try:
        with open(output_file, 'r') as f:
            header = f.readline().strip()
            # Look for typical RSparsePro output columns
            if any(col in header.lower() for col in ['pip', 'beta', 'z_estimated', 'cs']):
                return True
    except (IOError, UnicodeDecodeError):
        return False
    
    return False

def generate_commands(locus_files, output_dir, wrapper_script, varcol, beta_col, se_col, 
                     skip_existing=False, force=False):
    """
    Generate RSparsePro commands for each locus
    
    Returns:
        tuple: (list of command strings, list of skipped loci)
    """
    commands = []
    skipped_loci = []
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    for locus_name, gwas_file, ld_file in locus_files:
        # Create output file path
        output_file = os.path.join(output_dir, f"{locus_name}.rsparsepro.txt")
        
        # Check if output already exists
        if not force and skip_existing and check_output_exists(output_file):
            print(f"Skipping {locus_name}: output file already exists at {output_file}")
            skipped_loci.append(locus_name)
            continue
        elif not force and not skip_existing and check_output_exists(output_file):
            print(f"Warning: Output file already exists for {locus_name}: {output_file}")
            print("  Use --skip_existing to skip existing files, or --force to overwrite")
        
        # Create command with all arguments
        command = f"{wrapper_script} {gwas_file} {ld_file} {output_file} {varcol} {beta_col} {se_col}"
        commands.append(command)
    
    return commands, skipped_loci

def write_command_file(commands, command_file):
    """Write commands to file"""
    try:
        with open(command_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# RSparsePro fine-mapping commands\n")
            f.write(f"# Generated by {os.path.basename(__file__)}\n\n")
            
            for command in commands:
                f.write(command + '\n')
        
        # Make command file executable
        os.chmod(command_file, 0o755)
        
        print(f"Successfully wrote {len(commands)} commands to {command_file}")
        
    except IOError as e:
        print(f"Error writing command file: {e}")
        sys.exit(1)

def main():
    """Main function"""
    args = parse_arguments()
    
    # Validate mutually exclusive options
    if args.skip_existing and args.force:
        print("Error: --skip_existing and --force options are mutually exclusive")
        sys.exit(1)
    
    # Validate input directory
    if not os.path.isdir(args.input_dir):
        print(f"Error: Input directory does not exist: {args.input_dir}")
        sys.exit(1)
    
    # Validate wrapper script
    if not os.path.isfile(args.wrapper_script):
        print(f"Error: Wrapper script not found: {args.wrapper_script}")
        sys.exit(1)
    
    # Validate single GWAS file if provided
    if args.gwas_file and not os.path.isfile(args.gwas_file):
        print(f"Error: GWAS file not found: {args.gwas_file}")
        sys.exit(1)
    
    print(f"Searching for locus files in: {args.input_dir}")
    print(f"LD suffix: {args.ld_suffix}")
    if args.gwas_file:
        print(f"Using single GWAS file: {args.gwas_file}")
    else:
        print("Looking for locus-specific GWAS files")
    
    if args.skip_existing:
        print("Skipping loci with existing output files")
    elif args.force:
        print("Force mode: regenerating all commands regardless of existing files")
    
    # Find locus files
    locus_files = find_locus_files(args.input_dir, args.ld_suffix, args.gwas_file)
    
    if not locus_files:
        print("No matching locus files found!")
        print(f"Make sure files with suffix '{args.ld_suffix}' and corresponding .vars files exist in {args.input_dir}")
        if not args.gwas_file:
            print("Also ensure GWAS files exist for each locus, or provide --gwas_file for a single GWAS file")
        sys.exit(1)
    
    print(f"Found {len(locus_files)} loci:")
    for locus_name, gwas_file, ld_file in locus_files:
        print(f"  {locus_name}: {os.path.basename(gwas_file)}, {os.path.basename(ld_file)}")
    
    # Generate commands
    commands, skipped_loci = generate_commands(locus_files, args.output_dir, args.wrapper_script, 
                                             args.varcol, args.beta_col, args.se_col,
                                             args.skip_existing, args.force)
    
    # Write command file
    write_command_file(commands, args.command_file)
    
    # Summary
    print(f"Output directory: {args.output_dir}")
    print(f"Column settings: varcol={args.varcol}, beta_col={args.beta_col}, se_col={args.se_col}")
    print(f"Commands written to: {args.command_file}")
    
    if skipped_loci:
        print(f"Skipped {len(skipped_loci)} loci with existing output files:")
        for locus in skipped_loci:
            print(f"  {locus}")
    
    if commands:
        print(f"\nTo run {len(commands)} commands:")
        print(f"  bash {args.command_file}")
    else:
        print("\nNo commands to run (all loci already processed)")

if __name__ == "__main__":
    main()