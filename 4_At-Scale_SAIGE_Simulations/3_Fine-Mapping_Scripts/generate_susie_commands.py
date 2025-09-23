#!/usr/bin/env python3
"""
Script to generate SuSiE fine-mapping commands for all loci.

This script creates a command file where each line runs SuSiE fine-mapping
for a different locus across all available matrix types.
"""

import os
import glob
import argparse
from pathlib import Path

def create_susie_commands(gwas_dir, ld_base_dir, output_base_dir, wrapper_script, output_commands_file):
    """
    Generate SuSiE fine-mapping commands for all loci.

    Parameters:
    - gwas_dir: Directory containing GWAS summary statistics
    - ld_base_dir: Base directory containing LD matrices  
    - output_base_dir: Base directory for fine-mapping results
    - wrapper_script: Path to finemap_SuSiE.R script
    - output_commands_file: File to write commands to
    """

    # Define matrix subdirectories and their corresponding sub-folders
    matrix_configs = {
        "raw_in-sample": "weighted",
        "PC-corrected": "weighted", 
        "beta-inferred_2000": "Meta"
    }

    commands = []

    # Find all directories with "_causal" in the name
    print(f"Searching for causal directories in: {gwas_dir}")
    causal_dirs = glob.glob(os.path.join(gwas_dir, "*_causal*"))
    print(f"Found {len(causal_dirs)} causal directories")

    for causal_dir in causal_dirs:
        if not os.path.isdir(causal_dir):
            continue

        # Find the .meta.out file in this directory
        meta_files = glob.glob(os.path.join(causal_dir, "*.meta.out"))

        if len(meta_files) != 1:
            print(f"Warning: Expected 1 .meta.out file in {causal_dir}, found {len(meta_files)}")
            if len(meta_files) > 1:
                print(f"  Files found: {[os.path.basename(f) for f in meta_files]}")
            continue

        meta_file = meta_files[0]
        meta_filename = os.path.basename(meta_file)

        # Parse the filename to extract components
        # Format: Pheno#.#_causal.heritability_#.continuous/binary_0.25.chr#.meta.out
        filename_parts = meta_filename.replace('.meta.out', '').split('.')

        if len(filename_parts) < 4:
            print(f"Warning: Unexpected filename format: {meta_filename}")
            print(f"  Parts: {filename_parts}")
            continue

        try:
            pheno_causal = filename_parts[0] + "." + filename_parts[1]  # e.g., "Pheno1.2_causal"
            heritability = filename_parts[2]  # e.g., "heritability_0.01"
            trait_type = filename_parts[3]    # e.g., "continuous" or "binary_0.25"
            chromosome = filename_parts[4]    # e.g., "chr1"

            print(f"Processing: {pheno_causal} ({trait_type}, {chromosome})")

        except IndexError:
            print(f"Warning: Could not parse filename: {meta_filename}")
            continue

        # For each matrix configuration
        for matrix_type, subfolder in matrix_configs.items():
            # Construct matrix file path
            matrix_filename = f"{pheno_causal}.{matrix_type}.unphased.vcor1"
            matrix_path = os.path.join(ld_base_dir, matrix_type, subfolder, matrix_filename)

            # Check if matrix file exists (optional check)
            if not os.path.exists(matrix_path):
                print(f"Warning: Matrix file not found: {matrix_path}")
                # Still add the command - it might be created later

            # Construct output directory and file
            output_dir = os.path.join(output_base_dir, matrix_type)
            output_filename = f"{pheno_causal}.{subfolder}.txt"
            output_path = os.path.join(output_dir, output_filename)

            # Create output directory command (mkdir -p)
            mkdir_cmd = f"mkdir -p {output_dir}"
            if mkdir_cmd not in commands:  # Avoid duplicates
                commands.append(mkdir_cmd)

            # Create the SuSiE command
            # Usage: Rscript finemap_SuSiE.R summary_stats_file ld_file out_file [optional_params]
            susie_cmd = f"Rscript {wrapper_script} {meta_file} {matrix_path} {output_path}"

            commands.append(susie_cmd)

    # Write commands to file
    print(f"\nWriting {len(commands)} commands to: {output_commands_file}")
    with open(output_commands_file, 'w') as f:
        for cmd in commands:
            f.write(cmd + '\n')

    return commands

def main():
    """Main function with command line argument parsing."""

    parser = argparse.ArgumentParser(description='Generate SuSiE fine-mapping commands')

    # Define default paths
    default_wrapper = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/Meta_Fine-Mapping/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts/finemap_SuSiE.R"
    default_gwas = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results"
    default_ld = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/LD_matrices"
    default_output = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results/susie_results"
    default_commands = "susie_commands.txt"

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
    parser.add_argument('--dry-run', '-d', action='store_true',
                        help='Print first 10 commands without writing to file')

    args = parser.parse_args()

    # Check if directories exist
    for dir_path, name in [(args.gwas_dir, 'GWAS'), (args.ld_dir, 'LD matrices')]:
        if not os.path.exists(dir_path):
            print(f"Error: {name} directory does not exist: {dir_path}")
            return 1

    if not os.path.exists(args.wrapper_script):
        print(f"Error: Wrapper script does not exist: {args.wrapper_script}")
        return 1

    # Generate commands
    commands = create_susie_commands(
        gwas_dir=args.gwas_dir,
        ld_base_dir=args.ld_dir,
        output_base_dir=args.output_dir,
        wrapper_script=args.wrapper_script,
        output_commands_file=args.commands_file if not args.dry_run else None
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
