#!/usr/bin/env python3

"""
cleanup_multisusie_junk_files.py

Identify and remove junk files from MultiSuSiE fine-mapping output directories.

This script uses the same logic as generate_multisusie_commands.py to build a
whitelist of REAL output files based on the GWAS summary statistics and
matrix configurations. Any files in the output directory that aren't on
this whitelist are considered junk.

By default, runs in SAFE MODE: only reports junk files without deleting.
Use --delete-mode to actually delete identified junk files.

MultiSuSiE Output Format:
- Output files are located at: {output_base_dir}/{matrix_type}/{pheno_causal}.{trait_type}.{matrix_type}.multisusie.txt
- For multi-population analysis, each combination represents cross-population results

Usage:

python cleanup_multisusie_junk_files.py \
    --gwas-dir /path/to/gwas \
    --pheno-dir /path/to/phenotypes \
    --ld-dir /path/to/ld \
    --output-dir /path/to/output \
    [--populations "AFR,EUR"] \
    [--heritability 0.0003] \
    [--matrix-types "matrix1,matrix2"] \
    [--delete-mode]

"""

import argparse
import glob
import os
import sys
from pathlib import Path


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
    
    if not os.path.exists(ld_base_dir):
        print(f"Warning: LD base directory does not exist: {ld_base_dir}")
        return matrix_types
    
    for item in os.listdir(ld_base_dir):
        matrix_dir = os.path.join(ld_base_dir, item)
        if not os.path.isdir(matrix_dir):
            continue
        
        # Apply matrix filter if specified
        if matrix_filter and item not in matrix_filter:
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
        else:
            # Still include if directory exists but no population subdirs found
            # (in case populations are stored differently)
            matrix_types.append(item)
    
    return matrix_types


def normalize_heritability(heritability_str):
    """
    Normalize heritability string for comparison.
    
    Parameters:
    - heritability_str: Heritability string from filename
    
    Returns:
    - float: Normalized heritability value
    """
    try:
        if heritability_str.startswith("heritability_"):
            heritability_str = heritability_str.replace("heritability_", "")
        return float(heritability_str)
    except (ValueError, TypeError):
        return None


def parse_gwas_filename(meta_filename, target_heritability=None):
    """
    Parse GWAS filename to extract components.
    
    Expected format: Pheno#.#_causal.heritability_#.#.continuous/binary_0.25.chr#.step2.txt
    
    Parameters:
    - meta_filename: The GWAS filename
    - target_heritability: Target heritability value to filter for
    
    Returns:
    - dict: Dictionary with parsed components or None if parsing fails
    """
    # Remove the .step2.txt extension
    base_name = meta_filename.replace('.step2.txt', '')
    
    # Split by dots
    parts = base_name.split('.')
    
    if len(parts) < 5:
        return None
    
    try:
        # First two parts: Pheno# and #_causal
        pheno_causal = parts[0] + "." + parts[1]
        
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
        
        # Reconstruct heritability
        heritability = ".".join(heritability_parts)
        
        # Filter by heritability if specified
        if target_heritability is not None:
            parsed_heritability = normalize_heritability(heritability)
            if parsed_heritability is None:
                return None
            
            # Compare with tolerance for floating point
            if abs(parsed_heritability - target_heritability) > 1e-10:
                return None
        
        # Handle trait type
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


def build_valid_files_whitelist(gwas_base_dir, pheno_base_dir, matrix_types, output_base_dir, populations, target_heritability=None):
    """
    Build a set of all REAL output file paths based on GWAS files and matrix configurations.
    
    For each GWAS file and each matrix configuration, the real output file is:
    - {output_base_dir}/{matrix_type}/{pheno_causal}.{trait_type}.{matrix_type}.multisusie.txt
    
    Parameters:
    - gwas_base_dir: Base directory containing GWAS summary statistics (with population subdirs)
    - pheno_base_dir: Base directory containing phenotype files (with population subdirs)
    - matrix_types: List of matrix type names
    - output_base_dir: Base output directory
    - populations: List of population names
    - target_heritability: Target heritability value to filter for
    
    Returns:
    - set: Valid output file paths (absolute paths)
    """
    valid_files = set()
    
    # Check if GWAS base directory exists
    if not os.path.exists(gwas_base_dir):
        print(f"Warning: GWAS base directory does not exist: {gwas_base_dir}")
        return valid_files
    
    # Check if phenotype base directory exists
    if not os.path.exists(pheno_base_dir):
        print(f"Warning: Phenotype base directory does not exist: {pheno_base_dir}")
        return valid_files
    
    # Use first population as reference for finding unique phenotype patterns
    first_pop = populations[0] if populations else None
    
    if not first_pop:
        print(f"Warning: No populations specified")
        return valid_files
    
    pop_gwas_dir = os.path.join(gwas_base_dir, first_pop)
    pop_pheno_dir = os.path.join(pheno_base_dir, first_pop)
    
    if not os.path.exists(pop_gwas_dir):
        print(f"Warning: Population GWAS directory does not exist: {pop_gwas_dir}")
        return valid_files
    
    if not os.path.exists(pop_pheno_dir):
        print(f"Warning: Population phenotype directory does not exist: {pop_pheno_dir}")
        return valid_files
    
    # Find all causal GWAS files from the first population
    causal_files = glob.glob(os.path.join(pop_gwas_dir, "*_causal*.step2.txt"))
    
    # Filter out unwanted patterns
    causal_files = [f for f in causal_files
                    if 'drop_causal' not in os.path.basename(f)
                    and 'matched_to_' not in os.path.basename(f)]
    
    # Parse each GWAS file
    for gwas_file in causal_files:
        gwas_filename = os.path.basename(gwas_file)
        
        # Parse the filename
        parsed = parse_gwas_filename(gwas_filename, target_heritability)
        
        if parsed is None:
            continue
        
        pheno_causal = parsed['pheno_causal']
        trait_type = parsed['trait_type']
        heritability = parsed['heritability']
        
        # Verify that phenotype file exists for this pattern
        # (MultiSuSiE requires phenotype files, so we skip if missing)
        pheno_filename = f"{pheno_causal}.{heritability}.{trait_type}.txt"
        pheno_path = os.path.join(pop_pheno_dir, pheno_filename)
        
        if not os.path.exists(pheno_path):
            continue
        
        # For each matrix configuration
        for matrix_type in matrix_types:
            # Output filename format: Pheno#.#_causal.continuous/binary_0.25.matrix_type.multisusie.txt
            output_filename = f"{pheno_causal}.{trait_type}.{matrix_type}.multisusie.txt"
            output_path = os.path.join(output_base_dir, matrix_type, output_filename)
            
            # Add to whitelist
            valid_files.add(output_path)
    
    return valid_files


def find_all_junk_files_recursive(output_base_dir, valid_files_whitelist):
    """
    Recursively find all junk files in output directory structure.
    
    Parameters:
    - output_base_dir: Base output directory
    - valid_files_whitelist: Set of valid file paths to keep
    
    Returns:
    - list: All junk file paths found
    """
    all_junk_files = []
    
    for root, dirs, files in os.walk(output_base_dir):
        for filename in sorted(files):
            filepath = os.path.join(root, filename)
            
            if filepath not in valid_files_whitelist:
                all_junk_files.append(filepath)
    
    return sorted(all_junk_files)


def format_file_size(bytes_size):
    """
    Format bytes to human-readable size string.
    
    Parameters:
    - bytes_size: Size in bytes
    
    Returns:
    - str: Formatted string (e.g., "1.5 MB")
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_size < 1024.0:
            return f"{bytes_size:.1f} {unit}"
        bytes_size /= 1024.0
    return f"{bytes_size:.1f} PB"


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
    
    matrix_list = [name.strip() for name in matrix_string.split(',')]
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


def main():
    parser = argparse.ArgumentParser(
        description="Identify and remove junk files from MultiSuSiE fine-mapping output directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Safe mode: just list junk files (DEFAULT)
  python cleanup_multisusie_junk_files.py --gwas-dir /gwas --pheno-dir /pheno --ld-dir /ld --output-dir /output
  
  # Actually delete junk files (requires --delete-mode)
  python cleanup_multisusie_junk_files.py --gwas-dir /gwas --pheno-dir /pheno --ld-dir /ld --output-dir /output --delete-mode
  
  # With custom populations and heritability filter
  python cleanup_multisusie_junk_files.py --gwas-dir /gwas --pheno-dir /pheno --ld-dir /ld --output-dir /output \\
    --populations "AFR,EUR,EAS" --heritability 0.001
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--gwas-dir", 
        type=str, 
        required=True,
        help="Base directory containing GWAS summary statistics with population subdirectories"
    )
    
    parser.add_argument(
        "--pheno-dir", 
        type=str, 
        required=True,
        help="Base directory containing phenotype files with population subdirectories"
    )
    
    parser.add_argument(
        "--ld-dir", 
        type=str, 
        required=True,
        help="Base directory containing LD matrices"
    )
    
    parser.add_argument(
        "--output-dir", 
        type=str, 
        required=True,
        help="Base output directory to scan for junk files"
    )
    
    # Optional arguments
    parser.add_argument(
        "--populations",
        type=str,
        default="AFR,EUR",
        help="Comma-delimited list of population names (default: AFR,EUR)"
    )
    
    parser.add_argument(
        "--heritability",
        type=float,
        default=0.0003,
        help="Target heritability value to filter traits [default: 0.0003]"
    )
    
    parser.add_argument(
        "--matrix-types",
        type=str,
        default=None,
        help="Comma-delimited list of matrix subdirectory names to process (e.g., 'raw_in-sample,PC-corrected')"
    )
    
    parser.add_argument(
        "--delete-mode",
        action="store_true",
        default=False,
        help="ENABLE deletion of junk files. By default, only lists files without deleting (SAFE MODE)"
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    for dir_path, name in [(args.gwas_dir, 'GWAS'), (args.pheno_dir, 'Phenotype'), (args.ld_dir, 'LD'), (args.output_dir, 'Output')]:
        if not os.path.exists(dir_path):
            print(f"Error: {name} directory does not exist: {dir_path}", file=sys.stderr)
            sys.exit(1)
    
    # Validate heritability
    if args.heritability < 0:
        print(f"Error: Heritability must be non-negative, got {args.heritability}", file=sys.stderr)
        sys.exit(1)
    
    # Parse populations
    populations = parse_populations(args.populations)
    
    print(f"MultiSuSiE Fine-Mapping Junk File Cleanup")
    print(f"{'='*60}")
    print(f"GWAS directory: {args.gwas_dir}")
    print(f"Phenotype directory: {args.pheno_dir}")
    print(f"LD directory: {args.ld_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Populations: {', '.join(populations)}")
    print(f"Heritability filter: {args.heritability}")
    
    if args.delete_mode:
        print(f"Mode: DELETE MODE (files WILL be deleted)")
    else:
        print(f"Mode: SAFE MODE (no files will be deleted)")
    
    print(f"{'='*60}\n")
    
    # Parse matrix filter
    matrix_filter = parse_matrix_filter(args.matrix_types)
    
    # Detect matrix configurations
    print("Detecting matrix subdirectories...")
    matrix_types = detect_matrix_subdirs(args.ld_dir, matrix_filter)
    
    if not matrix_types:
        print("Error: No valid matrix configurations found!", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(matrix_types)} matrix type(s): {', '.join(matrix_types)}\n")
    
    # Build whitelist
    print("Building whitelist of valid output files...")
    valid_files = build_valid_files_whitelist(
        args.gwas_dir,
        args.pheno_dir,
        matrix_types,
        args.output_dir,
        populations,
        args.heritability
    )
    
    print(f"Whitelist contains {len(valid_files)} valid output file(s)\n")
    
    # Find junk files
    print("Scanning output directory for junk files...")
    all_junk_files = find_all_junk_files_recursive(args.output_dir, valid_files)
    
    if not all_junk_files:
        print("✓ No junk files found - output directory is clean!")
        sys.exit(0)
    
    print(f"\nFound {len(all_junk_files)} junk file(s):\n")
    
    total_size = 0
    for filepath in all_junk_files:
        try:
            file_size = os.path.getsize(filepath)
            total_size += file_size
            size_str = format_file_size(file_size)
            print(f"  {filepath} ({size_str})")
        except Exception as e:
            print(f"  {filepath} (error reading size: {e})")
    
    print(f"\nTotal size of junk files: {format_file_size(total_size)}")
    
    if not args.delete_mode:
        print(f"\n[SAFE MODE] Would delete {len(all_junk_files)} file(s)")
        print("Run with --delete-mode to actually delete these files")
        sys.exit(0)
    
    # Confirm deletion
    print(f"\n⚠️  DANGER: About to delete {len(all_junk_files)} file(s)!")
    response = input("Are you sure? Type 'yes' to confirm: ")
    
    if response.lower() != 'yes':
        print("Deletion cancelled - no files were deleted")
        sys.exit(0)
    
    # Delete files
    print(f"\nDeleting junk files...\n")
    
    deleted_count = 0
    failed_count = 0
    
    for filepath in all_junk_files:
        try:
            os.remove(filepath)
            print(f"✓ Deleted: {filepath}")
            deleted_count += 1
        except Exception as e:
            print(f"✗ Failed to delete {filepath}: {e}", file=sys.stderr)
            failed_count += 1
    
    print(f"\n{'='*60}")
    print(f"Deletion complete:")
    print(f"  Successfully deleted: {deleted_count} file(s)")
    if failed_count > 0:
        print(f"  Failed to delete: {failed_count} file(s)")
    print(f"  Freed space: {format_file_size(total_size)}")
    print(f"{'='*60}")
    
    if failed_count > 0:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
