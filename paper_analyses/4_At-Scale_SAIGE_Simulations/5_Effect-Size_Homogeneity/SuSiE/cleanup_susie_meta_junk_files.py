#!/usr/bin/env python3

"""
cleanup_susie_meta_junk_files.py

Identify and remove junk files from meta-analysis SuSiE fine-mapping output directories.

This script uses the same logic as generate_susie_commands.py to build a
whitelist of REAL output files based on the GWAS meta-analysis summary statistics
and matrix configurations. Any files in the output directory that aren't on
this whitelist are considered junk.

By default, runs in SAFE MODE: only reports junk files without deleting.
Use --delete-mode to actually delete identified junk files.

Meta-Analysis SuSiE Output Format:
- Output files are located at: {output_base_dir}/{matrix_type}/{pheno_causal}.{trait_type}.{subfolder}.txt
- Subfolders are "weighted" or "Meta" depending on matrix configuration
- Files come from meta-analysis of two populations with one population either:
  * Having the causal variant removed (drop_causal)
  * Having effect sizes matched to the other population (matched_to_*)

Usage:

python cleanup_susie_meta_junk_files.py \
    --gwas-dir /path/to/gwas/meta_analysis_results \
    --ld-dir /path/to/ld \
    --output-dir /path/to/output \
    [--populations "AFR,EUR"] \
    [--analysis-type "random"] \
    [--heritability 0.0003] \
    [--matrix-types "matrix1,matrix2"] \
    [--filter-keyword "drop_causal"] \
    [--delete-mode]

"""

import argparse
import glob
import os
import sys
from pathlib import Path


def detect_matrix_subdirs(ld_base_dir, matrix_filter=None):
    """
    Automatically detect matrix subdirectories and their subfolder structure.
    
    Parameters:
    - ld_base_dir: Base directory containing LD matrices
    - matrix_filter: Optional list of matrix subdirectory names to include
    
    Returns:
    - dict: mapping of matrix_type -> list of subfolder names (weighted, Meta, etc.)
    """
    matrix_configs = {}
    
    if not os.path.exists(ld_base_dir):
        print(f"Warning: LD base directory does not exist: {ld_base_dir}")
        return matrix_configs
    
    for item in os.listdir(ld_base_dir):
        matrix_dir = os.path.join(ld_base_dir, item)
        if not os.path.isdir(matrix_dir):
            continue
        
        # Apply matrix filter if specified
        if matrix_filter and item not in matrix_filter:
            continue
        
        # Check for 'weighted' or 'Meta' subdirectory
        weighted_path = os.path.join(matrix_dir, "weighted")
        meta_path = os.path.join(matrix_dir, "Meta")
        
        matrix_configs[item] = []
        
        if os.path.isdir(weighted_path):
            matrix_configs[item].append("weighted")
        
        if os.path.isdir(meta_path):
            matrix_configs[item].append("Meta")
        
        if len(matrix_configs[item]) > 0:
            print(f"Found matrix type: {item} -> {', '.join(matrix_configs[item])}")
        else:
            print(f"Warning: No 'weighted' or 'Meta' subdirectory found in {matrix_dir}")
    
    return matrix_configs


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


def parse_gwas_filename(meta_filename, analysis_type='random', target_heritability=None):
    """
    Parse meta-analysis GWAS filename to extract components.
    
    Expected format: Pheno#.#_causal.heritability_#.#.continuous/binary_0.25.chr#.random_meta.out or .fixed_meta.out
    
    Parameters:
    - meta_filename: The GWAS filename
    - analysis_type: Type of analysis ('random' or 'fixed')
    - target_heritability: Target heritability value to filter for
    
    Returns:
    - dict: Dictionary with parsed components or None if parsing fails
    """
    # Determine the expected extension based on analysis type
    if analysis_type == 'random':
        expected_extension = '.random_meta.out'
    elif analysis_type == 'fixed':
        expected_extension = '.fixed_meta.out'
    else:
        expected_extension = '.meta.out'
    
    # Check if the filename has the expected extension
    if not meta_filename.endswith(expected_extension):
        return None
    
    # Remove the extension
    base_name = meta_filename.replace(expected_extension, '')
    
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


def build_valid_files_whitelist(gwas_base_dir, matrix_configs, output_base_dir, analysis_type='random', target_heritability=None, filter_keyword=None):
    """
    Build a set of all REAL output file paths based on GWAS meta-analysis files and matrix configurations.
    
    For each GWAS meta-analysis file and each matrix configuration, the real output file is:
    - {output_base_dir}/{matrix_type}/{pheno_causal}.{trait_type}.{subfolder}.txt
    
    Parameters:
    - gwas_base_dir: Base directory containing GWAS meta-analysis results (with causal directories)
    - matrix_configs: Dict of matrix_type -> list of subfolders
    - output_base_dir: Base output directory
    - analysis_type: Type of analysis ('random' or 'fixed')
    - target_heritability: Target heritability value to filter for
    - filter_keyword: Keyword to include in filenames (None = exclude default keywords)
    
    Returns:
    - set: Valid output file paths (absolute paths)
    """
    valid_files = set()
    
    # Check if GWAS base directory exists
    if not os.path.exists(gwas_base_dir):
        print(f"Warning: GWAS base directory does not exist: {gwas_base_dir}")
        return valid_files
    
    # Find all directories with "_causal" in the name
    causal_dirs = glob.glob(os.path.join(gwas_base_dir, "*_causal*"))
    
    # Determine the file extension pattern based on analysis type
    if analysis_type == 'random':
        file_pattern = "*.random_meta.out"
    elif analysis_type == 'fixed':
        file_pattern = "*.fixed_meta.out"
    else:
        file_pattern = "*.meta.out"
    
    # Process each causal directory
    for causal_dir in causal_dirs:
        causal_dir_name = os.path.basename(causal_dir)
        
        # Find all meta-analysis files in this directory
        meta_files = glob.glob(os.path.join(causal_dir, file_pattern))
        
        # Apply filter keyword logic
        if filter_keyword is not None:
            # Include only files containing the filter keyword
            meta_files = [f for f in meta_files
                         if filter_keyword in os.path.basename(f)]
        else:
            # Default behavior: exclude files with unwanted patterns
            meta_files = [f for f in meta_files
                         if 'drop_causal' not in os.path.basename(f)
                         and 'matched_to_' not in os.path.basename(f)]
        
        # Process each meta file
        for meta_file in meta_files:
            meta_filename = os.path.basename(meta_file)
            
            # Parse the filename
            parsed = parse_gwas_filename(meta_filename, analysis_type, target_heritability)
            
            if parsed is None:
                continue
            
            pheno_causal = parsed['pheno_causal']
            trait_type = parsed['trait_type']
            
            # For each matrix configuration
            for matrix_type, subfolders in matrix_configs.items():
                for subfolder in subfolders:
                    # Output filename format: Pheno#.#_causal.continuous/binary_0.25.weighted/Meta.txt
                    output_filename = f"{pheno_causal}.{trait_type}.{subfolder}.txt"
                    
                    # Construct output directory path
                    if len(matrix_configs[matrix_type]) > 1:
                        # If multiple subfolders, include subfolder in dir name
                        output_dir = os.path.join(output_base_dir, f"{matrix_type}.{subfolder}")
                    else:
                        # Single subfolder case
                        output_dir = os.path.join(output_base_dir, matrix_type)
                    
                    output_path = os.path.join(output_dir, output_filename)
                    
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
        dirs[:] = [d for d in dirs if d != 'timing']
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


def main():
    parser = argparse.ArgumentParser(
        description="Identify and remove junk files from meta-analysis SuSiE fine-mapping output directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Safe mode: just list junk files (DEFAULT)
  python cleanup_susie_meta_junk_files.py --gwas-dir /gwas --ld-dir /ld --output-dir /output
  
  # Actually delete junk files (requires --delete-mode)
  python cleanup_susie_meta_junk_files.py --gwas-dir /gwas --ld-dir /ld --output-dir /output --delete-mode
  
  # With custom analysis type and heritability filter
  python cleanup_susie_meta_junk_files.py --gwas-dir /gwas --ld-dir /ld --output-dir /output \\
    --analysis-type fixed --heritability 0.0005
  
  # Filter for specific causal manipulation (e.g., drop_causal)
  python cleanup_susie_meta_junk_files.py --gwas-dir /gwas --ld-dir /ld --output-dir /output \\
    --filter-keyword "drop_causal"
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--gwas-dir", 
        type=str, 
        required=True,
        help="Base directory containing GWAS meta-analysis results with causal directories"
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
        "--analysis-type",
        type=str,
        choices=['random', 'fixed'],
        default='random',
        help="Type of meta-analysis: random or fixed effect (default: random)"
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
        help="Comma-delimited list of matrix subdirectory names to process (e.g., 'beta-inferred_2000,beta-inferred_3000')"
    )
    
    parser.add_argument(
        "--filter-keyword",
        type=str,
        default=None,
        help="Keyword to include in filenames (e.g., 'drop_causal', 'matched_to_'). If specified, only files containing this keyword are processed. If not specified, files containing 'drop_causal' or 'matched_to_' are EXCLUDED (default behavior)."
    )
    
    parser.add_argument(
        "--delete-mode",
        action="store_true",
        default=False,
        help="ENABLE deletion of junk files. By default, only lists files without deleting (SAFE MODE)"
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    for dir_path, name in [(args.gwas_dir, 'GWAS'), (args.ld_dir, 'LD'), (args.output_dir, 'Output')]:
        if not os.path.exists(dir_path):
            print(f"Error: {name} directory does not exist: {dir_path}", file=sys.stderr)
            sys.exit(1)
    
    # Validate heritability
    if args.heritability < 0:
        print(f"Error: Heritability must be non-negative, got {args.heritability}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Meta-Analysis SuSiE Fine-Mapping Junk File Cleanup")
    print(f"{'='*60}")
    print(f"GWAS directory: {args.gwas_dir}")
    print(f"LD directory: {args.ld_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Analysis type: {args.analysis_type.upper()}")
    print(f"Heritability filter: {args.heritability}")
    
    if args.filter_keyword:
        print(f"File filter: INCLUDING files containing '{args.filter_keyword}'")
    else:
        print(f"File filter: EXCLUDING files containing 'drop_causal' or 'matched_to_' (default)")
    
    if args.delete_mode:
        print(f"Mode: DELETE MODE (files WILL be deleted)")
    else:
        print(f"Mode: SAFE MODE (no files will be deleted)")
    
    print(f"{'='*60}\n")
    
    # Parse matrix filter
    matrix_filter = parse_matrix_filter(args.matrix_types)
    
    # Detect matrix configurations
    print("Detecting matrix subdirectories...")
    matrix_configs = detect_matrix_subdirs(args.ld_dir, matrix_filter)
    
    if not matrix_configs:
        print("Error: No valid matrix configurations found!", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(matrix_configs)} matrix type(s)\n")
    
    # Build whitelist
    print("Building whitelist of valid output files...")
    valid_files = build_valid_files_whitelist(
        args.gwas_dir,
        matrix_configs,
        args.output_dir,
        args.analysis_type,
        args.heritability,
        args.filter_keyword
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
