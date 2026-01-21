#!/usr/bin/env python3

"""
Script to identify and remove junk files from SuSiEx fine-mapping output directories.

This script scans SuSiEx output directories and removes files that were not
explicitly generated as final results.

For meta-analysis matched effect-size runs:
- Supports analysis-type filtering (fixed_matched_to_EUR, fixed_matched_to_AFR)
- Supports heritability filtering
- Supports filter keywords (matched_to_EUR, matched_to_AFR)

SAFE MODE: By default, the script operates in safe mode (no deletion).
Use --delete-mode flag to enable actual file deletion.

Expected output filename format:
{output_base_dir}/{pheno_causal}.{trait_type}.{matrix_type}.txt

Whitelist generation logic:
1. Parse GWAS directory to identify all causal files matching target heritability
2. Extract pheno_causal, trait_type combinations
3. For each trait_type × pheno_causal combo, build list of expected outputs
4. Expected output filenames: {pheno_causal}.{trait_type}.{matrix_type}.txt
5. Scan output directory for unexpected files (not in whitelist)
"""

import os
import glob
import argparse
import sys
from pathlib import Path
from collections import defaultdict


def detect_matrix_subdirs(ld_base_dir):
    """
    Automatically detect matrix subdirectories in LD directory.
    
    Parameters:
    - ld_base_dir: Base directory containing matrix subdirectories
    
    Returns:
    - list: List of matrix type directory names
    """
    matrix_types = []
    
    if not os.path.exists(ld_base_dir):
        print(f"Warning: LD base directory does not exist: {ld_base_dir}")
        return matrix_types
    
    for item in os.listdir(ld_base_dir):
        matrix_dir = os.path.join(ld_base_dir, item)
        if os.path.isdir(matrix_dir):
            matrix_types.append(item)
    
    return sorted(matrix_types)


def normalize_heritability(heritability_str):
    """
    Normalize heritability string for comparison.
    
    Parameters:
    - heritability_str: Heritability string (e.g., "heritability_0.0003")
    
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
    base_name = meta_filename.replace('.step2.txt', '')
    parts = base_name.split('.')
    
    if len(parts) < 5:
        return None
    
    try:
        # First two parts: Pheno# and #_causal
        pheno_causal = parts[0] + "." + parts[1]
        
        # Find heritability parts
        heritability_parts = []
        trait_type_index = None
        
        for i in range(2, len(parts)):
            part = parts[i]
            if part in ['continuous', 'binary_0'] or (part == '25' and i > 2 and parts[i-1] == 'binary_0'):
                trait_type_index = i
                break
            else:
                heritability_parts.append(part)
        
        if trait_type_index is None:
            return None
        
        heritability = ".".join(heritability_parts)
        
        # Filter by heritability if specified
        if target_heritability is not None:
            parsed_heritability = normalize_heritability(heritability)
            if parsed_heritability is None:
                return None
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


def build_whitelist(gwas_dir, matrix_types, target_heritability=None, filter_keyword=None):
    """
    Build whitelist of expected output files based on GWAS summary statistics.
    
    Parameters:
    - gwas_dir: Directory containing GWAS summary statistics
    - matrix_types: List of matrix type subdirectories
    - target_heritability: Filter for specific heritability value
    - filter_keyword: Filter for specific keyword in filename (e.g., "matched_to_EUR")
    
    Returns:
    - set: Set of expected output filenames
    """
    whitelist = set()
    
    if not os.path.exists(gwas_dir):
        print(f"Warning: GWAS directory does not exist: {gwas_dir}")
        return whitelist
    
    # Find all causal files in GWAS directory
    causal_pattern = os.path.join(gwas_dir, "*_causal*.step2.txt")
    causal_files = glob.glob(causal_pattern)
    
    # Parse each GWAS file and build expected output list
    for gwas_file in causal_files:
        filename = os.path.basename(gwas_file)
        
        # Apply filter keyword if specified
        if filter_keyword and filter_keyword not in filename:
            continue
        
        parsed = parse_gwas_filename(filename, target_heritability=target_heritability)
        
        if parsed is None:
            continue
        
        pheno_causal = parsed['pheno_causal']
        trait_type = parsed['trait_type']
        
        # For each matrix type, create expected output filename
        for matrix_type in matrix_types:
            output_filename = f"{pheno_causal}.{trait_type}.{matrix_type}.txt"
            whitelist.add(output_filename)
    
    return whitelist


def find_junk_files(output_dir, whitelist):
    """
    Find files in output directory that are not in the whitelist.
    
    Parameters:
    - output_dir: Output directory to scan
    - whitelist: Set of whitelisted filenames
    
    Returns:
    - list: List of junk file paths
    """
    junk_files = []
    
    if not os.path.exists(output_dir):
        print(f"Warning: Output directory does not exist: {output_dir}")
        return junk_files
    
    # Find all txt files in output directory
    for root, dirs, files in os.walk(output_dir):
        for filename in files:
            # Look for txt files (SuSiEx output format)
            if filename.endswith('.txt'):
                # Check if file is in whitelist
                if filename not in whitelist:
                    full_path = os.path.join(root, filename)
                    junk_files.append(full_path)
    
    return junk_files


def delete_files_with_confirmation(junk_files):
    """
    Delete files after user confirmation.
    
    Parameters:
    - junk_files: List of file paths to delete
    
    Returns:
    - int: Number of files successfully deleted
    """
    if not junk_files:
        print("No junk files to delete.")
        return 0
    
    print(f"\nFound {len(junk_files)} junk file(s) to delete:")
    for i, filepath in enumerate(junk_files[:10], 1):
        print(f"  {i}. {filepath}")
    
    if len(junk_files) > 10:
        print(f"  ... and {len(junk_files) - 10} more files")
        
    deleted_count = 0
    failed_count = 0
    
    for filepath in junk_files:
        try:
            if os.path.exists(filepath):
                os.remove(filepath)
                deleted_count += 1
                print(f"✓ Deleted: {filepath}")
            else:
                print(f"⚠ File not found: {filepath}")
        except Exception as e:
            failed_count += 1
            print(f"✗ Failed to delete {filepath}: {e}")
    
    print(f"\nDeletion summary:")
    print(f"  Successfully deleted: {deleted_count}")
    print(f"  Failed: {failed_count}")
    
    return deleted_count


def cleanup_susiex(gwas_dir, ld_dir, output_dir, heritability=None, analysis_type=None, 
                   filter_keyword=None, delete_mode=False):
    """
    Main cleanup function.
    
    Parameters:
    - gwas_dir: Directory containing GWAS summary statistics
    - ld_dir: Directory containing LD matrices
    - output_dir: SuSiEx output directory to clean
    - heritability: Target heritability value
    - analysis_type: Analysis type string (for logging)
    - filter_keyword: Filter for specific keyword in GWAS filenames
    - delete_mode: Whether to actually delete files
    """
    print("=" * 80)
    print("SUSIEX FINE-MAPPING JUNK FILE CLEANUP")
    print("=" * 80)
    print(f"GWAS directory: {gwas_dir}")
    print(f"LD directory: {ld_dir}")
    print(f"Output directory: {output_dir}")
    if heritability:
        print(f"Heritability filter: {heritability}")
    if analysis_type:
        print(f"Analysis type: {analysis_type}")
    if filter_keyword:
        print(f"Filter keyword: {filter_keyword}")
    print(f"Mode: {'DELETE' if delete_mode else 'SAFE (no deletion)'}")
    print("=" * 80)
    print()
    
    # Detect matrix types
    print(f"Detecting matrix subdirectories in {ld_dir}...")
    matrix_types = detect_matrix_subdirs(ld_dir)
    
    if not matrix_types:
        print("ERROR: No matrix subdirectories found!")
        return False
    
    print(f"Found {len(matrix_types)} matrix type(s):")
    for matrix_type in matrix_types[:5]:
        print(f"  - {matrix_type}")
    if len(matrix_types) > 5:
        print(f"  ... and {len(matrix_types) - 5} more")
    print()
    
    # Build whitelist
    print(f"Building whitelist of expected output files...")
    whitelist = build_whitelist(gwas_dir, matrix_types, target_heritability=heritability, 
                               filter_keyword=filter_keyword)
    
    print(f"Whitelist contains {len(whitelist)} expected output file(s)")
    if len(whitelist) == 0:
        print("WARNING: Whitelist is empty! This may indicate an issue with the GWAS directory.")
        print("Skipping cleanup to avoid accidentally deleting all files.")
        return False
    
    # Show sample whitelist entries
    print(f"Sample whitelist entries (first 5):")
    for filename in sorted(list(whitelist))[:5]:
        print(f"  - {filename}")
    print()
    
    # Find junk files
    print(f"Scanning output directory for junk files...")
    junk_files = find_junk_files(output_dir, whitelist)
    
    print(f"Found {len(junk_files)} potential junk file(s)")
    
    if len(junk_files) == 0:
        print("✓ No junk files found. Output directory is clean!")
        print()
        return True
    
    # Show summary of junk files
    print(f"\nJunk files summary:")
    
    # Group by directory
    junk_by_dir = defaultdict(list)
    for filepath in junk_files:
        dirpath = os.path.dirname(filepath)
        junk_by_dir[dirpath].append(filepath)
    
    for dirpath, files in sorted(junk_by_dir.items()):
        print(f"\n  Directory: {dirpath}")
        print(f"  Count: {len(files)} file(s)")
        for filepath in files[:3]:
            print(f"    - {os.path.basename(filepath)}")
        if len(files) > 3:
            print(f"    ... and {len(files) - 3} more")
    
    print()
    
    # Handle deletion
    if delete_mode:
        deleted_count = delete_files_with_confirmation(junk_files)
        return deleted_count > 0
    else:
        print("SAFE MODE: Files will NOT be deleted")
        print(f"To enable deletion, run with --delete-mode flag")
        print()
        return True


def main():
    """Main function with command line argument parsing."""
    parser = argparse.ArgumentParser(
        description='Clean up junk files from SuSiEx fine-mapping output directories',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES:

Safe mode scan (default, no deletion):
  python cleanup_susiex_meta_junk_files.py \\
    --gwas-dir /path/to/meta_analysis_results \\
    --ld-dir /path/to/LD_matrices \\
    --output-dir /path/to/susiex_results \\
    --heritability 0.0003

Delete mode with confirmation:
  python cleanup_susiex_meta_junk_files.py \\
    --gwas-dir /path/to/meta_analysis_results \\
    --ld-dir /path/to/LD_matrices \\
    --output-dir /path/to/susiex_results \\
    --analysis-type fixed_matched_to_EUR \\
    --filter-keyword matched_to_EUR \\
    --delete-mode

List detected matrix types:
  python cleanup_susiex_meta_junk_files.py \\
    --ld-dir /path/to/LD_matrices \\
    --list-matrices
"""
    )
    
    parser.add_argument('--gwas-dir', '-g', required=True,
                       help='Directory containing GWAS summary statistics')
    parser.add_argument('--ld-dir', '-l', required=True,
                       help='Base directory containing LD matrices')
    parser.add_argument('--output-dir', '-o', required=True,
                       help='SuSiEx output directory to clean')
    parser.add_argument('--heritability', type=float, default=0.0003,
                       help='Target heritability value for filtering (default: 0.0003)')
    parser.add_argument('--analysis-type', '-a', type=str, default=None,
                       help='Analysis type string (for logging purposes)')
    parser.add_argument('--filter-keyword', '-f', type=str, default=None,
                       help='Filter keyword for GWAS filenames (e.g., "matched_to_EUR")')
    parser.add_argument('--delete-mode', action='store_true',
                       help='Enable actual file deletion (default: safe mode, no deletion)')
    parser.add_argument('--list-matrices', action='store_true',
                       help='List detected matrix types and exit')
    
    args = parser.parse_args()
    
    # Handle list-matrices option
    if args.list_matrices:
        print(f"Scanning for matrix subdirectories in: {args.ld_dir}")
        matrix_types = detect_matrix_subdirs(args.ld_dir)
        
        if matrix_types:
            print(f"\nFound {len(matrix_types)} matrix type(s):")
            for matrix_type in sorted(matrix_types):
                print(f"  - {matrix_type}")
        else:
            print("No matrix subdirectories found.")
        
        return 0
    
    # Validate required arguments
    if not os.path.exists(args.gwas_dir):
        print(f"ERROR: GWAS directory does not exist: {args.gwas_dir}")
        return 1
    
    if not os.path.exists(args.ld_dir):
        print(f"ERROR: LD directory does not exist: {args.ld_dir}")
        return 1
    
    if not os.path.exists(args.output_dir):
        print(f"ERROR: Output directory does not exist: {args.output_dir}")
        return 1
    
    # Run cleanup
    success = cleanup_susiex(
        gwas_dir=args.gwas_dir,
        ld_dir=args.ld_dir,
        output_dir=args.output_dir,
        heritability=args.heritability,
        analysis_type=args.analysis_type,
        filter_keyword=args.filter_keyword,
        delete_mode=args.delete_mode
    )
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
