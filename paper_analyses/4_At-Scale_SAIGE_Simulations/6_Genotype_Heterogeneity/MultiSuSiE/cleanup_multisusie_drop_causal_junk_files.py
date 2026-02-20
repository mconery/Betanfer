#!/usr/bin/env python3

"""
Script to identify and remove junk files from MultiSuSiE fine-mapping output directories.

This script scans MultiSuSiE output directories and removes files that were not
explicitly generated as final results.

For meta-analysis dropped causal runs:
- GWAS files are organized in population subdirectories: {gwas_dir}/{population}/*_causal*.step2.txt
- LD matrices are in matrix_type/population/ structure
- Output directory contains all matrix_type subdirectories (flat structure)
- Supports heritability filtering
- Supports filter keyword: drop_causal

SAFE MODE: By default, the script operates in safe mode (no deletion).
Use --delete-mode flag to enable actual file deletion.

Expected output filename format:
{output_base_dir}/{matrix_type}/{pheno_causal}.{trait_type}.{matrix_type}.txt

Whitelist generation logic for DROP_CAUSAL analysis:
1. Central population uses files with ".drop_causal" in name
2. Other populations use standard files (no ".drop_causal")
3. Both types used to build expected outputs
4. Scan output directory for unexpected files (not in whitelist)
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
    
    For MultiSuSiE drop_causal analysis:
    Expected formats:
    - Central pop: Pheno#.#_causal.heritability_#.#.continuous/binary_0.25.chr#.drop_causal.step2.txt
    - Other pops: Pheno#.#_causal.heritability_#.#.continuous/binary_0.25.chr#.step2.txt
    
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
        
        # Find heritability parts and check for drop_causal keyword
        heritability_parts = []
        trait_type_index = None
        has_drop_causal = False
        
        for i in range(2, len(parts)):
            part = parts[i]
            # Check if this is trait type (continuous or binary_0.25)
            if part in ['continuous', 'binary_0'] or (part == '25' and i > 2 and parts[i-1] == 'binary_0'):
                trait_type_index = i
                break
            # Check for drop_causal keyword
            elif part == 'drop_causal':
                has_drop_causal = True
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
            'chromosome': chromosome,
            'has_drop_causal': has_drop_causal
        }
        
    except (IndexError, ValueError):
        return None


def build_whitelist(gwas_dir, matrix_types, central_pop, target_heritability=None):
    """
    Build whitelist of expected output files based on GWAS summary statistics.
    
    For drop_causal analysis:
    - Central population files have ".drop_causal" in the name
    - Other population files are standard (no ".drop_causal")
    
    Parameters:
    - gwas_dir: Base directory containing population subdirectories with GWAS stats
    - matrix_types: List of matrix type subdirectories
    - central_pop: Central population name (the one with drop_causal files)
    - target_heritability: Filter for specific heritability value
    
    Returns:
    - set: Set of expected output filenames
    """
    whitelist = set()
    
    if not os.path.exists(gwas_dir):
        print(f"Warning: GWAS directory does not exist: {gwas_dir}")
        return whitelist
    
    # For drop_causal analysis, scan all population subdirectories
    causal_files = []
    for item in os.listdir(gwas_dir):
        pop_dir = os.path.join(gwas_dir, item)
        if os.path.isdir(pop_dir):
            # Search for causal files in this population directory
            causal_pattern = os.path.join(pop_dir, "*_causal*.step2.txt")
            causal_files.extend(glob.glob(causal_pattern))
    
    if not causal_files:
        print(f"Warning: No causal GWAS files found in {gwas_dir}")
        return whitelist
    
    # Parse each GWAS file and build expected output list
    for gwas_file in causal_files:
        filename = os.path.basename(gwas_file)
        pop_dir = os.path.dirname(gwas_file)
        pop_name = os.path.basename(pop_dir)
        
        # For drop_causal analysis:
        # - Central pop should have .drop_causal in filename
        # - Other pops should NOT have .drop_causal in filename
        if pop_name == central_pop:
            # Central pop: only include files with drop_causal
            if 'drop_causal' not in filename:
                continue
        else:
            # Other pops: only include standard files (no drop_causal)
            if 'drop_causal' in filename:
                continue
        
        parsed = parse_gwas_filename(filename, target_heritability=target_heritability)
        if parsed is None:
            continue
        
        pheno_causal = parsed['pheno_causal']
        trait_type = parsed['trait_type']
        
        # For each matrix type, create expected output filename
        for matrix_type in matrix_types:
            output_filename = f"{pheno_causal}.{trait_type}.{matrix_type}.multisusie.txt"
            whitelist.add(output_filename)
    
    return whitelist


def find_junk_files(output_dir, whitelist):
    """
    Find files in output directory that are not in the whitelist.
    
    For MultiSuSiE, looks for unexpected txt files in matrix_type subdirectories.
    
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
    
    # Find all txt files in output directory and subdirectories
    for root, dirs, files in os.walk(output_dir):
        for filename in files:
            # Look for txt files (MultiSuSiE output format)
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
    
    # Ask for confirmation
    response = input(f"\nPermanently delete these {len(junk_files)} file(s)? (yes/no): ").strip().lower()
    
    if response != 'yes':
        print("Deletion cancelled.")
        return 0
    
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


def cleanup_multisusie_drop_causal(gwas_dir, ld_dir, output_dir, central_pop, heritability=None, delete_mode=False):
    """
    Main cleanup function for MultiSuSiE drop_causal fine-mapping results.
    
    Parameters:
    - gwas_dir: Base directory containing population subdirectories with GWAS stats
    - ld_dir: Directory containing LD matrices
    - output_dir: MultiSuSiE output directory to clean
    - central_pop: Central population (the one with drop_causal files)
    - heritability: Target heritability value
    - delete_mode: Whether to actually delete files
    """
    print("=" * 80)
    print("MULTISUSIE FINE-MAPPING JUNK FILE CLEANUP (DROP_CAUSAL ANALYSIS)")
    print("=" * 80)
    print(f"GWAS base directory: {gwas_dir}")
    print(f"  (contains population subdirectories with step2_results)")
    print(f"LD directory: {ld_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Central population (drop_causal): {central_pop}")
    if heritability:
        print(f"Heritability filter: {heritability}")
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
    print(f"  Central population ({central_pop}): scanning for files with .drop_causal")
    print(f"  Other populations: scanning for standard files (no .drop_causal)")
    whitelist = build_whitelist(gwas_dir, matrix_types, central_pop, target_heritability=heritability)
    
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
        description='Clean up junk files from MultiSuSiE drop_causal fine-mapping output directories',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
IMPORTANT: This script is specifically for DROP_CAUSAL analysis, where:
- Central population uses files with .drop_causal in the name
- Other populations use standard files (without .drop_causal)

For matched_effects analysis, use the standard cleanup_multisusie_meta_junk_files.py script.

EXAMPLES:

Safe mode scan (default, no deletion):
  python cleanup_multisusie_drop_causal_junk_files.py \\
    --gwas-dir /path/to/step2_results \\
    --ld-dir /path/to/LD_matrices \\
    --output-dir /path/to/multisusie_results \\
    --central-pop EUR \\
    --heritability 0.0003

Delete mode with confirmation:
  python cleanup_multisusie_drop_causal_junk_files.py \\
    --gwas-dir /path/to/step2_results \\
    --ld-dir /path/to/LD_matrices \\
    --output-dir /path/to/multisusie_results \\
    --central-pop AFR \\
    --delete-mode

List detected matrix types:
  python cleanup_multisusie_drop_causal_junk_files.py \\
    --ld-dir /path/to/LD_matrices \\
    --list-matrices
"""
    )
    
    parser.add_argument('--gwas-dir', '-g', required=True,
                       help='Base directory containing population subdirectories with GWAS stats (step2_results)')
    parser.add_argument('--ld-dir', '-l', required=True,
                       help='Base directory containing LD matrices')
    parser.add_argument('--output-dir', '-o', required=True,
                       help='MultiSuSiE output directory to clean')
    parser.add_argument('--central-pop', '-cp', required=True,
                       help='Central population (the one with drop_causal files)')
    parser.add_argument('--heritability', type=float, default=0.0003,
                       help='Target heritability value for filtering (default: 0.0003)')
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
    
    # Validate central population directory exists
    central_pop_dir = os.path.join(args.gwas_dir, args.central_pop)
    if not os.path.exists(central_pop_dir):
        print(f"ERROR: Central population directory does not exist: {central_pop_dir}")
        return 1
    
    # Run cleanup
    success = cleanup_multisusie_drop_causal(
        gwas_dir=args.gwas_dir,
        ld_dir=args.ld_dir,
        output_dir=args.output_dir,
        central_pop=args.central_pop,
        heritability=args.heritability,
        delete_mode=args.delete_mode
    )
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
