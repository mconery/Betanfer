#!/usr/bin/env python3

"""
Script to identify and clean up junk files in SuSiEx output directories.

This script repurposes the file-checking logic from generate_susiex_commands.py
to identify files in output directories that aren't legitimate SuSiEx output files.

Features:
- Scans output directories recursively
- Identifies legitimate SuSiEx output files by checking format and content
- Supports safe mode (default) for listing files without deleting
- Supports delete mode (--delete-mode) for actual file removal

Output files are considered valid if they:
1. Have the correct naming pattern (Pheno#.#_causal...step2.txt)
2. Contain proper TSV format with expected header
3. Have at least 2 lines (header + data)
4. Meet minimum file size requirements
"""

import os
import argparse
import sys
from pathlib import Path


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


def parse_gwas_filename(filename):
    """
    Parse GWAS filename to validate it matches expected SuSiEx output format.
    
    Expected format: Pheno#.#_causal.heritability_#.#.continuous/binary_0.25.chr#.step2.txt
    
    Parameters:
    - filename: The filename to validate
    
    Returns:
    - dict: Dictionary with parsed components if valid, None otherwise
    """
    # Remove the .step2.txt extension
    if not filename.endswith('.step2.txt'):
        return None
    
    base_name = filename.replace('.step2.txt', '')
    
    # Split by dots
    parts = base_name.split('.')
    
    if len(parts) < 5:  # Need at least: Pheno#, #_causal, heritability_#, #, trait_type, chr#
        return None
    
    try:
        # First two parts: Pheno# and #_causal
        if not parts[0].startswith('Pheno') or not parts[1].endswith('_causal'):
            return None
        
        pheno_causal = parts[0] + "." + parts[1]  # e.g., "Pheno1.2_causal"
        
        # Find heritability parts - they start with "heritability_" and continue until we hit trait type
        heritability_parts = []
        trait_type_index = None
        
        for i in range(2, len(parts)):
            part = parts[i]
            
            # Check if this is a trait type (continuous or binary_0)
            if part in ['continuous', 'binary_0'] or (part == '25' and i > 2 and parts[i-1] == 'binary_0'):
                trait_type_index = i
                break
            else:
                heritability_parts.append(part)
        
        if trait_type_index is None:
            return None
        
        # Reconstruct heritability (e.g., "heritability_0.01" from ["heritability_0", "01"])
        heritability = ".".join(heritability_parts)
        
        # Validate heritability format
        parsed_heritability = normalize_heritability(heritability)
        if parsed_heritability is None:
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
        
        # Validate chromosome format (chr# where # is a number or letter)
        if not chromosome.startswith('chr'):
            return None
        
        return {
            'pheno_causal': pheno_causal,
            'heritability': heritability,
            'trait_type': trait_type,
            'chromosome': chromosome
        }
    
    except (IndexError, ValueError):
        return None


def is_valid_output_file(file_path, min_file_size=100):
    """
    Check if a file is a valid SuSiEx output file.
    
    Parameters:
    - file_path: Path to the file
    - min_file_size: Minimum file size in bytes to consider complete
    
    Returns:
    - bool: True if file is valid SuSiEx output
    """
    filename = os.path.basename(file_path)
    
    # First check: filename format
    parsed = parse_gwas_filename(filename)
    if parsed is None:
        return False
    
    # Second check: file size
    try:
        file_size = os.path.getsize(file_path)
        if file_size < min_file_size:
            return False
    except (IOError, OSError):
        return False
    
    # Third check: file content
    try:
        with open(file_path, 'r') as f:
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
        return False


def find_junk_files(output_dir, recursive=True):
    """
    Find all junk files in output directory.
    
    Parameters:
    - output_dir: Base output directory to scan
    - recursive: Whether to scan subdirectories
    
    Returns:
    - list: Paths to files that appear to be junk
    """
    junk_files = []
    
    if not os.path.isdir(output_dir):
        print(f"Error: Output directory does not exist: {output_dir}")
        return junk_files
    
    # Walk through directory
    if recursive:
        for root, dirs, files in os.walk(output_dir):
            # Skip timing directories
            if 'timing' in root:
                continue
            
            for filename in files:
                file_path = os.path.join(root, filename)
                
                # Skip directories in the timing folder
                if 'timing' in root:
                    continue
                
                # Check if this looks like a SuSiEx output file
                if filename.endswith('.txt'):
                    if not is_valid_output_file(file_path):
                        junk_files.append(file_path)
    else:
        # Non-recursive: just check immediate files
        for filename in os.listdir(output_dir):
            file_path = os.path.join(output_dir, filename)
            
            if os.path.isfile(file_path) and filename.endswith('.txt'):
                if not is_valid_output_file(file_path):
                    junk_files.append(file_path)
    
    return junk_files


def print_summary(junk_files, total_size):
    """
    Print a summary of junk files found.
    
    Parameters:
    - junk_files: List of junk file paths
    - total_size: Total size of junk files in bytes
    """
    print(f"\n{'='*70}")
    print(f"JUNK FILE CLEANUP SUMMARY")
    print(f"{'='*70}")
    print(f"Total junk files found: {len(junk_files)}")
    print(f"Total size: {total_size / (1024*1024):.2f} MB ({total_size:,} bytes)")
    print(f"{'='*70}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Identify and optionally delete junk SuSiEx output files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List junk files (safe mode, default)
  python cleanup_susiex_meta_junk_files.py /path/to/output

  # Actually delete junk files
  python cleanup_susiex_meta_junk_files.py /path/to/output --delete-mode

  # Non-recursive (only check immediate directory)
  python cleanup_susiex_meta_junk_files.py /path/to/output --no-recursive
        """
    )
    
    parser.add_argument(
        'output_dir',
        type=str,
        help='Base output directory to scan for junk files'
    )
    
    parser.add_argument(
        '--delete-mode',
        action='store_true',
        default=False,
        help='Actually delete junk files (default is safe mode: list only)'
    )
    
    parser.add_argument(
        '--no-recursive',
        action='store_true',
        default=False,
        help='Only scan immediate directory, do not recurse into subdirectories'
    )
    
    parser.add_argument(
        '--min-file-size',
        type=int,
        default=100,
        help='Minimum file size in bytes to consider a valid output file (default: 100)'
    )
    
    args = parser.parse_args()
    
    output_dir = args.output_dir
    delete_mode = args.delete_mode
    recursive = not args.no_recursive
    min_file_size = args.min_file_size
    
    # Validate output directory
    if not os.path.isdir(output_dir):
        print(f"Error: Output directory does not exist: {output_dir}")
        sys.exit(1)
    
    # Find junk files
    print(f"Scanning directory: {output_dir}")
    print(f"Recursive scanning: {recursive}")
    print(f"Minimum file size: {min_file_size} bytes")
    print()
    
    junk_files = find_junk_files(output_dir, recursive=recursive)
    
    if not junk_files:
        print("No junk files found!")
        sys.exit(0)
    
    # Calculate total size
    total_size = 0
    try:
        for file_path in junk_files:
            if os.path.exists(file_path):
                total_size += os.path.getsize(file_path)
    except (IOError, OSError):
        pass
    
    # Print junk files
    print(f"Junk files identified ({len(junk_files)} total):\n")
    for i, file_path in enumerate(junk_files, 1):
        try:
            file_size = os.path.getsize(file_path)
            size_str = f"{file_size / 1024:.1f} KB" if file_size < 1024*1024 else f"{file_size / (1024*1024):.1f} MB"
            print(f"{i:4d}. {file_path:<80s} ({size_str})")
        except (IOError, OSError):
            print(f"{i:4d}. {file_path:<80s} (size unknown)")
    
    print_summary(junk_files, total_size)
    
    # Handle deletion
    if delete_mode:
        print("⚠️  DELETE MODE ACTIVE")
        print("Proceeding to delete junk files...")
        
        deleted_count = 0
        failed_count = 0
        failed_files = []
        
        for file_path in junk_files:
            try:
                os.remove(file_path)
                deleted_count += 1
                print(f"✓ Deleted: {file_path}")
            except (IOError, OSError) as e:
                failed_count += 1
                failed_files.append((file_path, str(e)))
                print(f"✗ Failed to delete: {file_path}")
                print(f"  Error: {e}")
        
        print(f"\n{'='*70}")
        print(f"DELETION SUMMARY")
        print(f"{'='*70}")
        print(f"Successfully deleted: {deleted_count}/{len(junk_files)} files")
        print(f"Failed to delete: {failed_count}/{len(junk_files)} files")
        print(f"{'='*70}\n")
        
        if failed_files:
            print("Failed deletions:")
            for file_path, error in failed_files:
                print(f"  - {file_path}")
                print(f"    {error}")
    
    else:
        print("SAFE MODE: Files listed above but NOT deleted.")
        print("To actually delete these files, run with --delete-mode flag:")
        print(f"\n  python cleanup_susiex_meta_junk_files.py {output_dir} --delete-mode\n")


if __name__ == '__main__':
    main()
