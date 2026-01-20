#!/usr/bin/env python3

"""
cleanup_junk_ld_files.py

Identify and remove junk files from LD matrix output directories.

This script uses the same logic as create_beta_inference_commands.py to build a
whitelist of REAL output files based on the corresponding input variant files.
Any files in the output directory that aren't on this whitelist are considered junk.

Usage:

python cleanup_junk_files.py \
    --variant_dir /path/to/variant/files \
    --output_dir /path/to/output \
    [--delete-mode] \
    [--n_traits 2000] \
    [--pattern "Pheno*raw_in-sample*.vars"]

Without --delete-mode: Lists junk files to console without deleting
With --delete-mode: Actually deletes identified junk files

"""

import argparse
import glob
import os
import sys


def find_variant_files(variant_dir, pattern="Pheno*raw_in-sample*.vars"):
    """
    Find all variant files matching the specified pattern.
    
    Args:
        variant_dir: Directory containing variant files
        pattern: Glob pattern to match variant files
    
    Returns:
        List of variant file paths
    """
    search_pattern = os.path.join(variant_dir, pattern)
    variant_files = glob.glob(search_pattern)
    return sorted(variant_files)


def generate_output_filename(variant_file, output_dir, n_traits):
    """
    Generate output filename based on input variant file.
    Mirrors the logic from create_beta_inference_commands.py
    
    Args:
        variant_file: Path to input variant file
        output_dir: Output directory
        n_traits: Number of traits
    
    Returns:
        Output file path for the matrix file (.unphased.vcor1)
    """
    # Extract base filename without path and extension
    base_name = os.path.basename(variant_file)
    
    # Remove .vars extension
    if base_name.endswith('.vars'):
        base_name = base_name[:-5]
    
    # Transform filename: remove raw_in-sample.unphased.vcor1 and add inferred_N-trait.unphased.vcor1
    if '.raw_in-sample.unphased.vcor1' in base_name:
        base_name = base_name.replace('.raw_in-sample.unphased.vcor1', f'.inferred_{n_traits}-trait.unphased.vcor1')
    elif '.raw_in-sample' in base_name:
        base_name = base_name.replace('.raw_in-sample', f'.inferred_{n_traits}-trait')
    elif not base_name.endswith('.unphased.vcor1'):
        base_name += '.unphased.vcor1'
    else:
        # Fallback: just append the suffix
        base_name += f'.inferred_{n_traits}-trait.unphased.vcor1'
    
    return os.path.join(output_dir, base_name)


def build_valid_files_whitelist(variant_files, output_dir, n_traits):
    """
    Build a set of all REAL output file paths based on variant files.
    
    For each variant file, the real output files are:
    - {output_filename} (the matrix file)
    - {output_filename}.vars (the variant copy)
    
    Args:
        variant_files: List of input variant file paths
        output_dir: Output directory
        n_traits: Number of traits
    
    Returns:
        Set of valid output file paths (absolute paths)
    """
    valid_files = set()
    
    for variant_file in variant_files:
        output_file = generate_output_filename(variant_file, output_dir, n_traits)
        
        # Add both the matrix file and the .vars copy to whitelist
        valid_files.add(output_file)  # e.g., /path/output/Pheno_X.inferred_2000-trait.unphased.vcor1
        valid_files.add(f"{output_file}.vars")  # e.g., /path/output/Pheno_X.inferred_2000-trait.unphased.vcor1.vars
    
    return valid_files


def find_junk_files(output_dir, valid_files_whitelist):
    """
    Find all files in output_dir that are NOT on the whitelist.
    
    Args:
        output_dir: Directory to scan for junk files
        valid_files_whitelist: Set of valid file paths to keep
    
    Returns:
        List of junk file paths
    """
    junk_files = []
    
    # List all files in the output directory (non-recursive)
    try:
        all_files = os.listdir(output_dir)
    except Exception as e:
        print(f"Error reading directory {output_dir}: {e}", file=sys.stderr)
        return junk_files
    
    for filename in all_files:
        filepath = os.path.join(output_dir, filename)
        
        # Skip subdirectories - only process files
        if os.path.isdir(filepath):
            continue
        
        # Check if this file is on the whitelist
        if filepath not in valid_files_whitelist:
            junk_files.append(filepath)
    
    return sorted(junk_files)


def main():
    parser = argparse.ArgumentParser(
        description="Identify and remove junk files from LD matrix output directories",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Safe mode: just list junk files
  python cleanup_junk_files.py --variant_dir /variants --output_dir /output --delete-mode
  
  # Actually delete junk files
  python cleanup_junk_files.py --variant_dir /variants --output_dir /output
  
  # Process with custom n_traits value (must match your pipeline)
  python cleanup_junk_files.py --variant_dir /variants --output_dir /output --n_traits 5000 --delete-mode
        """
    )
    
    # Required arguments
    parser.add_argument(
        "--variant_dir", 
        type=str, 
        required=True,
        help="Directory containing input variant files"
    )
    
    parser.add_argument(
        "--output_dir", 
        type=str, 
        required=True,
        help="Output directory to scan for junk files"
    )
    
    # Optional arguments
    parser.add_argument(
        "--n_traits",
        type=int,
        default=2000,
        help="Number of traits (must match pipeline configuration) [default: 2000]"
    )
    
    parser.add_argument(
        "--pattern",
        type=str,
        default="Pheno*raw_in-sample*.vars",
        help="Glob pattern for variant files [default: Pheno*raw_in-sample*.vars]"
    )
    
    parser.add_argument(
        "--delete-mode",
        action="store_true",
        default=False,
        help="Delete mode: Actually deletes files instead of listing"
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.variant_dir):
        print(f"Error: Variant directory does not exist: {args.variant_dir}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.exists(args.output_dir):
        print(f"Error: Output directory does not exist: {args.output_dir}", file=sys.stderr)
        sys.exit(1)
    
    if not os.path.isdir(args.output_dir):
        print(f"Error: Path is not a directory: {args.output_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Variant directory: {args.variant_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Traits: {args.n_traits}")
    
    if args.safe_mode:
        print("Mode: SAFE MODE (no files will be deleted)\n")
    else:
        print("Mode: DELETE MODE (files will be permanently deleted)\n")
    
    # Find variant files
    print(f"Searching for variant files matching: {args.pattern}")
    variant_files = find_variant_files(args.variant_dir, args.pattern)
    
    if not variant_files:
        print(f"Error: No variant files found matching pattern: {args.pattern}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(variant_files)} variant file(s)\n")
    
    # Build whitelist of valid output files
    print("Building whitelist of valid output files...")
    valid_files = build_valid_files_whitelist(variant_files, args.output_dir, args.n_traits)
    print(f"Whitelist contains {len(valid_files)} valid output file(s)\n")
    #Print a sample of white-listed files
    print("Sample of whitelisted files:")
    for i in range(min(len(valid_files), 5)):
        print("  " + list(valid_files)[i])
    
    # Find junk files
    print(f"Scanning output directory for junk files...")
    all_junk_files = find_junk_files(args.output_dir, valid_files)
    
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
        print(f"  Failed to delete: {failed_count} file(s)", file=sys.stderr)
    print(f"  Freed space: {format_file_size(total_size)}")
    print(f"{'='*60}")
    
    if failed_count > 0:
        sys.exit(1)
    sys.exit(0)


def format_file_size(bytes_size):
    """
    Format bytes to human-readable size string.
    
    Args:
        bytes_size: Size in bytes
    
    Returns:
        Formatted string (e.g., "1.5 MB")
    """
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_size < 1024.0:
            return f"{bytes_size:.1f} {unit}"
        bytes_size /= 1024.0
    return f"{bytes_size:.1f} PB"


if __name__ == "__main__":
    main()
