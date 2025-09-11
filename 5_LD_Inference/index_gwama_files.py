#!/usr/bin/env python3
"""
index_gwama_files.py

Parallelized utility to tabix-index GWAMA/BOLT summary statistic files,
extracting genomic coordinates directly from variant IDs using a flexible, user-specified format.

This script identifies all files matching a user-supplied pattern and processes them in parallel,
writing BGZF-compressed/Tabix-indexed versions to a specified scratch directory. Genomic
position columns (chromosome and position) are extracted from a variant ID field such as CHR:POS:REF:ALT,
enabling indexing even when no BIM is provided and coordinate mapping must be inferred.

Usage
-----
python index_gwama_files.py \
  --scratch_dir /path/to/scratch_dir \
  --pattern "/path/to/GWAS/pheno*chr22*.out" \
  [--variant_col rs_number] \
  [--variant_id_format CHR:POS:REF:ALT] \
  [--chrom_col CHR] \
  [--pos_col POS] \
  [--filename_regex "^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.meta\\.out$"]

Arguments
---------
--scratch_dir         Directory to write .gz and .tbi outputs [REQUIRED]
--pattern             Glob pattern to search for input summary statistics files [REQUIRED]
--variant_col         Column name holding variant IDs [default: rs_number]
--variant_id_format   Format string for variant IDs (fields delimited by non-alphanumerics) [default: CHR:POS:REF:ALT]
--chrom_col           Output column name for chromosome [default: CHR]
--pos_col             Output column name for position [default: POS]
--filename_regex      Regex pattern to filter filenames (applied to basename only) [optional]

This script automatically detects the number of available CPUs and parallelizes across them.

Example
-------
python index_gwama_files.py \
  --scratch_dir /lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats_testing \
  --pattern "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results/*/Pheno*.meta.out" \
  --filename_regex "^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.meta\\.out$"

"""

import os
import glob
import argparse
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

from extraction_utils import index_gwama_file_with_variant_ids

def find_valid_meta_files(pattern, filename_regex=None):
    """
    Find files matching the glob pattern and optionally filter by regex on basename.
    
    Args:
        pattern: Glob pattern to search for files
        filename_regex: Optional regex string to filter filenames (applied to basename only)
        
    Returns:
        List of file paths that match both pattern and regex filter
    """
    # Get all files matching the glob pattern
    all_files = glob.glob(pattern)
    
    # If no regex filter specified, return all files
    if not filename_regex:
        return all_files
    
    # Compile the regex pattern
    try:
        regex_pattern = re.compile(filename_regex)
    except re.error as e:
        raise ValueError(f"Invalid regex pattern '{filename_regex}': {e}")
    
    # Filter files by regex on basename
    filtered_files = [
        f for f in all_files 
        if regex_pattern.match(os.path.basename(f))
    ]
    
    return filtered_files

def main():
    parser = argparse.ArgumentParser(description="Tabix-index GWAMA summary statistics files in parallel, extracting positions from variant IDs.")
    parser.add_argument("--scratch_dir", type=str, required=True,
                        help="Directory where indexed files (.gz and .tbi) will be written.")
    parser.add_argument("--pattern", type=str, required=True,
                        help="Glob pattern to search for GWAMA summary stats files (e.g. '/path/to/Pheno*.meta.out').")
    parser.add_argument("--variant_col", type=str, default="rs_number",
                        help="Column name for variant IDs in the files. [default: rs_number]")
    parser.add_argument("--variant_id_format", type=str, default="CHR:POS:REF:ALT",
                        help="Format string describing variant ID structure. [default: 'CHR:POS:REF:ALT']")
    parser.add_argument("--chrom_col", type=str, default="CHR",
                        help="Output column name for chromosome [default: CHR]")
    parser.add_argument("--pos_col", type=str, default="POS",
                        help="Output column name for position [default: POS]")
    parser.add_argument("--filename_regex", type=str, default=None,
                        help="Optional regex pattern to filter filenames (applied to basename only). Example: '^Pheno\\\\d+\\\\.(binary|continuous)\\\\.chr\\\\d+\\\\.meta\\\\.out$'")
    
    args = parser.parse_args()

    # Find files matching pattern and optional regex filter
    gwama_files = find_valid_meta_files(args.pattern, args.filename_regex)
    
    if not gwama_files:
        print(f"No files found matching pattern: {args.pattern}")
        if args.filename_regex:
            print(f"with filename regex: {args.filename_regex}")
        return

    print(f"Found {len(gwama_files)} files matching criteria")
    if args.filename_regex:
        print(f"Pattern: {args.pattern}")
        print(f"Filename regex: {args.filename_regex}")
        print("Sample matched files:")
        for f in gwama_files[:5]:  # Show first 5 files
            print(f"  {f}")
        if len(gwama_files) > 5:
            print(f"  ... and {len(gwama_files) - 5} more files")

    os.makedirs(args.scratch_dir, exist_ok=True)

    task_args = [
        (fname, args.scratch_dir, args.variant_col, args.chrom_col, args.pos_col, args.variant_id_format)
        for fname in gwama_files
    ]
    max_workers = min(os.cpu_count() or 1, len(task_args))

    print(f"Indexing {len(task_args)} files using {max_workers} threads...")

    success_count = 0
    error_count = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {
            executor.submit(index_gwama_file_with_variant_ids, args_): args_[0]
            for args_ in task_args
        }
        for future in as_completed(future_to_file):
            fname = future_to_file[future]
            try:
                bgzf_path, tbi_path = future.result()
                print(f"Indexed {os.path.basename(fname)} -> {os.path.basename(bgzf_path)}, {os.path.basename(tbi_path)}")
                success_count += 1
            except Exception as exc:
                print(f"Error indexing {os.path.basename(fname)}: {exc}")
                error_count += 1

    print(f"\nIndexing complete: {success_count} successful, {error_count} failed")

if __name__ == "__main__":
    main()