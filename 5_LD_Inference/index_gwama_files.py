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
  [--pos_col POS]

Arguments
---------
--scratch_dir         Directory to write .gz and .tbi outputs [REQUIRED]
--pattern             Glob pattern to search for input summary statistics files [REQUIRED]
--variant_col         Column name holding variant IDs [default: rs_number]
--variant_id_format   Format string for variant IDs (fields delimited by non-alphanumerics) [default: CHR:POS:REF:ALT]
--chrom_col           Output column name for chromosome [default: CHR]
--pos_col             Output column name for position [default: POS]

This script automatically detects the number of available CPUs and parallelizes across them.

"""


import os
import glob
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

from extraction_utils import index_gwama_file_with_variant_ids

def main():
    parser = argparse.ArgumentParser(description="Tabix-index GWAMA summary statistics files in parallel, extracting positions from variant IDs.")
    parser.add_argument("--scratch_dir", type=str, required=True,
                        help="Directory where indexed files (.gz and .tbi) will be written.")
    parser.add_argument("--pattern", type=str, required=True,
                        help="Glob pattern to search for GWAMA summary stats files (e.g. '/path/to/Pheno*.*.chr22.meta.out').")
    parser.add_argument("--variant_col", type=str, default="rs_number",
                        help="Column name for variant IDs in the files. [default: rs_number]")
    parser.add_argument("--variant_id_format", type=str, default="CHR:POS:REF:ALT",
                        help="Format string describing variant ID structure. [default: 'CHR:POS:REF:ALT']")
    parser.add_argument("--chrom_col", type=str, default="CHR",
                        help="Output column name for chromosome [default: CHR]")
    parser.add_argument("--pos_col", type=str, default="POS",
                        help="Output column name for position [default: POS]")
    args = parser.parse_args()

    gwama_files = glob.glob(args.pattern)
    if not gwama_files:
        print(f"No files found matching pattern: {args.pattern}")
        return

    os.makedirs(args.scratch_dir, exist_ok=True)

    task_args = [
        (fname, args.scratch_dir, args.variant_col, args.chrom_col, args.pos_col, args.variant_id_format)
        for fname in gwama_files
    ]
    max_workers = min(os.cpu_count() or 1, len(task_args))

    print(f"Indexing {len(task_args)} files using {max_workers} threads...")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {
            executor.submit(index_gwama_file_with_variant_ids, args_): args_[0]
            for args_ in task_args
        }
        for future in as_completed(future_to_file):
            fname = future_to_file[future]
            try:
                bgzf_path, tbi_path = future.result()
                print(f"Indexed {fname} -> {bgzf_path}, {tbi_path}")
            except Exception as exc:
                print(f"Error indexing {fname}: {exc}")

if __name__ == "__main__":
    main()
