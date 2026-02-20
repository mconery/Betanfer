#!/usr/bin/env python3
"""
index_single_gwama_file.py

Single-file tabix indexing utility for GWAMA/BOLT summary statistic files.
Extracts genomic coordinates directly from variant IDs and creates a tabix-indexed file.

This script processes a single input file, extracting chromosome and position information
from variant IDs (e.g., CHR:POS:REF:ALT format) and writes a BGZF-compressed,
tabix-indexed version to the specified scratch directory.

Usage
-----
python index_single_gwama_file.py \
  --input_file /path/to/input/file.meta.out \
  --scratch_dir /path/to/scratch_dir \
  [--variant_col rs_number] \
  [--variant_id_format CHR:POS:REF:ALT] \
  [--chrom_col CHR] \
  [--pos_col POS]

Arguments
---------
--input_file          Path to the input summary statistics file [REQUIRED]
--scratch_dir         Directory to write .gz and .tbi outputs [REQUIRED]
--variant_col         Column name holding variant IDs [default: rs_number]
--variant_id_format   Format string for variant IDs [default: CHR:POS:REF:ALT]
--chrom_col           Output column name for chromosome [default: CHR]
--pos_col             Output column name for position [default: POS]

Example
-------
python index_single_gwama_file.py \
  --input_file /path/to/Pheno123.binary.chr7.meta.out \
  --scratch_dir /lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats \
  --variant_col rs_number \
  --variant_id_format CHR:POS:REF:ALT

"""

import os
import argparse
import sys
from pathlib import Path

from extraction_utils import index_gwama_file_with_variant_ids

def main():
    parser = argparse.ArgumentParser(description="Tabix-index a single GWAMA summary statistics file, extracting positions from variant IDs.")
    parser.add_argument("--input_file", type=str, required=True,
                        help="Path to the input summary statistics file.")
    parser.add_argument("--scratch_dir", type=str, required=True,
                        help="Directory where indexed files (.gz and .tbi) will be written.")
    parser.add_argument("--variant_col", type=str, default="rs_number",
                        help="Column name for variant IDs in the file. [default: rs_number]")
    parser.add_argument("--variant_id_format", type=str, default="CHR:POS:REF:ALT",
                        help="Format string describing variant ID structure. [default: 'CHR:POS:REF:ALT']")
    parser.add_argument("--chrom_col", type=str, default="CHR",
                        help="Output column name for chromosome [default: CHR]")
    parser.add_argument("--pos_col", type=str, default="POS",
                        help="Output column name for position [default: POS]")
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file does not exist: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    
    # Create scratch directory if it doesn't exist
    os.makedirs(args.scratch_dir, exist_ok=True)
    
    # Prepare arguments for the indexing function
    task_args = (
        args.input_file,
        args.scratch_dir,
        args.variant_col,
        args.chrom_col,
        args.pos_col,
        args.variant_id_format
    )
    
    try:
        print(f"Indexing file: {os.path.basename(args.input_file)}")
        bgzf_path, tbi_path = index_gwama_file_with_variant_ids(task_args)
        print(f"Successfully indexed -> {os.path.basename(bgzf_path)}, {os.path.basename(tbi_path)}")
        print(f"Output files:")
        print(f"  BGZF: {bgzf_path}")
        print(f"  Index: {tbi_path}")
        
    except Exception as exc:
        print(f"Error indexing {os.path.basename(args.input_file)}: {exc}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()