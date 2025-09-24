#!/usr/bin/env python3

"""
index_single_saige_file.py

Single-file tabix indexing utility for SAIGE summary statistic files.

This script processes SAIGE output files that already have CHR and POS columns,
creating a tabix-indexed version for fast genomic region queries.

Unlike GWAMA files that require extracting coordinates from variant IDs,
SAIGE files already contain separate chromosome and position columns.

Usage
-----
python index_single_saige_file.py \
--input_file /path/to/saige/output.txt \
--scratch_dir /path/to/scratch_dir \
[--chrom_col CHR] \
[--pos_col POS] \
[--zero_based]

Arguments
---------
--input_file        Path to the input SAIGE summary statistics file [REQUIRED]
--scratch_dir       Directory to write .gz and .tbi outputs [REQUIRED]  
--chrom_col         Column name for chromosome [default: CHR]
--pos_col           Column name for position [default: POS]
--variant_col         Column name holding variant IDs [default: MarkerID]
--variant_id_format   Format string for variant IDs [default: CHR:POS:REF:ALT]
--zero_based        Set this flag if positions are 0-based [default: False, assumes 1-based]
--force             Overwrite existing output files [default: True]

Example
-------
python index_single_saige_file.py \
--input_file /path/to/saige_output.step2.singleAssoc.txt \
--scratch_dir /lustre/orion/bif154/scratch/mconery/indexed_saige_sumstats \
--chrom_col CHR \
--variant_col MarkerID \
--pos_col POS \
--variant_id_format CHR:POS:REF:ALT

Notes
-----
- SAIGE files must be sorted by chromosome and position for tabix indexing
- The script assumes 1-based positions (typical for GWAS) unless --zero_based is specified
- Input file should have a header line with column names
"""

import os
import argparse
import sys
from pathlib import Path
from extraction_utils import bgzip_and_tbi_sumstats_with_variant_ids

def main():
    parser = argparse.ArgumentParser(
        description="Tabix-index a single SAIGE summary statistics file with existing CHR and POS columns."
    )
    parser.add_argument("--input_file", type=str, required=True,
                        help="Path to the input SAIGE summary statistics file.")
    parser.add_argument("--scratch_dir", type=str, required=True,
                        help="Directory where indexed files (.gz and .tbi) will be written.")
    parser.add_argument("--chrom_col", type=str, default="CHR",
                        help="Column name for chromosome in the file. [default: CHR]")
    parser.add_argument("--pos_col", type=str, default="POS",
                        help="Column name for position in the file. [default: POS]")
    parser.add_argument("--variant_col", type=str, default="MarkerID",
                        help="Column name for variant IDs in the file. [default: MarkerID]")
    parser.add_argument("--variant_id_format", type=str, default="CHR:POS:REF:ALT",
                        help="Format string describing variant ID structure. [default: 'CHR:POS:REF:ALT']")
    parser.add_argument("--zero_based", action="store_true",
                        help="Set this flag if positions are 0-based. [default: False, assumes 1-based]")
    parser.add_argument("--force", action="store_true", default=True,
                        help="Overwrite existing output files. [default: True]")

    args = parser.parse_args()

    # Validate input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file does not exist: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    # Create scratch directory if it doesn't exist
    os.makedirs(args.scratch_dir, exist_ok=True)

    try:
        print(f"Indexing SAIGE file: {os.path.basename(args.input_file)}")
        print(f"Using chromosome column: {args.chrom_col}")
        print(f"Using position column: {args.pos_col}")
        print(f"Position format: {'0-based' if args.zero_based else '1-based'}")
        
        # Call the appropriate function from extraction_utils       
        bgzf_path, tbi_path = bgzip_and_tbi_sumstats_with_variant_ids(
            input_path=args.input_file,
            output_dir=args.scratch_dir,
            chrom_col_label=args.chrom_col,
            pos_col_label=args.pos_col,
            variant_col=args.variant_col,
            variant_id_format = args.variant_id_format,
            zero_based=args.zero_based,
            force=args.force,
            variant_col="MarkerID"
        )
        
        print(f"Successfully indexed -> {os.path.basename(bgzf_path)}, {os.path.basename(tbi_path)}")
        print(f"Output files:")
        print(f"  BGZF: {bgzf_path}")
        print(f"  Index: {tbi_path}")

    except Exception as exc:
        print(f"Error indexing {os.path.basename(args.input_file)}: {exc}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()