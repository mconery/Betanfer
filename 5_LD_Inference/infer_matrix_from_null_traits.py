#!/usr/bin/env python3
"""
infer_matrix_from_null_traits.py

Create LD (linkage disequilibrium) matrices from tabix-indexed GWAS summary statistics
using beta correlation inference. This script samples a specified number of traits
from available indexed summary statistics files and computes correlations between
their effect sizes to infer LD structure.

This approach uses the correlation of effect sizes across multiple traits as a proxy
for linkage disequilibrium, which can be useful when reference genotype data is not
available or when working with summary statistics from diverse populations.

Usage
-----
python create_ld_matrix_from_indexed_sumstats.py \
  --indexed_dir /path/to/indexed/sumstats \
  --variant_file /path/to/variants.txt \
  --output_file /path/to/ld_matrix.tsv.gz \
  --n_traits 100 \
  [--random_seed 5] \
  [--trait_type all] \
  [--variant_col rs_number] \
  [--beta_col beta] \
  [--variant_id_format CHR:POS:REF:ALT] \
  [--device cuda] \
  [--save_binary_prefix /path/to/binary_output]

Arguments
---------
--indexed_dir         Directory containing tabix-indexed summary statistics [REQUIRED]
--variant_file        File containing variant IDs to extract (one per line) [REQUIRED]
--output_file         Output file for the LD matrix [REQUIRED]
--n_traits           Number of traits to sample for correlation inference [REQUIRED]
--random_seed         Random seed for trait sampling [default: 5]
--trait_type          Trait type filter: 'all', 'continuous', or 'binary' [default: all]
--variant_col         Column name for variant IDs in summary stats [default: rs_number]
--beta_col           Column name for beta values [default: beta]
--variant_id_format   Format string for variant IDs [default: CHR:POS:REF:ALT]
--device             Device for computation ('cuda' or 'cpu') [default: cuda]
--save_binary_prefix  Optional prefix for PLINK binary format output

Features
--------
- Automatically discovers tabix-indexed summary statistics files
- Filter traits by type (continuous, binary, or all)
- Random sampling of traits with reproducible seeding
- Parallel beta extraction using tabix indexing
- GPU-accelerated correlation computation
- Optional binary matrix output in PLINK format

"""

import os
import glob
import argparse
import random
import sys
from pathlib import Path
import pandas as pd

from extraction_utils import (
    load_variant_list,
    make_var_dict_for_positional_ids,
    extract_betas_batch_region,
    compute_beta_correlation_matrix,
    save_beta_matrix, sample_traits,
    find_indexed_sumstats_files
)

def filter_traits_by_type(indexed_files, trait_type):
    """
    Filter indexed files based on trait type (continuous or binary).

    This function identifies trait type by examining the filename for keywords:
    - 'binary' or 'case_control' -> binary traits
    - 'quantitative' or 'continuous' -> continuous traits

    Args:
        indexed_files: List of (gz_file, tbi_file) tuples
        trait_type: Filter type - 'all', 'continuous', or 'binary'

    Returns:
        Filtered list of (gz_file, tbi_file) tuples
    """
    if trait_type == 'all':
        return indexed_files

    filtered = []

    for gz_file, tbi_file in indexed_files:
        basename = os.path.basename(gz_file).lower()

        # Detect trait type from filename
        is_binary = any(keyword in basename for keyword in ['binary', 'case_control', 'cc'])
        is_continuous = any(keyword in basename for keyword in ['quantitative', 'continuous', 'quant'])

        # If trait type not explicitly in filename, try to infer (optional)
        # You may need to adjust this logic based on your naming convention

        if trait_type == 'binary' and is_binary:
            filtered.append((gz_file, tbi_file))
        elif trait_type == 'continuous' and is_continuous:
            filtered.append((gz_file, tbi_file))
        elif trait_type == 'continuous' and not is_binary and not is_continuous:
            # Default to continuous if type not clear (conservative approach)
            filtered.append((gz_file, tbi_file))

    return filtered

def main():
    parser = argparse.ArgumentParser(description="Create LD matrix from tabix-indexed GWAS summary statistics using beta correlation inference.")
    parser.add_argument("--indexed_dir", type=str, required=True,
                        help="Directory containing tabix-indexed summary statistics files.")
    parser.add_argument("--variant_file", type=str, required=True,
                        help="File containing variant IDs to extract (one per line).")
    parser.add_argument("--output_file", type=str, required=True,
                        help="Output file for the LD matrix (will be compressed).")
    parser.add_argument("--n_traits", type=int, required=True,
                        help="Number of traits to sample for correlation inference.")
    parser.add_argument("--random_seed", type=int, default=5,
                        help="Random seed for trait sampling [default: 5]")
    parser.add_argument("--trait_type", type=str, default="all", choices=["all", "continuous", "binary"],
                        help="Trait type filter: 'all', 'continuous', or 'binary' [default: all]")
    parser.add_argument("--variant_col", type=str, default="rs_number",
                        help="Column name for variant IDs in summary stats [default: rs_number]")
    parser.add_argument("--beta_col", type=str, default="beta",
                        help="Column name for beta values [default: beta]")
    parser.add_argument("--variant_id_format", type=str, default="CHR:POS:REF:ALT",
                        help="Format string for variant IDs [default: CHR:POS:REF:ALT]")
    parser.add_argument("--device", type=str, default="cuda", choices=["cuda", "cpu"],
                        help="Device for computation [default: cuda]")
    parser.add_argument("--save_binary_prefix", type=str, default=None,
                        help="Optional prefix for PLINK binary format output (.vcor1 and .vcor1.vars)")

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.indexed_dir):
        print(f"Error: Indexed directory does not exist: {args.indexed_dir}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.variant_file):
        print(f"Error: Variant file does not exist: {args.variant_file}", file=sys.stderr)
        sys.exit(1)

    if args.n_traits <= 0:
        print(f"Error: Number of traits must be positive, got: {args.n_traits}", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(os.path.abspath(args.output_file))
    os.makedirs(output_dir, exist_ok=True)

    print(f"Creating LD matrix from indexed summary statistics")
    print(f"Input directory: {args.indexed_dir}")
    print(f"Variant file: {args.variant_file}")
    print(f"Output file: {args.output_file}")
    print(f"Number of traits: {args.n_traits}")
    print(f"Trait type: {args.trait_type}")
    print(f"Random seed: {args.random_seed}")
    print(f"Device: {args.device}")

    # Step 1: Find indexed summary statistics files
    print("\nStep 1: Discovering indexed summary statistics files...")
    indexed_files = find_indexed_sumstats_files(args.indexed_dir)

    if not indexed_files:
        print(f"Error: No tabix-indexed files found in {args.indexed_dir}", file=sys.stderr)
        print("Expected to find pairs of .gz and .gz.tbi files.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(indexed_files)} indexed summary statistics files")

    # Show some examples
    print("Sample files:")
    for i, (gz_file, tbi_file) in enumerate(indexed_files[:5]):
        print(f"  {os.path.basename(gz_file)}")
    if len(indexed_files) > 5:
        print(f"  ... and {len(indexed_files) - 5} more files")

    # Step 1b: Filter traits by type
    if args.trait_type != "all":
        print(f"\nStep 1b: Filtering for {args.trait_type} traits...")
        indexed_files = filter_traits_by_type(indexed_files, args.trait_type)

        if not indexed_files:
            print(f"Error: No {args.trait_type} trait files found after filtering", file=sys.stderr)
            sys.exit(1)

        print(f"After filtering: {len(indexed_files)} {args.trait_type} trait files")
        print("Sample files:")
        for i, (gz_file, tbi_file) in enumerate(indexed_files[:5]):
            print(f"  {os.path.basename(gz_file)}")
        if len(indexed_files) > 5:
            print(f"  ... and {len(indexed_files) - 5} more files")

    # Step 2: Load variants to extract
    print("\nStep 2: Loading variants to extract...")
    try:
        # Assume the variant file has a single column with variant IDs and no header
        variants_to_extract = set(pd.read_csv(args.variant_file, header=None)[0])      
    except Exception as e:
        print(f"Error loading variant file: {e}", file=sys.stderr)
        print("Trying alternative approach...", file=sys.stderr)

        # Fallback: read file as plain text, one variant per line
        try:
            with open(args.variant_file, 'r') as f:
                lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
            variants_to_extract = set(lines)
            print(f"Loaded {len(variants_to_extract)} variants using fallback method")
        except Exception as e2:
            print(f"Error with fallback method: {e2}", file=sys.stderr)
            sys.exit(1)

    if not variants_to_extract:
        print(f"Error: No variants loaded from {args.variant_file}", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(variants_to_extract)} variants to extract")

    # Step 3: Create variant position mapping
    print("\nStep 3: Creating variant position mapping...")
    try:
        variant_map = make_var_dict_for_positional_ids(set(variants_to_extract), args.variant_id_format)
        print(f"Created position mapping for {len(variant_map)} variants")
        #Get active chromosome
        active_chromosome = list(set([x[0] for x in variant_map.values()]))
        if len(active_chromosome) > 1:
            print("Error: Variant list spans multiple chromosomes.")
            sys.exit(1)
        else:
            active_chromosome = active_chromosome[0]
    except Exception as e:
        print(f"Error creating variant position mapping: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 4: Sample traits
    print(f"\nStep 4: Sampling {args.n_traits} traits...")
    try:
        selected_files = sample_traits(indexed_files, args.n_traits, args.random_seed, active_chromosome=active_chromosome)
        print(f"Selected {len(selected_files)} files for analysis")
        print("Selected files:")
        for i, file_path in enumerate(selected_files[:10]):
            print(f"  {i+1:3d}: {os.path.basename(file_path)}")
        if len(selected_files) > 10:
            print(f"  ... and {len(selected_files) - 10} more files")

    except Exception as e:
        print(f"Error sampling traits: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 5: Extract betas
    print("\nStep 5: Extracting beta values...")
    try:
        beta_matrix = extract_betas_batch_region(
            file_list=selected_files,
            variant_col=args.variant_col,
            beta_col=args.beta_col,
            variants_to_extract=variants_to_extract,
            variant_map=variant_map
        )
        print(f"Extracted betas: {beta_matrix.shape[0]} variants × {beta_matrix.shape[1]} traits")

        # Check for missing data
        missing_fraction = beta_matrix.isna().sum().sum() / (beta_matrix.shape[0] * beta_matrix.shape[1])
        print(f"Missing data fraction: {missing_fraction:.3f}")

    except Exception as e:
        print(f"Error extracting beta values: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 6: Compute correlation matrix
    print("\nStep 6: Computing LD matrix using beta correlations...")
    try:
        ld_matrix = compute_beta_correlation_matrix(
            beta_matrix,
            device=args.device,
            save_prefix=args.save_binary_prefix
        )
        print(f"Computed LD matrix: {ld_matrix.shape[0]} × {ld_matrix.shape[1]}")

    except Exception as e:
        print(f"Error computing correlation matrix: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 7: Save matrix
    print("\nStep 7: Saving LD matrix...")
    try:
        # Convert numpy array back to DataFrame with proper row/column names
        import pandas as pd
        variants_list = sorted(variants_to_extract)
        ld_df = pd.DataFrame(ld_matrix, index=variants_list, columns=variants_list)

        save_beta_matrix(ld_df, args.output_file, compress=False)
        print(f"Saved LD matrix to: {args.output_file}")

        if args.save_binary_prefix:
            print("Also saved binary format:")
            print(f"  Matrix: {args.save_binary_prefix}.vcor1")
            print(f"  Variants: {args.save_binary_prefix}.vcor1.vars")

    except Exception as e:
        print(f"Error saving matrix: {e}", file=sys.stderr)
        sys.exit(1)

    # Summary
    print(f"\n{'='*60}")
    print("LD MATRIX CREATION COMPLETE")
    print(f"{'='*60}")
    print(f"Input directory: {args.indexed_dir}")
    print(f"Available files: {len(indexed_files)}")
    print(f"Trait type filter: {args.trait_type}")
    print(f"Sampled traits: {args.n_traits}")
    print(f"Variants analyzed: {len(variants_to_extract)}")
    print(f"Output file: {args.output_file}")
    print(f"Matrix dimensions: {ld_matrix.shape[0]} × {ld_matrix.shape[1]}")
    print(f"Random seed used: {args.random_seed}")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
