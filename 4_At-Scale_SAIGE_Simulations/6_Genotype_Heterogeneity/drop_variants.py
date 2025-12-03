#!/usr/bin/env python3
"""
Filter GWAS summary statistics by removing causal variants and/or random subsets.

This script processes GWAS summary statistics files, allowing removal of:
1. Causal variants (from a separate file)
2. A random percentage of variants
3. Both causal and random variants (with causal variants counted toward random %)

Output is sorted by genomic position (chromosome, then position).
"""

import argparse
import pandas as pd
import numpy as np
import random
import sys
from pathlib import Path


def natural_chr_sort_key(chr_val):
    """
    Generate sort key for chromosome values to handle X, Y, MT properly.
    Works with both "1" and "chr1" formats.
    """
    chr_str = str(chr_val).strip()
    
    # Remove "chr" prefix if present
    if chr_str.lower().startswith("chr"):
        chr_str = chr_str[3:]
    
    # Map special chromosomes
    chr_mapping = {
        "x": (23, 0),
        "y": (24, 0),
        "m": (25, 0),
        "mt": (25, 0)
    }
    
    chr_lower = chr_str.lower()
    if chr_lower in chr_mapping:
        return chr_mapping[chr_lower]
    
    try:
        return (int(chr_str), 0)
    except ValueError:
        # If conversion fails, sort alphabetically at the end
        return (999, chr_str)


def load_causal_variants(causal_file, variant_id_col):
    """
    Load causal variants from file.
    
    Parameters
    ----------
    causal_file : str
        Path to file containing causal variants
    variant_id_col : str
        Name of the column containing variant IDs
        
    Returns
    -------
    set
        Set of causal variant IDs
    """
    try:
        causal_df = pd.read_csv(causal_file, sep=r'\s+', engine='python')
        if variant_id_col not in causal_df.columns:
            raise ValueError(
                f"Column '{variant_id_col}' not found in {causal_file}. "
                f"Available columns: {list(causal_df.columns)}"
            )
        causal_variants = set(causal_df[variant_id_col].astype(str).values)
        print(f"Loaded {len(causal_variants)} causal variants from {causal_file}")
        return causal_variants
    except Exception as e:
        print(f"Error loading causal variants: {e}", file=sys.stderr)
        sys.exit(1)


def load_summary_stats(sumstats_file):
    """
    Load GWAS summary statistics file.
    
    Parameters
    ----------
    sumstats_file : str
        Path to GWAS summary statistics file
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing summary statistics
    """
    try:
        df = pd.read_csv(sumstats_file, sep=r'\s+', engine='python')
        print(f"Loaded {len(df)} variants from {sumstats_file}")
        return df
    except Exception as e:
        print(f"Error loading summary statistics: {e}", file=sys.stderr)
        sys.exit(1)


def filter_variants(df, drop_causal=False, causal_variants=None, 
                   drop_fraction=0.0, variant_id_col="MarkerID",
                   chr_col="CHR", pos_col="POS"):
    """
    Filter variants from summary statistics.
    
    Parameters
    ----------
    df : pd.DataFrame
        Summary statistics dataframe
    drop_causal : bool
        Whether to drop causal variants
    causal_variants : set or None
        Set of causal variant IDs
    drop_fraction : float
        Fraction of variants to randomly drop (0.0 to 1.0, exclusive of 1.0)
    variant_id_col : str
        Name of variant ID column
    chr_col : str
        Name of chromosome column
    pos_col : str
        Name of position column
        
    Returns
    -------
    pd.DataFrame
        Filtered and sorted summary statistics
    """
    df = df.copy()
    
    # Initialize variant mask (True = keep)
    keep_mask = pd.Series([True] * len(df), index=df.index)
    
    # Drop causal variants if requested
    if drop_causal and causal_variants:
        causal_mask = df[variant_id_col].astype(str).isin(causal_variants)
        n_causal_in_sumstats = causal_mask.sum()
        keep_mask = keep_mask & ~causal_mask
        print(f"Dropped {n_causal_in_sumstats} causal variants")
    
    # Handle random variant dropping
    if drop_fraction > 0:
        # Get indices of variants to potentially drop
        eligible_indices = df[keep_mask].index.tolist()
        n_to_drop = max(1, int(np.round(len(eligible_indices) * drop_fraction)))
        
        # Determine which variants to drop
        if drop_causal and causal_variants:
            # If causal variants are also being dropped, they count toward the percentage
            # We need to drop n_to_drop total from the original set
            n_to_drop_random = max(0, n_to_drop - n_causal_in_sumstats)
        else:
            # Only dropping random variants - must exclude causal variants
            n_to_drop_random = n_to_drop
            if causal_variants:
                # Remove causal variants from eligible pool
                eligible_indices = [
                    idx for idx in eligible_indices 
                    if df.loc[idx, variant_id_col] not in causal_variants
                ]
        
        # Random sample without replacement
        if n_to_drop_random > 0 and len(eligible_indices) > 0:
            indices_to_drop = random.sample(
                eligible_indices, 
                min(n_to_drop_random, len(eligible_indices))
            )
            keep_mask.loc[indices_to_drop] = False
            print(f"Randomly dropped {len(indices_to_drop)} variants "
                  f"({drop_fraction*100:.1f}% of eligible)")
    
    # Apply mask
    df_filtered = df[keep_mask].copy()
    
    # Sort by chromosome and position
    # Create sortable chromosome representation
    df_filtered[f'__{chr_col}_sort__'] = df_filtered[chr_col].apply(natural_chr_sort_key)
    df_filtered = df_filtered.sort_values(
        by=[f'__{chr_col}_sort__', pos_col],
        na_position='last'
    )
    df_filtered = df_filtered.drop(columns=[f'__{chr_col}_sort__'])
    
    print(f"Output: {len(df_filtered)} variants after filtering")
    print(f"Removed {len(df) - len(df_filtered)} total variants")
    
    return df_filtered


def main():
    parser = argparse.ArgumentParser(
        description="Filter GWAS summary statistics by removing causal variants "
                    "and/or random subsets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Remove all causal variants
  %(prog)s -i sumstats.txt -o filtered.txt -c -C causal.txt
  
  # Drop 10% of random variants
  %(prog)s -i sumstats.txt -o filtered.txt -r 0.1
  
  # Drop all causal variants, then drop 5% random from remaining
  %(prog)s -i sumstats.txt -o filtered.txt -c -C causal.txt -r 0.05
        """
    )
    
    # Required arguments
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input GWAS summary statistics file"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to output filtered summary statistics file"
    )
    
    # Causal variant filtering arguments
    parser.add_argument(
        "-c", "--drop-causal",
        action="store_true",
        help="Remove causal variants from summary statistics"
    )
    parser.add_argument(
        "-C", "--causal-file",
        help="Path to file containing causal variants (required if --drop-causal is used)"
    )
    
    # Column name arguments
    parser.add_argument(
        "-m", "--marker-col",
        default="MarkerID",
        help="Name of variant ID column (default: MarkerID)"
    )
    parser.add_argument(
        "--chr-col",
        default="CHR",
        help="Name of chromosome column (default: CHR)"
    )
    parser.add_argument(
        "--pos-col",
        default="POS",
        help="Name of position column (default: POS)"
    )
    
    # Random variant dropping argument
    parser.add_argument(
        "-r", "--drop-fraction",
        type=float,
        default=0.0,
        help="Fraction of variants to randomly drop (0.0 to 1.0, exclusive). "
             "Default: 0.0"
    )
    
    # Seed for reproducibility
    parser.add_argument(
        "-s", "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible variant selection (optional)"
    )
    
    args = parser.parse_args()
    
    # Validation
    if args.drop_fraction < 0 or args.drop_fraction >= 1:
        parser.error("--drop-fraction must be between 0 (inclusive) and 1 (exclusive)")
    
    if args.drop_causal and not args.causal_file:
        parser.error("--causal-file is required when --drop-causal is specified")
    
    if args.causal_file and not args.drop_causal:
        print("Warning: --causal-file specified but --drop-causal not set. "
              "Causal variants will not be removed.", file=sys.stderr)
    
    # Set random seed if provided
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)
        print(f"Using random seed: {args.seed}")
    
    # Load data
    df = load_summary_stats(args.input)
    
    causal_variants = None
    if args.drop_causal:
        causal_variants = load_causal_variants(args.causal_file, args.marker_col)
    
    # Filter variants
    df_filtered = filter_variants(
        df,
        drop_causal=args.drop_causal,
        causal_variants=causal_variants,
        drop_fraction=args.drop_fraction,
        variant_id_col=args.marker_col,
        chr_col=args.chr_col,
        pos_col=args.pos_col
    )
    
    # Write output
    try:
        df_filtered.to_csv(args.output, sep="\t", index=False)
        print(f"\nOutput written to {args.output}")
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
