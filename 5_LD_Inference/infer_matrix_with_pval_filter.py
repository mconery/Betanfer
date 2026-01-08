#!/usr/bin/env python3

"""
infer_matrix_with_pval_filter.py

Create LD (linkage disequilibrium) matrices from tabix-indexed GWAS summary statistics
using beta correlation inference, with p-value filtering to exclude traits with significant
associations in the region of interest.

This script extends the null-trait approach by:
1. Checking p-values for all variants in the LD matrix
2. Optionally including a genomic buffer region around the matrix variants
3. Dropping traits where ANY variant in the matrix (or buffer) is significant
4. Computing LD matrices only from null traits (no significant associations in region)

This filtering approach ensures that the LD matrix reflects true linkage disequilibrium
rather than confounding due to pleiotropy or multiple causal variants.

Usage
-----

python infer_ld_with_pval_filter.py \\
    --indexed_dir /path/to/indexed/sumstats \\
    --variant_file /path/to/variants.txt \\
    --output_file /path/to/ld_matrix.tsv.gz \\
    --pval_col pvalue \\
    --pval_threshold 0.05 \\
    [--random_seed 5] \\
    [--variant_col rs_number] \\
    [--beta_col beta] \\
    [--pos_col POS] \\
    [--buffer 100000] \\
    [--variant_id_format CHR:POS:REF:ALT] \\
    [--device cuda] \\
    [--save_binary_prefix /path/to/binary_output] \\
    [--stats_output /path/to/filtering_stats.txt]

Arguments
---------

--indexed_dir Directory containing tabix-indexed summary statistics [REQUIRED]

--variant_file File containing variant IDs to extract (one per line) [REQUIRED]

--output_file Output file for the LD matrix [REQUIRED]

--pval_col Column name for p-values in summary stats [REQUIRED]

--pval_threshold P-value threshold for significance [REQUIRED]

--random_seed Random seed for reproducibility [default: 5]

--variant_col Column name for variant IDs in summary stats [default: rs_number]

--beta_col Column name for beta values [default: beta]

--pos_col Column name for variant positions [default: POS]
             [OPTIONAL unless --buffer is specified]

--buffer Genomic distance (bp) to extend buffer region around matrix variants
         [default: None (no buffer)]. When specified, variants within this distance
         upstream/downstream of the matrix boundaries are also checked for significance.

--variant_id_format Format string for variant IDs [default: CHR:POS:REF:ALT]

--device Device for computation ('cuda' or 'cpu') [default: cuda]

--save_binary_prefix Optional prefix for PLINK binary format output

--stats_output Optional output file for filtering statistics and trait information

Features
--------

- Filters traits based on p-value significance thresholds
- Optional genomic buffer region to exclude nearby associations
- Uses ALL traits that pass p-value filtering (no sampling)
- Logs detailed statistics on trait filtering
- Automatically discovers tabix-indexed summary statistics files
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
from typing import List, Tuple, Dict, Set
import pandas as pd
import numpy as np
import pysam
import logging
import re

from extraction_utils import (
    load_variant_list,
    make_var_dict_for_positional_ids,
    extract_betas_batch,
    compute_beta_correlation_matrix,
    save_beta_matrix
)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def find_indexed_sumstats_files(indexed_dir):
    """
    Find all tabix-indexed summary statistics files in the specified directory.
    
    Args:
        indexed_dir: Directory containing indexed files
        
    Returns:
        List of file paths for valid .gz files with corresponding .tbi indices
    """
    indexed_files = []
    
    # Find all .gz files
    gz_pattern = os.path.join(indexed_dir, "*.gz")
    gz_files = glob.glob(gz_pattern)
    
    for gz_file in gz_files:
        tbi_file = gz_file + ".tbi"
        if os.path.exists(tbi_file):
            indexed_files.append(gz_file)
    
    return indexed_files


def get_buffer_region_variants(variant_map: Dict[str, Tuple], buffer_bp: int) -> Tuple[str, int, int]:
    """
    Define the buffer region boundaries based on the matrix variants.
    
    Args:
        variant_map: Dictionary mapping variant ID to (chromosome, position)
        buffer_bp: Buffer size in base pairs
        
    Returns:
        Tuple of (chromosome, start_pos, end_pos) for the buffer region
    """
    if not variant_map:
        raise ValueError("Variant map is empty")
    
    # Get chromosome and positions
    chroms = set()
    positions = []
    
    for variant, (chrom, pos) in variant_map.items():
        chroms.add(chrom)
        try:
            positions.append(int(pos))
        except (ValueError, TypeError):
            logger.warning(f"Skipping variant {variant} with invalid position {pos}")
    
    if len(chroms) > 1:
        raise ValueError(f"Variant map spans multiple chromosomes: {chroms}")
    
    if not positions:
        raise ValueError("No valid positions found in variant map")
    
    chrom = chroms.pop()
    min_pos = min(positions)
    max_pos = max(positions)
    
    # Extend buffer
    buffer_start = max(1, min_pos - buffer_bp)  # Ensure >= 1
    buffer_end = max_pos + buffer_bp
    
    logger.info(f"Buffer region: {chrom}:{buffer_start}-{buffer_end} "
                f"(original range: {min_pos}-{max_pos}, buffer: {buffer_bp} bp)")
    
    return chrom, buffer_start, buffer_end


def check_trait_pvalue_significance(trait_file: str, variants_to_check: Set[str],
                                   variant_map: Dict[str, Tuple],
                                   variant_col: str, pval_col: str,
                                   pval_threshold: float,
                                   buffer_region: Tuple = None) -> Tuple[bool, Dict]:
    """
    Check if a trait has any significant p-values in the variants to check.
    
    Args:
        trait_file: Path to tabix-indexed trait file
        variants_to_check: Set of variant IDs to check p-values for
        variant_map: Dictionary mapping variant ID to (chromosome, position)
        variant_col: Column name for variant IDs
        pval_col: Column name for p-values
        pval_threshold: P-value significance threshold
        buffer_region: Optional tuple (chrom, start, end) for genomic region to check
        
    Returns:
        Tuple of (is_significant, stats_dict):
            - is_significant: True if ANY variant is significant (trait should be excluded)
            - stats_dict: Dictionary with details on significant variants found
    """
    stats = {
        'trait_file': os.path.basename(trait_file),
        'variants_checked': 0,
        'variants_found': 0,
        'significant_variants': [],
        'significant_in_matrix': False,
        'significant_in_buffer': False,
        'min_pval': 1.0,
        'variants_with_data': []
    }
    
    try:
        with pysam.TabixFile(trait_file) as tabix_file:
            header = tabix_file.header[0].strip().split('\t') if tabix_file.header else []
            
            if not header:
                logger.warning(f"Could not read header from {trait_file}")
                return False, stats
            
            try:
                variant_idx = header.index(variant_col)
                pval_idx = header.index(pval_col)
            except ValueError as e:
                logger.warning(f"Could not find columns in {trait_file}: {e}")
                return False, stats
            
            # Check matrix variants
            for variant in variants_to_check:
                stats['variants_checked'] += 1
                chrom, pos = variant_map.get(variant, (None, None))
                
                if not chrom or not pos:
                    continue
                
                try:
                    records = tabix_file.fetch(reference=chrom, start=int(pos)-1, end=int(pos)+1)
                    
                    for record in records:
                        fields = record.split('\t')
                        if fields[variant_idx] == variant:
                            stats['variants_found'] += 1
                            stats['variants_with_data'].append(variant)
                            
                            try:
                                pval = float(fields[pval_idx])
                                stats['min_pval'] = min(stats['min_pval'], pval)
                                
                                if pval < pval_threshold:
                                    stats['significant_variants'].append((variant, pval))
                                    stats['significant_in_matrix'] = True
                            except (ValueError, IndexError):
                                logger.debug(f"Could not parse p-value for {variant} in {trait_file}")
                            
                            break
                
                except Exception as e:
                    logger.debug(f"Error fetching {variant} from {trait_file}: {e}")
            
            # Check buffer region if specified
            if buffer_region:
                chrom, buffer_start, buffer_end = buffer_region
                
                try:
                    records = tabix_file.fetch(reference=chrom, start=buffer_start, end=buffer_end)
                    
                    for record in records:
                        fields = record.split('\t')
                        try:
                            variant_id = fields[variant_idx]
                            pval = float(fields[pval_idx])
                            
                            # Skip if we already checked this variant
                            if variant_id in variants_to_check:
                                continue
                            
                            stats['min_pval'] = min(stats['min_pval'], pval)
                            
                            if pval < pval_threshold:
                                stats['significant_variants'].append((variant_id, pval, 'buffer'))
                                stats['significant_in_buffer'] = True
                        
                        except (ValueError, IndexError):
                            continue
                
                except Exception as e:
                    logger.debug(f"Error fetching buffer region from {trait_file}: {e}")
    
    except Exception as e:
        logger.warning(f"Error checking p-values in {trait_file}: {e}")
        return False, stats
    
    # Trait is significant if any variant in matrix or buffer is significant
    is_significant = stats['significant_in_matrix'] or stats['significant_in_buffer']
    
    return is_significant, stats


def filter_traits_by_pvalue(trait_files: List[str], variants_to_check: Set[str],
                           variant_map: Dict[str, Tuple],
                           variant_col: str, pval_col: str, pval_threshold: float,
                           buffer_region: Tuple = None) -> Tuple[List[str], Dict]:
    """
    Filter trait files to keep only those with no significant p-values in the region.
    Uses ALL traits that pass the filter (no sampling).
    
    Args:
        trait_files: List of tabix-indexed trait files
        variants_to_check: Set of variant IDs to check
        variant_map: Dictionary mapping variant ID to (chromosome, position)
        variant_col: Column name for variant IDs
        pval_col: Column name for p-values
        pval_threshold: P-value significance threshold
        buffer_region: Optional tuple (chrom, start, end) for genomic buffer region
        
    Returns:
        Tuple of (filtered_trait_files, filtering_stats)
    """
    filtered_files = []
    filtering_stats = {
        'total_traits': len(trait_files),
        'traits_checked': 0,
        'traits_retained': 0,
        'traits_filtered_out': 0,
        'min_pval_across_all': 1.0,
        'trait_details': []
    }
    
    for i, trait_file in enumerate(trait_files):
        filtering_stats['traits_checked'] += 1
        
        if (i + 1) % 100 == 0:
            logger.info(f"Checked {i + 1}/{len(trait_files)} traits...")
        
        # Check if trait has significant p-values
        is_significant, stats = check_trait_pvalue_significance(
            trait_file=trait_file,
            variants_to_check=variants_to_check,
            variant_map=variant_map,
            variant_col=variant_col,
            pval_col=pval_col,
            pval_threshold=pval_threshold,
            buffer_region=buffer_region
        )
        
        # Update global statistics
        filtering_stats['min_pval_across_all'] = min(
            filtering_stats['min_pval_across_all'],
            stats['min_pval']
        )
        filtering_stats['trait_details'].append(stats)
        
        # Retain trait if NOT significant
        if not is_significant:
            filtered_files.append(trait_file)
            filtering_stats['traits_retained'] += 1
        else:
            filtering_stats['traits_filtered_out'] += 1
    
    return filtered_files, filtering_stats


def write_filtering_stats(stats: Dict, output_path: str, pval_threshold: float, buffer_bp: int = None):
    """
    Write detailed filtering statistics to output file.
    
    Args:
        stats: Dictionary with filtering statistics
        output_path: Path to output statistics file
        pval_threshold: P-value threshold used for filtering
        buffer_bp: Optional buffer size used
    """
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("P-VALUE BASED TRAIT FILTERING STATISTICS\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("FILTERING PARAMETERS:\n")
        f.write(f"  P-value threshold: {pval_threshold}\n")
        if buffer_bp:
            f.write(f"  Genomic buffer: {buffer_bp} bp\n")
        else:
            f.write(f"  Genomic buffer: None (matrix variants only)\n")
        f.write("\n")
        
        f.write("FILTERING SUMMARY:\n")
        f.write(f"  Total traits evaluated: {stats['total_traits']}\n")
        f.write(f"  Traits checked: {stats['traits_checked']}\n")
        f.write(f"  Traits retained (null): {stats['traits_retained']}\n")
        f.write(f"  Traits filtered out: {stats['traits_filtered_out']}\n")
        
        if stats['traits_checked'] > 0:
            retention_rate = 100 * stats['traits_retained'] / stats['traits_checked']
            f.write(f"  Retention rate: {retention_rate:.1f}%\n")
        
        f.write(f"  Minimum p-value across all traits: {stats['min_pval_across_all']:.2e}\n")
        f.write("\n")
        
        # Write per-trait details
        f.write("TRAIT-LEVEL DETAILS:\n")
        f.write("-" * 80 + "\n")
        f.write("Trait\tChecked\tFound\tSignificant\tMin_Pval\tStatus\n")
        f.write("-" * 80 + "\n")
        
        for detail in stats['trait_details']:
            status = "FILTERED" if detail['significant_in_matrix'] or detail['significant_in_buffer'] else "RETAINED"
            f.write(f"{detail['trait_file']}\t"
                   f"{detail['variants_checked']}\t"
                   f"{detail['variants_found']}\t"
                   f"{len(detail['significant_variants'])}\t"
                   f"{detail['min_pval']:.2e}\t"
                   f"{status}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Create LD matrix from tabix-indexed GWAS summary statistics "
                   "with p-value-based filtering to exclude traits with significant associations."
    )
    
    # Required arguments
    parser.add_argument("--indexed_dir", type=str, required=True,
                       help="Directory containing tabix-indexed summary statistics files.")
    parser.add_argument("--variant_file", type=str, required=True,
                       help="File containing variant IDs to extract (one per line).")
    parser.add_argument("--output_file", type=str, required=True,
                       help="Output file for the LD matrix (will be compressed).")
    parser.add_argument("--pval_col", type=str, required=True,
                       help="Column name for p-values in summary stats.")
    parser.add_argument("--pval_threshold", type=float, required=True,
                       help="P-value threshold for significance.")
    
    # Optional arguments
    parser.add_argument("--random_seed", type=int, default=5,
                       help="Random seed for reproducibility [default: 5]")
    parser.add_argument("--variant_col", type=str, default="rs_number",
                       help="Column name for variant IDs in summary stats [default: rs_number]")
    parser.add_argument("--beta_col", type=str, default="beta",
                       help="Column name for beta values [default: beta]")
    parser.add_argument("--pos_col", type=str, default="POS",
                       help="Column name for variant positions [default: POS]. "
                            "Required if --buffer is specified.")
    parser.add_argument("--buffer", type=int, default=None,
                       help="Genomic distance (bp) for buffer region around matrix variants. "
                            "When specified, variants within this distance of the matrix boundaries "
                            "are also checked for significance. Default: None (no buffer)")
    parser.add_argument("--variant_id_format", type=str, default="CHR:POS:REF:ALT",
                       help="Format string for variant IDs [default: CHR:POS:REF:ALT]")
    parser.add_argument("--device", type=str, default="cuda", choices=["cuda", "cpu"],
                       help="Device for computation [default: cuda]")
    parser.add_argument("--save_binary_prefix", type=str, default=None,
                       help="Optional prefix for PLINK binary format output (.vcor1 and .vcor1.vars)")
    parser.add_argument("--stats_output", type=str, default=None,
                       help="Optional output file for filtering statistics and trait information")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.indexed_dir):
        logger.error(f"Error: Indexed directory does not exist: {args.indexed_dir}")
        sys.exit(1)
    
    if not os.path.exists(args.variant_file):
        logger.error(f"Error: Variant file does not exist: {args.variant_file}")
        sys.exit(1)
    
    if args.buffer and not args.pos_col:
        logger.error("Error: --pos_col must be specified when using --buffer")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(os.path.abspath(args.output_file))
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info(f"Creating LD matrix from indexed summary statistics with p-value filtering")
    logger.info(f"Input directory: {args.indexed_dir}")
    logger.info(f"Variant file: {args.variant_file}")
    logger.info(f"Output file: {args.output_file}")
    logger.info(f"P-value threshold: {args.pval_threshold}")
    logger.info(f"P-value column: {args.pval_col}")
    if args.buffer:
        logger.info(f"Genomic buffer: {args.buffer} bp")
    else:
        logger.info(f"Genomic buffer: None (matrix variants only)")
    logger.info(f"Device: {args.device}")
    
    # Step 1: Find indexed summary statistics files
    logger.info("\nStep 1: Discovering indexed summary statistics files...")
    indexed_files = find_indexed_sumstats_files(args.indexed_dir)
    
    if not indexed_files:
        logger.error(f"Error: No tabix-indexed files found in {args.indexed_dir}")
        logger.error("Expected to find pairs of .gz and .gz.tbi files.")
        sys.exit(1)
    
    logger.info(f"Found {len(indexed_files)} indexed summary statistics files")
    logger.info("Sample files:")
    for i, gz_file in enumerate(indexed_files[:5]):
        logger.info(f"  {os.path.basename(gz_file)}")
    if len(indexed_files) > 5:
        logger.info(f"  ... and {len(indexed_files) - 5} more files")
    
    # Step 2: Load variants to extract
    logger.info("\nStep 2: Loading variants to extract...")
    try:
        variants_to_extract = set(pd.read_csv(args.variant_file, header=None)[0])
    except Exception as e:
        logger.warning(f"Error loading variant file: {e}")
        try:
            with open(args.variant_file, 'r') as f:
                lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]
            variants_to_extract = set(lines)
            logger.info(f"Loaded {len(variants_to_extract)} variants using fallback method")
        except Exception as e2:
            logger.error(f"Error with fallback method: {e2}")
            sys.exit(1)
    
    if not variants_to_extract:
        logger.error(f"Error: No variants loaded from {args.variant_file}")
        sys.exit(1)
    
    logger.info(f"Loaded {len(variants_to_extract)} variants to extract")
    
    # Step 3: Create variant position mapping
    logger.info("\nStep 3: Creating variant position mapping...")
    try:
        variant_map = make_var_dict_for_positional_ids(variants_to_extract, args.variant_id_format)
        logger.info(f"Created position mapping for {len(variant_map)} variants")
        
        # Get active chromosome
        active_chromosome = list(set([x[0] for x in variant_map.values()]))
        
        if len(active_chromosome) > 1:
            logger.error("Error: Variant list spans multiple chromosomes.")
            sys.exit(1)
        
        active_chromosome = active_chromosome[0]
    except Exception as e:
        logger.error(f"Error creating variant position mapping: {e}")
        sys.exit(1)
    
    # Step 4: Define buffer region if specified
    buffer_region = None
    if args.buffer:
        logger.info("\nStep 4: Defining genomic buffer region...")
        try:
            buffer_region = get_buffer_region_variants(variant_map, args.buffer)
            logger.info(f"Buffer region defined for chromosome filtering")
        except Exception as e:
            logger.error(f"Error defining buffer region: {e}")
            sys.exit(1)
    
    # Step 5: Filter traits by p-value (use ALL traits that pass filtering)
    logger.info(f"\nStep {'5' if args.buffer else '4'}: Filtering traits by p-value...")
    logger.info(f"Checking for associations with p-value < {args.pval_threshold}")
    logger.info(f"All traits passing filter will be used for LD matrix computation")
    
    try:
        filtered_files, filtering_stats = filter_traits_by_pvalue(
            trait_files=indexed_files,
            variants_to_check=variants_to_extract,
            variant_map=variant_map,
            variant_col=args.variant_col,
            pval_col=args.pval_col,
            pval_threshold=args.pval_threshold,
            buffer_region=buffer_region
        )
        
        logger.info(f"\nFiltering summary:")
        logger.info(f"  Total traits evaluated: {filtering_stats['traits_checked']}")
        logger.info(f"  Traits retained (null): {filtering_stats['traits_retained']}")
        logger.info(f"  Traits filtered out: {filtering_stats['traits_filtered_out']}")
        
        if filtering_stats['traits_checked'] > 0:
            retention_rate = 100 * filtering_stats['traits_retained'] / filtering_stats['traits_checked']
            logger.info(f"  Retention rate: {retention_rate:.1f}%")
        
        logger.info(f"  Minimum p-value observed: {filtering_stats['min_pval_across_all']:.2e}")
        
        # Check if we have any null traits
        if len(filtered_files) == 0:
            logger.error("ERROR: No null traits available! All traits were filtered out.")
            sys.exit(1)
        
        traits_to_use = filtered_files
        logger.info(f"Using {len(traits_to_use)} null traits for LD matrix computation")
    
    except Exception as e:
        logger.error(f"Error filtering traits by p-value: {e}")
        sys.exit(1)
    
    # Step 6: Extract betas from filtered traits
    logger.info(f"\nStep {'6' if args.buffer else '5'}: Extracting beta values from null traits...")
    try:
        beta_matrix = extract_betas_batch(
            file_list=traits_to_use,
            variant_col=args.variant_col,
            beta_col=args.beta_col,
            variants_to_extract=variants_to_extract,
            variant_map=variant_map
        )
        
        logger.info(f"Extracted betas: {beta_matrix.shape[0]} variants × {beta_matrix.shape[1]} traits")
        
        # Check for missing data
        missing_fraction = beta_matrix.isna().sum().sum() / (beta_matrix.shape[0] * beta_matrix.shape[1])
        logger.info(f"Missing data fraction: {missing_fraction:.3f}")
    
    except Exception as e:
        logger.error(f"Error extracting beta values: {e}")
        sys.exit(1)
    
    # Step 7: Compute correlation matrix
    logger.info(f"\nStep {'7' if args.buffer else '6'}: Computing LD matrix using beta correlations...")
    try:
        ld_matrix = compute_beta_correlation_matrix(
            gwas_files=traits_to_use,
            variants_to_extract=variants_to_extract,
            variant_map=variant_map,
            variant_col=args.variant_col,
            beta_col=args.beta_col,
            device=args.device,
            save_prefix=args.save_binary_prefix
        )
        
        logger.info(f"Computed LD matrix: {ld_matrix.shape[0]} × {ld_matrix.shape[1]}")
    
    except Exception as e:
        logger.error(f"Error computing correlation matrix: {e}")
        sys.exit(1)
    
    # Step 8: Save matrix
    logger.info(f"\nStep {'8' if args.buffer else '7'}: Saving LD matrix...")
    try:
        variants_list = sorted(variants_to_extract)
        ld_df = pd.DataFrame(ld_matrix, index=variants_list, columns=variants_list)
        save_beta_matrix(ld_df, args.output_file, compress=False)
        logger.info(f"Saved LD matrix to: {args.output_file}")
        
        if args.save_binary_prefix:
            logger.info("Also saved binary format:")
            logger.info(f"  Matrix: {args.save_binary_prefix}.vcor1")
            logger.info(f"  Variants: {args.save_binary_prefix}.vcor1.vars")
    
    except Exception as e:
        logger.error(f"Error saving matrix: {e}")
        sys.exit(1)
    
    # Step 9: Save filtering statistics if requested
    if args.stats_output:
        logger.info(f"\nStep {'9' if args.buffer else '8'}: Saving filtering statistics...")
        try:
            write_filtering_stats(filtering_stats, args.stats_output, 
                                args.pval_threshold, args.buffer)
            logger.info(f"Saved filtering statistics to: {args.stats_output}")
        except Exception as e:
            logger.error(f"Error saving filtering statistics: {e}")
            # Don't exit, this is optional
    
    # Final summary
    logger.info(f"\n{'='*70}")
    logger.info("LD MATRIX CREATION WITH P-VALUE FILTERING - COMPLETE")
    logger.info(f"{'='*70}")
    logger.info(f"Input directory: {args.indexed_dir}")
    logger.info(f"Available files: {len(indexed_files)}")
    logger.info(f"Total traits checked: {filtering_stats['traits_checked']}")
    logger.info(f"Traits retained (null): {filtering_stats['traits_retained']}")
    logger.info(f"Traits used for LD matrix: {len(traits_to_use)}")
    logger.info(f"Variants analyzed: {len(variants_to_extract)}")
    logger.info(f"P-value threshold: {args.pval_threshold}")
    if args.buffer:
        logger.info(f"Genomic buffer: {args.buffer} bp")
    logger.info(f"Output file: {args.output_file}")
    logger.info(f"Matrix dimensions: {ld_matrix.shape[0]} × {ld_matrix.shape[1]}")
    logger.info(f"{'='*70}")


if __name__ == "__main__":
    main()
