#!/usr/bin/env python3

"""
Fine-mapping Results Summary Script

Summarizes collated fine-mapping results at four grouping levels, computing:
  - Credible-set coverage (fraction of causal variants in credible sets)
  - Credible-set precision (fraction of credible sets containing a causal variant)
  - Variant-level AUPRC (from pooled PIPs across replicates)
  - Mean credible set size (weighted by number of credible sets per replicate)
  - Mean execution time
  - Number of replicates

Processes files one at a time to avoid memory issues with large datasets.
Aggregates at the most granular grouping level per file, then derives
coarser grouping levels from the granular results.

Author: Generated for genetic fine-mapping analysis
Date: February 2026
"""

import os
import pandas as pd
import numpy as np
from collections import defaultdict
import logging
import argparse
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 100 thresholds from 0.01 to 1.00
THRESHOLDS = np.round(np.arange(0.01, 1.005, 0.01), 2)[:100]

# The most granular grouping columns (including experimental conditions from filename)
GRANULAR_COLS = [
    'meta_analysis_type', 'simulation_type', 'population',
    'snp_heritability', 'trait_type', 'causal_variant_count',
    'method', 'ld_matrix'
]

# Non-causal PIP column names
NC_PIP_COLS = [f'non_causal_pip_ge_{t:.2f}' for t in THRESHOLDS]

# Columns to aggregate by simple summation
SIMPLE_SUM_COLS = [
    'causal_variant_count', 'causal_in_credible_sets',
    'num_credible_sets', 'credible_sets_with_causal',
    'weighted_cs_size', 'exec_time_clean', 'exec_time_valid',
    'total_variants', "failed_locus"
]

#Create a function to identify if a locus failed based on the causal_variant_pips column
def is_failed_locus(pips_str: str) -> bool:
    """
    Return True if this locus failed, defined as causal_variant_pips being
    a comma-delimited list of 'nan' (e.g., 'nan', 'nan,nan', etc.).
    """
    if pd.isna(pips_str):
        return False  # treat missing as not-a-failure for safety
    
    parts = [p.strip().lower() for p in str(pips_str).split(',')]
    # Drop empty tokens
    parts = [p for p in parts if p != '']
    
    if not parts:
        return False
    
    # All entries must be 'nan' (optionally allow 'na' if you ever use that)
    return all(p in ('nan', 'na') for p in parts)


def parse_filename(filename):
    """
    Parse experimental conditions from results filename.

    Filename format examples:
    - random_homo_3e-4.results_file.txt
    - fixed_drop_causal_AFR_5e-5.results_file.txt
    - random_matched_to_EUR_3e-4.results_file.txt
    - sub-sample_LD.fixed_homo_3e-4.results_file.txt

    Returns:
        Tuple of (meta_analysis_type, simulation_type, population, snp_heritability)
    """
    # Remove .results_file.txt suffix and optional sub-sample_LD. prefix
    name = filename.replace('.results_file.txt', '').replace('sub-sample_LD.', '')

    parts = name.split('_')

    # First part is meta-analysis type: 'random' or 'fixed'
    meta_analysis_type = parts[0]

    # Second part starts the simulation type
    if parts[1] == 'homo':
        simulation_type = 'homo'
        population = 'NA'
        heritability = parts[2]  # e.g., '3e-4'
    elif parts[1] == 'matched' and parts[2] == 'to':
        simulation_type = 'matched_to'
        population = parts[3]  # e.g., 'AFR' or 'EUR'
        heritability = parts[4]  # e.g., '3e-4'
    elif parts[1] == 'drop' and parts[2] == 'causal':
        simulation_type = 'drop_causal'
        population = parts[3]  # e.g., 'AFR' or 'EUR'
        heritability = parts[4]  # e.g., '3e-4'
    else:
        logger.error(f"Could not parse filename: {filename}")
        return 'unknown', 'unknown', 'unknown', 'unknown'

    return meta_analysis_type, simulation_type, population, heritability


def make_default_accum():
    """Create a fresh accumulator dictionary for one group."""
    return {
        'sum_causal_variant_count': 0,
        'sum_causal_in_cs': 0,
        'sum_num_cs': 0,
        'sum_cs_with_causal': 0,
        'sum_weighted_cs_size': 0.0,
        'sum_exec_time': 0.0,
        'count_exec_time': 0,
        'n_replicates': 0,
        'sum_total_variants': 0,
        'total_causal_pips': 0,
        'causal_pip_ge': np.zeros(100, dtype=np.int64),
        'non_causal_pip_ge': np.zeros(100, dtype=np.int64),
        'sum_failed_loci': 0, 
    }


def process_files(input_dir):
    """
    Process all results files one at a time, accumulating aggregates
    at the most granular grouping level.

    Returns:
        Dictionary mapping group-key tuples to accumulator dicts.
    """
    global_acc = defaultdict(make_default_accum)

    files = sorted(f for f in os.listdir(input_dir) if f.endswith('.results_file.txt'))
    if not files:
        logger.error("No results files found!")
        sys.exit(1)

    for file_idx, filename in enumerate(files):
        filepath = os.path.join(input_dir, filename)
        logger.info(f"[{file_idx + 1}/{len(files)}] Processing {filename}...")

        # Parse experimental conditions from filename
        meta_type, sim_type, pop, heritability = parse_filename(filename)

        df = pd.read_csv(filepath, sep='\t')
        logger.info(f"  Loaded {len(df)} rows")

        # Add experimental condition columns
        df['meta_analysis_type'] = meta_type
        df['simulation_type'] = sim_type
        df['population'] = pop
        df['snp_heritability'] = heritability

        # Derived columns for aggregation
        df['weighted_cs_size'] = df['avg_variants_per_cs'] * df['num_credible_sets']
        df['exec_time_clean'] = df['execution_time_seconds'].fillna(0.0)
        df['exec_time_valid'] = df['execution_time_seconds'].notna().astype(int)
        df['failed_locus'] = df['causal_variant_pips'].apply(is_failed_locus).astype(int)

        # GroupBy at granular level: simple column sums
        grouped_sums = df.groupby(GRANULAR_COLS)[SIMPLE_SUM_COLS].sum()
        grouped_counts = df.groupby(GRANULAR_COLS).size()
        grouped_nc = df.groupby(GRANULAR_COLS)[NC_PIP_COLS].sum()

        # Process causal PIPs: explode to individual values with group labels
        pip_df = df[GRANULAR_COLS + ['causal_variant_pips']].copy()
        pip_df['pip'] = pip_df['causal_variant_pips'].str.split(',')
        pip_df = pip_df.explode('pip')
        pip_df['pip'] = pd.to_numeric(pip_df['pip'], errors='coerce')
        pip_df = pip_df.dropna(subset=['pip'])

        # Count total causal PIPs per group
        pip_counts = pip_df.groupby(GRANULAR_COLS)['pip'].count()

        # Compute causal PIP threshold counts per group
        logger.info("  Computing causal PIP threshold counts...")
        causal_pip_data = {}
        for group_key, group_pips_series in pip_df.groupby(GRANULAR_COLS)['pip']:
            if not isinstance(group_key, tuple):
                group_key = (group_key,)
            pips = group_pips_series.values
            # Vectorized comparison: pips[:, None] >= thresholds[None, :]
            counts = (pips[:, np.newaxis] >= THRESHOLDS[np.newaxis, :]).sum(axis=0).astype(np.int64)
            causal_pip_data[group_key] = counts

        # Merge partial results into global accumulator
        for idx in grouped_sums.index:
            key = idx if isinstance(idx, tuple) else (idx,)
            g = global_acc[key]

            row = grouped_sums.loc[idx]
            g['sum_causal_variant_count'] += int(row['causal_variant_count'])
            g['sum_causal_in_cs'] += int(row['causal_in_credible_sets'])
            g['sum_num_cs'] += int(row['num_credible_sets'])
            g['sum_cs_with_causal'] += int(row['credible_sets_with_causal'])
            g['sum_weighted_cs_size'] += float(row['weighted_cs_size'])
            g['sum_exec_time'] += float(row['exec_time_clean'])
            g['count_exec_time'] += int(row['exec_time_valid'])
            g['sum_total_variants'] += int(row['total_variants'])
            g['n_replicates'] += int(grouped_counts.loc[idx])
            g['sum_failed_loci'] += int(row['failed_locus'])

            nc_row = grouped_nc.loc[idx]
            g['non_causal_pip_ge'] += nc_row.values.astype(np.int64)

            if key in causal_pip_data:
                g['causal_pip_ge'] += causal_pip_data[key]
                g['total_causal_pips'] += int(pip_counts.loc[idx])

        logger.info(f"  Done. Running total: {len(global_acc)} unique granular groups")

    return dict(global_acc)


def aggregate_accumulators(granular_acc, group_col_indices):
    """
    Re-aggregate granular-level accumulators to a coarser grouping level.

    Args:
        granular_acc: Dict mapping granular group-key tuples to accumulators.
        group_col_indices: Indices into GRANULAR_COLS for the desired grouping.

    Returns:
        Dict mapping coarse group-key tuples to accumulators.
    """
    coarse_acc = defaultdict(make_default_accum)

    for granular_key, accum in granular_acc.items():
        coarse_key = tuple(granular_key[i] for i in group_col_indices)
        g = coarse_acc[coarse_key]

        g['sum_causal_variant_count'] += accum['sum_causal_variant_count']
        g['sum_causal_in_cs'] += accum['sum_causal_in_cs']
        g['sum_num_cs'] += accum['sum_num_cs']
        g['sum_cs_with_causal'] += accum['sum_cs_with_causal']
        g['sum_weighted_cs_size'] += accum['sum_weighted_cs_size']
        g['sum_exec_time'] += accum['sum_exec_time']
        g['count_exec_time'] += accum['count_exec_time']
        g['n_replicates'] += accum['n_replicates']
        g['sum_total_variants'] += accum['sum_total_variants']
        g['total_causal_pips'] += accum['total_causal_pips']
        g['causal_pip_ge'] += accum['causal_pip_ge']
        g['non_causal_pip_ge'] += accum['non_causal_pip_ge']
        g['sum_failed_loci'] += accum['sum_failed_loci']

    return dict(coarse_acc)


def compute_stats_from_accum(accum):
    """
    Compute final summary statistics from an accumulator.

    Returns:
        Dictionary with cs_power, cs_coverage, variant_auprc,
        mean_cs_size, mean_execution_time_seconds, n_replicates.
    """
    a = accum

    # Credible-set power
    power = (a['sum_causal_in_cs'] / a['sum_causal_variant_count']
                if a['sum_causal_variant_count'] > 0 else np.nan)

    # Credible-set precision
    coverage = (a['sum_cs_with_causal'] / a['sum_num_cs']
                 if a['sum_num_cs'] > 0 else np.nan)

    # Mean credible set size (weighted by number of CS per replicate)
    mean_cs_size = (a['sum_weighted_cs_size'] / a['sum_num_cs']
                    if a['sum_num_cs'] > 0 else np.nan)

    # Mean execution time
    mean_exec_time = (a['sum_exec_time'] / a['count_exec_time']
                      if a['count_exec_time'] > 0 else np.nan)

    # Variant-level AUPRC
    total_causal = a['total_causal_pips']

    if total_causal > 0:
        points = []

        for j in range(100):
            tp = int(a['causal_pip_ge'][j])
            fp = int(a['non_causal_pip_ge'][j])

            if tp + fp == 0:
                continue

            prec = tp / (tp + fp)
            rec = tp / total_causal
            points.append((rec, prec))

        if not points:
            auprc = 0.0
        else:
            # Sort by recall ascending (break ties by precision ascending)
            points.sort()
            sorted_recalls = [p[0] for p in points]
            sorted_precisions = [p[1] for p in points]

            # Prepend (recall=0, precision=1) if curve doesn't start at recall=0
            if sorted_recalls[0] > 0:
                sorted_recalls.insert(0, 0.0)
                sorted_precisions.insert(0, 1.0)

            # Compute AUC using trapezoidal rule
            auprc = float(np.trapz(sorted_precisions, sorted_recalls))
    else:
        auprc = np.nan

    return {
        'cs_power': round(power, 6) if pd.notna(power) else np.nan,
        'cs_coverage': round(coverage, 6) if pd.notna(coverage) else np.nan,
        'variant_auprc': round(auprc, 6) if pd.notna(auprc) else np.nan,
        'mean_cs_size': round(mean_cs_size, 3) if pd.notna(mean_cs_size) else np.nan,
        'mean_execution_time_seconds': round(mean_exec_time, 3) if pd.notna(mean_exec_time) else np.nan,
        'n_replicates': a['n_replicates'],
        'n_failed_loci': a['sum_failed_loci']
    }


def write_summary(acc, group_col_names, output_path):
    """Write summary table from accumulators to a TSV file."""
    rows = []
    for key in sorted(acc.keys()):
        row = dict(zip(group_col_names, key))
        row.update(compute_stats_from_accum(acc[key]))
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Wrote {len(df)} rows to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize fine-mapping results at multiple grouping levels",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python fine_mapping_summarizer.py \\
      --input-dir /path/to/fine-mapping \\
      --output-dir /path/to/fine-mapping_summary
        """
    )

    parser.add_argument(
        '--input-dir',
        required=True,
        help='Directory containing collated results files (*.results_file.txt)'
    )

    parser.add_argument(
        '--output-dir',
        required=True,
        help='Directory to write summary files'
    )

    args = parser.parse_args()

    # Create output directory if needed
    os.makedirs(args.output_dir, exist_ok=True)

    # Process all files incrementally at the most granular level
    logger.info("Phase 1: Processing files and accumulating granular-level aggregates...")
    granular_acc = process_files(args.input_dir)
    logger.info(f"Phase 1 complete. {len(granular_acc)} unique granular groups.")

    # Define grouping levels with indices into GRANULAR_COLS
    # GRANULAR_COLS = ['meta_analysis_type', 'simulation_type', 'population',
    #                   'snp_heritability', 'trait_type', 'causal_variant_count',
    #                   'method', 'ld_matrix']
    #                    index 0               index 1            index 2
    #                    index 3               index 4            index 5
    #                    index 6               index 7
    levels = [
        (
            [0, 1, 2, 3, 4, 5, 6, 7],
            ['meta_analysis_type', 'simulation_type', 'population', 'snp_heritability',
             'trait_type', 'causal_variant_count', 'method', 'ld_matrix'],
            'summary_by_trait_causal_method_ld.txt'
        ),
        (
            [0, 1, 2, 3, 4, 6, 7],
            ['meta_analysis_type', 'simulation_type', 'population', 'snp_heritability',
             'trait_type', 'method', 'ld_matrix'],
            'summary_by_trait_method_ld.txt'
        ),
        (
            [0, 1, 2, 3, 5, 6, 7],
            ['meta_analysis_type', 'simulation_type', 'population', 'snp_heritability',
             'causal_variant_count', 'method', 'ld_matrix'],
            'summary_by_causal_method_ld.txt'
        ),
        (
            [0, 1, 2, 3, 6, 7],
            ['meta_analysis_type', 'simulation_type', 'population', 'snp_heritability',
             'method', 'ld_matrix'],
            'summary_by_method_ld.txt'
        ),
    ]

    logger.info("Phase 2: Computing summaries at each grouping level...")
    for col_indices, col_names, filename in levels:
        output_path = os.path.join(args.output_dir, filename)

        if col_indices == [0, 1, 2, 3, 4, 5, 6, 7]:
            # Already at the granular level
            acc = granular_acc
        else:
            acc = aggregate_accumulators(granular_acc, col_indices)

        logger.info(f"Writing {filename} ({len(acc)} groups)...")
        write_summary(acc, col_names, output_path)

    logger.info("All summaries complete!")


if __name__ == "__main__":
    main()
