#!/usr/bin/env python3
"""
LD Matrix Combiner Script
Combines multiple LD matrices from plink2 --r-unphased square output using sample size weighting.
Designed for OLCF Frontier HPC system with multi-core and GPU support.

Usage:
    python combine_ld_matrices.py --ld-files file1.unphased.vcor1 file2.unphased.vcor1 ... \
                                  --sample-sizes N1 N2 ... \
                                  --output-prefix output_path \
                                  --threads 32
"""

import argparse
import numpy as np
import pandas as pd
import os
import sys
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial
import logging
from typing import List, Tuple, Dict, Optional
import time

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class LDMatrixCombiner:
    def __init__(self, ld_files: List[str], sample_sizes: List[int], 
                 output_prefix: str, n_threads: int = None):
        """
        Initialize the LD matrix combiner.

        Args:
            ld_files: List of paths to .unphased.vcor1 files
            sample_sizes: List of sample sizes corresponding to each LD file
            output_prefix: Output path prefix for combined matrix
            n_threads: Number of threads to use (default: all available cores)
        """
        self.ld_files = ld_files
        self.sample_sizes = np.array(sample_sizes)
        self.output_prefix = output_prefix
        self.n_threads = n_threads or cpu_count()

        # Validate inputs
        if len(ld_files) != len(sample_sizes):
            raise ValueError("Number of LD files must match number of sample sizes")
        if len(ld_files) < 2:
            raise ValueError("At least 2 LD matrices are required for combination")

        logger.info(f"Initialized combiner with {len(ld_files)} matrices")
        logger.info(f"Sample sizes: {sample_sizes}")
        logger.info(f"Using {self.n_threads} threads")

    def load_variants(self, vars_file: str) -> pd.DataFrame:
        """Load variant information from .unphased.vcor1.vars file."""
        try:
            # The .vars file typically has columns: CHROM, POS, ID, REF, ALT
            variants = pd.read_csv(vars_file, sep='\t', dtype=str, header=None)
            logger.info(f"Loaded {len(variants)} variants from {vars_file}")
            return variants
        except Exception as e:
            logger.error(f"Error loading variants from {vars_file}: {e}")
            raise

    def load_ld_matrix(self, ld_file: str) -> np.ndarray:
        """Load LD matrix from .unphased.vcor1 file."""
        try:
            # Load the square correlation matrix
            matrix = np.loadtxt(ld_file, dtype=np.float32)
            logger.info(f"Loaded {matrix.shape[0]}x{matrix.shape[1]} LD matrix from {ld_file}")

            # Verify it's square
            if matrix.shape[0] != matrix.shape[1]:
                raise ValueError(f"LD matrix is not square: {matrix.shape}")

            return matrix
        except Exception as e:
            logger.error(f"Error loading LD matrix from {ld_file}: {e}")
            raise

    def create_variant_union(self) -> Tuple[pd.DataFrame, Dict[str, int]]:
        """Create union of all variants across matrices."""
        logger.info("Creating union of variants across all matrices...")

        all_variants = []
        matrix_variant_sets = []

        for i, ld_file in enumerate(self.ld_files):
            vars_file = ld_file + '.vars'
            if not os.path.exists(vars_file):
                raise FileNotFoundError(f"Variants file not found: {vars_file}")

            variants = self.load_variants(vars_file)
            # Create unique identifier for each variant
            variants['variant_id'] = variants[:,0]
            all_variants.append(variants)
            matrix_variant_sets.append(set(variants['variant_id']))

        # Create union of all variant IDs
        union_variants = set()
        for variant_set in matrix_variant_sets:
            union_variants.update(variant_set)

        # Create master variant list maintaining order
        master_variants = []
        seen = set()
        for variants_df in all_variants:
            for _, variant in variants_df.iterrows():
                if variant['variant_id'] not in seen:
                    master_variants.append(variant)
                    seen.add(variant['variant_id'])

        # Add any remaining variants
        for variants_df in all_variants:
            for _, variant in variants_df.iterrows():
                if variant['variant_id'] in union_variants and variant['variant_id'] not in seen:
                    master_variants.append(variant)
                    seen.add(variant['variant_id'])

        master_variants_df = pd.DataFrame(master_variants)
        variant_to_index = {vid: i for i, vid in enumerate(master_variants_df['variant_id'])}

        logger.info(f"Created union with {len(master_variants_df)} variants")
        return master_variants_df, variant_to_index

    def process_matrix_chunk(self, args: Tuple) -> Tuple[int, int, float, int]:
        """Process a chunk of the combined matrix calculation."""
        i, j, master_variants, variant_to_index, matrix_data = args

        variant_i = master_variants.iloc[i]['variant_id']
        variant_j = master_variants.iloc[j]['variant_id']

        weighted_sum = 0.0
        total_weight = 0.0

        for mat_idx, (ld_matrix, variants_df, sample_size) in enumerate(matrix_data):
            variant_ids = variants_df['variant_id'].tolist()

            # Check if both variants exist in this matrix
            if variant_i in variant_ids and variant_j in variant_ids:
                idx_i = variant_ids.index(variant_i)
                idx_j = variant_ids.index(variant_j)

                correlation = ld_matrix[idx_i, idx_j]
                weight = sample_size

                weighted_sum += correlation * weight
                total_weight += weight

        if total_weight > 0:
            combined_correlation = weighted_sum / total_weight
            return i, j, combined_correlation, int(total_weight)
        else:
            return i, j, 0.0, 0

    def combine_matrices_parallel(self) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame]:
        """Combine LD matrices using parallel processing."""
        logger.info("Loading all matrices and variant information...")

        # Load all matrices and variants
        matrix_data = []
        for i, ld_file in enumerate(self.ld_files):
            vars_file = ld_file + '.vars'
            ld_matrix = self.load_ld_matrix(ld_file)
            variants_df = self.load_variants(vars_file)
            variants_df['variant_id'] = variants_df[:,0]

            matrix_data.append((ld_matrix, variants_df, self.sample_sizes[i]))

        # Create variant union
        master_variants, variant_to_index = self.create_variant_union()
        n_variants = len(master_variants)

        logger.info(f"Computing combined {n_variants}x{n_variants} LD matrix...")

        # Initialize output matrices
        combined_matrix = np.zeros((n_variants, n_variants), dtype=np.float32)
        weight_matrix = np.zeros((n_variants, n_variants), dtype=np.int32)

        # Create arguments for parallel processing
        # Only compute upper triangle (including diagonal) since LD matrix is symmetric
        args_list = []
        for i in range(n_variants):
            for j in range(i, n_variants):
                args_list.append((i, j, master_variants, variant_to_index, matrix_data))

        logger.info(f"Processing {len(args_list)} matrix elements using {self.n_threads} threads...")

        # Process in parallel
        start_time = time.time()
        with Pool(self.n_threads) as pool:
            results = pool.map(self.process_matrix_chunk, args_list)

        # Fill in results
        for i, j, correlation, weight in results:
            combined_matrix[i, j] = correlation
            combined_matrix[j, i] = correlation  # Symmetric
            weight_matrix[i, j] = weight
            weight_matrix[j, i] = weight

        elapsed_time = time.time() - start_time
        logger.info(f"Matrix combination completed in {elapsed_time:.2f} seconds")

        return combined_matrix, weight_matrix, master_variants

    def save_results(self, combined_matrix: np.ndarray, master_variants: pd.DataFrame):
        """Save the combined LD matrix and variants file."""
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(self.output_prefix)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Save LD matrix
        matrix_file = self.output_prefix + '.unphased.vcor1'
        logger.info(f"Saving combined LD matrix to {matrix_file}")
        np.savetxt(matrix_file, combined_matrix, fmt='%.6f', delimiter='\t')

        # Save variants file
        vars_file = matrix_file + '.vars'
        logger.info(f"Saving variants information to {vars_file}")
        # Remove the variant_id column before saving
        output_variants = master_variants.drop('variant_id', axis=1)
        output_variants.to_csv(vars_file, sep='\t', index=False)

        logger.info(f"Successfully saved combined matrix with {len(master_variants)} variants")

    def run(self):
        """Execute the complete LD matrix combination workflow."""
        logger.info("Starting LD matrix combination workflow...")
        start_time = time.time()

        try:
            # Combine matrices
            combined_matrix, weight_matrix, master_variants = self.combine_matrices_parallel()

            # Save results
            self.save_results(combined_matrix, master_variants)

            total_time = time.time() - start_time
            logger.info(f"Workflow completed successfully in {total_time:.2f} seconds")

            # Print summary statistics
            logger.info("\n=== SUMMARY ===")
            logger.info(f"Input matrices: {len(self.ld_files)}")
            logger.info(f"Sample sizes: {self.sample_sizes.tolist()}")
            logger.info(f"Total variants in union: {len(master_variants)}")
            logger.info(f"Output prefix: {self.output_prefix}")

        except Exception as e:
            logger.error(f"Workflow failed: {e}")
            raise

def main():
    parser = argparse.ArgumentParser(
        description='Combine multiple LD matrices using sample size weighting',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--ld-files', nargs='+', required=True,
                       help='List of .unphased.vcor1 LD matrix files')
    parser.add_argument('--sample-sizes', nargs='+', type=int, required=True,
                       help='Sample sizes corresponding to each LD matrix')
    parser.add_argument('--output-prefix', required=True,
                       help='Output path prefix for combined matrix')
    parser.add_argument('--threads', type=int, default=None,
                       help='Number of threads to use (default: all available)')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate inputs
    if len(args.ld_files) != len(args.sample_sizes):
        parser.error("Number of LD files must match number of sample sizes")

    # Check that all input files exist
    for ld_file in args.ld_files:
        if not os.path.exists(ld_file):
            parser.error(f"LD matrix file not found: {ld_file}")
        vars_file = ld_file + '.vars'
        if not os.path.exists(vars_file):
            parser.error(f"Variants file not found: {vars_file}")

    # Initialize and run combiner
    combiner = LDMatrixCombiner(
        ld_files=args.ld_files,
        sample_sizes=args.sample_sizes,
        output_prefix=args.output_prefix,
        n_threads=args.threads
    )

    combiner.run()

if __name__ == '__main__':
    main()
