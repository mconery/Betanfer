#!/usr/bin/env python3


"""
Fine-mapping Results Collation Script with Timing and Trait Type Extraction

This script collates results from genetic fine-mapping experiments using simulated datasets.
It processes results from multiple fine-mapping methods and generates a summary file with
key statistics for each experiment, including execution times and trait types.

Author: Generated for genetic fine-mapping analysis
Date: September 2025
"""

import os
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
from pathlib import Path
import re
from typing import Dict, List, Tuple, Optional
import logging
import argparse
import sys

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_filename_info(filename: str) -> Tuple[Optional[int], Optional[int], Optional[str]]:
    """
    Parse phenotype number, causal variant count, and trait type from filename.

    Args:
        filename: Name of the file (e.g., "Pheno1.10_causal.continuous_0.25.filler.txt")
    Returns:
        Tuple of (phenotype_number, causal_count, trait_type) or (None, None, None) if parsing fails

    """

    try:

        # Extract Pheno number and causal variant count using regex

        pattern = r'Pheno(\d+)\.(\d+)_causal'

        match = re.search(pattern, filename)

        if match:

            pheno_num = int(match.group(1))

            causal_count = int(match.group(2))

        else:

            logger.warning(f"Could not parse filename: {filename}")

            return None, None, None

        # Extract trait type (continuous or binary)

        trait_pattern = r'_causal\.(continuous|binary)'

        trait_match = re.search(trait_pattern, filename)

        if trait_match:

            trait_type = trait_match.group(1)

        else:

            logger.warning(f"Could not extract trait type from filename: {filename}")

            trait_type = 'unknown'

        return pheno_num, causal_count, trait_type

    except Exception as e:

        logger.error(f"Error parsing filename {filename}: {e}")

        return None, None, None

def load_causal_variants(causal_variants_dir: str, pheno_num: int, causal_count: int) -> set:

    """

    Load the set of causal variants for a given phenotype.

    Args:

        causal_variants_dir: Directory containing causal variant files

        pheno_num: Phenotype number

        causal_count: Number of causal variants

    Returns:

        Set of causal variant IDs

    """

    causal_file = f"Pheno{pheno_num}.{causal_count}_causal.variant_set.txt"

    causal_path = os.path.join(causal_variants_dir, causal_file)

    try:

        with open(causal_path, 'r') as f:

            causal_variants = set(line.strip() for line in f if line.strip())

        logger.debug(f"Loaded {len(causal_variants)} causal variants from {causal_file}")

        return causal_variants

    except FileNotFoundError:

        logger.error(f"Causal variants file not found: {causal_path}")

        return set()

    except Exception as e:

        logger.error(f"Error loading causal variants from {causal_path}: {e}")

        return set()

def extract_execution_time(result_file_path: str) -> Optional[float]:

    """

    Extract execution time from the corresponding timing file.

    Args:

        result_file_path: Path to the result file

    Returns:

        Execution time in seconds, or None if timing file not found or invalid

    """

    try:

        # Get the directory and filename of the result file

        result_dir = os.path.dirname(result_file_path)

        result_filename = os.path.basename(result_file_path)

        # Construct timing file path: same directory + timing subdirectory + same filename

        timing_dir = os.path.join(result_dir, "timing")

        timing_file_path = os.path.join(timing_dir, result_filename)

        # Check if timing file exists

        if not os.path.exists(timing_file_path):

            logger.warning(f"Timing file not found: {timing_file_path}")

            return None

        # Read the timing file

        with open(timing_file_path, 'r') as f:

            lines = f.readlines()

        # Check if we have at least 2 lines

        if len(lines) < 2:

            logger.warning(f"Timing file has fewer than 2 lines: {timing_file_path}")

            return None

        # Extract execution time from the second line

        execution_time_str = lines[1].strip()

        # Convert to float

        execution_time = float(execution_time_str)

        logger.debug(f"Extracted execution time: {execution_time}s from {timing_file_path}")

        return execution_time

    except FileNotFoundError:

        logger.warning(f"Timing file not found for result file: {result_file_path}")

        return None

    except ValueError as e:

        logger.error(f"Could not parse execution time from {timing_file_path}: {e}")

        return None

    except Exception as e:

        logger.error(f"Error extracting execution time for {result_file_path}: {e}")

        return None

def process_single_result_file(args: Tuple[str, str]) -> Optional[Dict]:

    """

    Process a single fine-mapping result file and extract summary statistics.

    Args:

        args: Tuple of (result_file_path, causal_variants_dir)

    Returns:

        Dictionary with summary statistics or None if processing fails

    """

    result_file_path, causal_variants_dir = args

    try:

        # Parse filename to get phenotype, causal variant count, and trait type

        filename = os.path.basename(result_file_path)

        pheno_num, causal_count, trait_type = parse_filename_info(filename)

        if pheno_num is None or causal_count is None:

            logger.warning(f"Skipping file with unparseable name: {filename}")

            return None

        # Load causal variants

        causal_variants = load_causal_variants(causal_variants_dir, pheno_num, causal_count)

        if not causal_variants:

            logger.warning(f"No causal variants found for {filename}")

            return None

        # Extract execution time

        execution_time = extract_execution_time(result_file_path)

        # Read the results file

        try:

            df = pd.read_csv(result_file_path, sep='\t')

        except pd.errors.EmptyDataError:

            logger.warning(f"Empty file (no credible sets found): {result_file_path}")

            # Handle empty file: treat as no credible sets, NA timing

            # Get method and LD matrix info from file path

            path_parts = result_file_path.split(os.sep)

            method = None

            ld_matrix = None

            for i, part in enumerate(path_parts):

                if part.startswith("fine-mapping"):

                    continue

                elif part.endswith('_results'):

                    method = part.replace('_results', '')

                elif method and i < len(path_parts) - 1:

                    if not part.endswith('.txt'):

                        ld_matrix = part

                        break

            # Create result dict for empty file (no credible sets found)

            non_causal_pip_counts = {}

            for threshold_int in range(1, 101):

                threshold = threshold_int / 100.0

                non_causal_pip_counts[f'non_causal_pip_ge_{threshold:.2f}'] = 0

            result_dict = {

                'phenotype_number': pheno_num,

                'trait_type': trait_type,

                'causal_variant_count': causal_count,

                'causal_in_credible_sets': 0,

                'num_credible_sets': 0,

                'credible_sets_with_causal': 0,

                'avg_variants_per_cs': 0,

                'causal_variant_pips': ','.join(['0.0'] * len(causal_variants)) if causal_variants else '0.0',

                'method': method or 'unknown',

                'ld_matrix': ld_matrix or 'unknown',

                'result_file': result_file_path,

                'total_variants': 0,

                'execution_time_seconds': execution_time

            }

            result_dict.update(non_causal_pip_counts)

            logger.info(f"Processed {filename}: empty file (no credible sets), execution time: {execution_time}s")

            return result_dict

        except Exception as e:

            logger.error(f"Error reading file {result_file_path}: {e}")

            return None

        # Check for required columns

        required_cols = ['VARIANT_ID', 'PIP', 'cs']

        missing_cols = [col for col in required_cols if col not in df.columns]

        if missing_cols:

            logger.error(f"Missing required columns {missing_cols} in {result_file_path}")

            logger.error(f"Available columns: {list(df.columns)}")

            return None

        # Filter for variants in credible sets (cs > 0)

        cs_variants = df[df['cs'] > 0].copy()

        # Count causal variants in credible sets

        causal_in_cs = sum(1 for variant in cs_variants['VARIANT_ID'] if variant in causal_variants)

        # Count number of non-zero credible sets mapped at this locus

        num_credible_sets = len(cs_variants['cs'].unique()) if len(cs_variants) > 0 else 0

        # Count number of credible sets that contain at least one true causal variant

        credible_sets_with_causal = 0

        if len(cs_variants) > 0:

            for cs_id in cs_variants['cs'].unique():

                cs_variants_in_set = cs_variants[cs_variants['cs'] == cs_id]['VARIANT_ID']

                if any(variant in causal_variants for variant in cs_variants_in_set):

                    credible_sets_with_causal += 1

        # Calculate average number of variants per credible set

        if len(cs_variants) > 0:

            cs_counts = cs_variants['cs'].value_counts()

            avg_variants_per_cs = cs_counts.mean()

        else:

            avg_variants_per_cs = 0

        # Get PIPs of all causal variants (including those not in credible sets)

        causal_pips = []

        for variant in causal_variants:

            variant_rows = df[df['VARIANT_ID'] == variant]

            if not variant_rows.empty:

                causal_pips.append(variant_rows['PIP'].iloc[0])

            else:

                # If causal variant not found in results, assume PIP = 0

                causal_pips.append(0.0)

        # Calculate total number of variants at this locus

        total_variants = len(df)

        # Get non-causal variants and their PIPs

        non_causal_variants = df[~df['VARIANT_ID'].isin(causal_variants)].copy()

        # Count non-causal variants with PIPs >= each threshold (0.01 to 1.00)

        non_causal_pip_counts = {}

        for threshold_int in range(1, 101):  # 1 to 100 (representing 0.01 to 1.00)

            threshold = threshold_int / 100.0

            count = sum(1 for pip in non_causal_variants['PIP'] if pip >= threshold)

            non_causal_pip_counts[f'non_causal_pip_ge_{threshold:.2f}'] = count

        # Get method and LD matrix info from file path

        path_parts = result_file_path.split(os.sep)

        method = None

        ld_matrix = None

        for i, part in enumerate(path_parts):

            if part.startswith("fine-mapping"):

                continue

            elif part.endswith('_results'):

                method = part.replace('_results', '')

            elif method and i < len(path_parts) - 1:  # LD matrix folder comes after method folder

                if not part.endswith('.txt'):

                    ld_matrix = part

                    break

        result_dict = {

            'phenotype_number': pheno_num,

            'trait_type': trait_type,

            'causal_variant_count': causal_count,

            'causal_in_credible_sets': causal_in_cs,

            'num_credible_sets': num_credible_sets,

            'credible_sets_with_causal': credible_sets_with_causal,

            'avg_variants_per_cs': round(avg_variants_per_cs, 3),

            'causal_variant_pips': ','.join(map(str, causal_pips)),

            'method': method or 'unknown',

            'ld_matrix': ld_matrix or 'unknown',

            'result_file': result_file_path,

            'total_variants': total_variants,

            'execution_time_seconds': execution_time

        }

        # Add the non-causal PIP count columns

        result_dict.update(non_causal_pip_counts)

        logger.debug(f"Processed {filename}: trait={trait_type}, {causal_in_cs}/{causal_count} causal variants in CS, "

        f"{num_credible_sets} total credible sets, {credible_sets_with_causal} contain causal variants, "

        f"{total_variants} total variants, execution time: {execution_time}s")

        return result_dict

    except Exception as e:

        logger.error(f"Error processing {result_file_path}: {e}")

        return None

def find_all_result_files(results_base_dir: str, methods: List[str] = None) -> List[str]:

    """

    Recursively find all result files in the directory structure.

    Args:

        results_base_dir: Base directory containing method result folders

        methods: List of method names to search for (default: predefined list)

    Returns:

        List of paths to result files

    """

    if methods is None:

        methods = ['susie', 'carma', 'rsparsepro', 'susiex', 'mesusie', 'multisusie']

    result_files = []

    for method in methods:

        method_dir = os.path.join(results_base_dir, f"{method}_results")

        if os.path.exists(method_dir):

            logger.info(f"Searching in method directory: {method_dir}")

            # Find all .txt files in subdirectories

            for root, dirs, files in os.walk(method_dir):

                if root.split("/")[-2] == f"{method}_results":

                    for file in files:

                        if file.endswith('.txt') and 'Pheno' in file and '_causal' in file:

                            result_files.append(os.path.join(root, file))

        else:

            logger.warning(f"Method directory not found: {method_dir}")

    logger.info(f"Found {len(result_files)} result files")

    return result_files

def validate_directory_structure(causal_variants_dir: str, results_base_dir: str) -> bool:

    """

    Validate that the required directory structure exists.

    Args:

        causal_variants_dir: Directory containing causal variant files

        results_base_dir: Base directory containing fine-mapping results

    Returns:

        True if structure is valid, False otherwise

    """

    valid = True

    # Check causal variants directory

    if not os.path.exists(causal_variants_dir):

        logger.error(f"Causal variants directory does not exist: {causal_variants_dir}")

        valid = False

    else:

        # Check for causal variant files

        causal_files = [f for f in os.listdir(causal_variants_dir)

        if f.endswith('.variant_set.txt') and 'Pheno' in f]

        if not causal_files:

            logger.warning(f"No causal variant files found in {causal_variants_dir}")

        else:

            logger.info(f"Found {len(causal_files)} causal variant files")

    # Check results base directory

    if not os.path.exists(results_base_dir):

        logger.error(f"Results base directory does not exist: {results_base_dir}")

        valid = False

    else:

        # Check for method directories

        methods = ['susie', 'carma', 'rsparsepro', 'susiex', 'mesusie', 'multisusie']

        found_methods = []

        for method in methods:

            method_dir = os.path.join(results_base_dir, f"{method}_results")

            if os.path.exists(method_dir):

                found_methods.append(method)

        if not found_methods:

            logger.error("No method result directories found")

            valid = False

        else:

            logger.info(f"Found method directories for: {found_methods}")

    return valid

def collate_fine_mapping_results(

    causal_variants_dir: str,

    results_base_dir: str,

    output_file: str,

    n_processes: Optional[int] = None,

    methods: Optional[List[str]] = None

) -> None:

    """

    Main function to collate fine-mapping results from multiple experiments.

    Args:

        causal_variants_dir: Directory containing causal variant files

        results_base_dir: Base directory containing fine-mapping results

        output_file: Path to output summary file

        n_processes: Number of processes for parallel processing (default: CPU count)

        methods: List of method names to process (default: all standard methods)

    """

    if n_processes is None:

        n_processes = cpu_count()

    logger.info(f"Starting collation with {n_processes} processes")

    # Validate directory structure

    if not validate_directory_structure(causal_variants_dir, results_base_dir):

        logger.error("Directory structure validation failed. Please check your paths.")

        return

    # Find all result files

    result_files = find_all_result_files(results_base_dir, methods)

    if not result_files:

        logger.error("No result files found!")

        return

    # Prepare arguments for parallel processing

    process_args = [(file, causal_variants_dir) for file in result_files]

    # Process files in parallel

    logger.info(f"Processing {len(result_files)} files...")

    with Pool(processes=n_processes) as pool:

        results = pool.map(process_single_result_file, process_args)

    # Filter out None results and convert to DataFrame

    valid_results = [r for r in results if r is not None]

    if not valid_results:

        logger.error("No valid results processed!")

        return

    logger.info(f"Successfully processed {len(valid_results)}/{len(result_files)} files")

    # Create summary DataFrame

    summary_df = pd.DataFrame(valid_results)

    # Reorder columns for the output

    base_columns = [

        'phenotype_number',

        'trait_type',

        'causal_variant_count',

        'causal_in_credible_sets',

        'num_credible_sets',

        'credible_sets_with_causal',

        'avg_variants_per_cs',

        'causal_variant_pips',

        'method',

        'ld_matrix',

        'total_variants',

        'execution_time_seconds'

    ]

    # Add the non-causal PIP count columns in order

    pip_columns = [f'non_causal_pip_ge_{threshold:.2f}' for threshold in np.arange(0.01, 1.01, 0.01)]

    output_columns = base_columns + pip_columns

    summary_df = summary_df[output_columns]

    # Sort by phenotype number, trait type, causal variant count, then by method

    summary_df = summary_df.sort_values(['phenotype_number', 'trait_type', 'causal_variant_count', 'method'])

    # Save to file

    summary_df.to_csv(output_file, sep='\t', index=False)

    logger.info(f"Summary results saved to: {output_file}")

    # Print some basic statistics

    logger.info("Summary statistics:")

    logger.info(f" Total experiments processed: {len(summary_df)}")

    logger.info(f" Unique phenotypes: {summary_df['phenotype_number'].nunique()}")

    logger.info(f" Trait types: {list(summary_df['trait_type'].unique())}")

    logger.info(f" Methods: {list(summary_df['method'].unique())}")

    logger.info(f" LD matrices: {list(summary_df['ld_matrix'].unique())}")

    # Print execution time statistics

    execution_times = summary_df['execution_time_seconds'].dropna()

    if len(execution_times) > 0:

        logger.info(" Execution time statistics:")

        logger.info(f" Mean: {execution_times.mean():.2f}s")

        logger.info(f" Median: {execution_times.median():.2f}s")

        logger.info(f" Min: {execution_times.min():.2f}s")

        logger.info(f" Max: {execution_times.max():.2f}s")

        logger.info(f" Missing timing data: {summary_df['execution_time_seconds'].isna().sum()} files")

    else:

        logger.warning("No execution timing data found")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(

        description="Collate fine-mapping results from multiple experiments"

    )

    parser.add_argument(

        "causal_variants_dir",

        help="Directory containing causal variant files"

    )

    parser.add_argument(

        "results_base_dir",

        help="Base directory containing fine-mapping results"

    )

    parser.add_argument(

        "output_file",

        help="Path to output summary file"

    )

    parser.add_argument(

        "--methods",

        nargs="+",

        default=['susie', 'carma', 'rsparsepro', 'susiex', 'mesusie', 'multisusie'],

        help="List of methods to process (default: all standard methods)"

    )

    parser.add_argument(

        "--processes",

        type=int,

        default=None,

        help="Number of processes for parallel processing (default: CPU count)"

    )

    args = parser.parse_args()

    collate_fine_mapping_results(

        args.causal_variants_dir,

        args.results_base_dir,

        args.output_file,

        n_processes=args.processes,

        methods=args.methods

    )
