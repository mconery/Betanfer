'''
collate_noisy_results.py

This script analyzes the results of genomic fine-mapping experiments using LD matrices
with varying levels of noise. It computes key metrics including coverage, power, 
credible set statistics, and posterior inclusion probabilities (PIPs) of causal variants.

***This needs to be run in the rsparsepro conda environment ***
'''

import argparse
import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
import multiprocessing as mp
from functools import partial
import glob
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze genomic fine-mapping experiment results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--loci-file', '-l', 
                       type=str, required=True,
                       help='Path to file containing list of loci names')
    
    parser.add_argument('--input-dir', '-i', 
                       type=str, required=True,
                       help='Directory containing subdirectories for each locus')
    
    parser.add_argument('--output-dir', '-o', 
                       type=str, required=True,
                       help='Directory containing fine-mapping output files')
    
    parser.add_argument('--results-file', '-r', 
                       type=str, required=True,
                       help='Path for output tab-delimited results file')
    
    parser.add_argument('--num-cores', '-n', 
                       type=int, default=mp.cpu_count(),
                       help='Number of CPU cores to use for parallel processing')
    
    return parser.parse_args()

def read_loci_list(loci_file):
    """Read the list of loci from the specified file."""
    try:
        with open(loci_file, 'r') as f:
            loci = [line.strip() for line in f if line.strip()]
        logger.info(f"Read {len(loci)} loci from {loci_file}")
        return loci
    except FileNotFoundError:
        logger.error(f"Loci file not found: {loci_file}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error reading loci file: {e}")
        sys.exit(1)

def read_causal_variants(locus_name, input_dir):
    """Read causal variants for a given locus."""
    variants_file = Path(input_dir) / locus_name / f"{locus_name}.variants.txt"
    
    try:
        with open(variants_file, 'r') as f:
            causal_variants = set(line.strip() for line in f if line.strip())
        logger.debug(f"Read {len(causal_variants)} causal variants for {locus_name}")
        return causal_variants
    except FileNotFoundError:
        logger.error(f"Causal variants file not found: {variants_file}")
        return set()
    except Exception as e:
        logger.error(f"Error reading causal variants for {locus_name}: {e}")
        return set()

def find_fine_mapping_files(locus_name, output_dir):
    """Find all fine-mapping output files for a given locus."""
    output_path = Path(output_dir)
    
    # Pattern for files with LD mismatch degree
    pattern1 = f"{locus_name}.*.rsparsepro.txt"
    # Pattern for files without LD mismatch degree (unmodified LD)
    pattern2 = f"{locus_name}.rsparsepro.txt"
    
    files = []
    files.extend(glob.glob(str(output_path / pattern1)))
    files.extend(glob.glob(str(output_path / pattern2)))
    
    # Remove duplicates and sort
    files = sorted(list(set(files)))
    
    logger.debug(f"Found {len(files)} fine-mapping files for {locus_name}")
    return files

def read_fine_mapping_results(file_path):
    """Read fine-mapping results from a file."""
    try:
        # Expected columns: RSID, Z, PIP, z_estimated, cs
        df = pd.read_csv(file_path, sep='\t')
        
        # Validate required columns
        required_cols = ['RSID', 'Z', 'PIP', 'z_estimated', 'cs']
        if not all(col in df.columns for col in required_cols):
            logger.error(f"Missing required columns in {file_path}")
            logger.error(f"Expected: {required_cols}, Found: {list(df.columns)}")
            return None
            
        return df
    except Exception as e:
        logger.error(f"Error reading fine-mapping results from {file_path}: {e}")
        return None

def extract_ld_mismatch_degree(file_path, locus_name):
    """Extract LD mismatch degree from filename."""
    filename = Path(file_path).name
    
    # Remove locus name and .rsparsepro.txt extension
    remaining = filename.replace(f"{locus_name}.", "").replace("rsparsepro.txt", "").replace(".","")
    
    # If nothing remains, it's the unmodified LD case
    if remaining == "":
        return "unmodified"
    else:
        return remaining

def compute_metrics(locus_name, causal_variants, fine_mapping_df, ld_mismatch_degree):
    """Compute fine-mapping metrics for a single analysis."""
    
    # Number of true causal variants
    num_causal = len(causal_variants)
    
    # Get credible sets (non-zero cs values)
    credible_sets = fine_mapping_df[fine_mapping_df['cs'] != 0]['cs'].unique()
    num_credible_sets = len(credible_sets)
    
    # Initialize metrics
    coverage_count = 0
    power_count = 0
    credible_set_sizes = []
    causal_pips = []
    
    # Analyze each credible set
    for cs_id in credible_sets:
        cs_variants = set(fine_mapping_df[fine_mapping_df['cs'] == cs_id]['RSID'])
        credible_set_sizes.append(len(cs_variants))
        
        # Check if this credible set contains any causal variant (coverage)
        if cs_variants & causal_variants:
            coverage_count += 1
    
    # Check power: how many causal variants are in any credible set
    all_cs_variants = set(fine_mapping_df[fine_mapping_df['cs'] != 0]['RSID'])
    power_count = len(causal_variants & all_cs_variants)
    
    # Get PIPs for causal variants
    for variant in causal_variants:
        variant_row = fine_mapping_df[fine_mapping_df['RSID'] == variant]
        if not variant_row.empty:
            causal_pips.append(variant_row['PIP'].iloc[0])
    
    # Calculate final metrics
    coverage = coverage_count / num_credible_sets if num_credible_sets > 0 else 0
    power = power_count / num_causal if num_causal > 0 else 0
    
    # Format multi-value fields with commas
    credible_set_sizes_str = ','.join(map(str, credible_set_sizes)) if credible_set_sizes else ""
    causal_pips_str = ','.join([f"{pip:.6f}" for pip in causal_pips]) if causal_pips else ""
    
    # New metrics calculation
    # Calculate top_power (causal variant with max |Z| in credible sets)
    top_power = 0
    if causal_variants:
        causal_df = fine_mapping_df[fine_mapping_df['RSID'].isin(causal_variants)]
        if not causal_df.empty:
            top_causal = causal_df.loc[causal_df['Z'].abs().idxmax()]
            top_power = 1 if top_causal['cs'] != 0 else 0
    
    # Calculate top_coverage (credible set with max |Z| variant contains causal)
    top_coverage = 0
    credible_df = fine_mapping_df[fine_mapping_df['cs'] != 0]
    if not credible_df.empty:
        top_variant = credible_df.loc[credible_df['Z'].abs().idxmax()]
        cs_id = top_variant['cs']
        cs_variants = set(fine_mapping_df[fine_mapping_df['cs'] == cs_id]['RSID'])
        top_coverage = 1 if cs_variants & causal_variants else 0
    
    # Return statement
    return {
        'locus': locus_name,
        'ld_mismatch_degree': ld_mismatch_degree,
        'num_causal_variants': num_causal,
        'num_credible_sets': num_credible_sets,
        'coverage': coverage,
        'power': power,
        'top_power': top_power,
        'top_coverage': top_coverage,
        'credible_set_sizes': credible_set_sizes_str,
        'causal_variant_pips': causal_pips_str
    }

def analyze_locus(locus_name, input_dir, output_dir):
    """Analyze a single locus and return metrics for all LD mismatch degrees."""
    logger.info(f"Analyzing locus: {locus_name}")
    
    # Read causal variants
    causal_variants = read_causal_variants(locus_name, input_dir)
    if not causal_variants:
        logger.warning(f"No causal variants found for {locus_name}")
        return []
    
    # Find fine-mapping output files
    fine_mapping_files = find_fine_mapping_files(locus_name, output_dir)
    if not fine_mapping_files:
        logger.warning(f"No fine-mapping files found for {locus_name}")
        return []
    
    results = []
    
    # Analyze each fine-mapping file
    for file_path in fine_mapping_files:
        ld_mismatch_degree = extract_ld_mismatch_degree(file_path, locus_name)
        
        # Read fine-mapping results
        df = read_fine_mapping_results(file_path)
        if df is None:
            continue
        
        # Compute metrics
        metrics = compute_metrics(locus_name, causal_variants, df, ld_mismatch_degree)
        results.append(metrics)
        
        logger.debug(f"Processed {locus_name} - {ld_mismatch_degree}")
    
    logger.info(f"Completed analysis for {locus_name}: {len(results)} conditions")
    return results

def main():
    """Main function to orchestrate the analysis."""
    args = parse_arguments()
    
    # Validate input directories and files
    if not Path(args.loci_file).exists():
        logger.error(f"Loci file does not exist: {args.loci_file}")
        sys.exit(1)
    
    if not Path(args.input_dir).exists():
        logger.error(f"Input directory does not exist: {args.input_dir}")
        sys.exit(1)
    
    if not Path(args.output_dir).exists():
        logger.error(f"Output directory does not exist: {args.output_dir}")
        sys.exit(1)
    
    # Create output directory for results file if needed
    Path(args.results_file).parent.mkdir(parents=True, exist_ok=True)
    
    # Read loci list
    loci = read_loci_list(args.loci_file)
    
    logger.info(f"Starting analysis with {args.num_cores} CPU cores")
    
    # Set up parallel processing
    partial_analyze = partial(analyze_locus, 
                             input_dir=args.input_dir, 
                             output_dir=args.output_dir)
    
    # Process loci in parallel
    with mp.Pool(processes=args.num_cores) as pool:
        locus_results = pool.map(partial_analyze, loci)
    
    # Flatten results list
    all_results = []
    for locus_result in locus_results:
        all_results.extend(locus_result)
    
    if not all_results:
        logger.error("No results were generated")
        sys.exit(1)
    
    # Create results DataFrame
    results_df = pd.DataFrame(all_results)
    
    # Reorder columns for better readability
    column_order = [
        'locus',
        'ld_mismatch_degree',
        'num_causal_variants',
        'num_credible_sets',
        'coverage',
        'power',
        'top_power',
        'top_coverage',
        'credible_set_sizes',
        'causal_variant_pips'
    ]
    results_df = results_df[column_order]
    
    # Write results to tab-delimited file
    results_df.to_csv(args.results_file, sep='\t', index=False, float_format='%.6f')
    
    logger.info(f"Analysis complete. Results written to: {args.results_file}")
    logger.info(f"Total records: {len(results_df)}")
    logger.info(f"Unique loci: {results_df['locus'].nunique()}")
    logger.info(f"Unique LD mismatch conditions: {results_df['ld_mismatch_degree'].nunique()}")

if __name__ == "__main__":
    main()