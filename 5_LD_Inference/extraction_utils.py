"""
extraction_utils.py

Utility functions for reading and processing GWAS summary statistics.
"""

import pandas as pd
import numpy as np
import os
import logging
from typing import List, Dict, Set, Optional, Union
import gzip
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_variant_list(file_path: str, variant_col: str = 'variant_id', sep: str = '\t') -> Set[str]:
    """
    Load the common set of variants to extract from a file.
    
    Args:
        file_path: Path to file containing variant IDs
        variant_col: Name of the column containing variant IDs
        sep: Column separator (default: tab)
        
    Returns:
        Set of variant IDs to extract
    """
    try:
        df = pd.read_csv(file_path, sep=sep)
        
        if variant_col not in df.columns:
            raise ValueError(f"Column '{variant_col}' not found in {file_path}")
        
        variant_set = set(df[variant_col].astype(str))
        logger.info(f"Loaded {len(variant_set)} variants from {file_path}")
        return variant_set
        
    except Exception as e:
        logger.error(f"Error loading variant list from {file_path}: {e}")
        raise


def validate_columns(df: pd.DataFrame, variant_col: str, beta_col: str, file_path: str = "") -> bool:
    """
    Verify that required columns exist in the dataframe.
    
    Args:
        df: DataFrame to validate
        variant_col: Name of variant ID column
        beta_col: Name of beta column
        file_path: File path for error reporting
        
    Returns:
        True if columns exist, raises ValueError otherwise
    """
    missing_cols = []
    if variant_col not in df.columns:
        missing_cols.append(variant_col)
    if beta_col not in df.columns:
        missing_cols.append(beta_col)
    
    if missing_cols:
        raise ValueError(f"Missing columns {missing_cols} in {file_path}. Available columns: {list(df.columns)}")
    
    return True


def read_gwas_file(file_path: str, variant_col: str, beta_col: str, 
                   variant_set: Set[str], chunksize: int = 50000, sep: str = '\t') -> pd.Series:
    """
    Read and filter a single GWAS file to extract beta values for specified variants.
    
    Args:
        file_path: Path to GWAS summary statistics file
        variant_col: Name of the variant ID column
        beta_col: Name of the beta column
        variant_set: Set of variants to extract
        chunksize: Number of rows to read at once for large files
        sep: Column separator (default: tab)
        
    Returns:
        Series with variant IDs as index and beta values
    """
    try:
        # Initialize result series
        result_betas = pd.Series(dtype=float, name=os.path.basename(file_path))
        
        # Read file in chunks for memory efficiency
        chunk_iter = pd.read_csv(file_path, sep=sep, chunksize=chunksize, low_memory=False)
        
        total_variants_found = 0
        
        for chunk_num, chunk in enumerate(chunk_iter):
            # Validate columns on first chunk
            if chunk_num == 0:
                validate_columns(chunk, variant_col, beta_col, file_path)
            
            # Convert variant IDs to string for consistent matching
            chunk[variant_col] = chunk[variant_col].astype(str)
            
            # Filter to variants of interest
            mask = chunk[variant_col].isin(variant_set)
            filtered_chunk = chunk[mask]
            
            if len(filtered_chunk) > 0:
                # Extract beta values, handling duplicates by taking the first occurrence
                chunk_betas = filtered_chunk.set_index(variant_col)[beta_col]
                chunk_betas = chunk_betas[~chunk_betas.index.duplicated(keep='first')]
                
                # Merge with result
                result_betas = pd.concat([result_betas, chunk_betas])
                total_variants_found += len(chunk_betas)
        
        # Remove any duplicates that might have occurred across chunks
        result_betas = result_betas[~result_betas.index.duplicated(keep='first')]
        
        logger.info(f"Extracted {len(result_betas)} variants from {file_path}")
        return result_betas
        
    except Exception as e:
        logger.error(f"Error reading {file_path}: {e}")
        raise


def extract_betas_batch(file_list: List[str], variant_col: str, beta_col: str, 
                       variant_set: Set[str], chunksize: int = 50000, sep: str = '\t') -> pd.DataFrame:
    """
    Process multiple GWAS files efficiently to extract beta values.
    
    Args:
        file_list: List of paths to GWAS summary statistics files
        variant_col: Name of the variant ID column
        beta_col: Name of the beta column
        variant_set: Set of variants to extract
        chunksize: Number of rows to read at once for large files
        sep: Column separator (default: tab)
        
    Returns:
        DataFrame with variants as rows and studies as columns
    """
    beta_dict = {}
    failed_files = []
    
    logger.info(f"Processing {len(file_list)} GWAS files...")
    
    for i, file_path in enumerate(file_list):
        try:
            logger.info(f"Processing file {i+1}/{len(file_list)}: {os.path.basename(file_path)}")
            
            # Extract betas for this file
            file_betas = read_gwas_file(file_path, variant_col, beta_col, variant_set, chunksize, sep)
            
            # Use filename as column name
            study_name = os.path.splitext(os.path.basename(file_path))[0]
            if study_name.endswith('.gz'):
                study_name = os.path.splitext(study_name)[0]
            
            beta_dict[study_name] = file_betas
            
        except Exception as e:
            logger.error(f"Failed to process {file_path}: {e}")
            failed_files.append(file_path)
            continue
    
    if failed_files:
        logger.warning(f"Failed to process {len(failed_files)} files: {failed_files}")
    
    # Combine all beta series into a DataFrame
    if beta_dict:
        beta_matrix = pd.DataFrame(beta_dict)
        logger.info(f"Created beta matrix with shape {beta_matrix.shape}")
        return beta_matrix
    else:
        logger.error("No files were successfully processed")
        return pd.DataFrame()


def memory_efficient_reader(file_path: str, variant_col: str, beta_col: str, 
                          variant_set: Set[str], chunksize: int = 10000, sep: str = '\t') -> pd.Series:
    """
    Memory-efficient reader for very large GWAS files.
    
    Args:
        file_path: Path to GWAS file
        variant_col: Name of variant ID column
        beta_col: Name of beta column
        variant_set: Set of variants to extract
        chunksize: Smaller chunk size for memory efficiency
        sep: Column separator (default: tab)
        
    Returns:
        Series with extracted beta values
    """
    return read_gwas_file(file_path, variant_col, beta_col, variant_set, chunksize, sep)


def validate_variant_format(variant_id: str) -> str:
    """
    Standardize variant ID formats (e.g., chr:pos:ref:alt).
    
    Args:
        variant_id: Raw variant ID
        
    Returns:
        Standardized variant ID
    """
    # Remove any whitespace
    variant_id = str(variant_id).strip()
    
    # Add more standardization rules as needed
    # For example, ensure chromosome has 'chr' prefix, standardize case, etc.
    
    return variant_id


def save_beta_matrix(beta_matrix: pd.DataFrame, output_path: str, 
                    compress: bool = True) -> None:
    """
    Save the beta matrix to file.
    
    Args:
        beta_matrix: DataFrame with beta values
        output_path: Path to save the file
        compress: Whether to compress the output
    """
    try:
        if compress and not output_path.endswith('.gz'):
            output_path += '.gz'
        
        if output_path.endswith('.gz'):
            beta_matrix.to_csv(output_path, sep='\t', compression='gzip')
        else:
            beta_matrix.to_csv(output_path, sep='\t')
        
        logger.info(f"Saved beta matrix to {output_path}")
        
    except Exception as e:
        logger.error(f"Error saving beta matrix: {e}")
        raise


def main():
    """
    Main function for testing utility functions.
    """
    # Example usage
    print("GWAS Beta Extraction Utilities")
    print("Available functions:")
    print("- load_variant_list()")
    print("- read_gwas_file()")
    print("- extract_betas_batch()")
    print("- save_beta_matrix()")


if __name__ == "__main__":
    main()
