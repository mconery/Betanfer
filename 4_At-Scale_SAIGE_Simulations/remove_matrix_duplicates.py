#!/usr/bin/env python3

"""
remove_matrix_duplicates.py

A helper script to remove duplicate variants from LD matrices and variant files.

This script reads paired .vcor1 (LD matrix) and .vcor1.vars (variant list) files,
identifies duplicate variants, removes the corresponding rows and columns from
the matrix, and rewrites the cleaned files back to their original locations.

Usage:
    python remove_matrix_duplicates.py <matrix_file> <vars_file>
    
    Where:
    - matrix_file: Path to the .vcor1 LD matrix file (square format)
    - vars_file: Path to the .vcor1.vars variant list file
    
    The script will overwrite the input files with cleaned versions.
    
Example:
    python remove_matrix_duplicates.py /path/to/locus.raw_in-sample.unphased.vcor1 /path/to/locus.raw_in-sample.unphased.vcor1.vars
"""

import sys
import numpy as np
from pathlib import Path
import tempfile
import shutil

def read_variants_file(vars_file):
    """
    Read variant IDs from .vars file and identify duplicates.
    
    Returns:
        variants (list): List of all variant IDs
        unique_indices (list): Indices of first occurrence of each unique variant
        duplicate_indices (list): Indices of duplicate occurrences
    """
    variants = []
    with open(vars_file, 'r') as f:
        for line in f:
            variant_id = line.strip()
            if variant_id:  # Skip empty lines
                variants.append(variant_id)
    
    # Find duplicates - keep track of first occurrence of each variant
    seen = {}
    unique_indices = []
    duplicate_indices = []
    
    for i, variant_id in enumerate(variants):
        if variant_id not in seen:
            # First occurrence - keep it
            seen[variant_id] = i
            unique_indices.append(i)
        else:
            # Duplicate occurrence - mark for removal
            duplicate_indices.append(i)
    
    return variants, unique_indices, duplicate_indices

def read_matrix_file(matrix_file):
    """
    Read LD matrix from .vcor1 file.
    
    PLINK --r square format is a square matrix with values separated by tabs/spaces.
    Each row represents correlations with all variants.
    
    Returns:
        matrix (numpy.ndarray): The LD matrix
    """
    try:
        # Read the matrix - handle both tab and space delimiters
        matrix = np.loadtxt(matrix_file, dtype=float)
        return matrix
    except Exception as e:
        print(f"Error reading matrix file {matrix_file}: {e}")
        raise

def write_matrix_file(matrix, matrix_file):
    """
    Write LD matrix to .vcor1 file in PLINK format.
    
    Args:
        matrix (numpy.ndarray): The LD matrix
        matrix_file (str): Path to output file
    """
    try:
        # Write in space-delimited format (PLINK standard)
        np.savetxt(matrix_file, matrix, fmt='%.6g', delimiter=' ')
    except Exception as e:
        print(f"Error writing matrix file {matrix_file}: {e}")
        raise

def write_variants_file(variants, vars_file):
    """
    Write variant list to .vars file.
    
    Args:
        variants (list): List of variant IDs
        vars_file (str): Path to output file
    """
    try:
        with open(vars_file, 'w') as f:
            for variant_id in variants:
                f.write(f"{variant_id}\n")
    except Exception as e:
        print(f"Error writing variants file {vars_file}: {e}")
        raise

def remove_duplicates(matrix_file, vars_file, verbose=True):
    """
    Remove duplicate variants from LD matrix and variant files.
    
    Args:
        matrix_file (str): Path to .vcor1 matrix file
        vars_file (str): Path to .vcor1.vars variant file
        verbose (bool): Print progress information
        
    Returns:
        tuple: (n_original, n_duplicates_removed, n_final)
    """
    
    if verbose:
        print(f"Processing: {matrix_file}")
        print(f"Variants file: {vars_file}")
    
    # Read variant list and identify duplicates
    variants, unique_indices, duplicate_indices = read_variants_file(vars_file)
    n_original = len(variants)
    n_duplicates = len(duplicate_indices)
    
    if verbose:
        print(f"Found {n_original} total variants, {n_duplicates} duplicates to remove")
    
    if n_duplicates == 0:
        if verbose:
            print("No duplicates found - files unchanged")
        return n_original, 0, n_original
    
    # Read LD matrix
    matrix = read_matrix_file(matrix_file)
    
    # Verify matrix dimensions match variant count
    if matrix.shape[0] != n_original or matrix.shape[1] != n_original:
        raise ValueError(f"Matrix dimensions {matrix.shape} don't match variant count {n_original}")
    
    # Remove duplicate rows and columns from matrix
    # Keep only the unique indices (first occurrence of each variant)
    cleaned_matrix = matrix[np.ix_(unique_indices, unique_indices)]
    
    # Create cleaned variant list
    cleaned_variants = [variants[i] for i in unique_indices]
    n_final = len(cleaned_variants)
    
    if verbose:
        print(f"Matrix reduced from {matrix.shape[0]}x{matrix.shape[1]} to {cleaned_matrix.shape[0]}x{cleaned_matrix.shape[1]}")
        print(f"Variant list reduced from {n_original} to {n_final} variants")
    
    # Create temporary files for safe writing
    matrix_temp = None
    vars_temp = None
    
    try:
        # Write to temporary files first
        with tempfile.NamedTemporaryFile(mode='w', delete=False, 
                                       dir=Path(matrix_file).parent,
                                       prefix='.tmp_matrix_') as f:
            matrix_temp = f.name
            
        with tempfile.NamedTemporaryFile(mode='w', delete=False,
                                       dir=Path(vars_file).parent, 
                                       prefix='.tmp_vars_') as f:
            vars_temp = f.name
        
        # Write cleaned data to temporary files
        write_matrix_file(cleaned_matrix, matrix_temp)
        write_variants_file(cleaned_variants, vars_temp)
        
        # Atomically replace original files
        shutil.move(matrix_temp, matrix_file)
        shutil.move(vars_temp, vars_file)
        
        if verbose:
            print(f"Successfully cleaned files written back to original locations")
        
        return n_original, n_duplicates, n_final
        
    except Exception as e:
        # Clean up temporary files if something went wrong
        if matrix_temp and Path(matrix_temp).exists():
            Path(matrix_temp).unlink()
        if vars_temp and Path(vars_temp).exists():
            Path(vars_temp).unlink()
        raise e

def main():
    if len(sys.argv) != 3:
        print("Usage: python remove_matrix_duplicates.py <matrix_file> <vars_file>")
        print()
        print("Example:")
        print("  python remove_matrix_duplicates.py locus.raw_in-sample.unphased.vcor1 locus.raw_in-sample.unphased.vcor1.vars")
        sys.exit(1)
    
    matrix_file = sys.argv[1]
    vars_file = sys.argv[2]
    
    # Verify input files exist
    if not Path(matrix_file).exists():
        print(f"Error: Matrix file does not exist: {matrix_file}")
        sys.exit(1)
        
    if not Path(vars_file).exists():
        print(f"Error: Variants file does not exist: {vars_file}")
        sys.exit(1)
    
    try:
        n_original, n_duplicates, n_final = remove_duplicates(matrix_file, vars_file, verbose=True)
        
        if n_duplicates > 0:
            print(f"✓ Successfully removed {n_duplicates} duplicate variants")
            print(f"✓ Matrix and variant files updated: {n_original} → {n_final} variants")
        else:
            print("✓ No duplicates found - files unchanged")
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()