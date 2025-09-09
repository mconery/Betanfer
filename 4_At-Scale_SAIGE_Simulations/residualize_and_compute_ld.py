#!/usr/bin/env python3

"""
residualize_and_compute_ld.py

Standalone script to residualize genotypes using principal components
and compute LD correlation matrices using pysnptools for efficient
binary PLINK file reading.

This script takes PLINK binary format genotype data (.bed/.bim/.fam) and 
PC eigenvec files, residualizes the genotypes by regressing out PCs, and 
computes correlation matrices for linkage disequilibrium analysis.

Usage
-----
python3 residualize_and_compute_ld.py \
    --geno_file /path/to/genotypes \
    --pc_file /path/to/pcs.eigenvec \
    --out_prefix /path/to/output_prefix \
    --num_pcs 10 \
    --locus_name locus1 \
    --population CEU \
    --chrom 1 \
    --start 1000000 \
    --end 2000000 \
    [--cleanup_temp]

Note: geno_file should be the prefix for .bed/.bim/.fam files (without extension)

Output files:
- {out_prefix}.ld: LD correlation matrix
- {out_prefix}.log: Analysis log file
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import sys

try:
    from pysnptools.snpreader import Bed
except ImportError:
    print("Error: pysnptools package is required. Install with: pip install pysnptools", file=sys.stderr)
    sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Residualize genotypes using PCs and compute LD correlation matrix."
    )
    
    parser.add_argument("--geno_file", required=True,
                       help="Path prefix to PLINK binary files (.bed/.bim/.fam) without extension.")
    parser.add_argument("--pc_file", required=True,
                       help="Path to PLINK eigenvec format PC file.")
    parser.add_argument("--out_prefix", required=True,
                       help="Output file prefix for LD matrix and log.")
    parser.add_argument("--num_pcs", type=int, default=10,
                       help="Number of PCs to regress out (default: 10).")
    parser.add_argument("--locus_name", required=True,
                       help="Name of the locus being analyzed.")
    parser.add_argument("--population", required=True,
                       help="Population name.")
    parser.add_argument("--chrom", required=True,
                       help="Chromosome number.")
    parser.add_argument("--start", type=int, required=True,
                       help="Start position of the region.")
    parser.add_argument("--end", type=int, required=True,
                       help="End position of the region.")
    parser.add_argument("--cleanup_temp", action="store_true",
                       help="Remove temporary genotype files after processing.")
    
    return parser.parse_args()

def load_genotype_data(geno_file_prefix):
    """Load PLINK binary format genotype data using pysnptools."""
    bed_file = f"{geno_file_prefix}.bed"
    if not Path(bed_file).exists():
        raise FileNotFoundError(f'Genotype file {bed_file} not found')
    
    print(f'Loading genotypes from {geno_file_prefix}')
    
    # Load using pysnptools - this is much faster than reading raw files
    bed_reader = Bed(geno_file_prefix, count_A1=False)  # count_A1=False for 0/1/2 encoding
    snp_data = bed_reader.read()
    
    # Debug: print the structure of the loaded data
    print(f"Debug - snp_data.iid shape: {snp_data.iid.shape}")
    print(f"Debug - snp_data.sid shape: {snp_data.sid.shape}")
    print(f"Debug - snp_data.pos shape: {snp_data.pos.shape}")
    
    # Extract sample IDs more robustly
    if snp_data.iid.ndim == 2:
        # Standard format with FID and IID columns
        sample_ids = pd.DataFrame({
            'FID': snp_data.iid[:, 0],
            'IID': snp_data.iid[:, 1]
        })
    elif snp_data.iid.ndim == 1:
        # Single column format - use same value for both FID and IID
        sample_ids = pd.DataFrame({
            'FID': snp_data.iid,
            'IID': snp_data.iid
        })
    else:
        raise ValueError(f"Unexpected iid structure: shape {snp_data.iid.shape}")
    
    # Extract SNP information more robustly
    n_snps = len(snp_data.sid)
    
    # Handle sid (SNP IDs) - can be 1D or 2D
    if snp_data.sid.ndim == 1:
        snp_ids = snp_data.sid
        chrom_info = [''] * n_snps  # No chromosome info available
    elif snp_data.sid.ndim == 2:
        if snp_data.sid.shape[1] >= 2:
            chrom_info = snp_data.sid[:, 0]
            snp_ids = snp_data.sid[:, 1]
        else:
            snp_ids = snp_data.sid[:, 0]
            chrom_info = [''] * n_snps
    else:
        raise ValueError(f"Unexpected sid structure: shape {snp_data.sid.shape}")
    
    # Handle pos (positions) - can be 1D or 2D
    if snp_data.pos.ndim == 1:
        positions = snp_data.pos
    elif snp_data.pos.ndim == 2:
        # Usually pos has columns: [chr, pos] or [chr, cM, pos]
        if snp_data.pos.shape[1] >= 2:
            positions = snp_data.pos[:, -1]  # Take the last column as BP position
        else:
            positions = snp_data.pos[:, 0]
    else:
        positions = np.zeros(n_snps)  # Fallback
    
    snp_info = pd.DataFrame({
        'CHR': chrom_info,
        'SNP': snp_ids,
        'BP': positions
    })
    
    genotypes = snp_data.val  # This is samples x SNPs
    n_samples, n_snps = genotypes.shape
    print(f'Loaded {n_samples} samples, {n_snps} SNPs')
    
    # Check for and handle missing values
    n_missing = np.isnan(genotypes).sum()
    if n_missing > 0:
        print(f'Found {n_missing} missing genotype calls')
    
    return sample_ids, genotypes, snp_info

def load_pc_data(pc_file, num_pcs):
    """Load PC data from PLINK eigenvec format."""
    print(f'Loading PCs from {pc_file}')
    pc_df = pd.read_csv(pc_file, sep=r'\s+', header=0,
                       names=['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)])
    
    # Select only the required number of PCs
    pc_cols = ['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)]
    pc_df = pc_df[pc_cols]
    
    return pc_df

def residualize_genotypes(sample_ids, genotypes, pc_df, num_pcs):
    """Residualize genotypes by regressing out PCs."""
    # Merge genotype and PC data on sample IDs
    merged_df = sample_ids.merge(pc_df, on=['FID', 'IID'], how='left')
    if merged_df.isnull().any().any():
        missing_samples = merged_df[merged_df.isnull().any(axis=1)][['FID', 'IID']]
        print(f"Warning: Missing PC data for {len(missing_samples)} samples")
        print("First few missing samples:")
        print(missing_samples.head())
        raise ValueError('Missing PC data for some samples')
    
    # Extract PC matrix
    pc_matrix = merged_df.iloc[:, 2:2+num_pcs].values
    print(f'Using {pc_matrix.shape[1]} PCs for residualization')
    
    n_samples, n_snps = genotypes.shape
    
    # Residualize each SNP
    print('Residualizing genotypes...')
    residual_genotypes = np.zeros_like(genotypes)
    
    n_processed = 0
    n_invariant = 0
    n_insufficient_data = 0
    
    for snp_idx in range(n_snps):
        if snp_idx % 1000 == 0 and snp_idx > 0:
            print(f'  Processed {snp_idx}/{n_snps} SNPs ({n_processed} residualized, {n_invariant} invariant, {n_insufficient_data} insufficient data)')
        
        y = genotypes[:, snp_idx]
        
        # Skip if all genotypes are missing or invariant
        valid_mask = ~np.isnan(y)
        if valid_mask.sum() < 10:
            residual_genotypes[:, snp_idx] = y
            n_insufficient_data += 1
            continue
            
        if np.var(y[valid_mask]) < 1e-8:
            residual_genotypes[:, snp_idx] = y
            n_invariant += 1
            continue
        
        # Get valid observations
        y_valid = y[valid_mask]
        X_valid = np.column_stack([np.ones(valid_mask.sum()), pc_matrix[valid_mask, :]])
        
        # Solve least squares: beta = (X'X)^-1 X'y
        try:
            XtX_inv = np.linalg.inv(X_valid.T @ X_valid)
            beta = XtX_inv @ X_valid.T @ y_valid
            
            # Compute residuals for all samples (including missing)
            X_all = np.column_stack([np.ones(n_samples), pc_matrix])
            y_pred = X_all @ beta
            residuals = y - y_pred
            
            # Keep original missing values as missing
            residuals[~valid_mask] = np.nan
            residual_genotypes[:, snp_idx] = residuals
            n_processed += 1
            
        except np.linalg.LinAlgError:
            # If matrix is singular, keep original genotypes
            residual_genotypes[:, snp_idx] = y
            n_insufficient_data += 1
    
    print(f'Residualization complete: {n_processed} SNPs residualized, {n_invariant} invariant, {n_insufficient_data} insufficient data')
    return residual_genotypes

def compute_correlation_matrix(residual_genotypes):
    """Compute correlation matrix from residualized genotypes using vectorized operations."""
    print('Computing correlation matrix from residuals...')
    
    n_samples, n_snps = residual_genotypes.shape
    
    # For large numbers of SNPs, use numpy's corrcoef which is much faster
    # but we need to handle missing values carefully
    if n_snps <= 10000:  # For reasonable sized matrices, compute full correlation
        # Transpose to get SNPs x samples for corrcoef
        genotypes_t = residual_genotypes.T
        
        # Handle missing values by computing correlations pairwise
        corr_matrix = np.full((n_snps, n_snps), np.nan)
        
        print(f'Computing {n_snps} x {n_snps} correlation matrix...')
        for i in range(n_snps):
            if i % 100 == 0 and i > 0:
                print(f'  Computing correlations for SNP {i}/{n_snps}')
            
            for j in range(i, n_snps):
                x, y = genotypes_t[i], genotypes_t[j]
                valid = ~(np.isnan(x) | np.isnan(y))
                
                if valid.sum() > 10:
                    if i == j:
                        corr_matrix[i, j] = 1.0
                    else:
                        corr = np.corrcoef(x[valid], y[valid])[0, 1]
                        corr_matrix[i, j] = corr_matrix[j, i] = corr if not np.isnan(corr) else 0.0
                else:
                    corr_matrix[i, j] = corr_matrix[j, i] = 0.0
    else:
        print(f"Large matrix ({n_snps} SNPs) - this may take some time...")
        # For very large matrices, we might want to implement chunked processing
        # or use sparse methods, but for now use the same approach
        corr_matrix = np.full((n_snps, n_snps), np.nan)
        
        for i in range(n_snps):
            if i % 500 == 0:
                print(f'  Computing correlations for SNP {i + 1}/{n_snps}')
            
            for j in range(i, n_snps):
                x, y = residual_genotypes[:, i], residual_genotypes[:, j]
                valid = ~(np.isnan(x) | np.isnan(y))
                
                if valid.sum() > 10:
                    if i == j:
                        corr_matrix[i, j] = 1.0
                    else:
                        corr = np.corrcoef(x[valid], y[valid])[0, 1]
                        corr_matrix[i, j] = corr_matrix[j, i] = corr if not np.isnan(corr) else 0.0
                else:
                    corr_matrix[i, j] = corr_matrix[j, i] = 0.0
    
    return corr_matrix

def write_output_files(corr_matrix, out_prefix, args, snp_info=None):
    """Write LD matrix and log files."""
    # Save correlation matrix in tab-delimited format
    ld_file = f"{out_prefix}.ld"
    print(f'Saving LD matrix to {ld_file}')
    np.savetxt(ld_file, corr_matrix, fmt='%.6f', delimiter='\t')
    
    # Create log file
    log_file = f"{out_prefix}.log"
    n_samples, n_snps = corr_matrix.shape[0], corr_matrix.shape[1]
    
    with open(log_file, 'w') as f:
        f.write(f'PC-residualized LD matrix for locus {args.locus_name}\n')
        f.write(f'Population: {args.population}\n')
        f.write(f'Chromosome: {args.chrom}\n')
        f.write(f'Region: {args.start}-{args.end}\n')
        f.write(f'Number of samples: {n_samples}\n')
        f.write(f'Number of SNPs: {n_snps}\n')
        f.write(f'Number of PCs regressed out: {args.num_pcs}\n')
        f.write(f'PC file: {args.pc_file}\n')
        f.write(f'Genotype file: {args.geno_file}\n')
        f.write(f'Method: Binary PLINK files with pysnptools\n')
        
        # Add some matrix statistics
        non_nan = ~np.isnan(corr_matrix)
        if non_nan.any():
            f.write(f'LD matrix mean: {np.nanmean(corr_matrix):.6f}\n')
            f.write(f'LD matrix std: {np.nanstd(corr_matrix):.6f}\n')
            f.write(f'Max LD (off-diagonal): {np.nanmax(corr_matrix[~np.eye(n_snps, dtype=bool)]):.6f}\n')
    
    print(f'Successfully created PC-residualized LD matrix: {ld_file}')

def cleanup_temp_files(geno_file_prefix):
    """Remove temporary genotype files."""
    temp_files = [
        f"{geno_file_prefix}.bed",
        f"{geno_file_prefix}.bim", 
        f"{geno_file_prefix}.fam",
        f"{geno_file_prefix}.log",
        f"{geno_file_prefix}.nosex"
    ]
    
    for temp_file in temp_files:
        try:
            Path(temp_file).unlink(missing_ok=True)
            if Path(temp_file).exists():
                print(f'Warning: Could not remove {temp_file}')
            else:
                print(f'Removed temporary file: {temp_file}')
        except Exception as e:
            print(f'Warning: Could not remove {temp_file}: {e}')

def main():
    args = parse_args()
    
    try:
        # Load data using pysnptools (much faster than raw format)
        sample_ids, genotypes, snp_info = load_genotype_data(args.geno_file)
        pc_df = load_pc_data(args.pc_file, args.num_pcs)
        
        # Residualize genotypes
        residual_genotypes = residualize_genotypes(sample_ids, genotypes, pc_df, args.num_pcs)
        
        # Compute correlation matrix
        corr_matrix = compute_correlation_matrix(residual_genotypes)
        
        # Write output files
        write_output_files(corr_matrix, args.out_prefix, args, snp_info)
        
        # Cleanup if requested
        if args.cleanup_temp:
            cleanup_temp_files(args.geno_file)
        
        print("Analysis completed successfully!")
        
    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()