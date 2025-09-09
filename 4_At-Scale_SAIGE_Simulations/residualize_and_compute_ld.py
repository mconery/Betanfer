#!/usr/bin/env python3

"""
residualize_and_compute_ld.py

Standalone script to residualize genotypes using principal components
and compute LD correlation matrices.

This script takes PLINK raw format genotype data and PC eigenvec files,
residualizes the genotypes by regressing out PCs, and computes correlation
matrices for linkage disequilibrium analysis.

Usage
-----
python3 residualize_and_compute_ld.py \
    --geno_file /path/to/genotypes.raw \
    --pc_file /path/to/pcs.eigenvec \
    --out_prefix /path/to/output_prefix \
    --num_pcs 10 \
    --locus_name locus1 \
    --population CEU \
    --chrom 1 \
    --start 1000000 \
    --end 2000000 \
    [--cleanup_temp]

Output files:
- {out_prefix}.ld: LD correlation matrix
- {out_prefix}.log: Analysis log file
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Residualize genotypes using PCs and compute LD correlation matrix."
    )
    
    parser.add_argument("--geno_file", required=True,
                       help="Path to PLINK raw format genotype file.")
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

def load_genotype_data(geno_file):
    """Load PLINK --export A format genotype data."""
    if not Path(geno_file).exists():
        raise FileNotFoundError(f'Genotype file {geno_file} not found')
    
    print(f'Loading genotypes from {geno_file}')
    geno_df = pd.read_csv(geno_file, sep='\t')
    
    # Extract sample IDs and genotype matrix
    sample_ids = geno_df[['FID', 'IID']].copy()
    genotypes = geno_df.iloc[:, 6:].values.astype(float)  # Skip FID, IID, PAT, MAT, SEX, PHENOTYPE
    
    n_samples, n_snps = genotypes.shape
    print(f'Loaded {n_samples} samples, {n_snps} SNPs')
    
    return sample_ids, genotypes

def load_pc_data(pc_file, num_pcs):
    """Load PC data from PLINK eigenvec format."""
    print(f'Loading PCs from {pc_file}')
    pc_df = pd.read_csv(pc_file, sep=r'\s+', header=None,
                       names=['FID', 'IID'] + [f'PC{i}' for i in range(1, 21)])
    
    # Select only the required number of PCs
    pc_cols = ['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)]
    pc_df = pc_df[pc_cols]
    
    return pc_df

def residualize_genotypes(sample_ids, genotypes, pc_df, num_pcs):
    """Residualize genotypes by regressing out PCs."""
    # Merge genotype and PC data on sample IDs
    merged_df = sample_ids.merge(pc_df, on=['FID', 'IID'], how='left')
    if merged_df.isnull().any().any():
        raise ValueError('Missing PC data for some samples')
    
    # Extract PC matrix
    pc_matrix = merged_df.iloc[:, 2:2+num_pcs].values
    print(f'Using {pc_matrix.shape[1]} PCs for residualization')
    
    n_samples, n_snps = genotypes.shape
    
    # Residualize each SNP
    print('Residualizing genotypes...')
    residual_genotypes = np.zeros_like(genotypes)
    
    for snp_idx in range(n_snps):
        if snp_idx % 1000 == 0:
            print(f'  Processing SNP {snp_idx + 1}/{n_snps}')
        
        y = genotypes[:, snp_idx]
        
        # Skip if all genotypes are missing or invariant
        valid_mask = ~np.isnan(y)
        if valid_mask.sum() < 10 or np.var(y[valid_mask]) < 1e-8:
            residual_genotypes[:, snp_idx] = y
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
            
        except np.linalg.LinAlgError:
            # If matrix is singular, keep original genotypes
            residual_genotypes[:, snp_idx] = y
    
    return residual_genotypes

def compute_correlation_matrix(residual_genotypes):
    """Compute correlation matrix from residualized genotypes."""
    print('Computing correlation matrix from residuals...')
    
    n_samples, n_snps = residual_genotypes.shape
    
    # Compute correlation matrix, handling missing values
    corr_matrix = np.full((n_snps, n_snps), np.nan)
    
    for i in range(n_snps):
        if i % 100 == 0:
            print(f'  Computing correlations for SNP {i + 1}/{n_snps}')
        
        for j in range(i, n_snps):
            x, y = residual_genotypes[:, i], residual_genotypes[:, j]
            valid = ~(np.isnan(x) | np.isnan(y))
            
            if valid.sum() > 10:
                corr = np.corrcoef(x[valid], y[valid])[0, 1]
                corr_matrix[i, j] = corr_matrix[j, i] = corr if not np.isnan(corr) else 0.0
            else:
                corr_matrix[i, j] = corr_matrix[j, i] = 0.0
    
    return corr_matrix

def write_output_files(corr_matrix, out_prefix, args):
    """Write LD matrix and log files."""
    # Save correlation matrix in PLINK format
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
    
    print(f'Successfully created PC-residualized LD matrix: {ld_file}')

def cleanup_temp_files(geno_file):
    """Remove temporary genotype files."""
    base_path = Path(geno_file).with_suffix('')
    temp_files = [
        f"{base_path}.raw",
        f"{base_path}.log", 
        f"{base_path}.nosex"
    ]
    
    for temp_file in temp_files:
        try:
            Path(temp_file).unlink(missing_ok=True)
            print(f'Removed temporary file: {temp_file}')
        except Exception as e:
            print(f'Warning: Could not remove {temp_file}: {e}')

def main():
    args = parse_args()
    
    try:
        # Load data
        sample_ids, genotypes = load_genotype_data(args.geno_file)
        pc_df = load_pc_data(args.pc_file, args.num_pcs)
        
        # Residualize genotypes
        residual_genotypes = residualize_genotypes(sample_ids, genotypes, pc_df, args.num_pcs)
        
        # Compute correlation matrix
        corr_matrix = compute_correlation_matrix(residual_genotypes)
        
        # Write output files
        write_output_files(corr_matrix, args.out_prefix, args)
        
        # Cleanup if requested
        if args.cleanup_temp:
            cleanup_temp_files(args.geno_file)
        
        print("Analysis completed successfully!")
        
    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()