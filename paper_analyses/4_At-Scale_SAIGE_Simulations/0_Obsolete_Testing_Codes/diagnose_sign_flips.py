#!/usr/bin/env python3

"""
diagnose_sign_flips.py

Diagnostic script to identify why correlations have opposite signs compared to plink2.
This will help us understand the real cause of the allele encoding issue.
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import sys

try:
    from pysnptools.snpreader import Bed
except ImportError:
    print("Error: pysnptools package is required", file=sys.stderr)
    sys.exit(1)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Diagnose correlation sign flips compared to plink2"
    )
    parser.add_argument("--geno_file", required=True,
                       help="Path prefix to PLINK binary files")
    parser.add_argument("--pc_file", required=True,
                       help="Path to PC file")
    parser.add_argument("--plink2_corr", required=True,
                       help="Path to plink2 correlation matrix (.unphased.vcor1)")
    parser.add_argument("--plink2_vars", required=True,
                       help="Path to plink2 variant list (.unphased.vcor1.vars)")
    parser.add_argument("--num_pcs", type=int, default=10)
    
    return parser.parse_args()

def load_genotype_data(geno_file_prefix):
    """Load genotype data and examine its properties."""
    print(f"Loading genotypes from {geno_file_prefix}")
    
    bed_reader = Bed(geno_file_prefix, count_A1=False)  # Default A2 counting
    snp_data = bed_reader.read()
    
    # Extract sample and SNP info
    if snp_data.iid.ndim == 2:
        sample_ids = pd.DataFrame({
            'FID': snp_data.iid[:, 0],
            'IID': snp_data.iid[:, 1]
        })
    else:
        sample_ids = pd.DataFrame({
            'FID': snp_data.iid,
            'IID': snp_data.iid
        })
    
    # Get SNP IDs
    if snp_data.sid.ndim == 1:
        snp_ids = snp_data.sid
    else:
        snp_ids = snp_data.sid[:, -1] if snp_data.sid.shape[1] > 1 else snp_data.sid[:, 0]
    
    genotypes = snp_data.val
    
    print(f"Loaded {genotypes.shape[0]} samples, {genotypes.shape[1]} SNPs")
    
    # Check allele frequencies to understand encoding
    print("\nAllele frequency analysis:")
    for i in range(min(5, len(snp_ids))):
        geno = genotypes[:, i]
        valid_geno = geno[~np.isnan(geno)]
        
        if len(valid_geno) > 0:
            freq_0 = np.mean(valid_geno == 0)
            freq_1 = np.mean(valid_geno == 1) 
            freq_2 = np.mean(valid_geno == 2)
            
            print(f"  {snp_ids[i]}: 0={freq_0:.3f}, 1={freq_1:.3f}, 2={freq_2:.3f}")
    
    return sample_ids, genotypes, snp_ids

def load_pc_data(pc_file, num_pcs):
    """Load PC data."""
    print(f"Loading PCs from {pc_file}")
    pc_df = pd.read_csv(pc_file, sep=r'\s+', header=0,
                       names=['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)])
    
    pc_cols = ['FID', 'IID'] + [f'PC{i}' for i in range(1, num_pcs + 1)]
    pc_df = pc_df[pc_cols]
    
    return pc_df

def residualize_genotypes_simple(sample_ids, genotypes, pc_df, num_pcs):
    """Simple CPU residualization for diagnostic purposes."""
    print("Performing residualization...")
    
    merged_df = sample_ids.merge(pc_df, on=['FID', 'IID'], how='left')
    if merged_df.isnull().any().any():
        raise ValueError('Missing PC data for some samples')
    
    pc_matrix = merged_df.iloc[:, 2:2+num_pcs].values
    n_samples, n_snps = genotypes.shape
    
    # Design matrix with intercept + PCs
    X = np.column_stack([np.ones(n_samples), pc_matrix])
    
    # Precompute for efficiency
    try:
        XtX_inv_Xt = np.linalg.inv(X.T @ X) @ X.T
    except np.linalg.LinAlgError:
        XtX_inv_Xt = np.linalg.pinv(X).T
    
    residual_genotypes = np.zeros_like(genotypes)
    
    for i in range(n_snps):
        if i % 100 == 0:
            print(f"  Processing SNP {i+1}/{n_snps}")
        
        y = genotypes[:, i]
        valid_mask = ~np.isnan(y)
        
        if valid_mask.sum() < 10 or np.var(y[valid_mask]) < 1e-8:
            residual_genotypes[:, i] = y
            continue
        
        try:
            beta = XtX_inv_Xt @ y
            residuals = y - X @ beta
            residuals[~valid_mask] = np.nan
            residual_genotypes[:, i] = residuals
        except:
            residual_genotypes[:, i] = y
    
    print("Residualization complete")
    return residual_genotypes

def compute_correlation_matrix_simple(residual_genotypes):
    """Simple correlation computation."""
    print("Computing correlation matrix...")
    
    n_samples, n_snps = residual_genotypes.shape
    corr_matrix = np.full((n_snps, n_snps), np.nan)
    
    for i in range(n_snps):
        if i % 50 == 0:
            print(f"  Computing correlations for SNP {i+1}/{n_snps}")
        
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

def load_plink2_results(plink2_corr, plink2_vars):
    """Load plink2 correlation matrix and variant list."""
    print(f"Loading plink2 results...")
    
    # Load variant list
    with open(plink2_vars, 'r') as f:
        plink2_snp_ids = [line.strip() for line in f]
    
    # Load correlation matrix
    plink2_corr_matrix = np.loadtxt(plink2_corr, delimiter='\t')
    
    print(f"Plink2: {len(plink2_snp_ids)} variants, {plink2_corr_matrix.shape} correlation matrix")
    
    return plink2_corr_matrix, plink2_snp_ids

def diagnose_differences(our_corr, our_snp_ids, plink2_corr, plink2_snp_ids):
    """Compare our results with plink2 and diagnose differences."""
    print("\n" + "=" * 60)
    print("DIAGNOSING CORRELATION DIFFERENCES")
    print("=" * 60)
    
    # Check if SNP lists match
    if len(our_snp_ids) != len(plink2_snp_ids):
        print(f"ERROR: Different number of SNPs!")
        print(f"  Our script: {len(our_snp_ids)}")
        print(f"  Plink2: {len(plink2_snp_ids)}")
        return
    
    # Check SNP ID order
    snp_order_matches = all(a == b for a, b in zip(our_snp_ids, plink2_snp_ids))
    print(f"SNP order matches: {snp_order_matches}")
    
    if not snp_order_matches:
        print("First few SNP ID differences:")
        for i in range(min(10, len(our_snp_ids))):
            if our_snp_ids[i] != plink2_snp_ids[i]:
                print(f"  Position {i}: Our='{our_snp_ids[i]}', Plink2='{plink2_snp_ids[i]}'")
    
    # Compare correlation matrices
    print(f"\nCorrelation matrix comparison:")
    print(f"  Our matrix shape: {our_corr.shape}")
    print(f"  Plink2 matrix shape: {plink2_corr.shape}")
    
    if our_corr.shape != plink2_corr.shape:
        print("ERROR: Matrix shapes don't match!")
        return
    
    # Calculate differences
    valid_mask = ~(np.isnan(our_corr) | np.isnan(plink2_corr))
    
    if np.sum(valid_mask) == 0:
        print("ERROR: No valid correlation pairs to compare!")
        return
    
    our_vals = our_corr[valid_mask]
    plink2_vals = plink2_corr[valid_mask]
    
    # Overall correlation between the two methods
    method_correlation = np.corrcoef(our_vals, plink2_vals)[0, 1]
    mean_abs_diff = np.mean(np.abs(our_vals - plink2_vals))
    
    print(f"\nOverall comparison:")
    print(f"  Correlation between methods: {method_correlation:.6f}")
    print(f"  Mean absolute difference: {mean_abs_diff:.6f}")
    
    # Check for sign flips
    sign_flips = np.sum((our_vals * plink2_vals) < 0)
    total_comparisons = len(our_vals)
    
    print(f"  Sign flips: {sign_flips} / {total_comparisons} ({100*sign_flips/total_comparisons:.1f}%)")
    
    # Find the worst mismatches
    abs_diffs = np.abs(our_vals - plink2_vals)
    worst_indices = np.argsort(abs_diffs)[-10:]  # 10 worst
    
    print(f"\nWorst mismatches (off-diagonal only):")
    count = 0
    for idx in reversed(worst_indices):
        if count >= 5:  # Show top 5
            break
        
        # Convert back to matrix indices
        i, j = np.unravel_index(np.where(valid_mask.flatten())[0][idx], our_corr.shape)
        
        if i != j:  # Skip diagonal
            print(f"  {our_snp_ids[i]} vs {our_snp_ids[j]}:")
            print(f"    Our: {our_vals[idx]:.6f}, Plink2: {plink2_vals[idx]:.6f}, Diff: {abs_diffs[idx]:.6f}")
            count += 1
    
    # Check if correlations are just negated
    negated_correlation = np.corrcoef(our_vals, -plink2_vals)[0, 1]
    print(f"\nCorrelation with negated plink2 values: {negated_correlation:.6f}")
    
    if abs(negated_correlation) > abs(method_correlation):
        print("*** POSSIBLE SYSTEMATIC SIGN FLIP DETECTED! ***")
        print("Your correlations might be systematically negated compared to plink2")

def main():
    args = parse_args()
    
    try:
        print("CORRELATION SIGN FLIP DIAGNOSTIC")
        print("=" * 60)
        
        # Load our data
        sample_ids, genotypes, our_snp_ids = load_genotype_data(args.geno_file)
        pc_df = load_pc_data(args.pc_file, args.num_pcs)
        
        # Compute our correlations
        residual_genotypes = residualize_genotypes_simple(sample_ids, genotypes, pc_df, args.num_pcs)
        our_corr = compute_correlation_matrix_simple(residual_genotypes)
        
        # Load plink2 results
        plink2_corr, plink2_snp_ids = load_plink2_results(args.plink2_corr, args.plink2_vars)
        
        # Diagnose differences
        diagnose_differences(our_corr, our_snp_ids, plink2_corr, plink2_snp_ids)
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()