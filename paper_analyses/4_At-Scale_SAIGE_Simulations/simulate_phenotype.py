#!/usr/bin/env python3

"""
simulate_phenotype.py

Simulate phenotypes from PLINK files with per-SNP heritability specification,
MAF-based effect‐size scaling, and inclusion of principal components.
Requires Python 3.9.
"""

import argparse
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed

def main():
    parser = argparse.ArgumentParser(
        description='Simulate phenotypes from PLINK files with per-SNP heritability'
    )
    parser.add_argument('--plink', required=True,
                        help='PLINK file prefix (without .bed/.bim/.fam extension)')
    parser.add_argument('--output', required=True,
                        help='Output phenotype file path')
    parser.add_argument('--variants', required=True,
                        help='File with causal variant IDs (one per line)')
    parser.add_argument('--trait-type', choices=['binary', 'continuous'], required=True,
                        help='Trait type')
    parser.add_argument('--per-snp-heritability', type=float, required=True,
                        help='Per-SNP heritability (0 < h2_snp < 1)')
    parser.add_argument('--prevalence', type=float, default=0.3,
                        help='Case prevalence for binary traits (0-1, default=0.3)')
    parser.add_argument('--pc-file', required=True,
                        help='File with principal component vectors (tab: FID IID PC1 ... PC10)')
    parser.add_argument('--seed', type=int, default=5,
                        help='Random seed (default=5)')
    args = parser.parse_args()

    # Validate per-SNP heritability
    M = 0  # will set after reading variants
    if args.per_snp_heritability <= 0 or args.per_snp_heritability >= 1:
        raise ValueError("per-SNP heritability must be between 0 and 1")

    # Set seed
    np.random.seed(args.seed)

    # Load causal variant IDs
    causal_variants = []
    with open(args.variants) as f:
        for line in f:
            vid = line.strip()
            if vid:
                causal_variants.append(vid)
    M = len(causal_variants)
    if M == 0:
        raise ValueError("No causal variant IDs provided")

    # Read PLINK genotype data
    bed = Bed(args.plink, count_A1=False)
    # Identify indices of causal SNPs
    causal_indices = []
    causal_sids = []
    for idx, sid in enumerate(bed.sid):
        if sid in causal_variants:
            causal_indices.append(idx)
            causal_sids.append(sid)
    if len(causal_indices) == 0:
        raise ValueError("No causal variants found in PLINK data")
    # Extract genotypes for causal SNPs
    G_mat = bed[:, causal_indices].read().val   # shape (n_samples, M)
    # Sample IDs
    iids = bed.iid

    # Compute allele frequencies per SNP
    p = np.mean(G_mat, axis=0) / 2.0  # length M

    # Compute effect sizes per SNP to achieve per-SNP heritability
    beta = np.zeros(M, dtype=float)
    for i in range(M):
        if p[i] > 0 and p[i] < 1:
            beta[i] = np.sqrt(args.per_snp_heritability / (2.0 * p[i] * (1.0 - p[i])))
        else:
            beta[i] = 0.0

    # Total trait heritability
    h2_total = M * args.per_snp_heritability
    if h2_total >= 1.0:
        raise ValueError(f"Total trait heritability ({h2_total:.3f}) must be < 1")

    # Normalize genetic score variance to 1
    raw_scores = G_mat @ beta
    var_raw = np.var(raw_scores)
    if var_raw <= 0:
        raise ValueError("Zero genetic variance: check variant list or MAFs")
    norm_factor = 1.0 / np.sqrt(var_raw)
    beta *= norm_factor
    # Recompute standardized genetic scores
    G_std = G_mat @ beta
    # Confirm Var(G_std) ~1.0
    # print(f"Genetic variance after normalization: {np.var(G_std):.6f}")

    # Simulate phenotype
    n = G_std.shape[0]
    if args.trait_type == 'continuous':
        # Noise variance = (1/h2_total) - 1
        noise_var = (1.0 / h2_total) - 1.0
        noise = np.random.normal(0, np.sqrt(noise_var), n)
        pheno = G_std + noise

    else:  # binary trait via liability threshold
        # Standardize G_std (mean 0, var 1)
        G_centered = G_std - np.mean(G_std)
        G_norm = G_centered / np.sqrt(np.var(G_centered)) if np.var(G_centered) > 0 else G_centered
        # Liability L = sqrt(h2_total)*G_norm + sqrt(1-h2_total)*ε
        eps = np.random.normal(0, 1, n)
        L = np.sqrt(h2_total) * G_norm + np.sqrt(1.0 - h2_total) * eps
        # Threshold for prevalence
        thresh = np.percentile(L, 100.0 * (1.0 - args.prevalence))
        pheno = (L > thresh).astype(int)

    # Build output DataFrame
    df_pheno = pd.DataFrame({
        '#FID': iids[:, 0],
        'IID':  iids[:, 1],
        'PHENO': pheno
    })
    # Load PCs and merge
    df_pcs = pd.read_csv(args.pc_file, sep='\t')
    df_out = pd.merge(df_pheno, df_pcs, on=['#FID', 'IID'])
    # Save
    df_out.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote simulated phenotypes with {M} causal SNPs and total h2={h2_total:.3f} to {args.output}")

if __name__ == '__main__':
    main()
