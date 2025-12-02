#!/usr/bin/env python3
"""

simulate_phenotype_common_effect_sizes.py

Simulate phenotypes from PLINK files with per-SNP heritability specification,
MAF-based effect‐size scaling, and inclusion of principal components.

"""

import argparse
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed

def main():
    parser = argparse.ArgumentParser(
        description='Simulate phenotypes from PLINK files with per-SNP heritability. '
                    'Effect sizes computed from reference population can be applied to target population.'
    )
    parser.add_argument('--plink-ref', required=True,
                        help='Reference PLINK file prefix for computing effect sizes (without .bed/.bim/.fam extension)')
    parser.add_argument('--plink-target', required=True,
                        help='Target PLINK file prefix for phenotype simulation (without .bed/.bim/.fam extension). '
                             'If omitted, defaults to --plink-ref.')
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

    # If target PLINK not specified, use reference
    plink_target = args.plink_target if args.plink_target else args.plink_ref
    # Validate per-SNP heritability
    M = 0 # will set after reading variants
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
    
    # Read PLINK genotype data from REFERENCE population
    bed_ref = Bed(args.plink_ref, count_A1=False)
    # Identify indices of causal SNPs in reference population
    causal_indices_ref = []
    causal_sids_ref = []
    for idx, sid in enumerate(bed_ref.sid):
        if sid in causal_variants:
            causal_indices_ref.append(idx)
            causal_sids_ref.append(sid)
    if len(causal_indices_ref) == 0:
        raise ValueError("No causal variants found in reference PLINK data")
    
    # Extract genotypes for causal SNPs from reference population
    G_ref = bed_ref[:, causal_indices_ref].read().val # shape (n_samples_ref, M)
    # Compute allele frequencies per SNP in REFERENCE population
    p_ref = np.mean(G_ref, axis=0) / 2.0 # length M
    
    # Compute effect sizes per SNP using REFERENCE population allele frequencies
    # to achieve target per-SNP heritability
    beta = np.zeros(M, dtype=float)
    for i in range(M):
        if p_ref[i] > 0 and p_ref[i] < 1:
            beta[i] = np.sqrt(args.per_snp_heritability / (2.0 * p_ref[i] * (1.0 - p_ref[i])))
        else:
            beta[i] = 0.0
    
    # Total trait heritability (same across populations)
    h2_total = M * args.per_snp_heritability
    if h2_total >= 1.0:
        raise ValueError(f"Total trait heritability ({h2_total:.3f}) must be < 1")
    
    # Normalize genetic score variance to 1 using reference population genotypes
    raw_scores_ref = G_ref @ beta
    var_raw_ref = np.var(raw_scores_ref)
    if var_raw_ref <= 0:
        raise ValueError("Zero genetic variance in reference: check variant list or MAFs")
    norm_factor = 1.0 / np.sqrt(var_raw_ref)
    beta *= norm_factor
    
    # Now read TARGET population PLINK data for phenotype simulation
    bed_target = Bed(plink_target, count_A1=False)
    # Identify indices of causal SNPs in target population
    causal_indices_target = []
    causal_sids_target = []
    for idx, sid in enumerate(bed_target.sid):
        if sid in causal_variants:
            causal_indices_target.append(idx)
            causal_sids_target.append(sid)
    if len(causal_indices_target) == 0:
        raise ValueError("No causal variants found in target PLINK data")
    
    # Extract genotypes for causal SNPs from target population
    G_target = bed_target[:, causal_indices_target].read().val # shape (n_samples_target, M)
    # Apply effect sizes computed from reference population to target population genotypes
    G_std_target = G_target @ beta

    # Sample IDs from target population
    iids_target = bed_target.iid
    # Simulate phenotype for TARGET population
    n = G_std_target.shape[0]
    if args.trait_type == 'continuous':
        # Noise variance = (1/h2_total) - 1
        noise_var = (1.0 / h2_total) - 1.0
        noise = np.random.normal(0, np.sqrt(noise_var), n)
        pheno = G_std_target + noise
    else: # binary trait via liability threshold
        # Standardize G_std_target (mean 0, var 1)
        G_centered = G_std_target - np.mean(G_std_target)
        G_norm = G_centered / np.sqrt(np.var(G_centered)) if np.var(G_centered) > 0 else G_centered
        
        # Liability L = sqrt(h2_total)*G_norm + sqrt(1-h2_total)*ε
        eps = np.random.normal(0, 1, n)
        L = np.sqrt(h2_total) * G_norm + np.sqrt(1.0 - h2_total) * eps
        
        # Threshold for prevalence
        thresh = np.percentile(L, 100.0 * (1.0 - args.prevalence))
        pheno = (L > thresh).astype(int)
    # Build output DataFrame
    df_pheno = pd.DataFrame({
        '#FID': iids_target[:, 0],
        'IID': iids_target[:, 1],
        'PHENO': pheno
    })
    
    # Load PCs and merge
    df_pcs = pd.read_csv(args.pc_file, sep='\t')
    df_out = pd.merge(df_pheno, df_pcs, on=['#FID', 'IID'])
    
    # Save
    df_out.to_csv(args.output, sep='\t', index=False)
    pop_label = f"target population ({plink_target})" if args.plink_ref != plink_target else "same population"
    print(f"Wrote simulated phenotypes with {M} causal SNPs and total h2={h2_total:.3f} to {args.output}")
    print(f"Effect sizes computed from reference population ({args.plink_ref})")
    print(f"Phenotypes simulated for {pop_label}")

if __name__ == '__main__':
    main()
