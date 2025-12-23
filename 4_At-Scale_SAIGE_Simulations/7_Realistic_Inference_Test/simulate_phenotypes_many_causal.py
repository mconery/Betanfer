#!/usr/bin/env python3

"""

simulate_phenotypes_many_causal.py

Simulate multiple phenotypes from PLINK files with per-SNP heritability specification,
random sampling of causal variants per phenotype, MAF-based effect-size scaling,
and inclusion of principal components. Uses parallel processing for acceleration.

Requires Python 3.9+, numpy, pandas, pysnptools, and multiprocessing

"""

import argparse
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
from multiprocessing import Pool, cpu_count
import os
import sys


def simulate_single_phenotype(args_dict):
    """
    Simulate a single phenotype with randomly sampled causal variants.
    
    Parameters
    ----------
    args_dict : dict
        Dictionary containing phenotype simulation parameters
        
    Returns
    -------
    tuple
        (pheno_id, phenotypes_array, h2_total, n_causal, sampled_variants)
    """
    pheno_id = args_dict['pheno_id']
    plink_file = args_dict['plink_file']
    all_causal_variants = args_dict['all_causal_variants']
    n_causal = args_dict['n_causal']
    h2_snp_min = args_dict['h2_snp_min']
    h2_snp_max = args_dict['h2_snp_max']
    trait_type = args_dict['trait_type']
    prevalence = args_dict['prevalence']
    seed = args_dict['seed']
    
    # Set seed for reproducibility within this process
    rng = np.random.RandomState(seed + pheno_id)
    
    # Load PLINK data
    bed = Bed(plink_file, count_A1=False)
    n_samples = bed.shape[0]
    
    # Randomly sample causal variants for this phenotype (unique within phenotype)
    sampled_variants = rng.choice(all_causal_variants, size=n_causal, replace=False)
    
    # Identify indices of causal SNPs in PLINK data
    causal_indices = []
    causal_sids = []
    for idx, sid in enumerate(bed.sid):
        if sid in sampled_variants:
            causal_indices.append(idx)
            causal_sids.append(sid)
    
    if len(causal_indices) == 0:
        raise ValueError(f"Phenotype {pheno_id}: No causal variants found in PLINK data")
    
    M = len(causal_indices)
    
    # Extract genotypes for causal SNPs
    G_mat = bed[:, causal_indices].read().val  # shape (n_samples, M)
    
    # Compute allele frequencies per SNP
    p = np.mean(G_mat, axis=0) / 2.0  # length M
    
    # Draw per-SNP heritability for each variant from uniform distribution
    h2_snp_per_variant = rng.uniform(h2_snp_min, h2_snp_max, size=M)
    
    # Compute effect sizes per SNP to achieve per-SNP heritability
    beta = np.zeros(M, dtype=float)
    for i in range(M):
        if p[i] > 0 and p[i] < 1:
            beta[i] = np.sqrt(h2_snp_per_variant[i] / (2.0 * p[i] * (1.0 - p[i])))
        else:
            beta[i] = 0.0
    
    # Total trait heritability
    h2_total = np.sum(h2_snp_per_variant)
    
    if h2_total >= 1.0:
        raise ValueError(f"Phenotype {pheno_id}: Total trait heritability ({h2_total:.3f}) must be < 1")
    
    # Normalize genetic score variance to 1
    raw_scores = G_mat @ beta
    var_raw = np.var(raw_scores)
    
    if var_raw <= 0:
        raise ValueError(f"Phenotype {pheno_id}: Zero genetic variance: check variant list or MAFs")
    
    norm_factor = 1.0 / np.sqrt(var_raw)
    beta *= norm_factor
    
    # Recompute standardized genetic scores
    G_std = G_mat @ beta
    
    # Simulate phenotype
    if trait_type == 'continuous':
        # Noise variance = (1/h2_total) - 1
        noise_var = (1.0 / h2_total) - 1.0
        noise = rng.normal(0, np.sqrt(noise_var), n_samples)
        pheno = G_std + noise
    else:  # binary trait via liability threshold
        # Standardize G_std (mean 0, var 1)
        G_centered = G_std - np.mean(G_std)
        G_norm = G_centered / np.sqrt(np.var(G_centered)) if np.var(G_centered) > 0 else G_centered
        
        # Liability L = sqrt(h2_total)*G_norm + sqrt(1-h2_total)*ε
        eps = rng.normal(0, 1, n_samples)
        L = np.sqrt(h2_total) * G_norm + np.sqrt(1.0 - h2_total) * eps
        
        # Threshold for prevalence
        thresh = np.percentile(L, 100.0 * (1.0 - prevalence))
        pheno = (L > thresh).astype(int)
    
    return pheno_id, pheno, h2_total, M, sampled_variants


def main():
    parser = argparse.ArgumentParser(
        description='Simulate multiple phenotypes from PLINK files with per-SNP heritability specification'
    )
    
    parser.add_argument('--plink', required=True,
                        help='PLINK file prefix (without .bed/.bim/.fam extension)')
    parser.add_argument('--output', required=True,
                        help='Output phenotype file path')
    parser.add_argument('--variants', required=True,
                        help='File with causal variant IDs (one per line)')
    parser.add_argument('--n-phenotypes', type=int, default=2000,
                        help='Number of phenotypes to simulate (default=2000)')
    parser.add_argument('--n-causal', type=int, default=500,
                        help='Number of causal variants per phenotype (default=500)')
    parser.add_argument('--h2-snp-min', type=float, required=True,
                        help='Minimum per-SNP heritability for uniform distribution')
    parser.add_argument('--h2-snp-max', type=float, required=True,
                        help='Maximum per-SNP heritability for uniform distribution')
    parser.add_argument('--trait-type', choices=['binary', 'continuous'], required=True,
                        help='Trait type')
    parser.add_argument('--prevalence', type=float, default=0.3,
                        help='Case prevalence for binary traits (0-1, default=0.3)')
    parser.add_argument('--pc-file', required=True,
                        help='File with principal component vectors (tab: FID IID PC1 ... PC10)')
    parser.add_argument('--n-jobs', type=int, default=None,
                        help='Number of parallel jobs (default: use all available cores)')
    parser.add_argument('--seed', type=int, default=5,
                        help='Random seed for reproducibility (default=5)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print detailed progress information')
    
    args = parser.parse_args()
    
    # Validate inputs
    if args.h2_snp_min <= 0 or args.h2_snp_max >= 1:
        raise ValueError("per-SNP heritability bounds must be in (0, 1)")
    
    if args.h2_snp_min > args.h2_snp_max:
        raise ValueError("h2-snp-min must be <= h2-snp-max")
    
    if args.n_phenotypes <= 0:
        raise ValueError("n-phenotypes must be > 0")
    
    if args.n_causal <= 0:
        raise ValueError("n-causal must be > 0")
    
    # Set seed
    np.random.seed(args.seed)
    
    # Load causal variant IDs
    if args.verbose:
        print(f"Loading causal variants from {args.variants}...")
    
    causal_variants = []
    with open(args.variants) as f:
        for line in f:
            vid = line.strip()
            if vid:
                causal_variants.append(vid)
    
    n_available = len(causal_variants)
    if n_available == 0:
        raise ValueError("No causal variant IDs provided")
    
    if args.n_causal > n_available:
        raise ValueError(
            f"Requested {args.n_causal} causal variants per phenotype, "
            f"but only {n_available} available in variant list"
        )
    
    if args.verbose:
        print(f"Loaded {n_available} causal variant IDs")
        print(f"Will simulate {args.n_phenotypes} phenotypes with {args.n_causal} causal variants each")
        print(f"Per-SNP heritability sampled from U({args.h2_snp_min}, {args.h2_snp_max})")
    
    # Load PLINK data once to get sample IDs
    if args.verbose:
        print(f"Loading PLINK data from {args.plink}...")
    
    bed = Bed(args.plink, count_A1=False)
    iids = bed.iid
    n_samples = bed.shape[0]
    
    if args.verbose:
        print(f"PLINK data has {n_samples} samples")
    
    # Load PCs
    if args.verbose:
        print(f"Loading PCs from {args.pc_file}...")
    
    df_pcs = pd.read_csv(args.pc_file, sep='\t')
    
    # Prepare arguments for parallel processing
    if args.verbose:
        print("Preparing simulation tasks...")
    
    n_jobs = args.n_jobs if args.n_jobs is not None else cpu_count()
    if args.verbose:
        print(f"Using {n_jobs} parallel jobs")
    
    simulation_args = []
    for pheno_id in range(args.n_phenotypes):
        sim_args = {
            'pheno_id': pheno_id,
            'plink_file': args.plink,
            'all_causal_variants': causal_variants,
            'n_causal': args.n_causal,
            'h2_snp_min': args.h2_snp_min,
            'h2_snp_max': args.h2_snp_max,
            'trait_type': args.trait_type,
            'prevalence': args.prevalence,
            'seed': args.seed
        }
        simulation_args.append(sim_args)
    
    # Run simulations in parallel
    if args.verbose:
        print(f"Starting parallel simulation of {args.n_phenotypes} phenotypes...")
    
    results = {}
    with Pool(n_jobs) as pool:
        for pheno_id, pheno_values, h2_total, n_causal_found, sampled_vars in pool.imap_unordered(
            simulate_single_phenotype, simulation_args
        ):
            results[pheno_id] = {
                'pheno': pheno_values,
                'h2_total': h2_total,
                'n_causal': n_causal_found,
                'variants': sampled_vars
            }
            
            if args.verbose and (pheno_id + 1) % 100 == 0:
                print(f"  Completed {pheno_id + 1}/{args.n_phenotypes} phenotypes")
    
    if args.verbose:
        print("Assembling output data...")
    
    # Sort results by pheno_id
    sorted_pheno_ids = sorted(results.keys())
    
    # Build output DataFrame with phenotypes
    pheno_array = np.column_stack([results[pid]['pheno'] for pid in sorted_pheno_ids])
    pheno_cols = [f'PHENO_{i}' for i in range(args.n_phenotypes)]
    df_pheno = pd.DataFrame(pheno_array, columns=pheno_cols)
    
    # Extract PC columns from df_pcs (columns 3 onwards, assuming format: #FID IID PC1 PC2 ... PC10)
    pc_col_names = [col for col in df_pcs.columns if col.startswith('PC')]
    
    # Build output with correct column order: #FID IID PC1 PC2 ... PC10 PHENO_0 PHENO_1 ...
    df_out = pd.DataFrame({'#FID': iids[:, 0], 'IID': iids[:, 1]})
    
    # Add PC columns
    for pc_col in pc_col_names:
        df_out[pc_col] = df_pcs[pc_col].values
    
    # Add phenotype columns
    for pheno_col in pheno_cols:
        df_out[pheno_col] = df_pheno[pheno_col].values
    
    # Save
    if args.verbose:
        print(f"Writing output to {args.output}...")
    
    df_out.to_csv(args.output, sep='\t', index=False)
    
    # Compute and report statistics
    h2_totals = [results[pid]['h2_total'] for pid in sorted_pheno_ids]
    h2_mean = np.mean(h2_totals)
    h2_std = np.std(h2_totals)
    h2_min = np.min(h2_totals)
    h2_max = np.max(h2_totals)
    
    print(f"\n{'='*70}")
    print(f"Simulation Complete!")
    print(f"{'='*70}")
    print(f"Output file: {args.output}")
    print(f"Phenotypes simulated: {args.n_phenotypes}")
    print(f"Causal variants per phenotype: {args.n_causal}")
    print(f"Total trait heritability (h2_total) statistics:")
    print(f"  Mean: {h2_mean:.6f}")
    print(f"  Std:  {h2_std:.6f}")
    print(f"  Min:  {h2_min:.6f}")
    print(f"  Max:  {h2_max:.6f}")
    print(f"Output dimensions: {df_out.shape[0]} samples x {df_out.shape[1]} columns")
    print(f"{'='*70}\n")
    
    # Optionally save metadata
    metadata_file = args.output.replace('.txt', '_metadata.txt')
    with open(metadata_file, 'w') as f:
        f.write("Phenotype Simulation Metadata\n")
        f.write(f"Number of phenotypes: {args.n_phenotypes}\n")
        f.write(f"Causal variants per phenotype: {args.n_causal}\n")
        f.write(f"Trait type: {args.trait_type}\n")
        f.write(f"Per-SNP heritability range: [{args.h2_snp_min}, {args.h2_snp_max}]\n")
        f.write(f"h2_total mean: {h2_mean:.6f}\n")
        f.write(f"h2_total std: {h2_std:.6f}\n")
        f.write(f"h2_total min: {h2_min:.6f}\n")
        f.write(f"h2_total max: {h2_max:.6f}\n")
        f.write(f"Random seed: {args.seed}\n")
        f.write(f"Prevalence (binary only): {args.prevalence}\n")
    
    if args.verbose:
        print(f"Metadata saved to {metadata_file}")


if __name__ == '__main__':
    main()
