#!/usr/bin/env python3
"""
simulate_phenotype_pval_target.py

Phenotype simulation script for fine-mapping studies
Based on the approach described in "Robust fine-mapping in the presence of linkage disequilibrium mismatch"
"""

import argparse
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed
import sys
import logging
from scipy.stats import norm
import statsmodels.api as sm

def setup_logging():
    """Set up logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Simulate phenotypes for fine-mapping studies using causal variants'
    )
    
    parser.add_argument('--variant-file', required=True,
                        help='Path to file containing causal variant IDs (one per line)')
    parser.add_argument('--plink-file', required=True,
                        help='Path to PLINK bed file (without .bed extension)')
    parser.add_argument('--pc-file', required=True,
                        help='Path to principal components file')
    parser.add_argument('--output-file', required=True,
                        help='Path to output phenotype file')
    parser.add_argument('--trait-type', choices=['binary', 'continuous'], required=True, 
                        help='Trait type [binary/continuous]')
    parser.add_argument('--target-pval', type=float, default=1e-15,
                       help='Target P-Value for Simulation')
    parser.add_argument('--random-seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')
    
    return parser.parse_args()

def read_causal_variants(variant_file):
    """Read causal variant IDs from file"""
    try:
        with open(variant_file, 'r') as f:
            variants = [line.strip().split("\t")[0] for line in f if line.strip()]
        logging.info(f"Loaded {len(variants)} causal variants from {variant_file}")
        return variants
    except FileNotFoundError:
        logging.error(f"Variant file {variant_file} not found")
        sys.exit(1)

def read_principal_components(pc_file):
    """Read principal components file"""
    try:
        # Assuming the file has headers: FID, IID, PC1, PC2, ..., PC10
        pc_data = pd.read_csv(pc_file, sep='\t')
        if pc_data.shape[1] != 12:
            logging.warning(f"Expected 12 columns in PC file, found {pc_data.shape[1]}")
            
        # Ensure we have FID, IID as first two columns
        if  pc_data.columns[0] == 'FID':
            pc_data.columns[0] = "#FID"
        if pc_data.columns[0] != '#FID' or pc_data.columns[1] != 'IID':
            logging.warning("First two columns should be #FID and IID")
        
        logging.info(f"Loaded PCs for {len(pc_data)} individuals")
        return pc_data
    except Exception as e:
        logging.error(f"Error reading PC file {pc_file}: {e}")
        sys.exit(1)

def perform_individual_regressions_with_pcs(null_phenotypes, causal_genotypes, pc_matrix):
    """
    Perform separate regressions of null_phenotypes onto each causal variant
    (one at a time) with PCs as covariates, using statsmodels OLS.
    Returns arrays of beta coefficients and p-values.
    """
    
    n_causal = causal_genotypes.shape[1]
    betas = []
    pvals = []
    
    for i in range(n_causal):
        # build design matrix for variant i
        variant = causal_genotypes[:, i]
        if pc_matrix is not None:
            X = np.column_stack([variant, pc_matrix])
        else:
            X = variant.reshape(-1, 1)
        X = sm.add_constant(X, has_constant='add')  # adds intercept
        
        # fit OLS model
        model = sm.OLS(null_phenotypes, X)
        result = model.fit()
        
        # variant coefficient is at index 1 (const→0, variant→1, then PCs from 2…)
        beta_i = result.params[1]
        pval_i = result.pvalues[1]
        
        betas.append(beta_i)
        pvals.append(pval_i)
        logging.info(f"Variant {i}: beta={beta_i:.4f}, p={pval_i:.3g}")
    return np.array(betas), np.array(pvals)

def simulate_phenotype(causal_genotypes, pc_matrix, target_pval, random_seed=5):
    """
    Simulate continuous phenotypes for GWAS fine-mapping, accounting
    for PCs and targeting a given p-value for each causal variant.
    
    Args:
        causal_genotypes (np.ndarray): n×m matrix (samples × causal variants).
        pc_matrix (np.ndarray): n×k matrix of principal components.
        target_pval (float): Desired two-sided p-value per variant.
        random_seed (int): RNG seed.
        n_jobs (int): Number of parallel jobs (−1 uses all CPUs).
    
    Returns:
        phenotypes (np.ndarray): Simulated trait for n samples.
        betas (np.ndarray): Length-m vector of per-variant effect sizes.
    """
    np.random.seed(random_seed)
    n, m = causal_genotypes.shape
    # Precompute PC projector: P = I - C(C'C)^{-1}C'
    C = pc_matrix
    CtC_inv = np.linalg.inv(C.T @ C)
    # Two-sided z-score target
    z = norm.isf(target_pval / 2)
    betas = np.zeros(m)
    # Compute residualized variance and per-variant beta
    for i in range(m):
        gi = causal_genotypes[:, i]
        temp = C.T @ gi
        temp2 = CtC_inv @ temp
        ri = gi - C @ temp2
        var_ri = np.var(ri, ddof=0)
        if var_ri <= 0:
            raise ValueError(f"Zero residual variance for variant {i}")
        betas[i] = z / np.sqrt(n * var_ri)
    # Genetic effect
    genetic_effect = causal_genotypes @ betas
    resid_var = max(0.0, 1.0 - np.var(genetic_effect, ddof=0))
    noise = np.random.normal(0, np.sqrt(resid_var), n)
    phenotypes = genetic_effect + noise
    return phenotypes, betas

def main():
    setup_logging()
    args = parse_arguments()

    # Set random seed
    np.random.seed(args.random_seed)

    logging.info(f"Starting phenotype simulation with seed {args.random_seed}")

    # Read causal variants
    causal_variants = read_causal_variants(args.variant_file)

    # Read principal components
    pc_data = read_principal_components(args.pc_file)
    #Extract PC matrix
    pc_columns = [col for col in pc_data.columns if col.startswith('PC')]
    if pc_columns:
        pc_matrix = pc_data[pc_columns].values
        logging.info(f"Using {pc_matrix.shape[1]} principal components as covariates")
    else:
        pc_matrix = None
        logging.warning("No principal components found in PC file")

    # Read genotype data using pysnptools
    logging.info(f"Loading genotype data from {args.plink_file}")
    try:
        bed_reader = Bed(args.plink_file, count_A1=False)

        # Get SNP and sample information
        snp_info = bed_reader.sid  # SNP IDs
        sample_info = bed_reader.iid  # Sample IDs (FID, IID pairs)

        logging.info(f"Loaded {len(sample_info)} samples and {len(snp_info)} SNPs")

        # Find indices of causal variants
        causal_indices = []
        missing_variants = []

        for variant in causal_variants:
            try:
                idx = list(snp_info).index(variant)
                causal_indices.append(idx)
            except ValueError:
                missing_variants.append(variant)

        if missing_variants:
            logging.warning(f"Could not find {len(missing_variants)} causal variants in genotype data: {missing_variants}")

        if not causal_indices:
            logging.error("No causal variants found in genotype data")
            sys.exit(1)

        # Read genotype data
        logging.info("Reading genotype data...")
        genotype_data = bed_reader.read().val
        causal_genotypes = genotype_data[:, causal_indices]
        
        # Simulate phenotype
        phenotype, effect_sizes = simulate_phenotype(causal_genotypes, causal_indices, args.per_variant_heritability, args.trait_type, args.target_pval, pc_matrix=pc_matrix)
        # Check regression results
        perform_individual_regressions_with_pcs(phenotype, causal_genotypes, pc_matrix)
        
        # Create output dataframe
        # Match samples between genotype and PC data
        sample_df = pd.DataFrame(sample_info, columns=['#FID', 'IID'])

        # Merge with PC data
        output_df = sample_df.merge(pc_data, on=['#FID', 'IID'], how='left')

        # Add phenotype
        output_df['PHENO'] = phenotype

        # Reorder columns: FID, IID, PHENO, PC1-PC10
        pc_cols = [col for col in output_df.columns if col.startswith('PC')]
        final_cols = ['#FID', 'IID', 'PHENO'] + pc_cols
        output_df = output_df[final_cols]

        # Write output
        output_df.to_csv(args.output_file, sep='\t', index=False, na_rep='NA')

        logging.info("Phenotype simulation completed successfully")
        logging.info(f"Output written to {args.output_file}")
        logging.info(f"Phenotype statistics - samples: {len(output_df)}, mean: {np.mean(phenotype):.4f}, std: {np.std(phenotype):.4f}")

        # Log effect sizes for causal variants
        for i, (var_idx, effect) in enumerate(zip(causal_indices, effect_sizes)):
            variant_id = snp_info[var_idx]
            logging.info(f"Causal variant {i+1}: {variant_id}, effect size: {effect:.6f}")

    except Exception as e:
        logging.error(f"Error processing genotype data: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()