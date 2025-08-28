'''
simulate_phenotype.py

The purpose of this script is to simulate phenotypes given a set of input files
and some simulation parameters.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***

Updated version that accounts for allele frequency when setting effect sizes
'''

#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from pysnptools.snpreader import Bed

def main():
    parser = argparse.ArgumentParser(description='Simulate phenotypes from PLINK files and include principal components')
    parser.add_argument('--plink', required=True, help='PLINK file prefix (without .bed/.bim/.fam extension)')
    parser.add_argument('--output', required=True, help='Output phenotype file path')
    parser.add_argument('--variants', required=True, help='File with causal variants and effect sizes (format: variant_id effect_size)')
    parser.add_argument('--trait-type', choices=['binary', 'continuous'], required=True, help='Trait type')
    parser.add_argument('--heritability', type=float, default=0.5, help='Heritability (0-1, default=0.5)')
    parser.add_argument('--prevalence', type=float, default=0.3, help='Case prevalence for binary traits (0-1, default=0.3)')
    parser.add_argument('--pc-file', required=True, help='File with principal component vectors (format: FID IID PC1 ... PC10)')
    parser.add_argument('--base-effect-size', type=float, default=0.1, help="Base effect size for common variants (MAF=0.5; default 0.1)")
    parser.add_argument('--maf-adjustment', choices=['none', 'equal_power', 'equal_variance'], default='equal_power', 
                       help='Method to adjust effect sizes based on MAF (default: equal_power)')
    parser.add_argument('--seed', type=int, default=5, help='Random seed (default=5)')
    
    args = parser.parse_args()
    
    # Set random seed for numpy
    np.random.seed(args.seed)
    
    # Load causal variants list
    causal_variants = []
    with open(args.variants) as f:
        for line in f:
            variant = line.strip()
            if variant:  # Skip empty lines
                causal_variants.append(variant)
    
    # Read PLINK data
    temp = Bed(args.plink, count_A1=False)
    
    # Find indices of causal variants in the PLINK data
    causal_indices = []
    causal_sids = []
    for i, sid in enumerate(temp.sid):
        if sid in causal_variants:
            causal_indices.append(i)
            causal_sids.append(sid)
    
    if len(causal_indices) == 0:
        raise ValueError("No causal variants found in PLINK data")
    
    # Read only causal SNPs into memory
    causal_genotypes = temp[:, causal_indices].read().val
    
    # Extract sample ids
    iids = temp.iid
    
    # Calculate allele frequencies for causal variants
    allele_freqs = np.mean(causal_genotypes, axis=0) / 2.0
    
    # Adjust effect sizes based on allele frequency
    causal_effects = {}
    
    if args.maf_adjustment == 'none':
        # Use same effect size for all variants
        for sid in causal_sids:
            causal_effects[sid] = args.base_effect_size
            
    elif args.maf_adjustment == 'equal_power':
        # Scale effect sizes to achieve equal statistical power
        # beta ∝ 1/sqrt(2*p*(1-p)) for equal power
        for i, sid in enumerate(causal_sids):
            p = allele_freqs[i]
            # Prevent division by zero for monomorphic variants
            if p == 0 or p == 1:
                causal_effects[sid] = 0
            else:
                # Scale relative to MAF=0.5 (maximum variance)
                scaling_factor = 1.0 / np.sqrt(2 * p * (1 - p))
                reference_scaling = 1.0 / np.sqrt(2 * 0.5 * (1 - 0.5))  # = 1.0
                causal_effects[sid] = args.base_effect_size * scaling_factor / reference_scaling
                
    elif args.maf_adjustment == 'equal_variance':
        # Scale effect sizes to achieve equal variance explained
        # Variance = 2*p*(1-p)*beta^2, so beta ∝ 1/sqrt(2*p*(1-p)) for equal variance
        for i, sid in enumerate(causal_sids):
            p = allele_freqs[i]
            if p == 0 or p == 1:
                causal_effects[sid] = 0
            else:
                # This is actually the same as equal_power scaling
                scaling_factor = 1.0 / np.sqrt(2 * p * (1 - p))
                reference_scaling = 1.0 / np.sqrt(2 * 0.5 * (1 - 0.5))
                causal_effects[sid] = args.base_effect_size * scaling_factor / reference_scaling
    
    # Create effect sizes array in the same order as genotypes
    effect_sizes = [causal_effects[sid] for sid in causal_sids]
    
    # Print summary of effect sizes and allele frequencies
    print(f"Summary of {len(causal_sids)} causal variants:")
    print(f"MAF adjustment method: {args.maf_adjustment}")
    print("Variant\tMAF\tEffect_Size\tVariance_Explained")
    for i, sid in enumerate(causal_sids):
        p = allele_freqs[i]
        beta = effect_sizes[i]
        var_exp = 2 * p * (1 - p) * beta**2
        print(f"{sid}\t{p:.4f}\t{beta:.4f}\t{var_exp:.6f}")
    
    # Calculate genetic scores
    genetic_score = causal_genotypes @ np.array(effect_sizes)
    
    # Simulate phenotype
    if args.trait_type == 'continuous':
        genetic_var = np.var(genetic_score)
        noise_var = genetic_var * (1/args.heritability - 1) if args.heritability > 0 else 1
        noise = np.random.normal(0, np.sqrt(noise_var), len(genetic_score))
        phenotype = genetic_score + noise
        
    elif args.trait_type == 'binary':
        genetic_score_std = (genetic_score - np.mean(genetic_score)) / np.std(genetic_score)
        liability = np.sqrt(args.heritability) * genetic_score_std + \
                   np.sqrt(1 - args.heritability) * np.random.normal(0, 1, len(genetic_score))
        threshold = np.percentile(liability, 100 * (1 - args.prevalence))
        phenotype = (liability > threshold).astype(int)
    
    # Save results
    df_pheno = pd.DataFrame({
        '#FID': iids[:, 0],
        'IID': iids[:, 1],
        'PHENO': phenotype
    })

    # Load principal components
    df_pcs = pd.read_csv(args.pc_file, sep='\t')
    
    # Merge phenotype and PCs by FID and IID
    df_out = pd.merge(df_pheno, df_pcs, on=['#FID', 'IID'])

    # Save final output
    df_out.to_csv(args.output, sep='\t', index=False)
    print(f"NOTE: Wrote {args.output} to file.")
    
    # Print final statistics
    print(f"\nSimulation Statistics:")
    print(f"Total genetic variance: {np.var(genetic_score):.6f}")
    print(f"Target heritability: {args.heritability}")
    print(f"Number of individuals: {len(iids)}")

if __name__ == '__main__':
    main()