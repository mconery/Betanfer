'''
simulate_phenotype.py

The purpose of this script is to simulate phenotypes given a set of input files
and some simulation parameters.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
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
    parser.add_argument('--effects', required=True, help='File with causal variants and effect sizes (format: variant_id effect_size)')
    parser.add_argument('--trait-type', choices=['binary', 'continuous'], required=True, help='Trait type')
    parser.add_argument('--heritability', type=float, default=0.5, help='Heritability (0-1, default=0.5)')
    parser.add_argument('--prevalence', type=float, default=0.3, help='Case prevalence for binary traits (0-1, default=0.3)')
    parser.add_argument('--pc-file', required=True, help='File with principal component vectors (format: FID IID PC1 ... PC10)')
    
    args = parser.parse_args()
        
    # Load causal effects
    causal_effects = {}
    with open(args.effects) as f:
        for line in f:
            if line.strip():
                variant, effect = line.split()
                causal_effects[variant] = float(effect)
    
    # Get list of causal variant IDs
    causal_sids = list(causal_effects.keys())
    # Read only causal SNPs into memory
    causal_genotypes = Bed(args.plink, count_A1=False)[:,causal_sids].read().val
    
    # Calculate genetic scores
    effect_sizes = [causal_effects[v] for v in causal_effects]
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

if __name__ == '__main__':
    main()
