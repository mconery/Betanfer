'''
simulate_null_phenotypes.py

This script will simulate null traits for testing SAIGE-GPU and estimating 
genotype correlations from beta correlations. 

***This needs to be run in the rsparsepro conda environment ***
'''

import argparse
import pandas as pd
import numpy as np

def simulate_phenotypes():
    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Simulate phenotypes for genomic data')
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    parser.add_argument('-n', '--num-pheno', type=int, default=1, 
                       help='Number of phenotypes to simulate (default: 1)')
    parser.add_argument('-t', '--type', choices=['binary', 'continuous'], default='continuous',
                       help='Phenotype type (default: continuous)')
    parser.add_argument('-p', '--param', type=float, default=None,
                       help='Simulation parameter (case prob for binary, SD for continuous)')
    
    args = parser.parse_args()
    
    # Set default parameters based on phenotype type
    if args.param is None:
        args.param = 0.5 if args.type == 'binary' else 1.0
    
    # Validate parameters
    if args.type == 'binary' and not (0 <= args.param <= 1):
        raise ValueError("Case probability must be between 0 and 1")
    if args.type == 'continuous' and args.param <= 0:
        raise ValueError("Standard deviation must be positive")
    
    # Read input data
    try:
        df = pd.read_csv(args.input, sep='\s+', header=0)
    except FileNotFoundError:
        raise SystemExit(f"Error: Input file {args.input} not found")
    
    # Simulate phenotypes
    pheno_cols = []
    for i in range(1, args.num_pheno + 1):
        col_name = f'Pheno{i}'
        if args.type == 'binary':
            pheno = np.random.binomial(n=1, p=args.param, size=len(df))
        else:
            pheno = np.random.normal(loc=0, scale=args.param, size=len(df))
        pheno_cols.append(pd.Series(pheno, name=col_name))
    
    # Combine with original data and save
    result_df = pd.concat([df] + pheno_cols, axis=1)
    result_df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')

if __name__ == '__main__':
    simulate_phenotypes()
