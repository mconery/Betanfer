'''
simulate_null_phenotypes.py

This script will simulate null traits for testing SAIGE-GPU and estimating 
genotype correlations from beta correlations. 

***This needs to be run in the rsparsepro conda environment ***
'''

import argparse
import pandas as pd
import numpy as np
from joblib import Parallel, delayed, cpu_count

def simulate_phenotypes_chunk(chunk, num_pheno, pheno_type, param, seed):
    rng = np.random.default_rng(seed)
    n = len(chunk)
    if pheno_type == 'binary':
        pheno_matrix = rng.binomial(1, param, (n, num_pheno))
    else:
        pheno_matrix = rng.normal(0, param, (n, num_pheno))
    return pd.concat([chunk.reset_index(drop=True), 
                     pd.DataFrame(pheno_matrix, columns=[f'Pheno{i+1}' for i in range(num_pheno)])], 
                     axis=1)

def write_phenotype(i, df, output_prefix):
    """Parallelized file writing function"""
    out_df = df[['FID', 'IID'] + [f'PC{j}' for j in range(1, 11)]].copy()
    out_df['y'] = df[f'Pheno{i}']
    out_file = f"{output_prefix}pheno{i}.txt"
    out_df.to_csv(out_file, sep='\t', index=False, float_format='%.6f')

def main():
    parser = argparse.ArgumentParser(description='Parallel phenotype simulator with distributed I/O')
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-n', '--num-pheno', type=int, default=1)
    parser.add_argument('-t', '--type', choices=['binary', 'continuous'], default='continuous')
    parser.add_argument('-p', '--param', type=float)
    parser.add_argument('--n-jobs', type=int, default=-1)
    args = parser.parse_args()

    # Validation and default parameters
    if args.param is None:
        args.param = 0.5 if args.type == 'binary' else 1.0
    if args.type == 'binary' and not (0 <= args.param <= 1):
        raise ValueError("Case probability must be between 0 and 1")
    if args.type == 'continuous' and args.param <= 0:
        raise ValueError("Standard deviation must be positive")

    # Read input data
    df = pd.read_csv(args.input, sep='\s+')
    n_jobs = cpu_count() if args.n_jobs == -1 else args.n_jobs

    # Parallel simulation
    chunks = np.array_split(df, n_jobs)
    seeds = [42 + i for i in range(n_jobs)]
    sim_results = Parallel(n_jobs=n_jobs)(
        delayed(simulate_phenotypes_chunk)(chunk, args.num_pheno, args.type, args.param, seed)
        for chunk, seed in zip(chunks, seeds)
    )
    final_df = pd.concat(sim_results, ignore_index=True)

    # Parallel writing using threading backend
    Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(write_phenotype)(i, final_df, args.output)
        for i in range(1, args.num_pheno + 1)
    )

if __name__ == '__main__':
    main()
