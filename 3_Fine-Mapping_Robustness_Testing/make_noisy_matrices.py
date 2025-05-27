'''
make_noisy_matrices.py

The purpose of this script is to add normally distributed noise to an unmodified
correlation matrix. The noises are drawn from the 1000G results and depend on 
both the number of traits simulated and the underlying genotype correlations

***This needs to be run in the rsparsepro conda environment ***
'''

import argparse
import os
import numpy as np
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Add simulated noise to genotype correlation matrices.")
    parser.add_argument('--matrix', required=True, help="Path to input .vcor1 file")
    parser.add_argument('--stddev_file', required=True, help="Path to standard deviation CSV file")
    parser.add_argument('--outdir', required=True, help="Directory for output files")
    parser.add_argument('--seed', type=int, default=5, help='Random seed for reproducibility (default: 5)')
    return parser.parse_args()

def load_stddev_table(stddev_file):
    df = pd.read_csv(stddev_file)
    df['bucket_tuple'] = df['bucket_range'].apply(lambda s: tuple(map(float, s.split('-'))))
    return df

def find_bucket(corr, buckets):
    for i, (low, high) in enumerate(buckets):
        if low <= abs(corr) < high or (high == 1.0 and abs(corr) <= high):
            return i
    raise ValueError(f"Correlation value {corr} not in any bucket")

def add_noise_to_matrix(matrix, stddevs, buckets):
    n = matrix.shape[0]
    noisy_matrix = matrix.copy()
    for i in range(n):
        for j in range(i+1, n):
            corr = matrix[i, j]
            bucket_idx = find_bucket(corr, buckets)
            std = stddevs[bucket_idx]
            noise = np.random.normal(loc=0.0, scale=std)
            new_corr = corr + noise
            new_corr = max(min(new_corr, 1.0), -1.0)
            noisy_matrix[i, j] = new_corr
            noisy_matrix[j, i] = new_corr
    return noisy_matrix

def main():
    args = parse_args()
    np.random.seed(args.seed)  # Seed set here after parsing arguments
    os.makedirs(args.outdir, exist_ok=True)

    stddev_df = load_stddev_table(args.stddev_file)
    trait_levels = sorted(stddev_df['n_traits'].unique())
    bucket_ranges = stddev_df[stddev_df['n_traits'] == trait_levels[0]]['bucket_tuple'].tolist()

    matrix = np.loadtxt(args.matrix)
    assert matrix.shape[0] == matrix.shape[1], "Matrix must be square"

    for n_traits in trait_levels:
        subset = stddev_df[stddev_df['n_traits'] == n_traits]
        stddevs = subset.sort_values('bucket_tuple')['stddev_error'].values
        noisy_matrix = add_noise_to_matrix(matrix, stddevs, bucket_ranges)

        base = os.path.basename(args.matrix)
        insert = f".{n_traits}_traits"
        if base.endswith('.unphased.vcor1'):
            outname = base.replace('.unphased.vcor1', f"{insert}.unphased.vcor1")
        else:
            outname = base + insert
        outpath = os.path.join(args.outdir, outname)

        np.savetxt(outpath, noisy_matrix, fmt="%.6f", delimiter='\t')
        print(f"NOTE: Wrote {n_traits}-trait file for {base}")
    print(f"SUCCESS: wrote all files for {base}")

if __name__ == "__main__":
    main()
