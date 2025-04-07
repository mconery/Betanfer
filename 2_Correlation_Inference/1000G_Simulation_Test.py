'''
1000G_Simulation_Test.py

The purpose of this script is to test the beta correlation approach on 1000G 
data using simulated traits. 

***Run in Conda finemap environment***
'''

import sys
import getopt
import os
import numpy as np
import pandas as pd
from pandas_plink import read_plink1_bin
import random
from concurrent.futures import ThreadPoolExecutor
from bisect import bisect_right

################################## TEST LOCATIONS #############################
'''
plink_path = "/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/1000G_Simulations/plink_files/chr21_10k_sample.bed"
n_pairs = 10000
n_traits = 1000
outdir = "/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/testing"
keep_file="/grand/projects/GeomicVar/mconery/resources/1000G/1000G.AFR.txt"
max_distance = 1e6
'''

def help(exit_num=0):
    print("""Usage: python script.py [OPTIONS]
Options:
  --plink      Path to PLINK bed file (required)
  --pairs      Number of variant pairs to test (required)
  --traits     Number of traits to simulate (required)
  --outdir     Output directory (required)
  --keep       File with individuals to keep (optional)
  --max_dist   Maximum distance of variant pairs (default: 1Mb; optional)
  -h, --help   Show this help message""")
    sys.exit(exit_num)


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND MULTI-PROCESS LOCI #############
###############################################################################

def main():
    
    #Fix a random seed
    np.random.seed(5)
    
    # Parse command-line arguments
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "plink=", "pairs=", "traits=", "keep=", "outdir=", "max_dist="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    plink_path = None
    n_pairs = None
    n_traits = None
    keep_file = None
    outdir = None
    max_distance = 1e6

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
            sys.exit()
        elif o == "--plink":
            plink_path = a
        elif o == "--pairs":
            n_pairs = int(a)
        elif o == "--traits":
            n_traits = int(a)
        elif o == "--keep":
            keep_file = a
        elif o == "--outdir":
            outdir = a
        elif  o == "--max_dist":
            max_distance = a
        else:
            assert False, "WARNING: Unhandled input option"
    
    # Validate required arguments
    if None in [plink_path, n_pairs, n_traits, outdir]:
        print("ERROR: Missing required arguments!")
        sys.exit(2)

    # Path Validation checks
    if not all(os.path.exists(f) for f in [plink_path, outdir]):
        print("ERROR: Given output directory or plink file not found")
        sys.exit(1)
    
    #Numeric Value Validation Checks
    def check_int(check_num, description):
        try:
           check_num = int(check_num)
        except ValueError:
            print(f"ERROR: Invalid {description} given.")
            sys.exit(1)
    check_int(max_distance, "maximum distance of variant pairs")
    check_int(max_distance, "number of traits")
    check_int(max_distance, "number of variant pairs")
    
    # Read and filter genotype data
    geno = read_plink1_bin(plink_path, verbose=False)
    if keep_file: #Checks for whether keep file given
        keep_ids = pd.read_csv(keep_file, sep='\t')['#IID'].values
        geno = geno.sel(sample=keep_ids)
        #Collect population name if filtering
        pop_name = os.path.basename(keep_file).lstrip('1000G.').rstrip('.txt')
    else:
        pop_name='ALL'
        
    # Precompute genotype matrix for faster access
    print("Precomputing genotypes...")
    geno_data = geno.values.astype(np.float32)
    variant_ids = geno['variant'].values
    positions = geno['pos'].values    
    
    # Sort variants by position
    sort_idx = np.argsort(positions)
    positions = positions[sort_idx]
    variant_ids = variant_ids[sort_idx]
    geno_data = geno_data[:, sort_idx]
    
    # Sort variants by position
    sort_idx = np.argsort(positions)
    positions = positions[sort_idx]
    variant_ids = variant_ids[sort_idx]
    geno_data = geno_data[:, sort_idx][:,:10000]  # Select first 10k
    
    # Multithreaded pair finding
    print(f"Finding pairs within {max_distance}bp...")
    def find_pairs_chunk(start, end):
        chunk_pairs = []
        pos_array = positions
        for i in range(start, end):
            max_pos = positions[i] + max_distance
            j_max = bisect_right(pos_array, max_pos, i+1)
            chunk_pairs.extend([(variant_ids[i], variant_ids[j]) for j in range(i+1, j_max)])
        return chunk_pairs

    n_vars = len(variant_ids)
    chunk_size = 500
    all_pairs = []
    with ThreadPoolExecutor() as executor:
        futures = []
        for start in range(0, n_vars, chunk_size):
            end = min(start + chunk_size, n_vars)
            futures.append(executor.submit(find_pairs_chunk, start, end))
        for future in futures:
            all_pairs.extend(future.result())
    print(f"Found {len(all_pairs)} valid pairs")
    
    #Sampling down to requested number of pairs
    all_pairs = random.sample(all_pairs, n_pairs)
    
    # Multithreaded genotype correlation calculation
    print("Calculating genotype correlations...")
    def process_correlation_chunk(pair_chunk):
        chunk_results = []
        for v1, v2 in pair_chunk:
            idx1 = np.where(variant_ids == v1)[0][0]
            idx2 = np.where(variant_ids == v2)[0][0]
            g1 = geno_data[:, idx1]
            g2 = geno_data[:, idx2]
            
            valid = ~(np.isnan(g1) | np.isnan(g2))
            if np.sum(valid) < 10:
                corr = np.nan
            else:
                corr = np.corrcoef(g1[valid], g2[valid])[0,1]
            
            chunk_results.append((v1, v2, corr))
        return chunk_results
    
    # Process correlations in parallel
    n_workers = os.cpu_count()
    pair_chunks = np.array_split(all_pairs, n_workers*4)
    geno_corrs = []
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(process_correlation_chunk, chunk) for chunk in pair_chunks]
        for future in futures:
            geno_corrs.extend(future.result())
    
    ################# ADDRESS TRAIT SIMULATIONS #####################
    # Simulate traits using filtered samples
    print("Simulating traits...")
    n_samples = len(geno.sample)
    traits = pd.DataFrame(
        np.random.normal(size=(n_samples, n_traits)),
        index=geno.sample.values,
        columns=[f"trait_{i}" for i in range(n_traits)]
    )
    
    print("Running associations...")
    Y = traits.values  # Convert to numpy array once
    n_workers = os.cpu_count()
    variant_indices = np.arange(len(variant_ids))
    chunk_size = max(100, len(variant_ids) // (n_workers * 4))
    #Function that runs regressions
    def process_regression_chunk(indices_chunk):
        chunk_betas = {}
        for idx in indices_chunk:
            #Extract Genotypes and center X/Y matrices
            X = geno_data[:, idx]
            X_centered = X - X.mean(axis=0, keepdims=True, )
            Y_centered = Y - Y.mean(axis=0, keepdims=True)
            mask = ~np.isnan(X_centered)
            X_valid = X_centered[mask]
            
            if len(X_valid) < 10:
                chunk_betas[variant_ids[idx]] = [np.nan] * n_traits
                continue
                
            Y_valid = Y_centered[mask]
            cov_xy = X_valid.T @ Y_valid
            var_x = np.sum(X_valid**2, axis=0)
            
            if var_x < 1e-6:  # Handle near-zero variance
                chunk_betas[variant_ids[idx]] = [np.nan] * n_traits
                continue
                
            betas = cov_xy / var_x
            chunk_betas[variant_ids[idx]] = betas
            
        return chunk_betas

    print("Calculating betas in parallel...")
    all_betas = {}
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        # Split into index-based chunks for zero lookups
        futures = []
        for i in range(0, len(variant_indices), chunk_size):
            chunk = variant_indices[i:i+chunk_size]
            futures.append(executor.submit(process_regression_chunk, chunk))
        for future in futures:
            all_betas.update(future.result())
    
    betas = pd.DataFrame.from_dict(all_betas, orient='index', columns=traits.columns)
    
    # Beta Correlation Section with Parallel Processing
    print("Calculating beta correlations...")
    # Precompute mapping and convert to numpy array
    var_to_idx = {var: idx for idx, var in enumerate(variant_ids)}
    beta_matrix = betas.values.astype(np.float32)  # traits × variants
    
    def calculate_beta_correlations(pair_chunk):
        """Process chunk of variant pairs for beta correlations"""
        chunk_results = []
        for var1, var2 in pair_chunk:
            try:
                i = var_to_idx[var1]
                j = var_to_idx[var2]
            except KeyError:
                continue  # Skip invalid pairs (shouldn't occur with proper filtering)
            b1 = beta_matrix[i]
            b2 = beta_matrix[j]
            # Vectorized NaN handling
            valid = ~(np.isnan(b1) | np.isnan(b2))
            n_valid = valid.sum()
            
            if n_valid < 2:
                corr = np.nan
            else:
                corr = np.corrcoef(b1[valid], b2[valid])[0, 1]
            chunk_results.append((var1, var2, corr))
        return chunk_results

    # Parallel processing setup
    pair_chunks = np.array_split(all_pairs, max(1, len(all_pairs)//chunk_size))
    beta_results = []
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(calculate_beta_correlations, chunk): chunk 
                  for chunk in pair_chunks}
        for future in futures:
            beta_results.extend(future.result())
    
    # Merge results with genotype correlations
    results_df = pd.DataFrame(all_pairs, columns=['var1', 'var2'])
    results_df = results_df.merge(
        pd.DataFrame(geno_corrs, columns=['var1', 'var2', 'geno_corr']),
        on=['var1', 'var2']
    ).merge(
        pd.DataFrame(beta_results, columns=['var1', 'var2', 'beta_corr']),
        on=['var1', 'var2']
    )
    
    # Write output
    print("Writing results...")
    out_file_name = f"{pop_name}.{n_traits}_traits.{n_pairs}_pairs.csv"
    output_path = os.path.join(outdir, out_file_name)
    results_df.to_csv(output_path, sep='\t', index=False)
    
    print("SUCCESS: Simulation complete!")

if __name__ == "__main__":
    main()
