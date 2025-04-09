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

#Define help function
def help(exit_num=0):
    print("""Usage: python script.py [OPTIONS]
Options:
  --plink         Path to PLINK bed file (required)
  --pairs         Number of variant pairs to test (required)
  --traits        Number of traits to simulate (required)
  --outdir        Output directory (required)
  --keep          File with individuals to keep (optional)
  --trait_type    Trait type: continuous/binary (default: continuous; optional)
  --case_prob     Case probability for binary traits (default: 0.5; optional)
  --max_dist      Maximum distance of variant pairs (default: 1Mb; optional)
  --n_causal      Number of causal variants (0-5, default: 0; optional)
  --min_distance  Causal variant minimum distance from tested pairs (default: 10Mb; optional)
  --effect_size   Per-allele effect size (default: 0.1; optional)
  --keep_conv     Keep file naming convention (default: 1000G.<POP>.txt; optional)
  --plink_conv    Plink file naming convention (default: <CHROMOSOME>_10k_sample.bed; optional)
  -h, --help      Show this help message""")
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
            ["help", "plink=", "pairs=", "traits=", "keep=", "outdir=", 
             "max_dist=", "trait_type=", "case_prob=", "keep_conv=", "plink_conv=",
             "n_causal=", "min_distance=", "effect_size="]
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
    case_prob = 0.5
    trait_type = 'continuous'
    n_causal = 0
    min_distance = 10000000  # 10Mb
    effect_size = 0.1
    keep_conv = "1000G.<POP>.txt"
    plink_conv = "<CHROMOSOME>_10k_sample.bed"

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
        elif o == "--trait_type":
            trait_type = a.lower()
        elif o == "--case_prob":
            case_prob = float(a)
            if not (0 < case_prob < 1):
                print("ERROR: --case_prob must be between 0 and 1")
                sys.exit(1)
        elif o == "--n_causal":
            n_causal = int(a)
        elif o == "--min_distance":
            min_distance = int(a)
        elif o == "--effect_size":
            effect_size = float(a)
            if not (effect_size > 0):
                print("ERROR: Specified effect size must be >0.")
                sys.exit(1)
        elif o == "--keep_conv":
            keep_conv = a
        elif o == "--plink_conv":
            plink_conv = a
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
           return check_num
        except ValueError:
            print(f"ERROR: Invalid {description} given.")
            sys.exit(1)
    max_distance = check_int(max_distance, "maximum distance of variant pairs")
    n_traits = check_int(n_traits, "number of traits")
    n_pairs = check_int(n_pairs, "number of variant pairs")
    n_causal = check_int(n_causal, "number of causal variants")
    min_distance = check_int(min_distance, "minimum distance of causal variants")
    
    # Validate trait type
    if trait_type not in ['continuous', 'binary']:
        print("ERROR: Invalid trait type give. Select continuous/binary.")
        sys.exit(2)
    
    #Check naming conventions for needed parameters
    if '<POP>' not in keep_conv and keep_file:
        print("ERROR: Invalid keep file naming convention given. Convention must include <POP>.")
        sys.exit(2)
    if '<CHROMOSOME>' not in plink_conv and keep_file:
        print("ERROR: Invalid plink file naming convention given. Convention must include <CHROMOSOME>.")
        sys.exit(2)
    
    #Set validation on number of causal variants
    if n_causal < 0 or n_causal > 5:
        print("ERROR: --n_causal must be 0-5")
        sys.exit(1)
    
    # Read and filter genotype data for samples
    geno = read_plink1_bin(plink_path, verbose=False)
    if keep_file: #Checks for whether keep file given
        #Collect population name if filtering
        keep_name_remnant = keep_conv.split("<POP>")
        pop_name = os.path.basename(keep_file)
        for remnant in keep_name_remnant:
            pop_name = pop_name.replace(remnant, "")
        #Filter genotype file
        keep_ids = pd.read_csv(keep_file, sep='\t')['#IID'].values
        geno = geno.sel(sample=keep_ids)
        print(f"Retained {len(geno.sample)} samples after filtering")
    else:
        pop_name='ALL'
    
    #Filter genotype data to remove monomorphic SNPs
    print("Removing monomorphic SNPs in filtered samples...")
    variances = np.var(geno.values, axis=0)
    monomorphic_mask = variances < 1e-6
    geno = geno.sel(variant=geno.variant.values[~monomorphic_mask])
    print(f"Retained {len(geno.variant)} polymorphic SNPs after filtering")
    
    #Get number of samples
    n_samples = len(geno)
    #Extract chromosome
    plink_name_remnant = plink_conv.split("<CHROMOSOME>")
    chromo = os.path.basename(plink_path)
    for remnant in plink_name_remnant:
        chromo = chromo.replace(remnant, "")
        
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
    
    #Create a function to select the causal region
    def select_causal_region(positions, variant_ids, n_causal, region_size, min_buffer):
        """Select causal variants from a random region with buffer space"""
        if n_causal == 0:
            return np.array([]), (0, 0)
        # Find viable regions considering buffer space
        chrom_length = positions[-1]
        viable_regions = []
        # Scan chromosome in 100kbp steps
        for start in range(0, chrom_length - region_size, 100_000):
            end = start + region_size
            if end + min_buffer > chrom_length:
                break  # Can't have buffer after region
            # Check for enough variants in region
            in_region = (positions >= start) & (positions <= end)
            if np.sum(in_region) >= n_causal:
                viable_regions.append((start, end))
        if not viable_regions:
            raise ValueError(f"No regions of {region_size/1e6}Mb with {min_buffer/1e6}Mb buffer fit on chromosome")
        # Randomly select a viable region
        region_start, region_end = random.choice(viable_regions)
        # Select causal variants within region
        region_indices = np.where((positions >= region_start) & (positions <= region_end))[0]
        causal_indices = np.random.choice(region_indices, n_causal, replace=False)
        return causal_indices, (region_start, region_end)

    #Create a function to filter out variant pairs too close to causal variants
    def filter_pairs(all_pairs, positions, variant_ids, region_coords, min_distance):
        """Filter pairs to stay outside causal region buffer"""
        region_start, region_end = region_coords
        buffer_start = region_start - min_distance
        buffer_end = region_end + min_distance
        # Precompute variant positions
        var_pos_map = {var: pos for var, pos in zip(variant_ids, positions)}
        valid_pairs = []
        for v1, v2 in all_pairs:
            pos1 = var_pos_map[v1]
            pos2 = var_pos_map[v2]
            # Check if both variants are outside buffer zone
            if (pos1 < buffer_start or pos1 > buffer_end) and \
               (pos2 < buffer_start or pos2 > buffer_end):
                valid_pairs.append((v1, v2))
        return valid_pairs    

    # Multithreaded pair finding
    print(f"Finding pairs within {int(max_distance)}bp...")
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
    print(f"NOTE: Found {len(all_pairs)} pairs in total")
    
    # After pair generation but before sampling:
    if n_causal > 0:
        # Select causal region with buffer space
        causal_indices, region_coords = select_causal_region(
            positions, variant_ids, n_causal, max_distance, min_distance
        )
        # Filter pairs outside causal buffer
        all_pairs = filter_pairs(all_pairs, positions, variant_ids, region_coords, min_distance)
        # Ensure sufficient pairs remain
        if len(all_pairs) < n_pairs:
            msg = (f"Only {len(all_pairs)} valid pairs remain after causal region filtering. "
                   f"Try reducing --min_distance (current: {min_distance/1e6}Mb)")
            raise ValueError(msg)
        print(f"NOTE: Found {len(all_pairs)} valid pairs after filtering for causal variants. ")
    
    #Sampling down to requested number of pairs
    if len(all_pairs) > n_pairs: #Make sure that enough pairs are left
        all_pairs = random.sample(all_pairs, n_pairs)
    else: 
        print("WARNING: Less than specified number of pairs remaining after filtering for causal region.")
    
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

    ################# TRAIT SIMULATIONS #####################
    if n_causal > 0:
        causal_geno = geno_data[:, causal_indices]  # (samples × n_causal)
        genetic_effects = causal_geno @ np.full(n_causal, effect_size)
    else:
        genetic_effects = np.zeros(n_samples)
    if trait_type == 'continuous':
        # Continuous: additive effects + noise
        traits = pd.DataFrame(
            genetic_effects[:, None] + np.random.normal(size=(n_samples, n_traits)),
            index=geno.sample.values,
            columns=[f"trait_{i}" for i in range(n_traits)]
        )
    elif trait_type == 'binary':
        # Binary: logistic effects
        intercept = np.log(case_prob / (1 - case_prob))
        log_odds = intercept + genetic_effects[:, None]
        probs = 1 / (1 + np.exp(-log_odds))
        
        traits = pd.DataFrame(
            np.random.binomial(1, probs, (n_samples, n_traits)),
            index=geno.sample.values,
            columns=[f"trait_{i}" for i in range(n_traits)]
        )
    
    print("Running associations...")
    Y = traits.values  # Convert to numpy array once
    n_workers = os.cpu_count()
    variant_indices = np.arange(len(variant_ids))
    chunk_size = max(100, len(variant_ids) // (n_workers * 4))
    #Create function that runs logistic regression quickly
    def calculate_logistic_betas_fast(X_valid, Y_valid, max_iter=1000, tol=1e-6, intercept_init = np.log(case_prob / (1 - case_prob))):
        """
        Vectorized Newton-Raphson implementation for logistic regression coefficients.
        Returns only the beta coefficients (slope parameters).
        """
        X = X_valid.ravel()  # Convert to 1D array
        Y = Y_valid.T        # Transpose for column-wise operations
        n_samples, n_targets = Y.shape
        # Initialize parameters (intercept and slope)
        theta = np.zeros((2, n_targets))
        theta[0] = intercept_init
        X_sq = X ** 2  # Precompute squared terms
        for _ in range(max_iter):
            # Compute linear predictor and probabilities
            logits = theta[0] + theta[1] * X[:, None]
            p = 1 / (1 + np.exp(-logits))
            # Compute weights and residuals
            W = p * (1 - p)
            residuals = Y - p
            # Compute gradient components
            grad0 = residuals.sum(axis=0)
            grad1 = (residuals * X[:, None]).sum(axis=0)
            # Compute Hessian components
            H0 = W.sum(axis=0)
            H1 = (W * X[:, None]).sum(axis=0)
            H2 = (W * X_sq[:, None]).sum(axis=0)
            # Compute determinants and inverse Hessian components
            det = H0 * H2 - H1 ** 2
            np.divide(1, det, out=det, where=det!=0)  # Safe division
            # Parameter updates
            delta0 = (H2 * grad0 - H1 * grad1) * det
            delta1 = (-H1 * grad0 + H0 * grad1) * det
            # Update parameters and check convergence
            theta += [delta0, delta1]
            if np.max(np.abs([delta0, delta1])) < tol:
                break
        return theta[1]  # Return only slope coefficients
    #Function that runs regressions
    def process_regression_chunk(indices_chunk, n_traits=Y.shape[1]):
        """Process regression for both linear and logistic models"""
        chunk_betas = {}
        for idx in indices_chunk:
            # Genotype processing
            X = geno_data[:, idx]
            mask = ~np.isnan(X)
            X_valid = X[mask]  # Keep original values
            # Skip invalid variants
            if len(X_valid) < 10 or np.sum(X_valid**2) < 1e-6:
                chunk_betas[variant_ids[idx]] = [np.nan] * n_traits
                continue
            # TRAIT TYPE BRANCHING
            if trait_type == 'continuous':
                # CENTER FOR LINEAR REGRESSION ONLY
                X_mean = np.nanmean(X_valid)
                X_centered = X_valid - X_mean
                Y_centered = Y[mask] - Y[mask].mean(axis=0, keepdims=True)
                cov_xy = X_centered.T @ Y_centered
                var_x = np.sum(X_centered**2)
                betas = cov_xy / var_x
            elif trait_type == 'binary':
                # LOGISTIC REGRESSION IMPLEMENTATION
                Y_valid = Y[mask]
                betas = calculate_logistic_betas_fast(X_valid, Y_valid.T)
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
    out_file_name = f"{chromo}.{pop_name}.{n_traits}_traits.{n_pairs}_pairs.{f'{trait_type}.{case_prob}' if trait_type == 'binary' else f'{trait_type}'}{f'.{n_causal}_causal_{effect_size}' if n_causal > 0 else ''}.csv"
    output_path = os.path.join(outdir, out_file_name)
    results_df.to_csv(output_path, sep='\t', index=False)
    
    print("SUCCESS: Simulation complete!")

if __name__ == "__main__":
    main()
