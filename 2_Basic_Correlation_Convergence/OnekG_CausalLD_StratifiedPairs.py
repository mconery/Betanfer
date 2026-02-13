#!/usr/bin/env python3

"""
OnekG_Simulation_StratifiedPairs_GPU_Version.py

Modified version that:
1. Calculates genotype correlations for ALL possible variant pairs on GPU
2. Stratifies pairs into 5 LD buckets: 0-0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1.0
3. Samples equal numbers of pairs from each bucket
4. Proceeds with trait simulation and beta correlation analysis

The purpose of this script is to test the beta correlation approach on 1000G
data using simulated traits with a SINGLE causal variant, stratified by LD.

***Run in Conda finemap environment***
"""

# For timing
import time
tic = time.perf_counter()

import sys
import getopt
import os
import numpy as np

# Define help function
def help(exit_num=0):
    print("""Usage: python OnekG_Simulation_StratifiedPairs_GPU_Version.py [OPTIONS]

Options:
  --plink                   Path to PLINK bed file (required)
  --pairs                   Number of variant pairs to test (required, must be divisible by 5)
  --traits                  Number of traits to simulate (required)
  --outdir                  Output directory (required)
  --keep                    File with individuals to keep (optional)
  --causal_maf_threshold    Minimum MAF for causal variant selection (default: 0.25; optional)
  --pair_maf_threshold      Minimum MAF for variants retained for pairs (default: 0.05; optional)
  --trait_proportion        Proportion of traits affected by causal variant (default: 1.0; optional)
  --effect_size             Per-allele effect size (default: 0.1; optional)
  --keep_conv               Keep file naming convention (default: 1000G.<POP>.txt; optional)
  --plink_conv              Plink file naming convention (default: <CHROMOSOME>_10k_sample.bed; optional)
  --random_seed             Random seed (default: 5; optional)
  --threads                 Number of threads (default: all available; optional)
  -h, --help                Show this help message""")
    sys.exit(exit_num)

###############################################################################
############### COLLECT / CHECK ARGUMENTS AND MULTI-PROCESS LOCI #############
###############################################################################

def main():
    # Parse command-line arguments
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "plink=", "pairs=", "traits=", "keep=", "outdir=",
             "causal_maf_threshold=", "pair_maf_threshold=", "trait_proportion=", 
             "effect_size=", "keep_conv=", "plink_conv=", "random_seed=", "threads="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    plink_path = None
    n_pairs = None
    n_traits = None
    keep_file = None
    outdir = None
    causal_maf_threshold = 0.25
    pair_maf_threshold = 0.05
    trait_proportion = 1.0
    effect_size = 0.1
    keep_conv = "1000G.<POP>.txt"
    plink_conv = "<CHROMOSOME>_10k_sample.bed"
    random_seed = 5
    n_workers = os.cpu_count()

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
        elif o == "--causal_maf_threshold":
            causal_maf_threshold = float(a)
        elif o == "--pair_maf_threshold":
            pair_maf_threshold = float(a)
        elif o == "--trait_proportion":
            trait_proportion = float(a)
        elif o == "--effect_size":
            effect_size = float(a)
        elif o == "--keep_conv":
            keep_conv = a
        elif o == "--plink_conv":
            plink_conv = a
        elif o == "--random_seed":
            random_seed = a
        elif o == "--threads":
            n_workers = a
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

    # Numeric Value Validation Checks
    def check_int(check_num, description):
        try:
            check_num = int(check_num)
            return check_num
        except ValueError:
            print(f"ERROR: Invalid {description} given.")
            sys.exit(1)

    n_traits = check_int(n_traits, "number of traits")
    n_pairs = check_int(n_pairs, "number of variant pairs")
    effect_size = float(effect_size)
    random_seed = check_int(random_seed, "random seed")
    n_workers = check_int(n_workers, "number of threads")
    causal_maf_threshold = float(causal_maf_threshold)
    pair_maf_threshold = float(pair_maf_threshold)
    trait_proportion = float(trait_proportion)

    # Validate that n_pairs is divisible by 5
    if n_pairs % 5 != 0:
        print(f"ERROR: --pairs must be divisible by 5 for stratified sampling (given: {n_pairs})")
        sys.exit(1)

    # Validate MAF thresholds
    if causal_maf_threshold < 0 or causal_maf_threshold > 0.5:
        print("ERROR: --causal_maf_threshold must be between 0 and 0.5")
        sys.exit(1)

    if pair_maf_threshold < 0 or pair_maf_threshold > 0.5:
        print("ERROR: --pair_maf_threshold must be between 0 and 0.5")
        sys.exit(1)

    # Validate trait proportion
    if trait_proportion <= 0 or trait_proportion > 1:
        print("ERROR: --trait_proportion must be between 0 and 1")
        sys.exit(1)

    # Validate effect size
    if not (effect_size > 0):
        print("ERROR: Specified effect size must be >0.")
        sys.exit(1)

    # Validate number of workers
    if n_workers > os.cpu_count():
        print(f"ERROR: Selected number of threads exceeds cores per node. Threads must be less than {os.cpu_count()}.")
        sys.exit(2)

    # Check naming conventions for needed parameters
    if '<POP>' not in keep_conv and keep_file:
        print("ERROR: Invalid keep file naming convention given. Convention must include .")
        sys.exit(2)

    if '<CHROMOSOME>' not in plink_conv and keep_file:
        print("ERROR: Invalid plink file naming convention given. Convention must include .")
        sys.exit(2)

    # Extract chromosome
    plink_name_remnant = plink_conv.split("<CHROMOSOME>")
    chromo = os.path.basename(plink_path)
    for remnant in plink_name_remnant:
        chromo = chromo.replace(remnant, "")

    # Extract population name
    if keep_file:  # Checks for whether keep file given
        # Collect population name if filtering
        keep_name_remnant = keep_conv.split("<POP>")
        pop_name = os.path.basename(keep_file)
        for remnant in keep_name_remnant:
            pop_name = pop_name.replace(remnant, "")
        if pop_name == 'ALL':
            keep_file = None
    else:
        pop_name = 'ALL'

    # Create output file name
    out_file_name = f"{chromo}.{pop_name}.{n_traits}_traits.{n_pairs}_pairs.continuous.maf_{causal_maf_threshold}.prop_{trait_proportion}.effect_{effect_size}.csv"
    output_path = os.path.join(outdir, out_file_name)

    # Check if output file exists and exit if it does
    if os.path.exists(output_path):
        print(f"NOTE: Output file already exists for {out_file_name}")
        sys.exit(0)

    # Import remainder of necessary packages
    import pandas as pd
    from pandas_plink import read_plink1_bin
    import random
    from concurrent.futures import ThreadPoolExecutor
    import dpnp
    import dpctl
    import warnings
    warnings.simplefilter(action='ignore', category=RuntimeWarning)

    # Find and select GPUs
    def select_gpu_devices():
        """Detect and select available GPU devices"""
        available_devices = [device for device in dpctl.get_devices() if device.is_gpu]
        if not available_devices:
            print("WARNING: No GPU devices found. Falling back to CPU.")
            return None, False
        print(f"NOTE: Found {len(available_devices)} GPU devices.")
        for i, device in enumerate(available_devices):
            print(f"  Device {i}: {device.name}")
        return available_devices, True

    available_gpus, gpu_enabled = select_gpu_devices()

    if not gpu_enabled:
        print("ERROR: This script requires GPU support for all-pairs correlation calculation")
        sys.exit(1)

    # Adjust random seed so that different regions are picked for each simulation in a given pop
    random_seed = random_seed + n_traits + round(effect_size * 10) + round(trait_proportion * 10)

    # Set random seeds
    np.random.seed(random_seed)
    random.seed(random_seed)

    ###############################################################################
    ################# STEP 1: READ AND FILTER GENOTYPE DATA ######################
    ###############################################################################

    print("NOTE: Reading genotype data from PLINK files...")
    geno = read_plink1_bin(plink_path, verbose=False)

    if keep_file:
        keep_ids = pd.read_csv(keep_file, sep='\t')['#IID'].values
        geno = geno.sel(sample=keep_ids)
        print(f"NOTE: Retained {len(geno.sample)} samples after filtering")

    # Filter genotype data to remove monomorphic SNPs
    print("NOTE: Removing monomorphic SNPs in filtered samples...")
    variances = np.var(geno.values, axis=0)
    monomorphic_mask = variances < 1e-6
    geno = geno.sel(variant=geno.variant.values[~monomorphic_mask])
    print(f"NOTE: Retained {len(geno.variant)} polymorphic SNPs after filtering")

    # Get number of samples
    n_samples = len(geno)

    # Precompute genotype matrix for faster access
    print("NOTE: Precomputing genotypes...")
    geno_data = geno.values.astype(np.float32)
    variant_ids = geno['variant'].values
    positions = geno['pos'].values

    # Sort variants by position
    sort_idx = np.argsort(positions)
    positions = positions[sort_idx]
    variant_ids = variant_ids[sort_idx]
    geno_data = geno_data[:, sort_idx]

    ###############################################################################
    ################# STEP 2: CALCULATE MAF FOR ALL VARIANTS #####################
    ###############################################################################

    def calculate_maf(geno_data):
        """Calculate minor allele frequencies for all variants"""
        allele_counts = np.nansum(geno_data, axis=0)
        total_alleles = 2 * np.sum(~np.isnan(geno_data), axis=0)
        maf = np.minimum(allele_counts / total_alleles,
                        1 - (allele_counts / total_alleles))
        return maf

    print("NOTE: Calculating MAF values...")
    maf_values = calculate_maf(geno_data)
    maf_mapping = dict(zip(variant_ids, maf_values))

    ###############################################################################
    ########## STEP 3: SELECT CAUSAL VARIANT BASED ON MAF THRESHOLD ##############
    ###############################################################################

    print(f"NOTE: Selecting causal variant with MAF >= {causal_maf_threshold}...")
    eligible_causal = np.where(maf_values >= causal_maf_threshold)[0]

    if len(eligible_causal) == 0:
        print(f"ERROR: No variants with MAF >= {causal_maf_threshold}")
        sys.exit(1)

    causal_idx = np.random.choice(eligible_causal)
    causal_pos = positions[causal_idx]
    causal_var = variant_ids[causal_idx]
    causal_maf = maf_values[causal_idx]

    print(f"NOTE: Selected causal variant {causal_var} at position {causal_pos} with MAF {causal_maf:.4f}")

    ###############################################################################
    ########## NEW STEP 4: GPU-ACCELERATED ALL-PAIRS GENOTYPE CORRELATION ########
    ###############################################################################

    print("NOTE: Calculating genotype correlations for all possible variant pairs on GPUs...")
    print("WARNING: This is computationally intensive and may take significant time")

    # Get indices of eligible variants (passing pair MAF threshold, excluding causal)
    eligible_indices = np.where((maf_values >= pair_maf_threshold) & (variant_ids != causal_var))[0]
    eligible_variant_ids = variant_ids[eligible_indices]
    n_eligible = len(eligible_indices)

    print(f"NOTE: Computing correlations for {n_eligible} eligible variants")
    print(f"NOTE: Total possible pairs: {n_eligible * (n_eligible - 1) // 2:,}")

    # Extract genotype data for eligible variants only
    eligible_geno_data = geno_data[:, eligible_indices].astype(np.float32)
    
    # Create SYCL queues for all GPUs
    context = dpctl.SyclContext(available_gpus)
    queues = [dpctl.SyclQueue(context, gpu) for gpu in available_gpus]
    
    def calculate_correlation_matrix_gpu(X, queues, chunk_size=500):
        """
        Calculate Pearson correlation matrix using GPU acceleration
        with pairwise deletion for missing values.
    
        X: genotype matrix (samples x variants), may contain NaNs
        Returns: correlation matrix (variants x variants)
        """
        n_samples, n_variants = X.shape
    
        # Precompute masks and masked data on CPU once
        # mask: 1 where non-missing, 0 where NaN
        mask = ~np.isnan(X)
        X_filled = np.where(mask, X, 0).astype(np.float32)
        mask = mask.astype(np.float32)
    
        corr_matrix = np.full((n_variants, n_variants), np.nan, dtype=np.float32)
        n_chunks = int(np.ceil(n_variants / chunk_size))
    
        print(f"NOTE: Computing correlation matrix in {n_chunks} row chunks across {len(queues)} GPUs")
        for chunk_idx in range(n_chunks):
            start_idx = chunk_idx * chunk_size
            end_idx = min((chunk_idx + 1) * chunk_size, n_variants)
    
            queue = queues[chunk_idx % len(queues)]
    
            # Move chunk and full matrices for this queue
            X_chunk = dpnp.asarray(X_filled[:, start_idx:end_idx], dtype=dpnp.float32, sycl_queue=queue)
            M_chunk = dpnp.asarray(mask[:, start_idx:end_idx], dtype=dpnp.float32, sycl_queue=queue)
            X_full  = dpnp.asarray(X_filled, dtype=dpnp.float32, sycl_queue=queue)
            M_full  = dpnp.asarray(mask, dtype=dpnp.float32, sycl_queue=queue)
    
            # Pairwise counts n_ij (number of samples with both non-missing)
            N = dpnp.matmul(M_chunk.T, M_full)  # (chunk_size x n_variants)
    
            # Sums of x and y over non-missing
            sum_x  = dpnp.matmul(X_chunk.T, M_full)      # (chunk_size x n_variants)
            sum_y  = dpnp.matmul(M_chunk.T, X_full)      # (chunk_size x n_variants)
    
            # Sums of squares
            sum_x2 = dpnp.matmul((X_chunk * X_chunk).T, M_full)
            sum_y2 = dpnp.matmul(M_chunk.T, (X_full * X_full))
    
            # Sum of products
            sum_xy = dpnp.matmul(X_chunk.T, X_full)
    
            # Means
            mean_x = sum_x / N
            mean_y = sum_y / N
    
            # Variances and covariance (pairwise)
            var_x = sum_x2 - (sum_x * sum_x) / N
            var_y = sum_y2 - (sum_y * sum_y) / N
            cov_xy = sum_xy - (sum_x * sum_y) / N
    
            # Pearson r with safeguards
            denom = dpnp.sqrt(var_x * var_y)
            with dpnp.errstate(divide='ignore', invalid='ignore'):
                r = cov_xy / denom
    
            # Enforce minimum n_ij and valid denom
            min_n = 10  # or whatever threshold you prefer
            r = dpnp.where((N < min_n) | (denom <= 0), dpnp.nan, r)
    
            corr_matrix[start_idx:end_idx, :] = r.asnumpy()
    
            if (chunk_idx + 1) % 10 == 0 or (chunk_idx + 1) == n_chunks:
                print(f"NOTE: Completed {chunk_idx + 1}/{n_chunks} chunks ({100*(chunk_idx+1)/n_chunks:.1f}%)")
    
        return corr_matrix
    
    # Calculate full correlation matrix on GPU
    all_corr_matrix = calculate_correlation_matrix_gpu(eligible_geno_data, queues, chunk_size=500)
    
    print("NOTE: Correlation matrix calculation complete")

    ###############################################################################
    ########## NEW STEP 5: CREATE PAIR INDEX AND EXTRACT CORRELATIONS ###########
    ###############################################################################

    print("NOTE: Extracting upper triangle correlations...")

    # Create list of all pairs with their correlations
    all_pairs_with_corr = []

    # Only iterate over upper triangle (i < j) to avoid duplicates
    for i in range(n_eligible):
        for j in range(i + 1, n_eligible):
            corr_val = all_corr_matrix[i, j]
            if not np.isnan(corr_val):  # Only keep valid correlations
                all_pairs_with_corr.append({
                    'var1': eligible_variant_ids[i],
                    'var2': eligible_variant_ids[j],
                    'geno_corr': corr_val,
                    'abs_geno_corr': abs(corr_val),
                    'idx1': eligible_indices[i],
                    'idx2': eligible_indices[j]
                })
        
        # Progress update every 1000 variants
        if (i + 1) % 1000 == 0:
            print(f"NOTE: Processed {i + 1}/{n_eligible} variants for pair extraction")

    print(f"NOTE: Extracted {len(all_pairs_with_corr):,} valid variant pairs")

    # Convert to DataFrame for easier manipulation
    import pandas as pd
    pairs_df = pd.DataFrame(all_pairs_with_corr)

    ###############################################################################
    ########## NEW STEP 6: STRATIFY PAIRS INTO LD BUCKETS ########################
    ###############################################################################

    print("NOTE: Stratifying pairs into LD buckets based on absolute genotype correlation...")

    # Define LD buckets
    increment = 0.01
    ld_buckets = [(x/100, (x + 100*increment)/100, f"{x/100}-{(x + 100*increment)/100}") for x in range(0, 100, int(100 * increment))]

    # Assign each pair to a bucket
    def assign_bucket(abs_corr):
        for lower, upper, label in ld_buckets:
            if lower <= abs_corr < upper:
                return label
        # Handle edge case of exactly 1.0
        if abs_corr >= 1.0:
            return "0.8-1.0"
        return None

    pairs_df['ld_bucket'] = pairs_df['abs_geno_corr'].apply(assign_bucket)

    # Display bucket distribution
    print("\nNOTE: Distribution of pairs across LD buckets:")
    bucket_counts = pairs_df['ld_bucket'].value_counts().sort_index()
    for bucket_label in [b[2] for b in ld_buckets]:
        count = bucket_counts.get(bucket_label, 0)
        pct = 100 * count / len(pairs_df) if len(pairs_df) > 0 else 0
        print(f"  {bucket_label}: {count:,} pairs ({pct:.2f}%)")

    ###############################################################################
    ########## NEW STEP 7: SAMPLE PAIRS FROM EACH BUCKET #########################
    ###############################################################################

    print(f"\nNOTE: Sampling {n_pairs} pairs total ({n_pairs // 5} from each bucket)...")

    pairs_per_bucket = n_pairs // len(ld_buckets)
    selected_pairs_list = []

    for lower, upper, label in ld_buckets:
        bucket_pairs = pairs_df[pairs_df['ld_bucket'] == label]
        
        if len(bucket_pairs) < pairs_per_bucket:
            print(f"WARNING: Bucket {label} has only {len(bucket_pairs)} pairs, but {pairs_per_bucket} requested")
            print(f"WARNING: Sampling all {len(bucket_pairs)} available pairs from this bucket")
            sampled = bucket_pairs
        else:
            sampled = bucket_pairs.sample(n=pairs_per_bucket, random_state=random_seed)
        
        selected_pairs_list.append(sampled)
        print(f"NOTE: Selected {len(sampled)} pairs from bucket {label}")

    # Combine all selected pairs
    selected_pairs_df = pd.concat(selected_pairs_list, ignore_index=True)

    print(f"\nNOTE: Total selected pairs: {len(selected_pairs_df)}")

    # Convert to format expected by downstream code
    selected_pairs = list(zip(selected_pairs_df['var1'], selected_pairs_df['var2']))

    ###############################################################################
    ########## STEP 8: CALCULATE LD WITH CAUSAL VARIANT ##########################
    ###############################################################################

    print("NOTE: Calculating LD between selected pairs and causal variant...")

    from concurrent.futures import ThreadPoolExecutor

    def calculate_ld_with_causal_for_pairs(pair_chunk):
        chunk_results = []
        causal_geno = geno_data[:, causal_idx]
        
        for v1, v2 in pair_chunk:
            idx1 = np.where(variant_ids == v1)[0][0]
            idx2 = np.where(variant_ids == v2)[0][0]
            
            g1 = geno_data[:, idx1]
            valid1 = ~(np.isnan(g1) | np.isnan(causal_geno))
            if np.sum(valid1) < 10:
                ld1 = np.nan
            else:
                ld1 = np.corrcoef(g1[valid1], causal_geno[valid1])[0, 1]
            
            g2 = geno_data[:, idx2]
            valid2 = ~(np.isnan(g2) | np.isnan(causal_geno))
            if np.sum(valid2) < 10:
                ld2 = np.nan
            else:
                ld2 = np.corrcoef(g2[valid2], causal_geno[valid2])[0, 1]
            
            chunk_results.append((v1, v2, ld1, ld2))
        
        return chunk_results

    # Process in parallel using CPUs
    print(f"NOTE: Using {n_workers} threads for LD calculations with causal variant")
    pair_chunks = np.array_split(selected_pairs, n_workers * 4)
    causal_lds = []

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(calculate_ld_with_causal_for_pairs, chunk) 
                   for chunk in pair_chunks]
        for future in futures:
            causal_lds.extend(future.result())

    ###############################################################################
    ########## STEP 9: SIMULATE TRAITS WITH ZERO-CENTERED NORMAL #################
    ###############################################################################

    print("NOTE: Simulating traits...")

    # Determine which traits are affected by causal variant
    n_affected_traits = max(1, int(n_traits * trait_proportion))
    affected_trait_indices = np.random.choice(n_traits, n_affected_traits, replace=False)

    print(f"NOTE: Causal variant affects {n_affected_traits} out of {n_traits} traits")

    # Initialize all traits with zero-centered normal distribution
    traits = np.random.normal(0, 1, size=(n_samples, n_traits))

    # Add causal effect to affected traits
    causal_geno = geno_data[:, causal_idx]

    # Orient genotypes to minor allele
    alt_pct = np.sum(causal_geno == 2) / len(causal_geno)
    ref_pct = np.sum(causal_geno == 0) / len(causal_geno)
    if alt_pct > ref_pct:
        causal_geno = 2 - causal_geno

    genetic_effects = causal_geno * effect_size

    for trait_idx in affected_trait_indices:
        traits[:, trait_idx] += genetic_effects

    traits_df = pd.DataFrame(
        traits,
        index=geno.sample.values,
        columns=[f"trait_{i}" for i in range(n_traits)]
    )

    ###############################################################################
    ########## STEP 10: GPU-ACCELERATED ASSOCIATION ANALYSIS #####################
    ###############################################################################

    print("NOTE: Running associations on GPUs...")

    import dpnp

    def calculate_linear_regression(X_valid, Y_valid, queue):
        """Calculate beta for linear regression"""
        X_mean = dpnp.nanmean(X_valid, axis=0)
        X_centered = X_valid - X_mean
        Y_centered = Y_valid - dpnp.nanmean(Y_valid, axis=0, keepdims=True)
        
        cov_matrix = dpnp.cov(X_centered.T, Y_centered.T)
        cov_xy = cov_matrix[dpnp.asarray(range(X_centered.shape[1]), dtype=int, sycl_queue=queue),:]
        cov_xy = cov_xy[:,dpnp.asarray(range(X_centered.shape[1], cov_matrix.shape[1]), sycl_queue=queue)]
        var_x = cov_matrix.diagonal()[dpnp.asarray(range(X_centered.shape[1]), dtype=int, sycl_queue=queue)]
        
        betas = cov_xy.T / var_x
        return betas.T

    def process_regression_chunk(X, Y, queue, min_var_calls=10):
        chunk_betas = dpnp.empty([X.shape[1], Y.shape[1]], dtype=dpnp.float32, sycl_queue=queue)
        mask = ~dpnp.isnan(X)
        n_valid = dpnp.sum(mask, axis=0, dtype=int)
        
        chunk_betas[dpnp.asarray(dpnp.where(n_valid < min_var_calls), sycl_queue=queue)[0,:],:] = dpnp.nan
        X_valid = X[:,dpnp.asarray(dpnp.where(n_valid >= min_var_calls), sycl_queue=queue)[0,:]]
        
        betas = calculate_linear_regression(X_valid, Y, queue)
        chunk_betas[dpnp.asarray(dpnp.where(n_valid >= min_var_calls), sycl_queue=queue)[0,:],:] = betas
        
        return chunk_betas

    print("NOTE: Calculating betas in parallel...")
    all_betas = []

    # Split variants across GPUs
    variant_indices = np.arange(len(variant_ids))
    chunk_count = 20
    var_chunks = np.array_split(variant_indices, chunk_count)
    var_chunks = [np.array_split(x, len(available_gpus)) for x in var_chunks]

    print(f"NOTE: Using {len(available_gpus)} GPUs for multi-thread regressions")

    for i in range(len(var_chunks)):
        sub_chunks = var_chunks[i]
        with ThreadPoolExecutor(max_workers=len(queues)) as executor:
            futures = []
            for q, chunk in zip(queues, sub_chunks):
                X = dpnp.asarray(geno_data[:,chunk], dtype=dpnp.float32, sycl_queue=q)
                Y = dpnp.asarray(traits_df.values, dtype=dpnp.float32, sycl_queue=q)
                futures.append((executor.submit(process_regression_chunk, X, Y, q), chunk))
            
            for future, chunk in futures:
                betas = pd.DataFrame(future.result().asnumpy(), 
                                    columns=traits_df.columns, 
                                    index=variant_ids[chunk])
                all_betas.append(betas)
        
        print(f"NOTE: Finished iteration {i + 1} of {chunk_count} for beta calculations.")

    betas = pd.concat(all_betas)

    ###############################################################################
    ########## STEP 11: CALCULATE BETA CORRELATIONS ##############################
    ###############################################################################

    print("NOTE: Calculating beta correlations...")

    def calculate_beta_correlations(pair_chunk):
        """Process chunk of variant pairs for beta correlations"""
        chunk_results = []
        for var1, var2 in pair_chunk:
            b1 = betas.loc[var1,:]
            b2 = betas.loc[var2,:]
            
            valid = ~(np.isnan(b1) | np.isnan(b2))
            n_valid = valid.sum()
            
            if n_valid < 2:
                corr = np.nan
            else:
                corr = np.corrcoef(b1[valid], b2[valid])[0, 1]
            
            chunk_results.append((var1, var2, corr, sum(valid)))
        
        return chunk_results

    chunk_size = 500
    pair_chunks = np.array_split(selected_pairs, max(1, len(selected_pairs)//chunk_size))
    beta_results = []

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(calculate_beta_correlations, chunk) for chunk in pair_chunks]
        for future in futures:
            beta_results.extend(future.result())

    ###############################################################################
    ########## STEP 12: COMBINE AND OUTPUT RESULTS ###############################
    ###############################################################################

    print("NOTE: Combining results...")

    # Start with selected pairs dataframe which already has genotype correlations and LD bucket
    results_df = selected_pairs_df[['var1', 'var2', 'geno_corr', 'ld_bucket']].copy()

    # Add MAF values
    results_df['maf1'] = results_df['var1'].map(maf_mapping)
    results_df['maf2'] = results_df['var2'].map(maf_mapping)

    # Add causal variant information
    results_df['causal_var'] = causal_var
    results_df['causal_maf'] = causal_maf
    results_df['n_affected_traits'] = n_affected_traits
    results_df['trait_proportion'] = trait_proportion

    # Merge LD with causal
    results_df = results_df.merge(
        pd.DataFrame(causal_lds, columns=['var1', 'var2', 'ld_var1_causal', 'ld_var2_causal']),
        on=['var1', 'var2']
    )

    # Merge beta correlations
    results_df = results_df.merge(
        pd.DataFrame(beta_results, columns=['var1', 'var2', 'beta_corr', 'beta_traits']),
        on=['var1', 'var2']
    )

    # Write output
    print("NOTE: Writing results...")
    results_df.to_csv(output_path, sep='\t', index=False)

    toc = time.perf_counter()
    print(f"SUCCESS: Simulation complete in {toc - tic:0.4f} seconds")

if __name__ == "__main__":
    main()
