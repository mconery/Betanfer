#!/usr/bin/env python3
"""
OnekG_Simulation_RandomPairs_GPU_Version.py

The purpose of this script is to test the beta correlation approach on 1000G
data using simulated traits with a SINGLE causal variant. This version:
1. Selects a causal variant based on MAF threshold
2. The causal variant affects only a specified proportion of traits
3. Randomly selects variant pairs from all remaining variants
4. Simulates continuous traits using zero-centered normal distribution
5. Calculates genotype correlations and beta correlations

It is designed to use CPUs and GPUs on Aurora.
CPUs are used for correlation calculations and GPUs for regressions.

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
    print("""Usage: python OnekG_Simulation_RandomPairs_GPU_Version.py [OPTIONS]

Options:
  --plink              Path to PLINK bed file (required)
  --pairs              Number of variant pairs to test (required)
  --traits             Number of traits to simulate (required)
  --outdir             Output directory (required)
  --keep               File with individuals to keep (optional)
  --causal_maf_threshold  Minimum MAF for causal variant selection (default: 0.25; optional)
  --trait_proportion   Proportion of traits affected by causal variant (default: 1.0; optional)
  --effect_size        Per-allele effect size (default: 0.1; optional)
  --keep_conv          Keep file naming convention (default: 1000G..txt; optional)
  --plink_conv         Plink file naming convention (default: _10k_sample.bed; optional)
  --random_seed        Random seed (default: 5; optional)
  --threads            Number of threads (default: all available; optional)
  -h, --help           Show this help message""")
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
             "causal_maf_threshold=", "trait_proportion=", "effect_size=",
             "keep_conv=", "plink_conv=", "random_seed=", "threads="]
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
    trait_proportion = 1.0
    effect_size = 0.1
    keep_conv = "1000G..txt"
    plink_conv = "_10k_sample.bed"
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
    trait_proportion = float(trait_proportion)

    # Validate MAF threshold
    if causal_maf_threshold < 0 or causal_maf_threshold > 0.5:
        print("ERROR: --causal_maf_threshold must be between 0 and 0.5")
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
    if '' not in keep_conv and keep_file:
        print("ERROR: Invalid keep file naming convention given. Convention must include .")
        sys.exit(2)

    if '' not in plink_conv and keep_file:
        print("ERROR: Invalid plink file naming convention given. Convention must include .")
        sys.exit(2)

    # Extract chromosome
    plink_name_remnant = plink_conv.split("")
    chromo = os.path.basename(plink_path)
    for remnant in plink_name_remnant:
        chromo = chromo.replace(remnant, "")

    # Extract population name
    if keep_file:  # Checks for whether keep file given
        # Collect population name if filtering
        keep_name_remnant = keep_conv.split("")
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

    # Adjust random seed so that different regions are picked for each simulation in a given pop
    random_seed = random_seed + n_traits + round(effect_size * 10)

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
    ########## STEP 4: SELECT RANDOM PAIRS FROM REMAINING VARIANTS ###############
    ###############################################################################

    print("NOTE: Identifying non-causal variants...")
    non_causal_indices = np.where(variant_ids != causal_var)[0]
    non_causal_vars = variant_ids[non_causal_indices]
    print(f"NOTE: Found {len(non_causal_vars)} non-causal variants")

    print("NOTE: Creating all possible pairs from non-causal variants...")
    all_possible_pairs = []
    for i in range(len(non_causal_vars)):
        for j in range(i+1, len(non_causal_vars)):
            all_possible_pairs.append((non_causal_vars[i], non_causal_vars[j]))

    print(f"NOTE: Created {len(all_possible_pairs)} total possible pairs")

    # Ensure sufficient pairs
    if len(all_possible_pairs) < n_pairs:
        print(f"ERROR: Only {len(all_possible_pairs)} valid pairs available, but {n_pairs} requested.")
        sys.exit(1)

    print(f"NOTE: Randomly sampling {n_pairs} pairs...")
    selected_pairs = random.sample(all_possible_pairs, n_pairs)

    ###############################################################################
    ########## STEP 5: CALCULATE GENOTYPE CORRELATIONS FOR PAIRS #################
    ###############################################################################

    print("NOTE: Calculating genotype correlations between paired variants using CPUs...")

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
                corr = np.corrcoef(g1[valid], g2[valid])[0, 1]
            chunk_results.append((v1, v2, corr))
        return chunk_results

    # Calculate LD with causal variant for paired variants
    print("NOTE: Calculating LD between paired variants and causal variant...")

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

    # Process correlations in parallel
    print(f"NOTE: Using {n_workers} threads for correlation calculations")
    pair_chunks = np.array_split(selected_pairs, n_workers*4)
    geno_corrs = []
    causal_lds = []

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        # Genotype pair correlations
        futures = [executor.submit(process_correlation_chunk, chunk) for chunk in pair_chunks]
        for future in futures:
            geno_corrs.extend(future.result())

        # Correlations with causal variant
        futures = [executor.submit(calculate_ld_with_causal_for_pairs, chunk) for chunk in pair_chunks]
        for future in futures:
            causal_lds.extend(future.result())

    ###############################################################################
    ########## STEP 6: SIMULATE TRAITS WITH ZERO-CENTERED NORMAL #################
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
    ########## STEP 7: GPU-ACCELERATED ASSOCIATION ANALYSIS ######################
    ###############################################################################

    print("NOTE: Running associations on GPUs...")

    # Create SYCL queues for all GPUs
    context = dpctl.SyclContext(available_gpus)
    queues = [dpctl.SyclQueue(context, gpu) for gpu in available_gpus]

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
            for q, chunk in zip(queues, sub_chunks):
                X = dpnp.asarray(geno_data[:,chunk], dtype=dpnp.float32, sycl_queue=q)
                Y = dpnp.asarray(traits_df.values, dtype=dpnp.float32, sycl_queue=q)
                betas = executor.submit(process_regression_chunk, X, Y, q)
                betas = pd.DataFrame(betas.result().asnumpy(), columns=traits_df.columns, index=variant_ids[chunk])
                all_betas.append(betas)
        print(f"NOTE: Finished iteration {i + 1} of {chunk_count} for beta calculations.")

    betas = pd.concat(all_betas)

    ###############################################################################
    ########## STEP 8: CALCULATE BETA CORRELATIONS ###############################
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
    ########## STEP 9: COMBINE AND OUTPUT RESULTS ################################
    ###############################################################################

    results_df = pd.DataFrame(selected_pairs, columns=['var1', 'var2'])

    # Add MAF values
    results_df['maf1'] = results_df['var1'].map(maf_mapping)
    results_df['maf2'] = results_df['var2'].map(maf_mapping)

    # Add causal variant information
    results_df['causal_var'] = causal_var
    results_df['causal_maf'] = causal_maf
    results_df['n_affected_traits'] = n_affected_traits
    results_df['trait_proportion'] = trait_proportion

    # Merge correlations
    results_df = results_df.merge(
        pd.DataFrame(geno_corrs, columns=['var1', 'var2', 'geno_corr']),
        on=['var1', 'var2']
    ).merge(
        pd.DataFrame(causal_lds, columns=['var1', 'var2', 'ld_var1_causal', 'ld_var2_causal']),
        on=['var1', 'var2']
    ).merge(
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
