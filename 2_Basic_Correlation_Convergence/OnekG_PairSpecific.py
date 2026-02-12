'''

OnekG_Simulation_PairSpecific_GPU_Version.py

The purpose of this script is to test beta correlation patterns on 1000G data
using simulated traits with PAIR-SPECIFIC causal effects. Each pair of variants
is tested with separate trait simulations where the first variant influences a 
user-specified number of traits and the second variant influences a separate 
user-specified number of traits.

It is designed to use the CPUs and GPUs on Aurora.
CPUs are used for correlation calculations and GPUs for the regressions.

***Run in Conda finemap environment***

'''

#For timing
import time
tic = time.perf_counter()
import sys
import getopt
import os
import numpy as np

#Define help function
def help(exit_num=0):
    print("""Usage: python OnekG_Simulation_PairSpecific_GPU_Version.py [OPTIONS]

Options:
  --plink             Path to PLINK bed file (required)
  --pairs             Number of variant pairs to test (required)
  --traits            Total number of traits to simulate (required)
  --traits_var1       Number of traits causally influenced by variant 1 (required)
  --traits_var2       Number of traits causally influenced by variant 2 (required)
  --outdir            Output directory (required)
  --keep              File with individuals to keep (optional)
  --effect_size       Per-allele effect size (default: 0.1; optional)
  --keep_conv         Keep file naming convention (default: 1000G..txt; optional)
  --plink_conv        Plink file naming convention (default: _10k_sample.bed; optional)
  --random_seed       Random seed (default: 5; optional)
  --threads           Number of threads (default: all available; optional)
  -h, --help          Show this help message""")
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
            ["help", "plink=", "pairs=", "traits=", "traits_var1=", "traits_var2=",
             "keep=", "outdir=", "keep_conv=", "plink_conv=",
             "effect_size=", "random_seed=", "threads="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    plink_path = None
    n_pairs = None
    n_traits = None
    n_traits_var1 = None
    n_traits_var2 = None
    keep_file = None
    outdir = None
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
        elif o == "--traits_var1":
            n_traits_var1 = int(a)
        elif o == "--traits_var2":
            n_traits_var2 = int(a)
        elif o == "--keep":
            keep_file = a
        elif o == "--outdir":
            outdir = a
        elif o == "--effect_size":
            effect_size = float(a)
            if not (effect_size > 0):
                print("ERROR: Specified effect size must be >0.")
                sys.exit(1)
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
    if None in [plink_path, n_pairs, n_traits, n_traits_var1, n_traits_var2, outdir]:
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

    n_traits = check_int(n_traits, "number of traits")
    n_traits_var1 = check_int(n_traits_var1, "number of traits for variant 1")
    n_traits_var2 = check_int(n_traits_var2, "number of traits for variant 2")
    n_pairs = check_int(n_pairs, "number of variant pairs")
    effect_size = float(effect_size)
    random_seed = check_int(random_seed, "random seed")
    n_workers = check_int(n_workers, "number of threads")

    # Validate trait number constraints
    if n_traits_var1 + n_traits_var2 > n_traits:
        print(f"ERROR: Sum of traits_var1 ({n_traits_var1}) and traits_var2 ({n_traits_var2}) cannot exceed total traits ({n_traits})")
        sys.exit(1)

    if n_traits_var1 < 0 or n_traits_var2 < 0:
        print("ERROR: Number of traits for each variant must be non-negative")
        sys.exit(1)

    #Validate number of workers
    if n_workers > os.cpu_count():
        print(f"ERROR: Selected number of threads exceeds cores per node. Threads must be less than {os.cpu_count()}.")
        sys.exit(2)

    #Check naming conventions for needed parameters
    if '' not in keep_conv and keep_file:
        print("ERROR: Invalid keep file naming convention given. Convention must include .")
        sys.exit(2)

    if '' not in plink_conv and keep_file:
        print("ERROR: Invalid plink file naming convention given. Convention must include .")
        sys.exit(2)

    #Extract chromosome
    plink_name_remnant = plink_conv.split("")
    chromo = os.path.basename(plink_path)
    for remnant in plink_name_remnant:
        chromo = chromo.replace(remnant, "")

    #Extract population name
    if keep_file: #Checks for whether keep file given
        #Collect population name if filtering
        keep_name_remnant = keep_conv.split("")
        pop_name = os.path.basename(keep_file)
        for remnant in keep_name_remnant:
            pop_name = pop_name.replace(remnant, "")
        if pop_name == 'ALL':
            keep_file = None
    else:
        pop_name='ALL'

    #Create output file name
    out_file_name = f"{chromo}.{pop_name}.{n_traits}_traits.{n_traits_var1}v1_{n_traits_var2}v2.{n_pairs}_pairs.continuous_{effect_size}.pair_specific.csv"
    output_path = os.path.join(outdir, out_file_name)

    #Check if output file exists and exit if it does
    if os.path.exists(output_path):
        print(f"NOTE: Output file already exists for {out_file_name}")
        sys.exit(0)

    #Import remainder of necessary packages (Done here to save on processing before file existence check)
    import pandas as pd
    from pandas_plink import read_plink1_bin
    import random
    from concurrent.futures import ThreadPoolExecutor
    import dpnp
    import dpctl
    import warnings
    warnings.simplefilter(action='ignore', category=RuntimeWarning) #Silence runtime warnings

    #Find and select GPUs
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

    # Call this after argument processing
    available_gpus, gpu_enabled = select_gpu_devices()

    #Adjust random seed so that different regions are picked for each simulation in a given pop
    random_seed = random_seed + n_traits + round(effect_size * 10)

    #Set random seeds
    np.random.seed(random_seed)
    random.seed(random_seed)

    ###############################################################################
    ################# STEP 1: READ AND FILTER GENOTYPE DATA ######################
    ###############################################################################

    print("NOTE: Reading genotype data from PLINK files...")
    geno = read_plink1_bin(plink_path, verbose=False)

    if keep_file: #Checks for whether keep file given
        #Filter genotype file
        keep_ids = pd.read_csv(keep_file, sep='\t')['#IID'].values
        geno = geno.sel(sample=keep_ids)
        print(f"NOTE: Retained {len(geno.sample)} samples after filtering")

    #Filter genotype data to remove monomorphic SNPs
    print("NOTE: Removing monomorphic SNPs in filtered samples...")
    variances = np.var(geno.values, axis=0)
    monomorphic_mask = variances < 1e-6
    geno = geno.sel(variant=geno.variant.values[~monomorphic_mask])
    print(f"NOTE: Retained {len(geno.variant)} polymorphic SNPs after filtering")

    #Get number of samples
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

    #Define function to calculate maf of all remaining snps
    def calculate_maf(geno_data):
        """Calculate minor allele frequencies for all variants"""
        allele_counts = np.nansum(geno_data, axis=0)
        total_alleles = 2 * np.sum(~np.isnan(geno_data), axis=0)
        maf = np.minimum(allele_counts / total_alleles,
                        1 - (allele_counts / total_alleles))
        return maf

    #Calculate MAF for all remaining variants
    print("NOTE: Calculating MAF values...")
    maf_values = calculate_maf(geno_data)
    maf_mapping = dict(zip(variant_ids, maf_values))

    ###############################################################################
    ######### STEP 2: RANDOMLY SELECT VARIANT PAIRS (NO LD FILTERING) ############
    ###############################################################################

    print(f"NOTE: Randomly selecting {n_pairs} unique variant pairs...")
    n_variants = len(variant_ids)

    # Calculate total possible pairs
    max_possible_pairs = n_variants * (n_variants - 1) // 2

    if n_pairs > max_possible_pairs:
        print(f"ERROR: Requested {n_pairs} pairs but only {max_possible_pairs} unique pairs possible with {n_variants} variants")
        sys.exit(1)

    # Randomly sample unique pairs
    all_pairs = set()
    while len(all_pairs) < n_pairs:
        idx1, idx2 = np.random.choice(n_variants, size=2, replace=False)
        # Ensure consistent ordering (smaller index first)
        pair = (variant_ids[min(idx1, idx2)], variant_ids[max(idx1, idx2)])
        all_pairs.add(pair)

    all_pairs = list(all_pairs)
    print(f"NOTE: Successfully selected {len(all_pairs)} unique variant pairs")

    ###############################################################################
    ########### STEP 3: CALCULATE GENOTYPE CORRELATIONS FOR PAIRS ################
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
                corr = np.corrcoef(g1[valid], g2[valid])[0,1]
            chunk_results.append((v1, v2, corr))
        return chunk_results

    # Process correlations in parallel
    print(f"NOTE: Using {n_workers} threads for correlation calculations")
    pair_chunks = np.array_split(all_pairs, n_workers*4)
    geno_corrs = []

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(process_correlation_chunk, chunk) for chunk in pair_chunks]
        for future in futures:
            geno_corrs.extend(future.result())

    ###############################################################################
    ############### STEP 4: SIMULATE BASE TRAIT VALUES ############################
    ###############################################################################

    print(f"NOTE: Simulating {n_traits} continuous traits with zero-centered normal distribution...")

    # Generate base trait values (no genetic effects yet)
    base_traits = np.random.normal(loc=0.0, scale=1.0, size=(n_samples, n_traits)).astype(np.float32)

    print(f"NOTE: Base traits shape: {base_traits.shape}")

    ###############################################################################
    ####### STEP 5: PROCESS EACH PAIR WITH PAIR-SPECIFIC TRAIT SIMULATION ########
    ###############################################################################

    print("NOTE: Processing pairs with pair-specific causal effects...")

    # Create SYCL queues for all GPUs
    context = dpctl.SyclContext(available_gpus)
    queues = [dpctl.SyclQueue(context, gpu) for gpu in available_gpus]

    #Create a function to run linear regression
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

    def process_pair_batch(pair_batch, queue):
        """Process a batch of pairs on a single GPU"""
        batch_results = []

        for v1, v2 in pair_batch:
            # Get indices for this pair
            idx1 = np.where(variant_ids == v1)[0][0]
            idx2 = np.where(variant_ids == v2)[0][0]

            # Extract genotypes for this pair
            g1 = geno_data[:, idx1].reshape(-1, 1)
            g2 = geno_data[:, idx2].reshape(-1, 1)

            # Orient genotypes to minor allele
            for i, g in enumerate([g1, g2]):
                alt_pct = np.sum(g == 2) / len(g)
                ref_pct = np.sum(g == 0) / len(g)
                if alt_pct > ref_pct:
                    if i == 0:
                        g1[:, 0] = 2 - g1[:, 0]
                    else:
                        g2[:, 0] = 2 - g2[:, 0]

            # Create pair-specific trait matrix
            pair_traits = base_traits.copy()

            # Add genetic effects from variant 1 to first n_traits_var1 traits
            if n_traits_var1 > 0:
                genetic_effects_v1 = g1[:, 0] * effect_size
                pair_traits[:, :n_traits_var1] += genetic_effects_v1.reshape(-1, 1)

            # Add genetic effects from variant 2 to next n_traits_var2 traits
            if n_traits_var2 > 0:
                genetic_effects_v2 = g2[:, 0] * effect_size
                pair_traits[:, n_traits_var1:n_traits_var1+n_traits_var2] += genetic_effects_v2.reshape(-1, 1)

            # Prepare data for GPU regression
            X = dpnp.asarray(np.column_stack([g1, g2]), dtype=dpnp.float32, sycl_queue=queue)
            Y = dpnp.asarray(pair_traits, dtype=dpnp.float32, sycl_queue=queue)

            # Check for sufficient valid samples
            mask = ~dpnp.isnan(X)
            n_valid = dpnp.sum(mask, axis=0, dtype=int)

            if dpnp.any(n_valid < 10):
                # Skip this pair if insufficient data
                batch_results.append((v1, v2, np.nan, n_traits_var1 + n_traits_var2))
                continue

            # Calculate betas for both variants across all traits
            betas = calculate_linear_regression(X, Y, queue)
            betas_np = betas.asnumpy()

            # Calculate beta correlation between the two variants
            b1 = betas_np[0, :]  # Betas for variant 1
            b2 = betas_np[1, :]  # Betas for variant 2

            # Only use traits with causal effects for beta correlation
            n_causal_traits = n_traits_var1 + n_traits_var2
            if n_causal_traits >= 2:
                valid = ~(np.isnan(b1[:n_causal_traits]) | np.isnan(b2[:n_causal_traits]))
                n_valid_corr = valid.sum()

                if n_valid_corr >= 2:
                    beta_corr = np.corrcoef(b1[:n_causal_traits][valid], 
                                           b2[:n_causal_traits][valid])[0, 1]
                else:
                    beta_corr = np.nan
            else:
                beta_corr = np.nan

            batch_results.append((v1, v2, beta_corr, n_causal_traits))

        return batch_results

    # Distribute pairs across GPUs
    print(f"NOTE: Using {len(available_gpus)} GPUs for pair-specific regressions")

    # Split pairs into batches for each GPU
    pairs_per_gpu = len(all_pairs) // len(available_gpus)
    pair_batches = []
    for i in range(len(available_gpus)):
        start_idx = i * pairs_per_gpu
        if i == len(available_gpus) - 1:
            end_idx = len(all_pairs)
        else:
            end_idx = (i + 1) * pairs_per_gpu
        pair_batches.append(all_pairs[start_idx:end_idx])

    # Process batches in parallel across GPUs
    beta_results = []
    with ThreadPoolExecutor(max_workers=len(queues)) as executor:
        futures = [executor.submit(process_pair_batch, batch, queue) 
                  for batch, queue in zip(pair_batches, queues)]
        for i, future in enumerate(futures):
            results = future.result()
            beta_results.extend(results)
            print(f"NOTE: Finished GPU {i+1} of {len(queues)} ({len(results)} pairs processed)")

    ###############################################################################
    ############### STEP 6: COMBINE RESULTS AND WRITE OUTPUT ######################
    ###############################################################################

    print("NOTE: Combining results...")

    # Create results dataframe
    results_df = pd.DataFrame(all_pairs, columns=['var1', 'var2'])

    # Add MAF values
    results_df['maf1'] = results_df['var1'].map(maf_mapping)
    results_df['maf2'] = results_df['var2'].map(maf_mapping)

    # Merge genotype correlations
    results_df = results_df.merge(
        pd.DataFrame(geno_corrs, columns=['var1', 'var2', 'geno_corr']),
        on=['var1', 'var2']
    )

    # Merge beta correlations
    results_df = results_df.merge(
        pd.DataFrame(beta_results, columns=['var1', 'var2', 'beta_corr', 'n_causal_traits']),
        on=['var1', 'var2']
    )

    # Add metadata columns
    results_df['total_traits'] = n_traits
    results_df['n_traits_var1'] = n_traits_var1
    results_df['n_traits_var2'] = n_traits_var2

    # Write output
    print("NOTE: Writing results...")
    results_df.to_csv(output_path, sep='\t', index=False)

    toc = time.perf_counter()
    print(f"SUCCESS: Simulation complete in {toc - tic:0.4f} seconds")

if __name__ == "__main__":
    main()
