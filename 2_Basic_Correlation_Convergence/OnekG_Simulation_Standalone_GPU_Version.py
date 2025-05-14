'''
OnekG_Simulation_Standalone_GPU_Version.py

The purpose of this script is to test the beta correlation approach on 1000G 
data using simulated traits. It is designed to use the CPUs and GPUs on Aurora.
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
    print("""Usage: python OnekG_Simulation_Standalone.py [OPTIONS]
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
  --random_seed   Random seed (default: 5; optional)
  --threads       Number of threads (default: all available; optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)

###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND MULTI-PROCESS LOCI #############
###############################################################################

def main():
        
    # Parse command-line arguments
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "plink=", "pairs=", "traits=", "keep=", "outdir=", 
             "max_dist=", "trait_type=", "case_prob=", "keep_conv=", "plink_conv=",
             "n_causal=", "min_distance=", "effect_size=", "random_seed=", "threads="]
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
    random_seed = check_int(random_seed, "random seed")
    n_workers = check_int(n_workers, "number of threads")
    
    #Validate number of workers
    if n_workers > os.cpu_count():
        print(f"ERROR: Selected number of threads exceeds cores per node. Threads must be less than {os.cpu_count()}.")
        sys.exit(2)
    
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
    
    #Extract chromosome
    plink_name_remnant = plink_conv.split("<CHROMOSOME>")
    chromo = os.path.basename(plink_path)
    for remnant in plink_name_remnant:
        chromo = chromo.replace(remnant, "")
    
    #Extract population name
    if keep_file: #Checks for whether keep file given
        #Collect population name if filtering
        keep_name_remnant = keep_conv.split("<POP>")
        pop_name = os.path.basename(keep_file)
        for remnant in keep_name_remnant:
            pop_name = pop_name.replace(remnant, "")
        if pop_name == 'ALL':
            keep_file = None
    else:
        pop_name='ALL'
       
    #Create output file name
    out_file_name = f"{chromo}.{pop_name}.{n_traits}_traits.{n_pairs}_pairs.{f'{trait_type}_{case_prob}' if trait_type == 'binary' else f'{trait_type}'}{f'.{n_causal}_causal_{effect_size}' if n_causal > 0 else f'.{n_causal}_causal'}.csv"
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
    from bisect import bisect_right
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
    random_seed = random_seed + n_traits + round(effect_size * 10) + n_causal
    if trait_type == "binary": 
        random_seed += 1
    #Set random seeds
    np.random.seed(random_seed)
    random.seed(random_seed)
    
    # Read and filter genotype data for samples
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
    print(f"NOTE: Finding pairs within {int(max_distance)}bp...")
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
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
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
    
    ###########################################################################################################    
    ################################# CPU-Based Variant Correlation Calculations ############################
    ###########################################################################################################
    
    # Multithreaded genotype correlation calculation
    print("NOTE: Calculating genotype correlations using CPUs...")
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
    print(f"NOTE: Using {len(available_gpus)} threads per node for correlation calculations")
    pair_chunks = np.array_split(all_pairs, n_workers*4)
    geno_corrs = []
    with ThreadPoolExecutor(max_workers=len(available_gpus)) as executor:
        futures = [executor.submit(process_correlation_chunk, chunk) for chunk in pair_chunks]
        for future in futures:
            geno_corrs.extend(future.result())
    
    ###########################################################################################################
    ############################################### Calculate Trait Values ####################################
    ###########################################################################################################
    
    ################# TRAIT SIMULATIONS #####################
    if n_causal > 0:
        causal_geno = geno_data[:, causal_indices]  # (samples × n_causal)
        alt_pct = np.sum(causal_geno == 2, axis=0)/len(causal_geno)
        ref_pct = np.sum(causal_geno == 0, axis=0)/len(causal_geno)
        #check minor allele in order to set genetic effects based on minor allele
        for i in range(n_causal):    
            if alt_pct[i] > ref_pct[i]: 
                causal_geno[:,i] = 2 - causal_geno[:,i]
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
        # Find intercept for each trait to achieve target case probability (They're the same for all traits)
        def f(intercept):
            probs = 1 / (1 + np.exp(-(intercept + genetic_effects)))
            return np.mean(probs) - case_prob
        # Newton-Raphson optimization
        intercept = 0.0  # Initial guess
        for _ in range(100):
            probs = 1 / (1 + np.exp(-(intercept + genetic_effects)))
            grad = np.mean(probs * (1 - probs))
            if abs(grad) < 1e-6: break
            intercept -= f(intercept)/grad
            if abs(f(intercept)) < 1e-4: break
        # Generate traits with per-trait intercepts
        log_odds = genetic_effects + intercept
        probs = 1 / (1 + np.exp(-log_odds))
        #Reshape probs to calculate trait values for all traits
        probs = np.tile(probs.reshape(-1, 1), (1, n_traits))
        traits = pd.DataFrame(
            np.random.binomial(1, probs, (n_samples, n_traits)),
            index=geno.sample.values,
            columns=[f"trait_{i}" for i in range(n_traits)]
        )
    
    ###########################################################################################################
    ######################### GPU-Enabled Calculation of Regression Coefficients ##############################
    ###########################################################################################################
        
    print("NOTE: Running associations on GPUs...")
    # Create SYCL queues for all GPUs
    context = dpctl.SyclContext(available_gpus)  # Shared across all 6 GPUs
    queues = [dpctl.SyclQueue(context, gpu) for gpu in available_gpus]
    #Create a function to run linear regression
    def calculate_linear_regression(X_valid, Y_valid, queue):
        """Calculate beta, SE, and z-scores for linear regression"""
        X_mean = dpnp.nanmean(X_valid, axis=0)
        X_centered = X_valid - X_mean
        Y_centered = Y_valid - dpnp.nanmean(Y_valid, axis=0, keepdims=True)
        cov_matrix = dpnp.cov(X_centered.T, Y_centered.T) #Calculate covariance matrix
        cov_xy = cov_matrix[dpnp.asarray(range(X_centered.shape[1]), dtype=int, sycl_queue=queue),]
        cov_xy = cov_xy[:,dpnp.asarray(range(X_centered.shape[1], cov_matrix.shape[1]), sycl_queue=queue)]
        var_x = cov_matrix.diagonal()[dpnp.asarray(range(X_centered.shape[1]), dtype=int, sycl_queue=queue)]
        betas = cov_xy.T / var_x
        return betas.T
    
    #Create function that runs logistic regression quickly
    def calculate_logistic_betas_fast(X, Y, queue, max_iter=500, tol=1e-6):
        """
        Vectorized Newton-Raphson implementation for logistic regression coefficients.
        Returns only the beta coefficients (slope parameters).
        """
        n_samples, n_variants = X.shape
        _, n_targets = Y.shape
        # Initialize parameters (intercept and slope)
        y_mean = Y.mean(axis=0, keepdims=True)
        intercept_init = dpnp.log(y_mean / (1 - y_mean))
        theta = dpnp.zeros((2, n_variants, n_targets), dtype = dpnp.float32, sycl_queue=queue)
        theta[0] = intercept_init
        del intercept_init, y_mean
        # Precompute reusable buffers
        X_expanded = X[:, :, None]  # (samples, variants, 1)
        X_sq = dpnp.square(X)  # In-place computation
        # Reusable buffers for intermediate calculations
        logits = dpnp.empty([n_samples, n_variants, n_traits], dtype=X.dtype, sycl_queue=queue)
        p = dpnp.empty_like(logits, sycl_queue=queue)
        W = dpnp.empty_like(logits, sycl_queue=queue)
        residuals = dpnp.empty_like(logits, sycl_queue=queue)
        #Redirect output to a null file
        for iteration in range(max_iter):
            # Vectorized logit/probability calculation
            dpnp.add(theta[0][None], X_expanded * theta[1][None], out=logits)
            dpnp.negative(logits, out=logits)
            dpnp.exp(logits, out=logits)
            dpnp.add(1, logits, out=logits)
            dpnp.reciprocal(logits, out=p)
            # Weights and residuals
            dpnp.subtract(1, p, out=W)
            dpnp.multiply(p, W, out=W)
            dpnp.subtract(Y[:, None, :], p, out=residuals)
            # Batch gradient calculation
            grad0 = residuals.sum(axis=0)
            grad1 = dpnp.sum(residuals * X_expanded, axis=0)
            # Parallel Hessian components
            H0 = W.sum(axis=0)
            H1 = dpnp.sum(W * X_expanded, axis=0)
            H2 = dpnp.sum(W * X_sq[:, :, None], axis=0)
            # Stabilized matrix inversion
            det = H0 * H2 - H1**2
            dpnp.divide(1, det, out=det)
            # Simultaneous parameter updates
            delta0 = (H2 * grad0 - H1 * grad1) * det
            delta1 = (-H1 * grad0 + H0 * grad1) * det
            theta[0, :, :] = theta[0, :, :] + delta0
            theta[1, :, :] = theta[1, :, :] + delta1
            # Early exit check using reduced memory
            max_delta = max(dpnp.max(dpnp.abs(delta0)), 
                           dpnp.max(dpnp.abs(delta1)))
            if max_delta < tol:
                break
        return theta[1]
    
    #Function that runs regressions
    def process_regression_chunk(X, Y, queue, min_var_calls=10):
        #Make matrix to hold outputs
        chunk_betas = dpnp.empty([X.shape[1], Y.shape[1]], dtype=dpnp.float32, sycl_queue=queue)
        #Drop variants with sufficient missingness
        mask = ~dpnp.isnan(X)
        n_valid = dpnp.sum(mask, axis=0, dtype=int)
        #Skip invalid variants
        chunk_betas[dpnp.asarray(dpnp.where(n_valid < min_var_calls), sycl_queue=queue)[0,:],:] = dpnp.nan
        X_valid = X[:,dpnp.asarray(dpnp.where(n_valid >= min_var_calls), sycl_queue=queue)[0,:]]
        # TRAIT TYPE BRANCHING
        if trait_type == 'continuous':
            betas = calculate_linear_regression(X_valid, Y, queue)
        elif trait_type == 'binary':
            betas = calculate_logistic_betas_fast(X_valid, Y, queue)
        #Fill the chunk-wide list for non-missing variants
        chunk_betas[dpnp.asarray(dpnp.where(n_valid >= min_var_calls), sycl_queue=queue)[0,:],:] = betas
        return chunk_betas
    
    print("NOTE: Calculating betas in parallel...")
    all_betas = []
    # Split pairs across 6 GPUs
    variant_indices = np.arange(len(variant_ids))
    chunk_count = 10 #Number of times to split up gpu loads
    var_chunks = np.array_split(variant_indices, chunk_count)
    var_chunks = [np.array_split(x, len(available_gpus)) for x in var_chunks]
    # Process correlations in parallel
    print(f"NOTE: Using {len(available_gpus)} GPUs for multi-thread regressions")
    #Loop over the chunks
    for i in range(len(var_chunks)):
        #Pull the current set of chunks
        sub_chunks = var_chunks[i]
        # Submit batches to GPU queues
        with ThreadPoolExecutor(max_workers=len(queues)) as executor:
            # Split into chunks for each gpu
            for q, chunk in zip(queues, sub_chunks):
                #Extract variant genotypes
                X = dpnp.asarray(geno_data[:,chunk], dtype=dpnp.float32, sycl_queue=q)
                Y = dpnp.asarray(traits.values, dtype=dpnp.float32, sycl_queue=q)
                betas = executor.submit(process_regression_chunk, X, Y, q)
                betas = pd.DataFrame(betas.result().asnumpy(), columns=traits.columns, index=variant_ids[chunk])
                all_betas.append(betas)
        print(f"NOTE: Finished iteration {i + 1} of {chunk_count} for beta calculations.")
    #Merge dataframes together       
    betas = pd.concat(all_betas)
    
    ###########################################################################################################
    ################################# CPU-Based Calculation of Beta Correlations ##############################
    ###########################################################################################################
    
    
    # Beta Correlation Section with Parallel Processing
    print("NOTE: Calculating beta correlations...")
    
    def calculate_beta_correlations(pair_chunk):
        """Process chunk of variant pairs for beta correlations"""
        chunk_results = []
        for var1, var2 in pair_chunk:
            b1 = betas.loc[var1,:]
            b2 = betas.loc[var2,:]
            # Vectorized NaN handling
            valid = ~(np.isnan(b1) | np.isnan(b2))
            n_valid = valid.sum()
            
            if n_valid < 2:
                corr = np.nan
            else:
                corr = np.corrcoef(b1[valid], b2[valid])[0, 1]
            chunk_results.append((var1, var2, corr, sum(valid)))
        return chunk_results
        
    # Parallel processing beta setup
    pair_chunks = np.array_split(all_pairs, max(1, len(all_pairs)//chunk_size))
    beta_results = []
    with ThreadPoolExecutor(max_workers=len(available_gpus)) as executor:
        futures = [executor.submit(calculate_beta_correlations, chunk) for chunk in pair_chunks]
        for future in futures:
            beta_results.extend(future.result())
    
    ###########################################################################################################
    #################################### Combine and Output Table of Results ##################################
    ###########################################################################################################    
    
    # Merge results with genotype correlations
    results_df = pd.DataFrame(all_pairs, columns=['var1', 'var2'])
    # Add MAF values
    results_df['maf1'] = results_df['var1'].map(maf_mapping)
    results_df['maf2'] = results_df['var2'].map(maf_mapping)
    # Merge correlations
    results_df = results_df.merge(
        pd.DataFrame(geno_corrs, columns=['var1', 'var2', 'geno_corr']),
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
