'''

OnekG_Simulation_CausalLD_GPU_Version.py

The purpose of this script is to test the beta correlation approach on 1000G
data using simulated traits with a SINGLE causal variant. This version selects
variant pairs based on LD with the causal variant: pairs consist of one variant
with low LD to the causal variant and one variant with high LD to the causal variant.
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
    print("""Usage: python OnekG_Simulation_CausalLD_GPU_Version.py [OPTIONS]

Options:
  --plink             Path to PLINK bed file (required)
  --pairs             Number of variant pairs to test (required)
  --traits            Number of traits to simulate (required)
  --outdir            Output directory (required)
  --keep              File with individuals to keep (optional)
  --trait_type        Trait type: continuous/binary (default: continuous; optional)
  --case_prob         Case probability for binary traits (default: 0.5; optional)
  --low_ld_threshold  Maximum absolute LD for 'not in LD' variants (default: 0.01; optional)
  --high_ld_threshold Minimum absolute LD for 'in LD' variants (default: 0.1; optional)
  --effect_size       Per-allele effect size (default: 0.1; optional)
  --trait_subset      Subset of traits to have causal effects (default 0; optional)
  --keep_conv         Keep file naming convention (default: 1000G.<POP>.txt; optional)
  --plink_conv        Plink file naming convention (default: <CHROMOSOME>_10k_sample.bed; optional)
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
            ["help", "plink=", "pairs=", "traits=", "keep=", "outdir=",
             "trait_type=", "case_prob=", "keep_conv=", "plink_conv=",
             "low_ld_threshold=", "high_ld_threshold=", "effect_size=", 
             "random_seed=", "trait_subset=", "variable_effects=", "threads="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    
    plink_path = None
    n_pairs = None
    n_traits = None
    keep_file = None
    outdir = None
    low_ld_threshold = 0.01
    high_ld_threshold = 0.1
    case_prob = 0.5
    trait_type = 'continuous'
    effect_size = 0.1
    keep_conv = "1000G.<POP>.txt"
    plink_conv = "<CHROMOSOME>_10k_sample.bed"
    random_seed = 5
    variable_effects = None
    trait_subset = 0
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
        elif o == "--low_ld_threshold":
            low_ld_threshold = float(a)
        elif o == "--high_ld_threshold":
            high_ld_threshold = float(a)
        elif o == "--trait_type":
            trait_type = a.lower()
        elif o == "--case_prob":
            case_prob = float(a)
            if not (0 < case_prob < 1):
                print("ERROR: --case_prob must be between 0 and 1")
                sys.exit(1)
        elif o == "--effect_size":
            effect_size = float(a)
            if not (effect_size > 0):
                print("ERROR: Specified effect size must be >0.")
                sys.exit(1)
        elif o == "--trait_subset":
            trait_subset = float(a)
            if not (0 < trait_subset <= 1):
                print("ERROR: --trait_subset must be between 0 and 1")
                sys.exit(1)
        elif o == "--variable_effects":
            variable_effects = a.lower() in ['true', '1', 'yes']
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
    
    n_traits = check_int(n_traits, "number of traits")
    n_pairs = check_int(n_pairs, "number of variant pairs")
    effect_size = float(effect_size)
    random_seed = check_int(random_seed, "random seed")
    n_workers = check_int(n_workers, "number of threads")
    low_ld_threshold = float(low_ld_threshold)
    high_ld_threshold = float(high_ld_threshold)
    
    # Validate LD thresholds
    if low_ld_threshold < 0 or low_ld_threshold > 1:
        print("ERROR: --low_ld_threshold must be between 0 and 1")
        sys.exit(1)
    if high_ld_threshold < 0 or high_ld_threshold > 1:
        print("ERROR: --high_ld_threshold must be between 0 and 1")
        sys.exit(1)
    if low_ld_threshold >= high_ld_threshold:
        print("ERROR: --low_ld_threshold must be less than --high_ld_threshold")
        sys.exit(1)
        
    #Validate trait subset
    if trait_subset < 0 or trait_subset > 1:
        print("ERROR: --trait_subset must be between 0 and 1")
        sys.exit(1)
        
    #Validate number of workers
    if n_workers > os.cpu_count():
        print(f"ERROR: Selected number of threads exceeds cores per node. Threads must be less than {os.cpu_count()}.")
        sys.exit(2)
    
    # Validate trait type
    if trait_type not in ['continuous', 'binary']:
        print("ERROR: Invalid trait type given. Select continuous/binary.")
        sys.exit(2)
    
    #Check naming conventions for needed parameters
    if '<POP>' not in keep_conv and keep_file:
        print("ERROR: Invalid keep file naming convention given. Convention must include <POP>.")
        sys.exit(2)
    
    if '<CHROMOSOME>' not in plink_conv and keep_file:
        print("ERROR: Invalid plink file naming convention given. Convention must include <CHROMOSOME>.")
        sys.exit(2)
    
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
    out_file_name = f"{chromo}.{pop_name}.{n_traits}_traits.{n_pairs}_pairs.{trait_type}{'_' + str(case_prob) if trait_type == 'binary' else ''}.subset_{trait_subset}.{'variable' if variable_effects else 'fixed'}_{effect_size}.csv"
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
    if trait_type == "binary":
        random_seed += 1
    
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
    ################# STEP 2: SELECT CAUSAL VARIANT ###############################
    ###############################################################################
    
    print("NOTE: Randomly selecting causal variant...")
    causal_idx = np.random.randint(0, len(variant_ids))
    causal_pos = positions[causal_idx]
    causal_var = variant_ids[causal_idx]
    print(f"NOTE: Selected causal variant {causal_var} at position {causal_pos}")
    
    ###############################################################################
    ######### STEP 3: CALCULATE LD WITH CAUSAL VARIANT FOR ALL VARIANTS ##########
    ###############################################################################
    
    print("NOTE: Calculating genotype correlations between all variants and causal variant...")
    
    def calculate_ld_with_causal_chunk(variant_chunk, causal_idx):
        """Calculate LD between variants in chunk and the causal variant"""
        chunk_results = []
        causal_geno = geno_data[:, causal_idx]
        
        for var_idx in variant_chunk:
            var_id = variant_ids[var_idx]
            
            # Skip the causal variant itself
            if var_idx == causal_idx:
                chunk_results.append((var_id, np.nan))
                continue
            
            var_geno = geno_data[:, var_idx]
            
            # Calculate LD with causal variant
            valid = ~(np.isnan(var_geno) | np.isnan(causal_geno))
            
            if np.sum(valid) < 10:
                ld_causal = np.nan
            else:
                ld_causal = np.corrcoef(var_geno[valid], causal_geno[valid])[0, 1]
            
            chunk_results.append((var_id, ld_causal))
        
        return chunk_results
    
    # Process in parallel
    n_vars = len(variant_ids)
    chunk_size = 500
    variant_chunks = [list(range(start, min(start + chunk_size, n_vars))) 
                     for start in range(0, n_vars, chunk_size)]
    
    all_ld_results = []
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(calculate_ld_with_causal_chunk, chunk, causal_idx) 
                  for chunk in variant_chunks]
        for future in futures:
            all_ld_results.extend(future.result())
    
    # Create dictionary of variant -> LD with causal
    ld_with_causal = {var_id: ld for var_id, ld in all_ld_results}
    
    ###############################################################################
    ######### STEP 4: FILTER VARIANTS INTO "NOT IN LD" AND "IN LD" SETS ##########
    ###############################################################################
    
    print(f"NOTE: Filtering variants by LD thresholds (low <= {low_ld_threshold}, high >= {high_ld_threshold})...")
    
    not_in_ld_vars = []
    in_ld_vars = []
    
    for var_id, ld in ld_with_causal.items():
        if var_id == causal_var:  # Skip causal variant itself
            continue
        
        if np.isnan(ld):  # Skip variants with missing LD
            continue
        
        abs_ld = abs(ld)
        
        if abs_ld <= low_ld_threshold:
            not_in_ld_vars.append(var_id)
        elif abs_ld >= high_ld_threshold:
            in_ld_vars.append(var_id)
    
    print(f"NOTE: Found {len(not_in_ld_vars)} variants NOT in LD (|r| <= {low_ld_threshold})")
    print(f"NOTE: Found {len(in_ld_vars)} variants IN LD (|r| >= {high_ld_threshold})")
    
    ###############################################################################
    ######### STEP 5: CREATE ALL POSSIBLE PAIRS BETWEEN THE TWO SETS #############
    ###############################################################################
    
    print("NOTE: Creating all possible pairs between NOT in LD and IN LD variants...")
    
    all_pairs = []
    for not_ld_var in not_in_ld_vars:
        for in_ld_var in in_ld_vars:
            all_pairs.append((not_ld_var, in_ld_var))
    
    print(f"NOTE: Created {len(all_pairs)} total possible pairs")
    
    # Ensure sufficient pairs remain
    if len(all_pairs) < n_pairs:
        msg = (f"ERROR: Only {len(all_pairs)} valid pairs available, but {n_pairs} requested. "
               f"Try adjusting --low_ld_threshold or --high_ld_threshold to get more variants in each set.")
        print(msg)
        sys.exit(1)
    
    ###############################################################################
    ######### STEP 6: SAMPLE DOWN TO REQUESTED NUMBER OF PAIRS ####################
    ###############################################################################
    
    print(f"NOTE: Sampling {n_pairs} pairs from {len(all_pairs)} available pairs...")
    if len(all_pairs) > n_pairs:
        all_pairs = random.sample(all_pairs, n_pairs)
    else:
        print(f"WARNING: Using all {len(all_pairs)} available pairs (equal to requested {n_pairs}).")
    
    ###########################################################################################################
    ################################# CPU-Based Variant Correlation Calculations ############################
    ###########################################################################################################
    
    # Multithreaded genotype correlation calculation
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
    
    # Calculate correlations with causal variant for pairs
    print("NOTE: Extracting LD values with causal variant for paired variants...")
    
    def extract_causal_ld_for_pairs(pair_chunk):
        """Extract already-calculated LD between each variant in pairs and the causal variant"""
        chunk_results = []
        
        for v1, v2 in pair_chunk:
            ld_causal_1 = ld_with_causal.get(v1, np.nan)
            ld_causal_2 = ld_with_causal.get(v2, np.nan)
            
            chunk_results.append((v1, v2, ld_causal_1, ld_causal_2))
        
        return chunk_results
    
    # Process correlations in parallel
    print(f"NOTE: Using {n_workers} threads for correlation calculations")
    pair_chunks = np.array_split(all_pairs, n_workers*4)
    geno_corrs = []
    causal_lds = []
    
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        # Genotype pair correlations
        futures = [executor.submit(process_correlation_chunk, chunk) for chunk in pair_chunks]
        for future in futures:
            geno_corrs.extend(future.result())
        
        # Correlations with causal variant
        futures = [executor.submit(extract_causal_ld_for_pairs, chunk) for chunk in pair_chunks]
        for future in futures:
            causal_lds.extend(future.result())
    
    ###########################################################################################################
    ############################################### Calculate Trait Values ####################################
    ###########################################################################################################
    
    print("NOTE: Simulating traits with causal variant effects...")
    
        ################# TRAIT SIMULATIONS #####################
    causal_geno = geno_data[:, causal_idx].reshape(-1, 1)

    # Orient genotypes to minor allele
    alt_pct = np.sum(causal_geno == 2, axis=0) / len(causal_geno)
    ref_pct = np.sum(causal_geno == 0, axis=0) / len(causal_geno)

    if alt_pct[0] > ref_pct[0]:
        causal_geno[:, 0] = 2 - causal_geno[:, 0]

    # Determine which traits are affected by the causal variant
    n_affected_traits = int(np.ceil(n_traits * trait_subset))
    affected_trait_indices = np.random.choice(n_traits, size=n_affected_traits, replace=False)

    # Generate effect sizes for each trait
    if variable_effects:
        # Draw from uniform distribution: [0, 2*effect_size] centered at effect_size
        trait_effect_sizes = np.random.uniform(0, 2 * effect_size, size=n_traits)
    else:
        # Use fixed effect size for all traits
        trait_effect_sizes = np.full(n_traits, effect_size)

    # Set effect size to 0 for non-affected traits
    trait_effect_sizes_masked = np.zeros(n_traits)
    trait_effect_sizes_masked[affected_trait_indices] = trait_effect_sizes[affected_trait_indices]

    # Calculate genetic effects for each trait (n_samples x n_traits matrix)
    genetic_effects_matrix = causal_geno[:, 0][:, None] * trait_effect_sizes_masked[None, :]

    if trait_type == 'continuous':
        # Continuous: additive effects + noise
        traits = pd.DataFrame(
            genetic_effects_matrix + np.random.normal(size=(n_samples, n_traits)),
            index=geno.sample.values,
            columns=[f"trait_{i}" for i in range(n_traits)]
        )

    elif trait_type == 'binary':
        # Generate binary traits with trait-specific effects
        traits_list = []

        for trait_idx in range(n_traits):
            trait_genetic_effects = genetic_effects_matrix[:, trait_idx]

            # Find intercept to achieve target case probability for this trait
            def f(intercept):
                probs = 1 / (1 + np.exp(-(intercept + trait_genetic_effects)))
                return np.mean(probs) - case_prob

            # Newton-Raphson optimization
            intercept = 0.0
            for _ in range(100):
                probs = 1 / (1 + np.exp(-(intercept + trait_genetic_effects)))
                grad = np.mean(probs * (1 - probs))
                if abs(grad) < 1e-6:
                    break
                intercept -= f(intercept)/grad
                if abs(f(intercept)) < 1e-4:
                    break

            # Generate binary trait
            log_odds = trait_genetic_effects + intercept
            probs = 1 / (1 + np.exp(-log_odds))
            trait_values = np.random.binomial(1, probs)
            traits_list.append(trait_values)

        # Combine all traits into DataFrame
        traits = pd.DataFrame(
            np.column_stack(traits_list),
            index=geno.sample.values,
            columns=[f"trait_{i}" for i in range(n_traits)]
        )

    ###########################################################################################################
    ######################### GPU-Enabled Calculation of Regression Coefficients ##############################
    ###########################################################################################################
    
    print("NOTE: Running associations on GPUs...")
    
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
    
    #Create function that runs logistic regression quickly
    def calculate_logistic_betas_fast(X, Y, queue, max_iter=100, tol=1e-6):
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
        X_expanded = X[:, :, None]
        X_sq = dpnp.square(X)
        
        # Reusable buffers for intermediate calculations
        logits = dpnp.empty([n_samples, n_variants, n_targets], dtype=X.dtype, sycl_queue=queue)
        p = dpnp.empty_like(logits, sycl_queue=queue)
        W = dpnp.empty_like(logits, sycl_queue=queue)
        residuals = dpnp.empty_like(logits, sycl_queue=queue)
        
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
    
    # Split variants across GPUs
    variant_indices = np.arange(len(variant_ids))
    chunk_count = 20
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
    
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
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
    
    # Add causal variant information
    results_df['causal_var'] = causal_var
    results_df['causal_maf'] = maf_mapping[causal_var]
    
    # Merge correlations between pair variants
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
