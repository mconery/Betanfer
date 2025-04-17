'''
OnekG_Simulation_Wrapper.py

The purpose of this script is to act as a wrapper for 1000G_Simulation_Test.py
and manage the multi-node distribution of jobs using mpi4py

***Run in Conda finemap environment***
'''

#For timing
import time
tic = time.perf_counter()
#Remainder of package imports
import sys
import getopt
import os
from mpi4py.futures import MPIPoolExecutor
from functools import partial

#Define help function
def help(exit_num=0):
    print("""Usage: python OnekG_Simulation_Wrapper.py [OPTIONS]
Options:
  --command_file        File of simulation commands to run (required)
  --plink_dir           Directory of 1000G Plink Files (required)
  --keep_dir            Directory of Population Keep Files for Plink Filtering (required)
  --out_dir             Output Directory (required)
  --chromosomes         Chromosomes to test (comma-delimted; optional)
  --populations         Populations to test (comma-delimted; optional)
  --trait_types         Trait types to test (comma-delimted; optional)
  --n_traits            Number of traits to test (comma-delimted; optional)
  --n_pairs             Number of SNP pairs to test (comma-delimited; optional)
  --n_causals           Number of causal SNPs to test (comma-delimited; optional)
  --effect_sizes        Effect Sizes to test (comma-delimited; optional)
  --case_probs_binary   Case proportions to test (comma-delimited; optional)
  -h, --help            Show this help message""")
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
            ["help", "plink_dir=", "keep_dir=", "out_dir=",
             "chromosomes=", "populations=", "trait_types=", "n_traits=", 
             "n_pairs=", "n_causals=", "effect_sizes=", "case_probs_binary="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
        
    #Set parameter default values
    command_file = None
    plink_dir = None
    keep_dir = None
    out_dir = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
            sys.exit()
        elif o == "--plink_dir":
            plink_dir = a
        elif o == "--keep_dir":
            keep_dir = a
        elif o == "--out_dir":
            out_dir = a
        else:
            print(o)
            assert False, "WARNING: Unhandled input option"
    
    # Validate required arguments
    if None in [plink_dir, keep_dir, out_dir]:
        print("ERROR: Missing required command file!")
        sys.exit(2)

    # Path Validation checks
    if not all(os.path.exists(f) for f in [command_file, plink_dir, keep_dir, out_dir]):
        print("ERROR: Input command file does not exist.")
        sys.exit(1)
        
    # Process list functions
    def process_float_list(flag, float_type, default_values, options_dict = dict(opts)):
        try:
            pvals = list(map(float, options_dict.get(flag, '').split(','))) if flag in options_dict else None
            if pvals is None or not pvals:
                PVAL_THRESHOLDS = default_values
            else:
                PVAL_THRESHOLDS = sorted(pvals, reverse=True)
        except ValueError:
            print(f"ERROR: Invalid {float_type} provided. Must be comma-separated numbers")
            sys.exit(1)
        return PVAL_THRESHOLDS
    def process_int_list(flag, int_type, default_values, options_dict = dict(opts)):
        try:
            pvals = list(map(int, options_dict.get(flag, '').split(','))) if flag in options_dict else None
            if pvals is None or not pvals:
                PVAL_THRESHOLDS = default_values
            else:
                PVAL_THRESHOLDS = sorted(pvals, reverse=True)
        except ValueError:
            print(f"ERROR: Invalid {int_type} provided. Must be comma-separated numbers")
            sys.exit(1)
        return PVAL_THRESHOLDS
    def process_string_list(flag, string_type, default_values, options_dict = dict(opts)):
        try:
            pvals = list(map(str, options_dict.get(flag, '').split(','))) if flag in options_dict else None
            if pvals is None or not pvals:
                PVAL_THRESHOLDS = default_values
            else:
                PVAL_THRESHOLDS = sorted(pvals, reverse=True)
        except ValueError:
            print(f"ERROR: Invalid {string_type} provided. Must be comma-separated numbers")
            sys.exit(1)
        return PVAL_THRESHOLDS
    #Call processing functions
    chromosomes = process_int_list("--chromosomes", "chromosomes", [2, 3, 4])
    populations = process_string_list("--populations", "populations", ["EUR", "AFR", "ALL"])
    trait_types = process_string_list("--trait_types" , "trait types", ["continuous", "binary"])
    n_traits = process_int_list("--n_traits", "number of traits", [50, 100, 200, 300, 400, 500, 1000])
    n_pairs = process_int_list("--n_traits", "number of pairs", [10000])
    n_causals = process_int_list("--n_causals", "number of causal SNPs", [0, 1, 2, 3, 4, 5])
    effect_sizes = process_float_list("--effect_sizes", "effect sizes", [0.1, 0.5, 1, 2, 5])
    case_probs_binary = process_float_list("--case_probs_binary", "binary trait case proportions", [0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    
    # Get location of simulation script
    script_path = os.path.realpath(__file__)
    onekg_script = os.path.join(os.path.dirname(script_path), 'OnekG_Simulation_Test.py')
    
    #Create list of testing parameters
    test_parameters = []
    #Loop over inputs and create testing inputs
    for n_pair in n_pairs:
        for chromosome in chromosomes:
            for pop in populations:
                for trait_type in trait_types:
                    for n_trait in n_traits:
                        for n_causal in n_causals:
                            temp = (
                                f"python {onekg_script} "
                                f"--plink {plink_dir}/chr{chromosome}_10k_sample.bed "
                                f"--pairs={n_pair} --traits={n_trait} --outdir={out_dir} "
                                f"--keep={keep_dir}/1000G.{pop}.txt "
                                f"--trait_type={trait_type} --n_causal={n_causal} --min_distance 30000000"
                            )
                            if n_causal != 0:
                                for effect_size in effect_sizes:
                                    temp_e = f"{temp} --effect_size={effect_size}"
                                    if trait_type == "binary":
                                        for case_prob in case_probs_binary:
                                            temp_eb = f"{temp_e} --case_prob={case_prob}"
                                            test_parameters.append(temp_eb)
                                    else:
                                        test_parameters.append(temp_e)
                            else:
                                if trait_type == "binary":
                                    for case_prob in case_probs_binary:
                                        temp_b = f"{temp} --case_prob={case_prob}"
                                        test_parameters.append(temp_b)
                                else:
                                    test_parameters(temp)
    
    #Create running function
    def run_command(command):
        os.system(command)
    # Create worker partial
    worker = partial(run_command)
    #Now loop over the test_parameters and submit the jobs
    with MPIPoolExecutor() as executor:
        results = executor.map(worker, test_parameters)
    del results
    
    toc = time.perf_counter()
    print(f"SUCCESS: Simulations complete in {(toc - tic)/60:0.4f} minutes")

if __name__ == "__main__":
    main()
