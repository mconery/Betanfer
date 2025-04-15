'''
define_sig_loci.py

The purpose of this script is to read GWAS results and define significant loci bounds
using a tile-based approach. Modified to support configurable column names and output base name.

***Run in Conda finemap environment***
'''

import os
import glob
import subprocess
import getopt
import sys
from multiprocessing import Pool, cpu_count
import pandas as pd
from intervaltree import Interval, IntervalTree

# Get Trait-level mapping script location
script_path = os.path.realpath(__file__)
MAP_SCRIPT = os.path.join(os.path.dirname(script_path), 'define_sig_loci_for_trait.py')

#Get number of threads
N_PROCESSES = min(32, cpu_count() // 2) # Adjusting down because we seem not to be able to use all cores without hitting memory issues

###############################################################################
#########################  ARGUMENT PARSING FUNCTION ##########################
###############################################################################

#Help function
def print_help():
    print("""-----------------------------------------------------------------
ARGUMENTS
    -i => Input directory with GWAS summary statistics (REQUIRED)
    -o => Output directory for loci files (REQUIRED)
    -h => Help""")
    sys.exit()

#Argument parser
def parse_arguments(argv):
    input_dir = None
    output_dir = None
    
    try:
        opts, args = getopt.getopt(argv, "i:o:h", ["help"])
    except getopt.GetoptError:
        print("ERROR: Invalid arguments")
        print_help()
    
    for opt, arg in opts:
        if opt == '-i':
            input_dir = arg.rstrip('/')
            if not os.path.isdir(input_dir):
                print(f"ERROR: Input directory {input_dir} not found")
                sys.exit(1)
        elif opt == '-o':
            output_dir = arg.rstrip('/')
            os.makedirs(output_dir, exist_ok=True)
        elif opt in ('-h', '--help'):
            print_help()
    
    if not input_dir or not output_dir:
        print("ERROR: Missing required arguments")
        print_help()
    
    # Extract population from input directory name
    population = os.path.basename(input_dir)
    return input_dir, output_dir, population

###############################################################################
#############################  PROCESSING FUNCTIONS ###########################
###############################################################################

#Trait processing function
def process_trait(file_path, output_dir, population):
    """Process a single trait file"""
    base_name = os.path.basename(file_path)
    trait = base_name.split(".")[0]
    output_file = os.path.join(output_dir, f"{trait}.{population}.sig_loci.csv")
    
    if os.path.exists(output_file):
        print(f"NOTE: Skipping {trait} - output exists")
        return

    cmd = [
        "python", MAP_SCRIPT,
        "-g", file_path,
        "-o", output_dir,
        "-p", "5e-8",
        "-s", "1e-6",
        "-t", "100000",
        "-b", f"{trait}.{population}",
        "--chr-col", "CHR",
        "--pos-col", "POS",
        "--p-col", "P"
    ]

    try:
        subprocess.run(cmd, check=True)
        print(f"SUCCESS: Completed {trait}")
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Error processing {trait}: {str(e)}")
        
#Create a function to read in each output file
def read_output(file_path):
    try:
        loci_list = pd.read_csv(file_path, sep=",")
        trait = os.path.basename(file_path).split('.')[0]
        loci_list['TRAIT'] = [trait for _ in range(len(loci_list))]
    except pd.errors.EmptyDataError:
        loci_list = None
    return loci_list

###############################################################################
##############################  MAIN FUNCTION #################################
###############################################################################

if __name__ == "__main__":
    # Parse command-line arguments
    input_dir, output_dir, population = parse_arguments(sys.argv[1:])
    
    # Get list of input files using population tag
    input_files = glob.glob(os.path.join(input_dir, f"*.{population}.mrmega.sumstats.txt.gz"))
    print(f"NOTE: Found {len(input_files)} trait files for population {population}")
    
    # Get number of threads
    N_PROCESSES = min(8, cpu_count() // 2)
    print(f"NOTE: Using {N_PROCESSES} parallel processes")

    # Process files in parallel through locus finding script
    with Pool(processes=N_PROCESSES) as pool:
        pool.starmap(process_trait, [(f, output_dir, population) for f in input_files])

    # Get output file names with population tag
    output_files = [os.path.join(output_dir, os.path.basename(f).replace(
        f'.{population}.mrmega.sumstats.txt.gz', 
        f'.{population}.sig_loci.csv'
    )) for f in input_files]

    # Read and combine output files
    with Pool(processes=N_PROCESSES) as pool:
        temp = pool.map(read_output, output_files)
    
    # Process combined results
    temp = [x for x in temp if isinstance(x, pd.DataFrame)]
    if temp:
        combined_df = pd.concat(temp, ignore_index=True)
        
        # Calculate overlap counts
        combined_df['OVERLAP_COUNT'] = 0
        for chr_name, chr_group in combined_df.groupby('CHR'):
            tree = IntervalTree()
            for idx, row in chr_group.iterrows():
                tree.add(Interval(row['START'], row['END'], idx))
            for idx, row in chr_group.iterrows():
                overlaps = tree.overlap(row['START'], row['END'])
                combined_df.loc[idx, 'OVERLAP_COUNT'] = len(overlaps) - 1
        
        # Write master file with population tag
        master_file = os.path.join(output_dir, f'sig_loci.{population}.master.csv')
        combined_df.to_csv(master_file, index=False)
        print(f"SUCCESS: Created master file {master_file}")
    else:
        print("NOTE: No significant loci found across all traits")
