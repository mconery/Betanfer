'''
define_sig_loci_for_bmd.py

The purpose of this script is to read GWAS results and define significant loci bounds
using a tile-based approach. Modified to support configurable column names and output base name.

***Run in Conda finemap environment***
'''

import os
import glob
import subprocess
from multiprocessing import Pool, cpu_count
import pandas as pd
from intervaltree import Interval, IntervalTree

# Set key directories and files
INPUT_DIR = "/grand/projects/GeomicVar/proj-mvp-ukb-meta/output_gia/metaanalysis_for_export/ALL"
OUTPUT_DIR = "/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci"
MAP_SCRIPT = "/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/1_Locus_Defining/define_sig_loci_for_trait.py"

#Get number of threads
N_PROCESSES = min(8, cpu_count() // 2) # Adjusting down because we seem not to be able to use all cores without hitting memory issues

#Create function to find loci for each trait
def process_trait(file_path):
    """Process a single trait file"""
    base_name = os.path.basename(file_path)
    trait = base_name.split(".")[0]
    output_file = os.path.join(OUTPUT_DIR, f"{trait}.sig_loci.csv")
    
    if os.path.exists(output_file):
        print(f"NOTE: Skipping {trait} - output exists")
        return
    
    cmd = [
        "python", MAP_SCRIPT,
        "-g", file_path,
        "-o", OUTPUT_DIR,
        "-p", "5e-8",
        "-s", "1e-6",
        "-t", "250000",
        "-b", trait,
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
        loci_list['TRAIT'] = [trait for x in range(len(loci_list))]
    except pd.errors.EmptyDataError:
        loci_list = None
    return loci_list

if __name__ == "__main__":
    # Get list of input files
    input_files = glob.glob(os.path.join(INPUT_DIR, "*.ALL.mrmega.sumstats.txt.gz"))
    
    print(f"NOTE: Found {len(input_files)} trait files")
    print(f"NOTE: Using {N_PROCESSES} parallel processes")
    
    # Create output directory if needed
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Process files in parallel through locus finding script
    with Pool(processes=N_PROCESSES) as pool:
        temp = pool.map(process_trait, input_files)
    
    #Get output file names
    output_files = [os.path.join(OUTPUT_DIR,os.path.basename(x).replace('.ALL.mrmega.sumstats.txt.gz', '.sig_loci.csv')) for x in input_files]
     
    #Read in all the output files in order to combine them together
    with Pool(processes=N_PROCESSES) as pool:
        temp = pool.map(read_output, output_files)
        
    #Remove empty entries
    temp = [x for x in temp if type(x) == pd.core.frame.DataFrame]
    #Combine lists and write to file
    combined_df = pd.concat(temp, ignore_index=True)
    
    #Add a column to the dataframe indicating the number of other loci overlapped by the locus
    # Initialize new column
    combined_df['OVERLAP_COUNT'] = 0
    
    # Process each chromosome separately
    for chr_name, chr_group in combined_df.groupby('CHR'):
        # Create interval tree for this chromosome
        tree = IntervalTree()
        for idx, row in chr_group.iterrows():
            tree.add(Interval(row['START'], row['END'], idx))
        
        # Calculate overlaps within chromosome
        for idx, row in chr_group.iterrows():
            overlaps = tree.overlap(row['START'], row['END'])
            combined_df.loc[idx, 'OVERLAP_COUNT'] = len(overlaps) - 1  # Exclude self
    
    #Write to file
    combined_df.to_csv(os.path.join(OUTPUT_DIR, 'sig_loci.master.csv'), index=False)

