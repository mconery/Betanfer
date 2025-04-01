'''
define_sig_loci_for_bmd.py

The purpose of this script is to read GWAS results and define significant loci bounds
using a tile-based approach. Modified to support configurable column names and output base name.

***Run on default Respulbica Python/3.9 due to version incompatibilities***
'''

import os
import glob
import subprocess
from multiprocessing import Pool, cpu_count

# Set key directories and files
INPUT_DIR = "/grand/projects/GeomicVar/proj-mvp-ukb-meta/output_gia/metaanalysis_for_export/ALL"
OUTPUT_DIR = "/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci"
MAP_SCRIPT = "/grand/projects/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/Meta_Fine-Mapping/1_Locus_Defining/define_sig_loci_for_trait.py"

#Get number of threads
N_PROCESSES = min(8, cpu_count() // 2) # Adjusting down because we seem not to be able to use all cores without hitting memory issues

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

if __name__ == "__main__":
    # Get list of input files
    input_files = glob.glob(os.path.join(INPUT_DIR, "*.ALL.mrmega.sumstats.txt.gz"))
    
    print(f"NOTE: Found {len(input_files)} trait files")
    print(f"NOTE: Using {N_PROCESSES} parallel processes")
    
    # Create output directory if needed
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
# Process files in parallel
with Pool(processes=N_PROCESSES) as pool:
    temp = pool.map(process_trait, input_files)
