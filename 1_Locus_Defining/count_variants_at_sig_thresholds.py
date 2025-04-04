'''
count_variants_at_sig_thresholds.py

The purpose of this script is to count how many significant variants would 
remain in long loci after pruning out variants below different significance 
thresholds. This script can also be used to count the number of variants in
loci in general by setting the significance threshold to 1 and the minimum 
length to 0. 

***Run in Conda finemap environment***
'''


import pandas as pd
import numpy as np
import getopt
import sys
import os
from multiprocessing import Pool, cpu_count
from functools import partial

# Configure number of processes
N_PROCESSES = min(8, cpu_count() // 2)

def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
-l => Input Loci CSV file REQUIRED
-o => Output CSV file REQUIRED
-g => GWAS directory (default: /grand/projects/GeomicVar/proj-mvp-ukb-meta/output_gia/metaanalysis_for_export/ALL) OPTIONAL
-m => Minimum locus length in Mb (default: 5) OPTIONAL
-n => GWAS file naming convention (default: <TRAIT>.ALL.mrmega.sumstats.txt.gz) OPTIONAL
-p => P-value thresholds (comma-separated, default: 5e-8,1e-10,1e-20,1e-30,1e-40,1e-50) OPTIONAL
--chr-col => Chromosome column name (default: CHR) OPTIONAL
--pos-col => Position column name (default: POS) OPTIONAL
--p-col => P-value column name (default: P) OPTIONAL
-h => help""")
    sys.exit(exit_num)

###############################################################################
###########################  PROCESS INDIVIDUAL LOCUS  ########################
###############################################################################

def process_trait(trait_group, gwas_folder, naming_conv, thresholds, chr_col, pos_col, p_col):
    """Process all loci for a single trait"""
    trait, group = trait_group
    try:
        # Load GWAS data once per trait
        gwas_file = os.path.join(gwas_folder, f"{naming_conv}{naming_conv}".replace("<TRAIT>", trait))
        df = pd.read_csv(gwas_file, sep='\t')
        df = df.rename(columns={chr_col: "CHR", pos_col: "POS", p_col: "P"})
        
        results = []
        for _, row in group.iterrows():
            geno_region = f"{row['CHR']}:{row['START']}-{row['END']}"
            chrom, region = geno_region.split(":")
            start, end = map(int, region.split("-"))
            
            # Filter using existing loaded DataFrame
            region_df = df[
                (df["CHR"].astype(str) == chrom.lstrip("chr")) &
                (df["POS"].between(start, end))
            ].copy()
            
            # Count variants at thresholds
            counts = {f"PV_{thresh:.0e}": (region_df["P"] < thresh).sum() for thresh in thresholds}
            results.append({**row.to_dict(), **counts})
            
            #Print Update
            print(f"NOTE: Completed processing {trait}.{geno_region}.")
        
        print(f"Processed {len(group)} loci for trait: {trait}")
        return results
        
    except Exception as e:
        print(f"ERROR processing {trait}: {str(e)}")
        return []

###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND MULTI-PROCESS LOCI #############
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "l:o:g:m:n:p:nhb:", 
                                 ["help", "chr-col=", "pos-col=", "p-col="])
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)
    
    if '-h' in options_dict or '--help' in options_dict:
        help()

    # Required arguments
    try:
        loci_file = options_dict['-l']
        output_file = options_dict['-o']
    except KeyError:
        print("ERROR: Missing required arguments (-g or -o)")
        help()

    # Configurable parameters with defaults
    gwas_dir = options_dict.get('-g', '/grand/projects/GeomicVar/proj-mvp-ukb-meta/output_gia/metaanalysis_for_export/ALL')
    name_conv = options_dict.get('-n', '<TRAIT>.ALL.mrmega.sumstats.txt.gz')
    chr_col = options_dict.get('--chr-col', 'CHR')
    pos_col = options_dict.get('--pos-col', 'POS')
    p_col = options_dict.get('--p-col', 'P')
    
    # P-value threshold handling
    try:
        pvals = list(map(float, options_dict.get('-p', '').split(','))) if '-p' in options_dict else None
        if pvals is None or not pvals:
            PVAL_THRESHOLDS = [5e-8, 1e-10, 1e-20, 1e-30, 1e-40, 1e-50]
        else:
            PVAL_THRESHOLDS = sorted(pvals, reverse=True)
    except ValueError:
        print("ERROR: Invalid p-value thresholds provided. Must be comma-separated numbers")
        sys.exit(1)
    
    # Numeric parameter
    try:
       min_len = int(str(options_dict.get('-m', 5)).lstrip('Mb').lstrip('mb').lstrip('MB'))
    except ValueError:
        print("ERROR: Invalid minimum length given")
        sys.exit(1)
    
    # Validation checks
    if not all(os.path.exists(f) for f in [gwas_dir, os.path.dirname(output_file), loci_file]):
        print("ERROR: Input directories/files not found")
        sys.exit(1)
    
    #Check name convention for <TRAIT>
    if '<TRAIT>' not in name_conv:
        print("ERROR: Invalid naming convention given. Convention must include <TRAIT>.")
    
    # Modified driver section
    loci_df = pd.read_csv(loci_file)
    long_df = loci_df[loci_df["LENGTH"].astype(int) >= min_len * 1e6].copy()
    trait_groups = long_df.groupby('TRAIT')
    
    
    # Create worker partial
    worker = partial(process_trait,
                    gwas_folder=gwas_dir,
                    naming_conv=name_conv,
                    thresholds=PVAL_THRESHOLDS,
                    chr_col=chr_col,
                    pos_col=pos_col,
                    p_col=p_col)
    
    # Process loci in parallel
    with Pool(processes=N_PROCESSES) as pool:
        results = pool.map(worker, trait_groups)
    
    # Combine results and save
    result_df = pd.DataFrame([item for sublist in results for item in sublist])
    result_df.to_csv(output_file, index=False)
    print(f"SUCCESS: Successfully processed {len(result_df)} loci")

if __name__ == "__main__":
    main(sys.argv[1:])
