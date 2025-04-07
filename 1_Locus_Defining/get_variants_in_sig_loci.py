'''
get_variants_in_sig_loci.py

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
-n => GWAS file naming convention (default: <TRAIT>.ALL.mrmega.sumstats.txt.gz) OPTIONAL
--chr-col => Chromosome column name (default: CHR) OPTIONAL
--pos-col => Position column name (default: POS) OPTIONAL
-h => help""")
    sys.exit(exit_num)

###############################################################################
###########################  PROCESS INDIVIDUAL TRAIT  ########################
###############################################################################

#Function to process individual trait and identify all variants in sig loci
def process_trait(trait_group, gwas_folder, naming_conv, chr_col, pos_col):
    """Collect all variants in significant loci for a trait"""
    trait, group = trait_group
    try:
        gwas_file = os.path.join(gwas_folder, f"{naming_conv}".replace("<TRAIT>", trait))
        df = pd.read_csv(gwas_file, sep='\t')
        df = df.rename(columns={chr_col: "CHR", pos_col: "POS"})
        df['VARIANT_ID'] = df['CHR'].astype(str) + ':' + df['POS'].astype(str)
        
        trait_variants = set()
        
        for _, row in group.iterrows():
            chrom = str(row['CHR'])
            start = row['START']
            end = row['END']
            
            region_df = df[
                (df["CHR"].astype(str) == chrom.lstrip("chr")) &
                (df["POS"].between(start, end))
            ]
            
            trait_variants.update(region_df['VARIANT_ID'].tolist())
            
        print(f"Processed {trait} with {len(trait_variants)} variants")
        return {trait: trait_variants}
        
    except Exception as e:
        print(f"ERROR processing {trait}: {str(e)}")
        return {}

###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND MULTI-PROCESS LOCI #############
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "l:o:g:n:nhb:", 
                                 ["help", "chr-col=", "pos-col="])
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
        
    # Validation checks
    if not all(os.path.exists(f) for f in [gwas_dir, os.path.dirname(output_file), loci_file]):
        print("ERROR: Input directories/files not found")
        sys.exit(1)
    
    #Check name convention for <TRAIT>
    if '<TRAIT>' not in name_conv:
        print("ERROR: Invalid naming convention given. Convention must include <TRAIT>.")
    
    # Modified driver section
    loci_df = pd.read_csv(loci_file)
    trait_groups = loci_df.groupby('TRAIT')
    
    # Create worker partial
    worker = partial(process_trait,
                    gwas_folder=gwas_dir,
                    naming_conv=name_conv,
                    chr_col=chr_col,
                    pos_col=pos_col)
    
    # Process loci in parallel
    with Pool(processes=N_PROCESSES) as pool:
        results = pool.map(worker, trait_groups)
    
    # Combine results
    variant_sets = {}
    all_variants = set()
    
    for result in results:
        if result:
            trait = list(result.keys())[0]
            variants = result[trait]
            variant_sets[trait] = variants
            all_variants.update(variants)
            
    # Create matrix
    matrix = pd.DataFrame(
        index=pd.Index(sorted(all_variants), name='VARIANT_ID'),
        columns=pd.Index(sorted(variant_sets.keys()), name='TRAIT'),
        data=0
    )
    
    for trait, variants in variant_sets.items():
        matrix.loc[matrix.index.isin(variants), trait] = 1
        
    matrix.to_csv(output_file, index=True)
    print("SUCCESS: Successfully processed all traits")

if __name__ == "__main__":
    main(sys.argv[1:])
