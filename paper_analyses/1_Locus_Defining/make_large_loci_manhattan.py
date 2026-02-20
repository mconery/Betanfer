'''
make_large_loci_manhattan.py

The purpose of this script is to make Manhattan plots for a given list of input
loci. The input file should be in the same format as all of my loci files. 

***Run in Conda finemap environment***
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import getopt
import sys
import os
from multiprocessing import Pool, cpu_count
from functools import partial

#Get number of threads
N_PROCESSES = min(8, cpu_count() // 2) # Adjusting down because we seem not to be able to use all cores without hitting memory issues

def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
-l => File of Input Loci REQUIRED
-o => Output Directory (default 5e-8) REQUIRED
-g => GWAS file directory (default is /grand/projects/GeomicVar/proj-mvp-ukb-meta/output_gia/metaanalysis_for_export/ALL)
-m => Minimum locus range for mapping (in Mb; default is 5Mb) OPTIONAL 
-n => GWAS Sumstats Naming Convention (default is <TRAIT>.ALL.mrmega.sumstats.txt.gz)
--chr-col => chromosome column name (default 'CHR') OPTIONAL  
--pos-col => position column name (default 'POS') OPTIONAL
--p-col => p-value column name (default 'P') OPTIONAL
-h => help""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "l:o:g:m:n:nhb:", 
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
        output_dir = options_dict['-o']
    except KeyError:
        print("ERROR: Missing required arguments (-g or -o)")
        help()

    # Configurable parameters with defaults
    gwas_dir = options_dict.get('-g', '/grand/projects/GeomicVar/proj-mvp-ukb-meta/output_gia/metaanalysis_for_export/ALL')
    name_conv = options_dict.get('-n', '<TRAIT>.ALL.mrmega.sumstats.txt.gz')
    chr_col = options_dict.get('--chr-col', 'CHR')
    pos_col = options_dict.get('--pos-col', 'POS')
    p_col = options_dict.get('--p-col', 'P')
    
    # Numeric parameter
    try:
       min_len = int(str(options_dict.get('-m', 5)).lstrip('Mb').lstrip('mb').lstrip('MB'))
    except ValueError:
        print("ERROR: Invalid minimum length given")
        sys.exit(1)
    
    # Validation checks
    if not all(os.path.exists(f) for f in [gwas_dir, output_dir, loci_file]):
        print("ERROR: Input directories/files not found")
        sys.exit(1)
    
    #Check name convention for <TRAIT>
    if '<TRAIT>' not in name_conv:
        print("ERROR: Invalid naming convention given. Convention must include <TRAIT>.")

    driver(loci_file, output_dir, gwas_dir, min_len, name_conv, chr_col, pos_col, p_col)

###############################################################################
#############################  PLOT MAKING FUNCTIONS ##########################
###############################################################################

def generate_manhattan(summary_stats_path, genomic_region, chr_col="CHR", pos_col="POS", p_col="P",
                      sig_threshold=5e-8, point_size=8, figsize=(14, 6), save_path=None):
    """
    Generate a Manhattan plot for a specific genomic region
    
    Parameters:
    summary_stats_path (str): Path to summary statistics file (CSV)
    genomic_region (str): Genomic region in "chr:start-end" format
    chr_col (str): Chromosome column name (default: "CHR")
    pos_col (str): Position column name (default: "POS")
    p_col (str): P-value column name (default: "P")
    sig_threshold (float): Genome-wide significance threshold (default: 5e-8)
    point_size (int): Marker size for points (default: 8)
    figsize (tuple): Figure dimensions (default: (14, 6))
    save_path (str): Optional path to save plot (default: None)
    
    Returns:
    matplotlib.figure.Figure: Generated figure object
    """
    # Read and prepare data
    df = pd.read_csv(summary_stats_path, sep="\t")
    df = df.rename(columns={chr_col: "CHR", pos_col: "POS", p_col: "P"})
    
    # Parse genomic region
    chrom, region = genomic_region.split(":")
    start, end = map(int, region.split("-"))
    
    # Filter data for specified region
    region_df = df[
        (df["CHR"].astype(str) == chrom.lstrip("chr")) & 
        (df["POS"].between(start, end))
    ].copy()
    
    # Calculate transformed p-values
    region_df["P"] = np.maximum(5e-324, region_df["P"])
    region_df["-log10P"] = -np.log10(region_df["P"])
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Color chromosomes alternately
    chroms = region_df["CHR"].unique()
    colors_cycle = ["#4E79A7", "#F28E2B"]  # Alternating colors
    
    for i, chrom in enumerate(chroms):
        chrom_df = region_df[region_df["CHR"] == chrom]
        ax.scatter(
            x=chrom_df["POS"],
            y=chrom_df["-log10P"],
            color=colors_cycle[i % len(colors_cycle)],
            s=point_size,
            alpha=0.7
        )
    
    # Add significance threshold line
    ax.axhline(-np.log10(sig_threshold), color="red", linestyle="--")
    
    # Format plot
    ax.set_xlabel(f"Chromosome {chrom} Position", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title(f"Manhattan Plot for {genomic_region}", fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    # Save if requested
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        
    return fig

#This function processes a single locus
def process_locus(row, gwas_folder, naming_conv, out_dir, chr_col, pos_col, p_col):
    """Process a single locus (modified from nested function)"""
    geno_region = f"{row['CHR']}:{row['START']}-{row['END']}"
    trait = str(row['TRAIT'])
    
    # Fix GWAS file path construction
    gwas_file = os.path.join(gwas_folder, f"{naming_conv}").replace('<TRAIT>', trait)
    out_file = os.path.join(out_dir, f"{trait}.{geno_region.replace(':', '_')}.jpg")

    if os.path.exists(out_file):
        print(f"WARNING: Skipping already finished {trait}.{geno_region}")
        return

    generate_manhattan(gwas_file, geno_region, 
                       chr_col=chr_col, pos_col=pos_col, p_col=p_col, 
                       save_path=out_file)
    print(f"NOTE: Finished {trait}.{geno_region}")

###############################################################################
#################################  DRIVER  ####################################
###############################################################################

def driver(loci_file, output_dir, gwas_dir, min_len, name_conv, chr_col, pos_col, p_col):
    
    #Add trailing slashes to directories
    output_dir = output_dir.rstrip('/') + '/'
    gwas_dir = gwas_dir.rstrip('/') + '/'
    
    #Read in loci_file
    loci_df = pd.read_csv(loci_file)
    #Get long loci
    long_df= loci_df[(loci_df["LENGTH"].astype(int) >= min_len * 1e6)].copy()
    
    #Verify that at least one locus remains
    if len(long_df) < 1:
        print("WARNING: No loci of minimum length detected. Exiting.")
        sys.exit(0)
    
   # Create partial function with fixed parameters
    worker = partial(process_locus,
                     gwas_folder=gwas_dir,
                     naming_conv=name_conv,
                     out_dir=output_dir,
                     chr_col=chr_col,
                     pos_col=pos_col,
                     p_col=p_col)

    # Use starmap with iterrows
    with Pool(processes=N_PROCESSES) as pool:
        pool.map(worker, [row for _, row in long_df.iterrows()])

###############################################################################
#############################  MAIN CALL ######################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv[1:])
