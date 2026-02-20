'''
define_sig_loci_for_trait.py

The purpose of this script is to read GWAS results and define significant loci bounds
using a tile-based approach. Modified to support configurable column names and output base name.

***Run in finemap Conda Environment***
'''

import sys
import getopt
import os
import pandas as pd
import numpy as np
from os.path import exists
import csv

def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
-g => GWAS summary statistics location REQUIRED
-o => output file directory REQUIRED
-p => primary p-value threshold (default 5e-8) OPTIONAL
-s => secondary p-value threshold (default 1e-6) OPTIONAL
-t => tile size (default 250kb; must be multiples of 1kb) OPTIONAL
-c => chromosome ranges file (default grch38_chr_bounds_mhc_only.bed) OPTIONAL
-b => base name for output files (default 'bmd') OPTIONAL
--chr-col => chromosome column name (default 'CHR') OPTIONAL  
--pos-col => position column name (default 'BP') OPTIONAL
--p-col => p-value column name (default 'P') OPTIONAL
--id-col => number of column for unique SNP id (zero-indexed; default 0) OPTIONAL
-h => help""")
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "g:o:p:s:t:c:nhb:", 
                                 ["help", "chr-col=", "pos-col=", "p-col=", "id-col="])
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)
    
    if '-h' in options_dict or '--help' in options_dict:
        help()

    # Required arguments
    try:
        gwas_file = options_dict['-g']
        output_loc = options_dict['-o']
    except KeyError:
        print("ERROR: Missing required arguments (-g or -o)")
        help()

    # Configurable parameters with defaults
    base_name = options_dict.get('-b', 'bmd')
    chr_col = options_dict.get('--chr-col', 'CHR')
    pos_col = options_dict.get('--pos-col', 'BP')
    p_col = options_dict.get('--p-col', 'P')

    # Numeric parameters
    try:
        prime_p_thresh = float(options_dict.get('-p', 5e-8))
        sec_p_thresh = float(options_dict.get('-s', 1e-6))
        stagger_bp = int(options_dict.get('-t', 250000))
        id_col = int(options_dict.get('--id-col', 0))
    except ValueError:
        print("ERROR: Invalid numeric parameter format")
        sys.exit(1)

    # Chromosome ranges file
    script_path = os.path.realpath(__file__)
    default_chr_file = os.path.join(os.path.dirname(script_path), 
                                  'grch38_chr_bounds_mhc_only.bed')
    chromo_lengths = options_dict.get('-c', default_chr_file)

    # Validation checks
    if not all(exists(f) for f in [gwas_file, chromo_lengths]):
        print("ERROR: Input file(s) not found")
        sys.exit(1)

    driver(gwas_file, output_loc, prime_p_thresh, sec_p_thresh,
           stagger_bp, chromo_lengths, base_name, chr_col, pos_col, p_col, id_col)

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

def driver(gwas_file, output_loc, prime_p_thresh, sec_p_thresh, stagger_bp,
          chromo_lengths, base_name, chr_col, pos_col, p_col, id_col):
    
    #Add slash to output dir and make it if necessary
    output_loc = output_loc.rstrip('/') + '/'
    os.makedirs(output_loc, exist_ok=True)
    # Final output with configurable base name
    output_loci_file = f"{output_loc}{base_name}.sig_loci.csv"

    # Read GWAS data with column flexibility
    cur_file = pd.read_csv(gwas_file, index_col=id_col, sep="\t")
    if cur_file.shape[1] <= 1:  # Fallback to space-separated
        cur_file = pd.read_csv(gwas_file, index_col=id_col, sep=" ")
    
    # Filter using configurable column names
    filt_file = cur_file[cur_file[p_col] <= sec_p_thresh]
    
    # Tile processing (unchanged logic)
    temp = filt_file.index.str.split(":")
    try:
        sig_tiles = [int(filt_file[chr_col][x]) + (int(filt_file[pos_col][x].lstrip('chr'))//stagger_bp)/1e9 for x in range(len(filt_file))]
    except AttributeError:
        sig_tiles = [int(filt_file[chr_col][x]) + (int(filt_file[pos_col][x])//stagger_bp)/1e9 for x in range(len(filt_file))]
    sig_tiles_no_dups = np.sort(list(set(sig_tiles)))

    if not sig_tiles_no_dups.size:
        print("NOTE: No significant loci identified for " + base_name)
        # Creates an empty file
        with open(output_loci_file, 'w') as out_path:
            pass
        sys.exit(0)

    #Make a list to hold the output data
    out_ranges = []
    #Loop over the significant tiles and append information to the out_ranges
    cur_locus = []
    for i in range(len(sig_tiles_no_dups)):
        cur_tile = sig_tiles_no_dups[i]
        chromo = int(cur_tile//1)
        start_bp = int(round(1000000000*(cur_tile%1)) * stagger_bp) #This gets the start position of the current tile (ending in 0)
        if cur_locus == []: #Handle the first interval
            cur_locus = [chromo, start_bp+1, start_bp+stagger_bp]
        elif start_bp == cur_locus[2] and chromo == cur_locus[0]: #Extend the current locus
            cur_locus[2] = cur_locus[2] + stagger_bp
        else: #We need to start a new locus in this case
            #Add a length to the cur_locus
            cur_locus.append(cur_locus[2] - cur_locus[1] + 1)
            out_ranges.append(cur_locus) #Append the current locus to the output list
            cur_locus = [chromo, start_bp+1, start_bp+stagger_bp]
    #Handle the final sig tile
    #Add a length to the cur_locus
    cur_locus.append(cur_locus[2] - cur_locus[1] + 1)
    out_ranges.append(cur_locus)
    #Sort the out ranges
    out_ranges = sorted(out_ranges, key = lambda x: (x[0], x[1], x[2]))
    
    #Filter down the filt_file again for the variants significant at the primary threshold
    temp = filt_file[filt_file[p_col] <= prime_p_thresh].index.str.split(":")
    prime_sig_vars = [int(x[0]) + float(x[1])/1e9 for x in temp]
    #Check that there are at least some significant variants at primary threshold
    if len(prime_sig_vars) == 0:
        # Creates an empty file
        with open(output_loci_file, 'w') as out_path:
            pass
        print("NOTE: No significant loci identified for " + base_name)
        sys.exit(0)

    #Read in chromosome lengths file
    chr_lens = [] #This file should already be ordered by chromosome 
    with open(chromo_lengths) as chromo_len_read:
        for line in chromo_len_read.readlines():
            temp = line.strip().split("\t") 
            temp[0] = temp[0][3:] #This will remove the "chr"s in the bed file
            chr_lens.append(temp)
        #Make a dictionary out of the chromosome lengths
        chr_lens_dict = {}
        for i in range(len(chr_lens)):
            if int(chr_lens[i][0]) in chr_lens_dict.keys():
                chr_lens_dict[int(chr_lens[i][0])].append([int(chr_lens[i][1]), int(chr_lens[i][2])])
            else: 
                chr_lens_dict[int(chr_lens[i][0])] = [[int(chr_lens[i][1]), int(chr_lens[i][2])]]

    #Filter for loci containing a primary significant variant and adjust range boundaries where needed
    semi_final_ranges = ["temp"]
    for i in range(len(out_ranges)):
        #Adjust the range boundaries if need be
        cur_chromo = out_ranges[i][0]
        for cur_chromo_range in chr_lens_dict[cur_chromo]:
            if (out_ranges[i][1] >= cur_chromo_range[0] and out_ranges[i][1] < cur_chromo_range[1]) or (out_ranges[i][2] > cur_chromo_range[0] and out_ranges[i][2] <= cur_chromo_range[1]):
                range_start = max(out_ranges[i][1] - stagger_bp, cur_chromo_range[0]) #Added padding equal to the stagger bp
                range_end = min(out_ranges[i][2] + stagger_bp, cur_chromo_range[1]) #Added padding equal to the stagger bp
                range_len = range_end - range_start + 1
                range_temp = [cur_chromo, range_start, range_end]
                #Verify that the locus contains at least one primary significant variant before appending
                if len([x for x in prime_sig_vars if x >= range_temp[0] + float(range_temp[1])/1e9 and x <= range_temp[0] + float(range_temp[2])/1e9]) > 0:
                    #Append new locus to output file list
                    range_temp.extend([range_len])
                    semi_final_ranges.append(range_temp)
    del semi_final_ranges[0]
    #Check that there are at least some significant variants after filtering out forbidden regions
    if len(semi_final_ranges) == 0:
        # Creates an empty file
        with open(output_loci_file, 'w') as out_path:
            pass
        print("NOTE: No significant loci identified for " + base_name)
        sys.exit(0)
    
    #Collapse loci if overlapping!
    semi_final_ranges = sorted(semi_final_ranges, key = lambda x: (x[0], x[1], x[2]))
    final_ranges = [semi_final_ranges[0]]
    for i in range(1, len(semi_final_ranges)):
        if semi_final_ranges[i][0] == final_ranges[-1][0] and semi_final_ranges[i][1] < final_ranges[-1][2]:
            #Collapse locus with last locus
            final_ranges[-1][2] = semi_final_ranges[i][2]
            final_ranges[-1][3] = final_ranges[-1][2] - final_ranges[-1][1] + 1
        else: 
            #No need to collapse loci
            final_ranges.append(semi_final_ranges[i])
    
    #Output message if no rows are left in final_ranges
    if len(final_ranges) == 0:
        # Creates an empty file
        with open(output_loci_file, 'w') as out_path:
            pass
        print("NOTE: No significant loci identified for " + base_name)
        sys.exit(0)
    
    #Write to file
    with open(output_loci_file, "w") as out_path:
        writer = csv.writer(out_path)
        writer.writerow(['CHR', 'START', 'END', 'LENGTH'])  # Header
        for locus in final_ranges:
            writer.writerow(locus)

if __name__ == "__main__":
    main(sys.argv[1:])
