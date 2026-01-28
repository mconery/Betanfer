'''

extract_matrices.py

This script extracts specific variants from plink bed/bim/fam files and creates
LD matrices. It handles variant ID matching with allele flipping and reversal.

Input variants are in format: chr:pos:ref:alt
Plink variants may have reversed alleles, which are detected and corrected.

'''

import sys
import getopt
from os.path import exists
import os
import pandas as pd
import numpy as np
import subprocess

def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS

-b => Base name of plink bed/bim/fam files (without extension) REQUIRED

-v => Input file of variants (format: chr:pos:ref:alt) REQUIRED

-c => Chromosome number to process REQUIRED

-o => Output directory for matrices REQUIRED

-p => Output file prefix REQUIRED

""")
    sys.exit(exit_num)

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "b:v:c:o:p:h")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    # Check for help flag
    if '-h' in options_dict:
        help()

    # Required arguments
    try:
        plink_base = options_dict['-b']
        variant_file = options_dict['-v']
        chr_num = int(options_dict['-c'])
        output_dir = options_dict['-o']
        output_prefix = options_dict['-p']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
    except ValueError:
        print("ERROR: Chromosome number must be an integer.")
        help()

    # Verify files exist
    if not exists(variant_file):
        print(f"ERROR: Variant file does not exist: {variant_file}")
        sys.exit(1)

    if not exists(plink_base + '.bim'):
        print(f"ERROR: Plink .bim file does not exist: {plink_base}.bim")
        sys.exit(1)

    if not exists(output_dir):
        print(f"ERROR: Output directory does not exist: {output_dir}")
        sys.exit(1)

    # Call driver
    driver(plink_base, variant_file, chr_num, output_dir, output_prefix)

def reverse_complement(allele):
    """Get reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(allele))

def parse_variant_id(variant_id):
    """Parse variant ID in format chr:pos:ref:alt"""
    parts = variant_id.split(':')
    if len(parts) != 4:
        return None
    try:
        chr_num = int(parts[0])
        pos = int(parts[1])
        ref = parts[2].upper()
        alt = parts[3].upper()
        return (chr_num, pos, ref, alt)
    except:
        return None

def driver(plink_base, variant_file, chr_num, output_dir, output_prefix):
    """Main driver function"""
    
    # Read in variant file and filter for current chromosome
    variants_raw = []
    with open(variant_file) as var_read:
        for line in var_read.readlines():
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parsed = parse_variant_id(line)
            if parsed and parsed[0] == chr_num:
                variants_raw.append(parsed)

    if not variants_raw:
        print(f"NOTE: No variants found for chromosome {chr_num}")
        return

    print(f"Found {len(variants_raw)} variants for chromosome {chr_num}")

    # Read plink bim file
    bim = pd.read_csv(plink_base + '.bim', sep='\t', header=None,
                      names=['chr', 'snp_id', 'cm', 'bp', 'a1', 'a2'])
    bim = bim[bim['chr'] == chr_num].reset_index(drop=True)

    if bim.empty:
        print(f"WARNING: No variants in plink file for chromosome {chr_num}")
        return

    # Create dictionaries for variant matching
    # Key: (chr, pos, ref, alt) -> bim index
    plink_variants = {}
    for idx, row in bim.iterrows():
        chr_val = int(row['chr'])
        pos = int(row['bp'])
        a1 = row['a1'].upper()
        a2 = row['a2'].upper()
        
        # Store both orientations
        plink_variants[(chr_val, pos, a1, a2)] = idx
        plink_variants[(chr_val, pos, a2, a1)] = idx

    # Map input variants to plink indices and track flipping
    variant_to_plink_idx = {}
    variant_flipped = {}  # Track if alleles were flipped
    
    for var in variants_raw:
        chr_val, pos, ref, alt = var
        
        # Try direct match
        if (chr_val, pos, ref, alt) in plink_variants:
            plink_idx = plink_variants[(chr_val, pos, ref, alt)]
            variant_to_plink_idx[var] = plink_idx
            variant_flipped[var] = False
            continue
        
        # Try flipped alleles
        if (chr_val, pos, alt, ref) in plink_variants:
            plink_idx = plink_variants[(chr_val, pos, alt, ref)]
            variant_to_plink_idx[var] = plink_idx
            variant_flipped[var] = True
            continue
        
        print(f"WARNING: Variant not found in plink file: {chr_val}:{pos}:{ref}:{alt}")

    if not variant_to_plink_idx:
        print(f"ERROR: No variants matched between input file and plink file")
        return

    print(f"Matched {len(variant_to_plink_idx)} of {len(variants_raw)} variants to plink file")

    # Get indices of matched variants in plink file
    matched_plink_indices = sorted(set(variant_to_plink_idx.values()))
    
    # Create extract file for plink
    extract_file = plink_base + '.extract.txt'
    with open(extract_file, 'w') as f:
        for idx in matched_plink_indices:
            f.write(bim.iloc[idx]['snp_id'] + '\n')

    # Extract variants from plink using plink1.9
    # Use plink1.9 for correlation matrix
    temp_prefix = plink_base + '.extracted'
    
    cmd = f'plink --bfile {plink_base} --extract {extract_file} --make-bed --out {temp_prefix} --memory 200000'
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running plink extraction: {result.stderr}")
        return

    # Make correlation matrix using plink
    matrix_prefix = os.path.join(output_dir, output_prefix)
    cmd = f'plink --bfile {temp_prefix} --r square --keep-allele-order --out {matrix_prefix} --memory 200000'
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR running plink matrix: {result.stderr}")
        return

    # Read correlation matrix
    matrix_file = matrix_prefix + '.ld'
    if not exists(matrix_file):
        print(f"ERROR: Matrix file not created: {matrix_file}")
        return

    # Read the matrix (space-separated, square format)
    with open(matrix_file, 'r') as f:
        lines = f.readlines()
    
    n_vars = len(lines)
    matrix = np.zeros((n_vars, n_vars))
    
    for i, line in enumerate(lines):
        values = [float(x) for x in line.strip().split()]
        matrix[i, :] = values

    # Apply allele flipping corrections
    # Create mapping from plink indices back to original variant order
    plink_to_original_idx = {}
    for orig_var, plink_idx in variant_to_plink_idx.items():
        # Find position in matched list
        pos_in_matched = sorted(set(variant_to_plink_idx.values())).index(plink_idx)
        plink_to_original_idx[pos_in_matched] = (orig_var, variant_flipped[orig_var])

    # Multiply rows and columns by -1 if alleles are flipped
    for i in range(n_vars):
        if plink_to_original_idx[i][1]:  # If flipped
            matrix[i, :] *= -1
            matrix[:, i] *= -1

    # Write corrected matrix
    with open(matrix_file, 'w') as f:
        for row in matrix:
            f.write(' '.join(f'{val:.6f}' for val in row) + '\n')

    # Create variant order file
    var_order_file = matrix_prefix + '.variants.txt'
    with open(var_order_file, 'w') as f:
        for i in range(n_vars):
            orig_var, _ = plink_to_original_idx[i]
            f.write(f'{orig_var[0]}:{orig_var[1]}:{orig_var[2]}:{orig_var[3]}\n')

    print(f"Successfully created matrix: {matrix_file}")
    print(f"Variant order: {var_order_file}")

    # Clean up temporary files
    for ext in ['.bed', '.bim', '.fam', '.log']:
        temp_file = temp_prefix + ext
        if exists(temp_file):
            os.remove(temp_file)
    
    if exists(extract_file):
        os.remove(extract_file)

    print(f"Completed matrix creation for {output_prefix}")

if __name__ == "__main__":
    main(sys.argv)
