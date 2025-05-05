'''
Select_Causal_Variants.py

The purpose of this script is to select the causal regions and variants to be 
mapped for my RSparsePro robustness test. 

***Run in Conda finemap environment***
'''

import argparse
import random
import os
from collections import defaultdict

###############################################################################
####################  PROCESS INDIVIDUAL SIMULATION FILE ######################
###############################################################################

def read_bed_regions(bed_file):
    """Read BED file into dictionary of {chrom: [(start, end)]}."""
    regions = defaultdict(list)
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            regions[chrom].append((start, end))
    return regions

def read_pruned_variants(prune_file, bed_regions=None):
    """Read and filter variants with optional BED region filtering."""
    variants = []
    with open(prune_file) as f:
        for line in f:
            vid = line.strip()
            try:
                chrom, pos_str, *_ = vid.split(':')
                pos = int(pos_str)
            except (ValueError, IndexError):
                raise ValueError(f"Invalid variant ID format: {vid}")

            # Apply BED filtering if provided
            if bed_regions is not None:
                in_region = False
                if chrom in bed_regions:
                    pos_0based = pos - 1  # Convert to BED's 0-based coordinate
                    for start, end in bed_regions[chrom]:
                        if start <= pos_0based < end:
                            in_region = True
                            break
                if not in_region:
                    continue  # Skip variants outside BED regions

            variants.append((chrom, pos, vid))
    return variants

def create_loci(variant_list, max_size):
    """Group variants into loci smaller than specified maximum size."""
    variant_list.sort(key=lambda x: (x[0], x[1]))
    loci = []
    current_locus = []
    
    for chrom, pos, vid in variant_list:
        if not current_locus or chrom != current_locus[0][0]:
            current_locus = [(chrom, pos, vid)]
        else:
            locus_size = pos - current_locus[0][1]
            if locus_size <= max_size:
                current_locus.append((chrom, pos, vid))
            else:
                loci.append(current_locus)
                current_locus = [(chrom, pos, vid)]
    if current_locus:
        loci.append(current_locus)
    return loci

def select_causal_variants(loci, num_loci, causal_per_locus=1):
    """Select specified number of causal variants per locus."""
    random.shuffle(loci)
    selected = []
    for locus in loci[:num_loci]:
        if len(locus) < causal_per_locus:
            raise ValueError(f"Locus has only {len(locus)} variants, can't select {causal_per_locus}")
        causals = random.sample(locus, causal_per_locus)
        for causal in causals:
            selected.append({
                'variant_id': causal[2],
                'chrom': causal[0],
                'pos': causal[1],
                'locus_size': locus[-1][1] - locus[0][1]
            })
    return selected

def main():
    parser = argparse.ArgumentParser(
        description='Select causal variants from LD-pruned loci',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--prune-in', required=True,
                        help='Input file with pruned variants (format: CHR:BP:REF:ALT)')
    parser.add_argument('--num-loci', type=int, required=True,
                        help='Number of loci to select')
    parser.add_argument('--max-size', type=int, required=True,
                        help='Maximum locus size in base pairs')
    parser.add_argument('--var_per_locus', type=int, required=True,
                        help='Number of variants to select per locus')
    parser.add_argument('--filter-file', 
                        help='BED file defining allowed genomic regions')
    parser.add_argument('--output', required=True,
                        help='Output file path')
    args = parser.parse_args()
    
    # Read BED regions if provided
    bed_regions = read_bed_regions(args.filter_file) if args.filter_file else None
    
    # Read and filter variants
    variants = read_pruned_variants(args.prune_in, bed_regions) 
    # Create loci
    loci = create_loci(variants, args.max_size)
    
    # Validate requested number of loci
    if len(loci) < args.num_loci:
        raise ValueError(f"ERROR: Only {len(loci)} loci available, requested {args.num_loci}")
    
    # Select causal variants
    causal_vars = select_causal_variants(loci, args.num_loci, args.var_per_locus)
    
    # Write output
    with open(args.output, 'w') as f:
        for var in causal_vars:
            f.write(f"{var['variant_id']}\n")

if __name__ == "__main__":
    main()
