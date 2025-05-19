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
            chrom = int(parts[0].lstrip("chr"))
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
                chrom = int(chrom)
                pos = int(pos_str)
            except (ValueError, IndexError):
                raise ValueError(f"ERROR: Invalid variant ID format: {vid}")
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
        chrom=int(chrom)
        pos=int(pos)
        if not current_locus or chrom != current_locus[0][0]:
            current_locus = [[chrom, pos, vid]]
        else:
            locus_size = pos - current_locus[0][1]
            if locus_size <= max_size:
                current_locus.append([chrom, pos, vid])
            else:
                loci.append(current_locus)
                current_locus = [[chrom, pos, vid]]
    if current_locus:
        loci.append(current_locus)
    return loci

def select_causal_variants(loci, num_loci, causal_per_locus=1, min_distance=1000000, padding=500000):
    """Select specified number of causal variants per locus."""
    random.shuffle(loci)
    selected = []
    selected_loci = []
    bed_entries = []
    for locus in loci:
        if len(locus) < causal_per_locus:
            raise ValueError(f"ERROR: Locus has only {len(locus)} variants, can't select {causal_per_locus}")
        #Get bounds of locus
        locus_chr = locus[0][0]
        locus_min = min([var[1] for var in locus])
        locus_max = max([var[1] for var in locus])
        #Verify that new locus is above min distance from current loci
        if any([True if selected_locus[0] == locus_chr and (abs(locus_min - selected_locus[2]) < min_distance or abs(selected_locus[1] - locus_max) < min_distance) else False for selected_locus in selected_loci]):
            continue #Skipping current locus because it's too close to previously selected locus
        else:
            selected_loci.append([locus_chr, locus_min, locus_max])
            bed_entries.append({'chrom':locus_chr, 'start':max(locus_min - padding, 0), 'end':(locus_max + padding)})
        causals = random.sample(locus, causal_per_locus)
        for causal in causals:
            selected.append({
                'variant_id': causal[2],
                'chrom': causal[0],
                'pos': causal[1],
                'locus_size': locus[-1][1] - locus[0][1]
                })
        if len(selected_loci) >= num_loci:
            return selected, bed_entries
    raise ValueError("ERROR: Insufficient number of loci found")

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
    parser.add_argument('--min_distance', type=int, required=True,
                        help='Minimum distance between loci in base pairs')
    parser.add_argument('--var_per_locus', type=int, required=True,
                        help='Number of variants to select per locus')
    parser.add_argument('--filter-file', 
                        help='BED file defining allowed genomic regions')
    parser.add_argument('--out_dir', required=True,
                        help='Output directory path')
    parser.add_argument('--out_prefix', required=True,
                        help='Prefix for output files')
    parser.add_argument('--padding', type=int, default=250000,
                        help='Padding distance for locus boundaries (bp)')
    args = parser.parse_args()
    
    #Create ouptut file names
    variant_out_file = os.path.join(args.out_dir, args.out_prefix + '.variants.txt')
    loci_out_file = os.path.join(args.out_dir, args.out_prefix + '.loci.bed')
    
    # Read BED regions if provided
    bed_regions = read_bed_regions(args.filter_file) if args.filter_file else None
    
    # Read and filter variants
    variants = read_pruned_variants(args.prune_in, bed_regions) 
    # Create loci
    loci = create_loci(variants, args.max_size)
    
    #Filter down loci for those with sufficient number of causal variants
    acceptable_loci = [x for x in loci if len(x) >= args.var_per_locus]
    
    # Validate requested number of loci
    if len(acceptable_loci) < args.num_loci:
        raise ValueError(f"ERROR: Only {len(loci)} loci with {args.var_per_locus} available, requested {args.num_loci}")
    
    # Select causal variants
    causal_vars, causal_loci = select_causal_variants(acceptable_loci, args.num_loci, args.var_per_locus, args.min_distance, args.padding)
    
    # Write variant output
    with open(variant_out_file, 'w') as f:
        for var in causal_vars:
            f.write(f"{var['variant_id']}\n")
    
    #Write locus output
    with open(loci_out_file, 'w') as f:
        for chrom, start, end in causal_loci:
            f.write(f"{chrom}\t{start}\t{end}\n")

if __name__ == "__main__":
    main()
