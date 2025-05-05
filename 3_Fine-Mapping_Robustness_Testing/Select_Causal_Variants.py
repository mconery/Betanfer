'''
Select_Causal_Variants.py

The purpose of this script is to select the causal regions and variants to be 
mapped for my RSparsePro robustness test. 

***Run in Conda finemap environment***
'''

import argparse
import random

###############################################################################
####################  PROCESS INDIVIDUAL SIMULATION FILE ######################
###############################################################################

def read_pruned_variants(prune_file):
    """Read pruned variants and parse genomic positions from variant IDs."""
    variants = []
    with open(prune_file) as f:
        for line in f:
            vid = line.strip()
            try:
                chrom, pos, *_ = vid.split(':')
                variants.append((chrom, int(pos), vid))
            except (ValueError, IndexError):
                raise ValueError(f"Invalid variant ID format: {vid}")
    return variants

def create_loci(variant_list, max_size, var_per_locus):
    """Group variants into loci smaller than specified maximum size."""
    variant_list.sort(key=lambda x: (x[0], x[1]))
    variant_list = [x for x in variant_list if len(x) >= var_per_locus]
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

def select_causal_variants(loci, num_loci, var_per_locus):
    """Randomly select loci and choose one causal variant per locus."""
    random.shuffle(loci)
    selected = []
    for locus in loci[:num_loci]:
        causal = random.sample(locus, var_per_locus)
        for variant in causal:
            selected.append({
                'variant_id': variant[2],
                'chrom': variant[0],
                'pos': variant[1],
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
    parser.add_argument('--output', required=True,
                        help='Output file path')
    args = parser.parse_args()

    # Read and process variants
    variants = read_pruned_variants(args.prune_in)
    loci = create_loci(variants, args.max_size, args.var_per_locus)
    
    # Validate requested number of loci
    if len(loci) < args.num_loci:
        raise ValueError(f"Only {len(loci)} loci available, requested {args.num_loci}")
    
    # Select causal variants
    causal_vars = select_causal_variants(loci, args.num_loci)
    
    # Write output
    with open(args.output, 'w') as f:
        for var in causal_vars:
            f.write(f"{var['variant_id']}\n")

if __name__ == "__main__":
    main()
