'''
count_causal_snps_in_sig_loci.py

The purpose of this script is to count how many of the causal SNPs are
contained within one of the significant loci from the simulated GWAS.

***This should be run in the rsparsepro conda environment***
'''

import sys
import getopt
import csv

def help(exit_code=0):
    print("""
-----------------------------------------------------------------
Counts how many causal SNPs are found within significant loci.

USAGE:
  python count_causal_snps_in_loci.py -l <loci_csv> -c <causal_snps_file>

ARGUMENTS:
  -l, --loci         Path to significant loci CSV file (output from define_sig_loci_for_trait.py) [REQUIRED]
  -c, --causal       Path to file with causal SNP IDs, one per line [REQUIRED]
  -h, --help         Show this help message and exit

The loci CSV must have columns: CHR, START, END.
The causal SNP file must have IDs in the format: chromosome:position:ref:alt
Only chromosome and position are used for matching.
-----------------------------------------------------------------
""")
    sys.exit(exit_code)

def read_loci(loci_csv):
    loci = []
    with open(loci_csv, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            chrom = str(row['CHR'])
            start = int(row['START'])
            end = int(row['END'])
            loci.append((chrom, start, end))
    return loci

def read_causal_snps(causal_file):
    snps = []
    with open(causal_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(':')
            if len(parts) < 2:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            snps.append((chrom, pos))
    return snps

def count_snps_in_loci(loci, snps):
    count = 0
    for chrom, pos in snps:
        for locus_chrom, start, end in loci:
            if chrom == locus_chrom and start <= pos <= end:
                count += 1
                break
    return count

def main(argv):
    loci_csv = None
    causal_file = None

    try:
        opts, args = getopt.getopt(argv, "hl:c:", ["help", "loci=", "causal="])
    except getopt.GetoptError as err:
        print("ERROR:", err)
        help(1)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            help(0)
        elif opt in ("-l", "--loci"):
            loci_csv = arg
        elif opt in ("-c", "--causal"):
            causal_file = arg

    if not loci_csv or not causal_file:
        print("ERROR: Both -l/--loci and -c/--causal arguments are required.")
        help(1)

    loci = read_loci(loci_csv)
    snps = read_causal_snps(causal_file)
    count = count_snps_in_loci(loci, snps)
    print(f"Out of {len(snps)} causal snps {count} are contained in a significant locus.")

if __name__ == "__main__":
    main(sys.argv[1:])
