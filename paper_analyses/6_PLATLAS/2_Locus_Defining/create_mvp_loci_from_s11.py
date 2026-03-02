'''
create_mvp_loci_from_s11.py

The purpose of this script is to read predefined locus boundaries from
Supplementary Table S11 and write per-population master locus CSV files
for the two MVP populations (AFR and EUR).

S11 contains one row per fine-mapping signal within each locus; this script
deduplicates to one row per unique (phenotype, chr, start, end) combination.
The same locus set is applied to both MVP populations.

Output master file format matches that produced by locus_parsl_allocator.py:
  CHR, START, END, LENGTH, TRAIT, OVERLAP_COUNT

Where:
  CHR          - chromosome number (integer, no "chr" prefix)
  START        - locus start position (integer)
  END          - locus end position (integer)
  LENGTH       - end - start + 1 (integer)
  TRAIT        - phenotype name from S11
  OVERLAP_COUNT - number of other loci (from any trait) overlapping this locus

Usage:
  python create_mvp_loci_from_s11.py --s11_file <path> --output_dir <path>

***Run in Conda finemap_old environment***
'''

import sys
import getopt
import os
import pandas as pd
from intervaltree import Interval, IntervalTree


def help(exit_num=0):
    print("""Usage: python create_mvp_loci_from_s11.py [OPTIONS]
Options:
  --s11_file    Path to Supplementary_Table-S11.txt (tab-delimited, required)
  --output_dir  Directory where master CSV files will be written (required)
  -h, --help    Show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  HELPER: COMPUTE OVERLAP_COUNT VIA INTERVALTREE  ##################
###############################################################################

def compute_overlap_counts(df):
    """Add OVERLAP_COUNT column: number of OTHER loci overlapping each locus.

    Replicates the logic from locus_parsl_allocator.py:combine_population.
    Groups by CHR, builds an IntervalTree, then counts overlaps minus self.
    """
    df = df.copy()
    df['OVERLAP_COUNT'] = 0
    for chr_name, chr_group in df.groupby('CHR'):
        tree = IntervalTree()
        for idx, row in chr_group.iterrows():
            tree.add(Interval(row['START'], row['END'], idx))
        for idx, row in chr_group.iterrows():
            overlaps = tree.overlap(row['START'], row['END'])
            df.loc[idx, 'OVERLAP_COUNT'] = len(overlaps) - 1  # subtract self
    return df


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND BUILD MASTER FILES  ############
###############################################################################

def main():

    # ── Parse arguments ────────────────────────────────────────────────────────
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "s11_file=", "output_dir="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    s11_file   = None
    output_dir = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--s11_file":
            s11_file = a
        elif o == "--output_dir":
            output_dir = a.rstrip('/')
        else:
            assert False, f"WARNING: Unhandled input option: {o}"

    # ── Validate required arguments ────────────────────────────────────────────
    if s11_file is None:
        print("ERROR: Missing required --s11_file")
        sys.exit(2)
    if output_dir is None:
        print("ERROR: Missing required --output_dir")
        sys.exit(2)

    for path, label in [(s11_file, "S11 file"), (output_dir, "output directory")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Read S11 and deduplicate loci ──────────────────────────────────────────
    print(f"NOTE: Reading S11 file: {s11_file}")
    s11 = pd.read_csv(s11_file, sep='\t', usecols=['phenotype', 'chr', 'start', 'end'])
    s11 = s11.drop_duplicates()

    # Strip "chr" prefix and cast to int
    s11['chr'] = s11['chr'].str.replace('^chr', '', regex=True).astype(int)
    s11['start'] = s11['start'].astype(int)
    s11['end']   = s11['end'].astype(int)

    # Rename to match master file column conventions
    s11 = s11.rename(columns={'phenotype': 'TRAIT', 'chr': 'CHR',
                               'start': 'START', 'end': 'END'})
    s11['LENGTH'] = s11['END'] - s11['START'] + 1

    n_traits = s11['TRAIT'].nunique()
    n_loci   = len(s11)
    print(f"NOTE: {n_loci} unique loci across {n_traits} traits found in S11")

    # ── Compute OVERLAP_COUNT ──────────────────────────────────────────────────
    s11 = compute_overlap_counts(s11)

    # ── Reorder columns to match master file format ────────────────────────────
    master = s11[['CHR', 'START', 'END', 'LENGTH', 'TRAIT', 'OVERLAP_COUNT']]

    # ── Write master files for both MVP populations ────────────────────────────
    for pop in ['MVP_AFR', 'MVP_EUR']:
        out_path = os.path.join(output_dir, f"sig_loci.{pop}.master.csv")
        master.to_csv(out_path, index=False)
        print(f"SUCCESS: Wrote {out_path} ({n_loci} loci from {n_traits} traits)")

    print("NOTE: create_mvp_loci_from_s11.py completed successfully.")


if __name__ == "__main__":
    main()
