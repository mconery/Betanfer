'''

sample_loci.py

The purpose of this script is to randomly sample 100 unique 750 kb locus
coordinates from the PLATLAS MVP master locus file and write them to a CSV
for use in the Frontier LD matrix pipeline.

Both MVP_EUR and MVP_AFR master files contain identical locus sets (both are
derived from Supplementary Table S11).  We therefore read only the EUR file
and produce a single coordinate list that is valid for both populations.

Output CSV columns: CHR, START, END  (hg19 coordinates)

Usage:
  python sample_loci.py [-i <master_csv>] [-o <output_csv>] [-n <n>] [-s <seed>]

***Run in Conda finemap_old environment***

'''

import sys
import getopt
import os
import pandas as pd


# ── Defaults ────────────────────────────────────────────────────────────────
DEFAULT_INPUT  = ("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/"
                  "mapping_loci/sig_loci.MVP_EUR.master.csv")
DEFAULT_OUTPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "sampled_loci_750kb.csv")
DEFAULT_N      = 100
DEFAULT_SEED   = 42
TARGET_LENGTH  = 750000


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
-i => Input master locus CSV (default: sig_loci.MVP_EUR.master.csv)
-o => Output sampled locus CSV (default: sampled_loci_750kb.csv next to this script)
-n => Number of loci to sample (default: 100)
-s => Random seed (default: 42)
-h => Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:n:s:h")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help(1)

    options_dict = dict(opts)

    if '-h' in options_dict:
        help()

    input_csv  = options_dict.get('-i', DEFAULT_INPUT)
    output_csv = options_dict.get('-o', DEFAULT_OUTPUT)
    n_sample   = int(options_dict.get('-n', DEFAULT_N))
    seed       = int(options_dict.get('-s', DEFAULT_SEED))

    if not os.path.exists(input_csv):
        print(f"ERROR: Input master CSV not found: {input_csv}")
        sys.exit(1)

    if n_sample < 1:
        print("ERROR: -n must be a positive integer.")
        sys.exit(1)

    print("Acceptable Inputs Given")

    driver(input_csv, output_csv, n_sample, seed)


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(input_csv, output_csv, n_sample, seed):

    # ── Read master locus file ─────────────────────────────────────────────
    print(f"NOTE: Reading master locus file: {input_csv}")
    master = pd.read_csv(input_csv)

    required_cols = {'CHR', 'START', 'END', 'LENGTH'}
    if not required_cols.issubset(master.columns):
        print(f"ERROR: Input CSV missing required columns. Expected: {required_cols}")
        sys.exit(1)

    n_total = len(master)
    print(f"NOTE: {n_total} total loci in master file")

    # ── Filter to 750 kb loci ─────────────────────────────────────────────
    master_750 = master[master['LENGTH'] == TARGET_LENGTH].copy()
    print(f"NOTE: {len(master_750)} loci with LENGTH == {TARGET_LENGTH}")

    if len(master_750) == 0:
        print(f"ERROR: No loci found with LENGTH == {TARGET_LENGTH}")
        sys.exit(1)

    # ── Deduplicate to unique (CHR, START, END) coordinates ───────────────
    unique_loci = (master_750[['CHR', 'START', 'END']]
                   .drop_duplicates()
                   .reset_index(drop=True))
    n_unique = len(unique_loci)
    print(f"NOTE: {n_unique} unique (CHR, START, END) coordinates after deduplication")

    if n_unique < n_sample:
        print(f"ERROR: Requested {n_sample} loci but only {n_unique} unique coordinates available.")
        sys.exit(1)

    # ── Random sample ─────────────────────────────────────────────────────
    print(f"NOTE: Setting random seed to {seed}")
    sampled = unique_loci.sample(n=n_sample, random_state=seed)
    sampled = sampled.sort_values(['CHR', 'START']).reset_index(drop=True)

    print(f"NOTE: Sampled {len(sampled)} loci")

    # ── Write output ──────────────────────────────────────────────────────
    out_dir = os.path.dirname(output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    sampled.to_csv(output_csv, index=False)
    print(f"SUCCESS: Wrote {len(sampled)} sampled loci to {output_csv}")

    # ── Sanity check ──────────────────────────────────────────────────────
    n_dup = sampled.duplicated(subset=['CHR', 'START', 'END']).sum()
    if n_dup > 0:
        print(f"WARNING: {n_dup} duplicate (CHR, START, END) entries in output — investigate!")
    else:
        print("NOTE: All sampled loci have unique (CHR, START, END) coordinates.")


###############################################################################
################################ CALL MAIN ###################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv)
