'''

build_rsid_lookup.py

The purpose of this script is to build an rsID -> genomic position lookup
table from the MVP hg19 summary statistics files for the 100 sampled PLATLAS
loci.

Betanfer LD matrices store variant IDs as rsIDs (taken from the RSID column of
the MVP sumstats), while the raw plink-based MVP matrices use chr:pos:a1:a2
format (from the plink BIM file).  To compare the two, we need a table that
maps each rsID to its genomic position and alleles.

The MVP hg19 sumstats are bgzip-compressed and tabix-indexed, so each locus
region is queried efficiently via tabix.  To maximise variant coverage, a small
number of sumstats files per population is scanned (default: 3).  The rsID->
position mapping is consistent across files for the same variant.

Output TSV columns: RSID, CHR, POS, REF, ALT

Usage:
  python build_rsid_lookup.py [-l <locus_csv>] [-s <sumstats_base>]
                              [-o <output_tsv>] [-n <n_files>]

***Run in Conda finemap_old environment on Aurora***

'''

import sys
import getopt
import os
import gzip
import glob
import subprocess
import csv
import pandas as pd

# ── Defaults ──────────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

DEFAULT_LOCUS_CSV     = os.path.join(_SCRIPT_DIR, "sampled_loci_750kb.csv")
DEFAULT_SUMSTATS_BASE = "/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19"
DEFAULT_OUTPUT        = os.path.join(_SCRIPT_DIR, "rsid_pos_lookup.tsv")
DEFAULT_N_FILES       = 3

POPULATIONS = ["EUR", "AFR"]


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
-l => Sampled locus CSV with columns CHR,START,END
      (default: sampled_loci_750kb.csv next to this script)
-s => MVP hg19 sumstats base directory containing EUR/ and AFR/ subdirs
      (default: /flare/.../sumstats/MVP_hg19)
-o => Output lookup TSV file path
      (default: rsid_pos_lookup.tsv next to this script)
-n => Number of sumstats files to scan per population (default: 3)
-h => Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "l:s:o:n:h")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help(1)

    options_dict = dict(opts)

    if '-h' in options_dict:
        help()

    locus_csv     = options_dict.get('-l', DEFAULT_LOCUS_CSV)
    sumstats_base = options_dict.get('-s', DEFAULT_SUMSTATS_BASE)
    output_tsv    = options_dict.get('-o', DEFAULT_OUTPUT)
    n_files       = int(options_dict.get('-n', DEFAULT_N_FILES))

    for path, label in [
        (locus_csv,     "locus CSV"),
        (sumstats_base, "sumstats base directory"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    print("Acceptable Inputs Given")

    driver(locus_csv, sumstats_base, output_tsv, n_files)


###############################################################################
############################# HELPFUL FUNCTIONS ##############################
###############################################################################

def get_header_cols(sumstats_file):
    """
    Read the first line of a gzip file and return a dict mapping
    column name -> column index.  Leading '#' is stripped.
    """
    with gzip.open(sumstats_file, 'rt') as fh:
        header_line = fh.readline().strip()
    cols = header_line.lstrip('#').split('\t')
    return {col.strip(): i for i, col in enumerate(cols)}


def tabix_query(sumstats_file, chr_num, start, end):
    """
    Run tabix to fetch lines in region chr_num:start-end.
    Returns a list of tab-split field lists (data lines only).
    """
    region = f"{chr_num}:{start}-{end}"
    result = subprocess.run(
        ['tabix', sumstats_file, region],
        capture_output=True, text=True
    )
    if result.returncode != 0 or not result.stdout.strip():
        return []
    lines = []
    for line in result.stdout.strip().split('\n'):
        if line and not line.startswith('#'):
            lines.append(line.split('\t'))
    return lines


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(locus_csv, sumstats_base, output_tsv, n_files):

    # ── Read sampled loci ──────────────────────────────────────────────────
    print(f"NOTE: Reading locus file: {locus_csv}")
    loci = pd.read_csv(locus_csv)
    print(f"NOTE: {len(loci)} loci to process")

    # ── Build lookup ───────────────────────────────────────────────────────
    # rsid → (chr, pos, ref, alt)
    lookup = {}

    for pop in POPULATIONS:
        pop_dir = os.path.join(sumstats_base, pop)
        if not os.path.exists(pop_dir):
            print(f"WARNING: Population directory not found: {pop_dir} — skipping {pop}")
            continue

        # Find sumstats files for this population (exclude hm3 files)
        all_gz = glob.glob(
            os.path.join(pop_dir, "MVP_R4.1000G_AGR.*.sumstats.hg19.txt.gz")
        )
        gz_files = sorted(
            [f for f in all_gz if not f.endswith("hm3_sumstats.txt.gz")]
        )[:n_files]

        if not gz_files:
            print(f"WARNING: No sumstats files found in {pop_dir} — skipping {pop}")
            continue

        print(f"NOTE: {pop}: scanning {len(gz_files)} sumstats file(s)")

        for sumstats_file in gz_files:
            print(f"  Reading: {os.path.basename(sumstats_file)}")

            # Check tabix index exists
            if not os.path.exists(sumstats_file + ".tbi"):
                print(f"    WARNING: No .tbi index found — skipping this file")
                continue

            # Get column positions from header
            col_idx = get_header_cols(sumstats_file)
            rsid_col = col_idx.get('RSID')
            chr_col  = col_idx.get('CHR')
            pos_col  = col_idx.get('POS')
            ref_col  = col_idx.get('REF')
            alt_col  = col_idx.get('ALT')

            if rsid_col is None:
                print(f"    WARNING: No RSID column in header — skipping this file")
                continue

            if any(c is None for c in [chr_col, pos_col, ref_col, alt_col]):
                print(f"    WARNING: Missing required column(s) in header — skipping")
                continue

            file_added = 0
            for _, locus in loci.iterrows():
                chr_num = int(locus['CHR'])
                start   = int(locus['START'])
                end     = int(locus['END'])

                rows = tabix_query(sumstats_file, chr_num, start, end)
                for fields in rows:
                    try:
                        rsid = fields[rsid_col].strip()
                        pos  = int(fields[pos_col])
                        ref  = fields[ref_col].upper().strip()
                        alt  = fields[alt_col].upper().strip()
                    except (IndexError, ValueError):
                        continue

                    # Skip missing rsIDs
                    if not rsid or rsid == '.' or rsid == 'NA':
                        continue

                    # Only add if not already in lookup (first file wins)
                    if rsid not in lookup:
                        lookup[rsid] = (chr_num, pos, ref, alt)
                        file_added += 1

            print(f"    Added {file_added} new rsID entries from this file "
                  f"(lookup total: {len(lookup)})")

    print(f"\nNOTE: Final lookup size: {len(lookup)} unique rsIDs")

    if len(lookup) == 0:
        print("ERROR: Lookup is empty — check sumstats paths and tabix availability")
        sys.exit(1)

    # ── Write output TSV ───────────────────────────────────────────────────
    out_dir = os.path.dirname(output_tsv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_tsv, 'w', newline='') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['RSID', 'CHR', 'POS', 'REF', 'ALT'])
        for rsid, (chrom, pos, ref, alt) in sorted(lookup.items()):
            writer.writerow([rsid, chrom, pos, ref, alt])

    print(f"SUCCESS: Wrote {len(lookup)} rsID entries to {output_tsv}")


###############################################################################
################################ CALL MAIN ###################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv)
