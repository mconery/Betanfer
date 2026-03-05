'''
convert_mvp_to_hg19.py

The purpose of this script is to convert a single MVP summary statistics file
from GRCh38 to hg19 coordinates.

The annotation file (/flare/GeomicVar/multi-biobank-meta/Annotation_gwPheWAS_dbGaP_b38.noAF_R2.txt.gz)
contains the mapping from GRCh38 positions (CHROM_b38, POS_b38) to hg19 positions
embedded in the MVPID column (format: CHR:POS:REF:ALT).  It also contains RSID.

For each variant in the sumstats, the script looks up the matching annotation
row by (CHR, POS, REF, ALT) in GRCh38 space.  Matched variants have their #ID,
CHR, and POS replaced with hg19 values from MVPID, and an RSID column is
inserted.  Variants with no annotation match are dropped.

Output column order:
  #ID  RSID  CHR  POS  REF  ALT  BETA  SE  P  LOG10P  ...  (all remaining cols)

Usage:
  python convert_mvp_to_hg19.py [OPTIONS]

  --input_file       Path to input .sumstats.txt.gz file (required)
  --output_file      Path for output .sumstats.hg19.txt file (required)
  --annotation_file  Path to annotation .txt.gz file (required)
  -h, --help         Show this help message

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os
import pandas as pd


def help(exit_num=0):
    print("""Usage: python convert_mvp_to_hg19.py [OPTIONS]
Options:
  --input_file       Path to input .sumstats.txt.gz file (required)
  --output_file      Path for output .sumstats.hg19.txt file (required)
  --annotation_file  Path to annotation .txt.gz mapping b38 to hg19 (required)
  -h, --help         Show this help message""")
    sys.exit(exit_num)


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND CONVERT FILE  #################
###############################################################################

def main():

    # ── Parse arguments ────────────────────────────────────────────────────────
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "input_file=", "output_file=", "annotation_file="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    input_file      = None
    output_file     = None
    annotation_file = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--input_file":
            input_file = a
        elif o == "--output_file":
            output_file = a
        elif o == "--annotation_file":
            annotation_file = a
        else:
            assert False, f"WARNING: Unhandled input option: {o}"

    # ── Validate required arguments ────────────────────────────────────────────
    if input_file is None:
        print("ERROR: Missing required --input_file")
        sys.exit(2)
    if output_file is None:
        print("ERROR: Missing required --output_file")
        sys.exit(2)
    if annotation_file is None:
        print("ERROR: Missing required --annotation_file")
        sys.exit(2)

    for path, label in [(input_file, "input file"), (annotation_file, "annotation file")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # Create output directory if needed
    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # ── Load annotation file ───────────────────────────────────────────────────
    # Only load columns needed for the join and coordinate replacement
    print(f"NOTE: Loading annotation file: {annotation_file}")
    annot = pd.read_csv(
        annotation_file,
        sep=' ',
        usecols=['MVPID', 'RSID', 'CHROM_b38', 'POS_b38'],
        dtype={'MVPID': str, 'RSID': str, 'CHROM_b38': str, 'POS_b38': int}
    )
    print(f"NOTE: Annotation loaded — {len(annot):,} rows")

    # Parse MVPID (hg19 CHR:POS:REF:ALT) into component columns
    mvpid_parts = annot['MVPID'].str.split(':', expand=True)
    annot['hg19_CHR'] = mvpid_parts[0].astype(int)
    annot['hg19_POS'] = mvpid_parts[1].astype(int)
    annot['mvp_REF']  = mvpid_parts[2]
    annot['mvp_ALT']  = mvpid_parts[3]

    # Normalize CHROM_b38: strip "chr" prefix if present, cast to int
    annot['b38_CHR'] = annot['CHROM_b38'].str.replace('^chr', '', regex=True).astype(int)

    # Build join key and create lookup DataFrame indexed on (b38_CHR, POS_b38, REF, ALT)
    annot_lookup = annot[['b38_CHR', 'POS_b38', 'mvp_REF', 'mvp_ALT',
                           'hg19_CHR', 'hg19_POS', 'RSID', 'MVPID']].copy()
    annot_lookup = annot_lookup.rename(columns={'mvp_REF': 'REF', 'mvp_ALT': 'ALT'})

    # ── Load sumstats ──────────────────────────────────────────────────────────
    print(f"NOTE: Loading sumstats: {input_file}")
    sumstats = pd.read_csv(
        input_file,
        sep='\t',
        dtype={'#ID': str, 'CHR': int, 'POS': int, 'REF': str, 'ALT': str}
    )
    n_input = len(sumstats)
    print(f"NOTE: Sumstats loaded — {n_input:,} variants")

    # ── Join sumstats to annotation on (CHR, POS, REF, ALT) in b38 space ──────
    merged = sumstats.merge(
        annot_lookup,
        left_on=['CHR', 'POS', 'REF', 'ALT'],
        right_on=['b38_CHR', 'POS_b38', 'REF', 'ALT'],
        how='left'
    )

    n_matched = merged['MVPID'].notna().sum()
    n_dropped = n_input - n_matched
    print(f"NOTE: {n_matched:,} variants matched in annotation")
    if n_dropped > 0:
        print(f"WARNING: {n_dropped:,} variants had no annotation match and will be dropped")

    # Drop unmatched rows
    merged = merged[merged['MVPID'].notna()].copy()

    # ── Replace b38 coordinates with hg19 ─────────────────────────────────────
    merged['#ID'] = merged['MVPID']
    merged['CHR'] = merged['hg19_CHR'].astype(int)
    merged['POS'] = merged['hg19_POS'].astype(int)

    # ── Assemble output column order: #ID RSID CHR POS REF ALT ... rest ────────
    # Identify original remaining columns (after REF ALT)
    original_cols = list(sumstats.columns)
    fixed_cols    = ['#ID', 'CHR', 'POS', 'REF', 'ALT']
    extra_cols    = [c for c in original_cols if c not in fixed_cols]

    output_cols = ['#ID', 'RSID', 'CHR', 'POS', 'REF', 'ALT'] + extra_cols
    out = merged[output_cols]

    # ── Write output (atomic: write to .tmp then rename) ──────────────────────
    tmp_file = output_file + ".tmp"
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    print(f"NOTE: Writing output: {output_file}")
    out.to_csv(tmp_file, sep='\t', index=False)
    os.rename(tmp_file, output_file)
    print(f"SUCCESS: Wrote {len(out):,} variants to {output_file}")


if __name__ == "__main__":
    main()
