'''

build_maf_lookup.py

Build an rsID -> minor allele frequency (MAF) lookup from MVP hg19 summary
statistics, separately for EUR and AFR populations.

The rsid_pos_lookup.json (built by build_rsid_lookup.py) defines the universe
of rsIDs that appear in the Betanfer LD matrices.  For each rsID in that
universe, this script finds the AAF value in a representative MVP sumstats
file per population and computes MAF = min(AAF, 1-AAF).

A single representative sumstats file is used per population (default:
A1C_Max_INT).  Variant AAF is derived from genotype data and is
trait-independent, so any large MVP trait file gives reliable MAF estimates.
rsIDs absent from the chosen file receive NA.

Output: TSV with columns rsid, EUR_maf, AFR_maf.

Usage:
  python build_maf_lookup.py [-r <rsid_lookup>] [-e <eur_sumstats>]
                             [-a <afr_sumstats>] [-o <output_tsv>]

***Run in Conda finemap_old environment on Aurora***

'''

import sys
import getopt
import os
import json
import pandas as pd

# ── Defaults ──────────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_SUMSTATS_BASE = "/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19"
_REP_TRAIT = "MVP_R4.1000G_AGR.A1C_Max_INT.{pop}.GIA.dbGaP.sumstats.hg19.txt"

DEFAULT_RSID_LOOKUP = os.path.join(_SCRIPT_DIR, "rsid_pos_lookup.json")
DEFAULT_EUR_SUMSTATS = os.path.join(_SUMSTATS_BASE, "EUR", _REP_TRAIT.format(pop="EUR"))
DEFAULT_AFR_SUMSTATS = os.path.join(_SUMSTATS_BASE, "AFR", _REP_TRAIT.format(pop="AFR"))
DEFAULT_OUTPUT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545/maf_lookup.tsv"


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
-r => rsID position lookup JSON (rsid_pos_lookup.json)
      (default: rsid_pos_lookup.json next to this script)
-e => EUR MVP hg19 sumstats .txt file (uncompressed)
      (default: A1C_Max_INT EUR file in standard sumstats base)
-a => AFR MVP hg19 sumstats .txt file (uncompressed)
      (default: A1C_Max_INT AFR file in standard sumstats base)
-o => Output TSV path (rsid, EUR_maf, AFR_maf)
      (default: mvp_results/RITM0322545/maf_lookup.tsv)
-h => Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "r:e:a:o:h")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help(1)

    options_dict = dict(opts)

    if '-h' in options_dict:
        help()

    rsid_lookup   = options_dict.get('-r', DEFAULT_RSID_LOOKUP)
    eur_sumstats  = options_dict.get('-e', DEFAULT_EUR_SUMSTATS)
    afr_sumstats  = options_dict.get('-a', DEFAULT_AFR_SUMSTATS)
    output_tsv    = options_dict.get('-o', DEFAULT_OUTPUT)

    for path, label in [
        (rsid_lookup,  "rsID lookup JSON"),
        (eur_sumstats, "EUR sumstats file"),
        (afr_sumstats, "AFR sumstats file"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    print("Acceptable Inputs Given")

    driver(rsid_lookup, eur_sumstats, afr_sumstats, output_tsv)


###############################################################################
############################# HELPFUL FUNCTIONS ##############################
###############################################################################

def load_rsids(rsid_lookup):
    """Return the set of rsIDs from rsid_pos_lookup.json."""
    print(f"NOTE: Loading rsID lookup: {rsid_lookup}")
    with open(rsid_lookup) as fh:
        lookup = json.load(fh)
    rsids = set(lookup.keys())
    print(f"NOTE: {len(rsids):,} rsIDs in lookup")
    return rsids


def read_aaf(sumstats_file, known_rsids, pop_label):
    """
    Read RSID and AAF columns from a sumstats file, filter to known_rsids,
    and return a Series indexed by rsID with MAF = min(AAF, 1-AAF) values.
    The header's first column may be '#ID'; pandas reads it as-is, so RSID
    and AAF column lookups are unaffected.
    """
    print(f"NOTE: Reading {pop_label} sumstats: {os.path.basename(sumstats_file)}")
    df = pd.read_csv(sumstats_file, sep='\t', usecols=['RSID', 'AAF'],
                     dtype={'RSID': str, 'AAF': float})
    print(f"NOTE: {len(df):,} rows read from {pop_label} sumstats")

    df = df[df['RSID'].isin(known_rsids)].copy()
    print(f"NOTE: {len(df):,} rows matching known rsIDs")

    df['maf'] = df['AAF'].clip(0, 1)
    df['maf'] = df['maf'].where(df['maf'] <= 0.5, 1 - df['maf'])

    # Keep first occurrence per rsID in case of duplicates
    df = df.drop_duplicates(subset='RSID', keep='first')
    return df.set_index('RSID')['maf']


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(rsid_lookup, eur_sumstats, afr_sumstats, output_tsv):

    known_rsids = load_rsids(rsid_lookup)

    eur_maf = read_aaf(eur_sumstats, known_rsids, "EUR")
    afr_maf = read_aaf(afr_sumstats, known_rsids, "AFR")

    # Build output DataFrame with all rsIDs from the lookup as rows
    out_df = pd.DataFrame({'rsid': sorted(known_rsids)})
    out_df['EUR_maf'] = out_df['rsid'].map(eur_maf)
    out_df['AFR_maf'] = out_df['rsid'].map(afr_maf)

    n_eur = out_df['EUR_maf'].notna().sum()
    n_afr = out_df['AFR_maf'].notna().sum()
    print(f"NOTE: EUR MAF coverage: {n_eur:,} / {len(out_df):,} rsIDs")
    print(f"NOTE: AFR MAF coverage: {n_afr:,} / {len(out_df):,} rsIDs")

    out_dir = os.path.dirname(output_tsv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    out_df.to_csv(output_tsv, sep='\t', index=False, na_rep='NA')
    print(f"SUCCESS: Wrote {len(out_df):,} rows to {output_tsv}")


###############################################################################
################################ CALL MAIN ###################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv)
