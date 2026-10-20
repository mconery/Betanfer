#!/usr/bin/env python3

'''

extract_mvp_ld_error_by_locus_bin.py

The purpose of this script is to make one streaming pass over the large
betanfer-vs-MVP LD comparison table and compute the mean absolute error of
betanfer LD against real MVP genotype LD, per locus, population and
|genotype correlation| bin, over every pair (not a sample).

Input:  mvp_matrix_comparison.tsv (locus, population, variant_1, variant_2,
        betanfer_r, mvp_r, delta_r), produced by compare_matrices.py.
Output: <out_dir>/mvp_ld_error_by_locus_bin.tsv
          Columns: locus, population, bin, n_pairs, mae
          bin is one of the five |mvp_r| bins of width 0.2: 0.0-0.2, 0.2-0.4,
          0.4-0.6, 0.6-0.8, 0.8-1.0. mae is mean(|delta_r|) within the bin.

Usage:
  python extract_mvp_ld_error_by_locus_bin.py [-c <comparison_tsv>] [-o <out_dir>]

***Run in Conda finemap_old environment on Aurora***

'''

import sys
import getopt
import os

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pv

# ── Defaults ──────────────────────────────────────────────────────────────────
_RITM_DIR = ("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/"
             "mvp_results/RITM0322545")
DEFAULT_COMPARISON = os.path.join(_RITM_DIR, "mvp_matrix_comparison.tsv")
DEFAULT_OUT_DIR     = _RITM_DIR
BLOCK_SIZE          = 1 << 27
BIN_EDGES           = np.array([0.2, 0.4, 0.6, 0.8])   # |mvp_r| bins of width 0.2
BIN_LABELS          = ["0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"]
COL_TYPES           = {"locus": pa.string(), "population": pa.string(),
                       "mvp_r": pa.float64(), "delta_r": pa.float64()}


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -c  Path to mvp_matrix_comparison.tsv
        (default: RITM0322545/mvp_matrix_comparison.tsv)
    -o  Output directory for mvp_ld_error_by_locus_bin.tsv
        (default: RITM0322545/)
    -h  Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(comparison, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    ro = pv.ReadOptions(block_size=BLOCK_SIZE, use_threads=True)
    po = pv.ParseOptions(delimiter="\t")
    co = pv.ConvertOptions(include_columns=list(COL_TYPES), column_types=COL_TYPES)

    print(f"NOTE: Streaming {comparison}...", flush=True)
    parts, n_rows = [], 0
    with pv.open_csv(comparison, read_options=ro, parse_options=po, convert_options=co) as reader:
        for batch in reader:
            df = batch.to_pandas()
            n_rows += len(df)
            df = df[np.isfinite(df.mvp_r) & np.isfinite(df.delta_r)]
            # Round before binning so values that land exactly on an edge
            # (e.g. 0.2) are not misassigned by floating-point error.
            df["bin"] = np.searchsorted(BIN_EDGES, np.round(df.mvp_r.abs(), 6), side="right")
            df["abs_delta"] = df.delta_r.abs()
            parts.append(df.groupby(["locus", "population", "bin"]).abs_delta.agg(["size", "sum"]))

    total = pd.concat(parts).groupby(level=[0, 1, 2]).sum().reset_index()
    total["mae"] = total["sum"] / total["size"]
    total["bin"] = total.bin.map(dict(enumerate(BIN_LABELS)))
    out = total.rename(columns={"size": "n_pairs"})[["locus", "population", "bin", "n_pairs", "mae"]]

    out_path = os.path.join(out_dir, "mvp_ld_error_by_locus_bin.tsv")
    out.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    print(f"NOTE: rows {n_rows:,}; kept {int(total['size'].sum()):,}; groups {len(out)}")
    print(f"NOTE: Written to {out_path}")


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, _ = getopt.getopt(argv[1:], "c:o:h")
    except getopt.GetoptError as e:
        print(f"ERROR: {e}")
        help(1)

    options_dict = dict(opts)
    if '-h' in options_dict:
        help()

    comparison = options_dict.get('-c', DEFAULT_COMPARISON)
    out_dir    = options_dict.get('-o', DEFAULT_OUT_DIR)

    if not os.path.exists(comparison):
        print(f"ERROR: comparison TSV not found: {comparison}")
        sys.exit(1)

    print("Acceptable Inputs Given")
    driver(comparison, out_dir)


if __name__ == "__main__":
    main(sys.argv)
