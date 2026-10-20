#!/usr/bin/env python3

'''

extract_mvp_ld_error_sample.py

The purpose of this script is to make one streaming pass over the large
betanfer-vs-MVP LD comparison table and produce the two inputs needed for a
scatter plot of betanfer beta correlation against real MVP genotype
correlation: a random sample of pairs to plot, and the exact Pearson
correlation (and other difference moments) computed over every pair.

Input:  mvp_matrix_comparison.tsv (locus, population, variant_1, variant_2,
        betanfer_r, mvp_r, delta_r), produced by compare_matrices.py.
Output: <out_dir>/mvp_ld_error_sample.tsv.gz
          Columns: locus, population, betanfer_r, mvp_r, delta_r
          A fixed-size random sample of pairs (sample size set by -n).
        <out_dir>/mvp_ld_error_moments.tsv
          Columns: population, n_pairs, pearson_r, mean_delta, mean_abs_delta, rmse_delta
          One row per population plus a row for "ALL", computed over every
          pair in the input (not just the sample).

Sampling: each row is kept independently with probability -k (Bernoulli), and
the output sample is then an exact random draw of -n rows from whatever was
kept, so the sample size does not depend on the exact row count of the input.

Usage:
  python extract_mvp_ld_error_sample.py [-c <comparison_tsv>] [-o <out_dir>]
                                        [-n <sample_size>] [-k <keep_prob>] [-s <seed>]

***Run in Conda finemap_old environment on Aurora***

'''

import sys
import getopt
import os
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pv

# ── Defaults ──────────────────────────────────────────────────────────────────
_RITM_DIR = ("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/"
             "mvp_results/RITM0322545")
DEFAULT_COMPARISON  = os.path.join(_RITM_DIR, "mvp_matrix_comparison.tsv")
DEFAULT_OUT_DIR      = _RITM_DIR
DEFAULT_SAMPLE_SIZE  = 100_000
DEFAULT_KEEP_PROB    = 1.7e-4
DEFAULT_SEED         = 20261002
BLOCK_SIZE           = 1 << 27
COLS                 = ["locus", "population", "betanfer_r", "mvp_r", "delta_r"]


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -c  Path to mvp_matrix_comparison.tsv
        (default: RITM0322545/mvp_matrix_comparison.tsv)
    -o  Output directory for the sample and moments tables
        (default: RITM0322545/)
    -n  Exact number of pairs to keep in the output sample (default: 100000)
    -k  Per-row Bernoulli keep probability used before the exact draw
        (default: 1.7e-4; must be large enough that at least -n rows survive)
    -s  Random seed (default: 20261002)
    -h  Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
##########  STREAMING READER  #################################################
###############################################################################

def batches(tsv_path):
    """Yield pyarrow record batches of the comparison TSV, five columns only."""
    ro = pv.ReadOptions(block_size=BLOCK_SIZE, use_threads=True)
    po = pv.ParseOptions(delimiter="\t")
    co = pv.ConvertOptions(include_columns=COLS, column_types={
        "locus": pa.string(), "population": pa.string(),
        "betanfer_r": pa.float64(), "mvp_r": pa.float64(), "delta_r": pa.float64()})
    with pv.open_csv(tsv_path, read_options=ro, parse_options=po, convert_options=co) as reader:
        for batch in reader:
            yield batch


###############################################################################
##########  MOMENT ACCUMULATOR  ###############################################
###############################################################################

class MomentAcc:
    """
    Per-population running sums needed for an exact Pearson r between mvp_r
    and betanfer_r, and for the mean / mean-abs / RMSE of delta_r, without
    holding all pairs in memory.
    """

    def __init__(self):
        self.sums = {}   # population -> [n, Sx, Sy, Sxx, Syy, Sxy, Sd, S|d|, Sdd]

    def add(self, population, betanfer_r, mvp_r, delta_r):
        ok = np.isfinite(betanfer_r) & np.isfinite(mvp_r) & np.isfinite(delta_r)
        b, m, d = betanfer_r[ok], mvp_r[ok], delta_r[ok]
        a = self.sums.setdefault(population, np.zeros(9))
        a += [len(b), m.sum(), b.sum(), (m * m).sum(), (b * b).sum(), (m * b).sum(),
              d.sum(), np.abs(d).sum(), (d * d).sum()]

    def table(self):
        """Parameters: none. Returns: DataFrame, one row per population plus 'ALL'."""
        rows = []
        total = sum(self.sums.values())
        for population, a in sorted(self.sums.items()) + [("ALL", total)]:
            n, sx, sy, sxx, syy, sxy, sd, sad, sdd = a
            cov = sxy / n - (sx / n) * (sy / n)
            var_x = sxx / n - (sx / n) ** 2
            var_y = syy / n - (sy / n) ** 2
            r = cov / np.sqrt(var_x * var_y)
            rows.append(dict(population=population, n_pairs=int(n), pearson_r=r,
                              mean_delta=sd / n, mean_abs_delta=sad / n, rmse_delta=np.sqrt(sdd / n)))
        return pd.DataFrame(rows)


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(comparison, out_dir, sample_size, keep_prob, seed):
    os.makedirs(out_dir, exist_ok=True)
    rng = np.random.default_rng(seed)
    acc = MomentAcc()
    kept_chunks = []
    n_rows, t0 = 0, time.time()

    print(f"NOTE: Streaming {comparison}...", flush=True)
    for batch in batches(comparison):
        df = batch.to_pandas()
        kept_chunks.append(df[rng.random(len(df)) < keep_prob])
        for population, g in df.groupby("population", sort=False):
            acc.add(population, g["betanfer_r"].to_numpy(), g["mvp_r"].to_numpy(), g["delta_r"].to_numpy())
        n_rows += len(df)
        if n_rows // 50_000_000 != (n_rows - len(df)) // 50_000_000:
            print(f"NOTE: {n_rows / 1e6:,.0f}M rows, {time.time() - t0:.0f}s", flush=True)

    sample = pd.concat(kept_chunks, ignore_index=True)
    sample = sample[np.isfinite(sample[["betanfer_r", "mvp_r", "delta_r"]]).all(axis=1)]
    if len(sample) < sample_size:
        print(f"ERROR: Only {len(sample)} rows survived Bernoulli keep (-k {keep_prob}); "
              f"need at least {sample_size}. Raise -k and retry.")
        sys.exit(1)
    pick = np.sort(rng.choice(len(sample), sample_size, replace=False))

    sample_path = os.path.join(out_dir, "mvp_ld_error_sample.tsv.gz")
    moments_path = os.path.join(out_dir, "mvp_ld_error_moments.tsv")
    sample.iloc[pick].to_csv(sample_path, sep="\t", index=False)
    moments = acc.table()
    moments.to_csv(moments_path, sep="\t", index=False, float_format="%.6g")

    print(moments.to_string(index=False))
    print(f"NOTE: rows {n_rows:,}; Bernoulli kept {len(sample):,}; sample {sample_size:,}; "
          f"{time.time() - t0:.0f}s")
    print(f"NOTE: Sample written to {sample_path}")
    print(f"NOTE: Moments written to {moments_path}")


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, _ = getopt.getopt(argv[1:], "c:o:n:k:s:h")
    except getopt.GetoptError as e:
        print(f"ERROR: {e}")
        help(1)

    options_dict = dict(opts)
    if '-h' in options_dict:
        help()

    comparison  = options_dict.get('-c', DEFAULT_COMPARISON)
    out_dir     = options_dict.get('-o', DEFAULT_OUT_DIR)
    sample_size = int(options_dict.get('-n', DEFAULT_SAMPLE_SIZE))
    keep_prob   = float(options_dict.get('-k', DEFAULT_KEEP_PROB))
    seed        = int(options_dict.get('-s', DEFAULT_SEED))

    if not os.path.exists(comparison):
        print(f"ERROR: comparison TSV not found: {comparison}")
        sys.exit(1)

    print("Acceptable Inputs Given")
    driver(comparison, out_dir, sample_size, keep_prob, seed)


if __name__ == "__main__":
    main(sys.argv)
