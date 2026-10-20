#!/usr/bin/env python3

'''

decompose_beta_vs_z_ld_error.py

The purpose of this script is to split the noise-corrected T_eff (see
collate_effective_trait_counts.py) into a weight term and a correlation
term, for the beta, z, and normalized-z estimators:

  T_eff = tr(C)^2 / ||C||_F^2 = K / (1 + omega), with
  K     = tr(C)^2 / sum_t C_tt^2               (driven by the spread of
                                                 per-trait weights/variances)
  omega = sum_{t != u} C_tu^2 / sum_t C_tt^2   (driven by cross-trait
                                                 correlation)

Squares are replaced by products over the odd/even chromosome halves, so the
estimate is not inflated by sampling noise in the single covariance matrix.

This duplicates the load()/cov()/shape() helpers from
collate_effective_trait_counts.py rather than importing them, so the two
scripts can be run independently.

Input:  <in_npz> (mvp_trait_correlation_matrix.npz, from
        build_mvp_trait_correlation_matrix.py): beta, se, aaf, chrom.
Output: <out_dir>/mvp_teff_decomposition.tsv
          Columns: estimator, T, K, omega, teff

Usage:
  python decompose_beta_vs_z_ld_error.py [-i <in_npz>] [-o <out_dir>]

***Run in the betainfer conda environment on Aurora***

'''

import sys
import getopt
import os

import numpy as np
import pandas as pd

# ── Defaults ──────────────────────────────────────────────────────────────────
_EFF_DIR = ("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/"
            "mvp_results/effective_traits")
DEFAULT_IN_NPZ  = os.path.join(_EFF_DIR, "mvp_trait_correlation_matrix.npz")
DEFAULT_OUT_DIR = _EFF_DIR


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -i  Path to mvp_trait_correlation_matrix.npz
        (default: mvp_results/effective_traits/mvp_trait_correlation_matrix.npz)
    -o  Output directory for mvp_teff_decomposition.tsv
        (default: mvp_results/effective_traits/)
    -h  Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
##########  LOADING AND COVARIANCE (duplicated from collate_effective_trait_counts.py) ####
###############################################################################

def load(in_npz):
    """z = BETA / SE; phenotypic scale s_t = median(SE * sqrt(2pq))."""
    d = np.load(in_npz)
    B, S = d["beta"].astype(np.float64), d["se"].astype(np.float64)
    pq = 2 * d["aaf"] * (1 - d["aaf"])
    s = np.nanmedian(S * np.sqrt(pq)[:, None], axis=0)
    return B / S, s, d["chrom"]


def cov(Z):
    """Trait x trait covariance across variants; missing values handled pairwise."""
    M = np.isfinite(Z)
    X = np.where(M, Z - np.nanmean(Z, axis=0), 0.0)
    m = M.astype(np.float64)
    return (X.T @ X) / np.maximum(m.T @ m - 1, 1)


def shape(C, kind, s):
    """Covariance of the vector a given Pearson estimator actually correlates."""
    w = {"z_pearson": np.ones(len(C)), "z_norm_pearson": 1 / np.sqrt(np.diag(C)), "beta_pearson": s}[kind]
    H = np.eye(len(C)) - 1.0 / len(C)
    return H @ (C * np.outer(w, w)) @ H


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(in_npz, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    Z, s, chrom = load(in_npz)
    odd = np.isin(chrom, np.arange(1, 23, 2))
    C, A, B = cov(Z), cov(Z[odd]), cov(Z[~odd])
    rows = []
    for k in ["beta_pearson", "z_pearson", "z_norm_pearson"]:
        F, Fa, Fb = (shape(M, k, s) for M in (C, A, B))
        dg = (np.diag(Fa) * np.diag(Fb)).sum()
        off = (Fa * Fb).sum() - dg  # off-diagonal contribution: cross-trait correlation
        tr2 = np.trace(F) ** 2
        rows.append({"estimator": k, "T": len(C), "K": tr2 / dg, "omega": off / dg, "teff": tr2 / (dg + off)})

    t = pd.DataFrame(rows)
    out_path = os.path.join(out_dir, "mvp_teff_decomposition.tsv")
    t.to_csv(out_path, sep="\t", index=False, float_format="%.5g")
    print(t.to_string(index=False))
    print(f"NOTE: Written to {out_path}")


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, _ = getopt.getopt(argv[1:], "i:o:h")
    except getopt.GetoptError as e:
        print(f"ERROR: {e}")
        help(1)

    options_dict = dict(opts)
    if '-h' in options_dict:
        help()

    in_npz  = options_dict.get('-i', DEFAULT_IN_NPZ)
    out_dir = options_dict.get('-o', DEFAULT_OUT_DIR)

    if not os.path.exists(in_npz):
        print(f"ERROR: input npz not found: {in_npz}")
        sys.exit(1)

    print("Acceptable Inputs Given")
    driver(in_npz, out_dir)


if __name__ == "__main__":
    main(sys.argv)
