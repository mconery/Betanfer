#!/usr/bin/env python3

'''

collate_effective_trait_counts.py

The purpose of this script is to estimate the effective number of
independent traits (T_eff) behind betanfer's effect-size (beta) correlation
estimator on the MVP data, and to compare it against the same estimator
applied to z-scores and to normalized z-scores. This quantifies how much of
the apparent cross-trait correlation structure that beta-based LD inference
relies on is real, versus how it would look under alternative scalings.

T_eff = tr(C)^2 / ||C||_F^2, where C is the trait x trait covariance of
whichever per-variant vector a given Pearson correlation estimator
implicitly correlates (beta, z, or normalized z). The naive ||C||_F^2 is
inflated by sampling noise in the estimated C; the noise-corrected version
replaces the square with the product of two independent halves of the data
(disjoint chromosome splits), which is unbiased because the two halves share
no variants.

Input:  <in_npz> (mvp_trait_correlation_matrix.npz, from
        build_mvp_trait_correlation_matrix.py): beta, se, aaf, chrom, trait, n.
Output: <out_dir>/mvp_teff_summary.tsv       - T_eff per estimator, overall
                                                and restricted to variants with
                                                no genome-wide-significant hit
        <out_dir>/mvp_teff_vs_t.tsv          - T_eff as a function of the
                                                number of traits included
        <out_dir>/mvp_trait_scale.tsv        - per-trait phenotypic scale and
                                                z-variance, for reference

Usage:
  python collate_effective_trait_counts.py [-i <in_npz>] [-o <out_dir>]

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
N_SPLITS, SEED  = 20, 20261001
GWS_Z           = 5.45    # |z| threshold used to flag a genome-wide-significant hit
SIZES           = [50, 100, 200, 250, 400, 600, 800, 1000, 1200]


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -i  Path to mvp_trait_correlation_matrix.npz
        (default: mvp_results/effective_traits/mvp_trait_correlation_matrix.npz)
    -o  Output directory for mvp_teff_summary.tsv, mvp_teff_vs_t.tsv,
        mvp_trait_scale.tsv (default: mvp_results/effective_traits/)
    -h  Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
##########  LOADING AND COVARIANCE  ###########################################
###############################################################################

def load(in_npz):
    """
    z = BETA / SE. The phenotypic scale s_t = median(SE * sqrt(2pq)) is the
    beta scale with the variance-of-genotype factor removed, leaving the
    phenotypic-variance factor alone.
    """
    d = np.load(in_npz)
    B, S = d["beta"].astype(np.float64), d["se"].astype(np.float64)
    pq = 2 * d["aaf"] * (1 - d["aaf"])
    s = np.nanmedian(S * np.sqrt(pq)[:, None], axis=0)
    return B / S, s, d["chrom"], d["trait"], d["n"]


def cov(Z):
    """Trait x trait covariance across variants; missing values handled pairwise."""
    M = np.isfinite(Z)
    X = np.where(M, Z - np.nanmean(Z, axis=0), 0.0)
    m = M.astype(np.float64)
    return (X.T @ X) / np.maximum(m.T @ m - 1, 1)


###############################################################################
##########  PER-ESTIMATOR COVARIANCE SHAPE  ###################################
###############################################################################

def shape(C, kind, s):
    """
    The covariance of the vector a given Pearson estimator actually
    correlates: weight each trait by w_t for that estimator, then
    double-center across traits (H).
    """
    w = {"z_pearson": np.ones(len(C)), "z_norm_pearson": 1 / np.sqrt(np.diag(C)), "beta_pearson": s}[kind]
    H = np.eye(len(C)) - 1.0 / len(C)
    return H @ (C * np.outer(w, w)) @ H


def teff(C_full, C_A, C_B):
    """Returns (raw, noise-corrected) T_eff."""
    tr2 = np.trace(C_full) ** 2
    return tr2 / (C_full ** 2).sum(), tr2 / (C_A * C_B).sum()


def kish(C, s):
    """T_eff under power differences alone, assuming independent traits (the beta Kish ratio)."""
    d = np.diag(C) * s ** 2
    return d.sum() ** 2 / (d ** 2).sum()


###############################################################################
##########  CHROMOSOME-HALF SPLITS  ############################################
###############################################################################

def splits(chrom):
    """Split 0 is the odd/even chromosome split; the rest are random with roughly balanced variant counts."""
    rng = np.random.default_rng(SEED)
    out = [np.isin(chrom, np.arange(1, 23, 2))]
    counts = np.bincount(chrom, minlength=23)[1:]
    while len(out) < N_SPLITS:
        pick = rng.permutation(22)[: rng.integers(9, 14)] + 1
        frac = counts[pick - 1].sum() / counts.sum()
        if 0.45 <= frac <= 0.55:
            out.append(np.isin(chrom, pick))
    return out


def summarize(Z, s, chrom, label):
    C = cov(Z)
    halves = [(cov(Z[a]), cov(Z[~a])) for a in splits(chrom)]
    rows = []
    for kind in ["z_pearson", "z_norm_pearson", "beta_pearson"]:
        F = shape(C, kind, s)
        est = [teff(F, shape(A, kind, s), shape(Bh, kind, s)) for A, Bh in halves]
        corr = np.array([e[1] for e in est])
        rows.append({"variants": label, "n_variants": len(Z), "estimator": kind, "T": len(C),
                     "raw": est[0][0], "corrected_parity": corr[0], "corrected_mean": corr.mean(),
                     "corrected_sd": corr.std(ddof=1)})
    rows.append({"variants": label, "n_variants": len(Z), "estimator": "beta_kish_only", "T": len(C),
                 "raw": kish(C, s), "corrected_parity": np.nan, "corrected_mean": np.nan,
                 "corrected_sd": np.nan})
    return rows, C, halves[0]


###############################################################################
##########  T_eff AS A FUNCTION OF THE NUMBER OF TRAITS  #######################
###############################################################################

def curve(C, half, s, names):
    rng = np.random.default_rng(SEED)
    rows = []
    for T in SIZES + [len(C)]:
        for rep in range(1 if T == len(C) else 10):
            k = np.sort(rng.choice(len(C), T, replace=False)) if T < len(C) else np.arange(T)
            sub = lambda M: M[np.ix_(k, k)]
            for kind in ["z_pearson", "z_norm_pearson", "beta_pearson"]:
                raw, cor = teff(shape(sub(C), kind, s[k]), shape(sub(half[0]), kind, s[k]),
                                shape(sub(half[1]), kind, s[k]))
                rows.append({"T": T, "rep": rep, "estimator": kind, "raw": raw, "corrected": cor})
    return pd.DataFrame(rows)


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(in_npz, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    Z, s, chrom, names, n = load(in_npz)
    rows, C, half = summarize(Z, s, chrom, "all_common")
    null = ~(np.abs(np.nan_to_num(Z)) > GWS_Z).any(axis=1)
    rows += summarize(Z[null], s, chrom[null], "no_gws_variant")[0]
    summ = pd.DataFrame(rows)
    summ_path = os.path.join(out_dir, "mvp_teff_summary.tsv")
    summ.to_csv(summ_path, sep="\t", index=False, float_format="%.4g")
    print(summ.to_string(index=False))

    cur = curve(C, half, s, names)
    cur_path = os.path.join(out_dir, "mvp_teff_vs_t.tsv")
    cur.to_csv(cur_path, sep="\t", index=False, float_format="%.4g")
    print(cur.groupby(["estimator", "T"])[["raw", "corrected"]].mean().round(1).to_string())

    d = np.sqrt(np.diag(C))
    R = C / np.outer(d, d)
    off = R[~np.eye(len(R), dtype=bool)]
    ev = np.linalg.eigvalsh(R)[::-1]
    print(f"NOTE: mean r^2 off-diag (raw) {np.mean(off ** 2):.4f}; top-5 eigen share "
          f"{ev[:5].sum() / ev.sum():.3f}; top eigen {np.round(ev[:5], 1)}")
    scale_path = os.path.join(out_dir, "mvp_trait_scale.tsv")
    pd.DataFrame({"trait": names, "scale_s": s, "median_n": n, "var_z": np.diag(C)}).to_csv(
        scale_path, sep="\t", index=False, float_format="%.4g")

    print(f"NOTE: Written to {summ_path}, {cur_path}, {scale_path}")


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
