#!/usr/bin/env python3

'''

build_mvp_trait_correlation_matrix.py

The purpose of this script is to pull a genome-wide random sample of common
variants and their effect size (BETA) and standard error (SE) from every MVP
EUR trait, so that downstream scripts can estimate how correlated the traits
are and how many effectively-independent traits betanfer's LD estimate rests
on (Supplementary Note on beta- vs z-based LD inference).

Input:  MVP_hg19/EUR/*.sumstats.hg19.txt.gz (tabix-indexed, read-only).
Output: <out_dir>/mvp_trait_correlation_matrix.npz
          beta, se  - (n_variants, n_traits) float32 arrays (NaN where a trait
                      has no record for that variant)
          n         - (n_traits,) per-trait median sample size
          trait     - (n_traits,) trait names
          id        - (n_variants,) variant IDs (as given in the sumstats file)
          chrom     - (n_variants,) chromosome
          pos       - (n_variants,) position
          aaf       - (n_variants,) alternate allele frequency in the anchor trait

Variants are drawn from -w random windows of width -b, spread over the
autosomes in proportion to chromosome length, with the MHC region excluded.
A variant is kept if it is common (0.05 <= AAF <= 0.95) in the first trait
file, which defines the anchor variant set that every other trait is aligned
to.

Usage:
  python build_mvp_trait_correlation_matrix.py [-d <mvp_eur_sumstats_dir>] [-o <out_dir>]
                                                [-w <n_windows>] [-b <window_bp>]
                                                [-s <seed>] [-n <n_workers>]

***Run in the betainfer conda environment on Aurora (needs pysam, which is
   not installed in the finemap / finemap_old environments)***

'''

import sys
import getopt
import glob
import os
import time
from multiprocessing import get_context

import numpy as np
import pysam

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_MVP_DIR  = "/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19/EUR"
DEFAULT_OUT_DIR  = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/effective_traits"
DEFAULT_N_WINDOW = 4000
DEFAULT_WINDOW_BP = 5000
DEFAULT_SEED      = 20261001
DEFAULT_N_WORKERS = 48
MHC = (6, 25_000_000, 34_000_000)
# hg19 autosome lengths, chr1..chr22
HG19_CHROM_LEN = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067,
                  159138663, 146364022, 141213431, 135534747, 135006516, 133851895,
                  115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
                  59128983, 63025520, 48129895, 51304566]


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -d  MVP EUR sumstats directory containing *.sumstats.hg19.txt.gz
        (default: /flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19/EUR)
    -o  Output directory for mvp_trait_correlation_matrix.npz
        (default: mvp_results/effective_traits/)
    -w  Number of random windows to sample (default: 4000)
    -b  Window width in bp (default: 5000)
    -s  Random seed (default: 20261001)
    -n  Number of parallel worker processes (default: 48)
    -h  Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
##########  WINDOWS AND ANCHOR VARIANT SET  ###################################
###############################################################################

def sample_windows(n_windows, window_bp, seed):
    """Spread windows over the autosomes in proportion to length, excluding the MHC."""
    rng = np.random.default_rng(seed)
    lengths = np.array(HG19_CHROM_LEN, dtype=float)
    chrs = rng.choice(np.arange(1, 23), size=n_windows * 2, p=lengths / lengths.sum())
    out = set()
    for c in chrs:
        start = int(rng.integers(1_000_000, HG19_CHROM_LEN[c - 1] - 1_000_000)) // window_bp * window_bp
        if not (c == MHC[0] and MHC[1] - window_bp <= start <= MHC[2]):
            out.add((int(c), start))
        if len(out) == n_windows:
            break
    return sorted(out)


def anchor_variants(anchor_trait_path, windows, window_bp):
    """
    Define the common-variant set from the first trait file (IDs include the
    alleles, so multi-allelic records are naturally distinct; an indel that
    spans two adjacent windows is kept only once).
    """
    ids, chrom, pos, aaf = [], [], [], []
    seen = set()
    with pysam.TabixFile(anchor_trait_path) as tf:
        for c, s in windows:
            for rec in tf.fetch(str(c), s, s + window_bp):
                f = rec.split("\t")
                try:
                    a = float(f[13])
                except ValueError:
                    continue
                if 0.05 <= a <= 0.95 and f[0] not in seen:
                    seen.add(f[0])
                    ids.append(f[0])
                    chrom.append(c)
                    pos.append(int(f[3]))
                    aaf.append(a)
    return ids, np.array(chrom), np.array(pos), np.array(aaf)


###############################################################################
##########  PER-TRAIT WORKER  #################################################
###############################################################################

def pull_trait(args):
    """
    Parameters: (path, windows, window_bp, index) where index maps a variant
    ID to its row in the anchor variant set.
    Returns: (beta, se, median_n) aligned to the anchor variant set, with NaN
    where this trait has no record for that variant.
    """
    path, windows, window_bp, index = args
    beta = np.full(len(index), np.nan, dtype=np.float32)
    se = np.full(len(index), np.nan, dtype=np.float32)
    ns = []
    with pysam.TabixFile(path) as tf:
        for c, s in windows:
            for rec in tf.fetch(str(c), s, s + window_bp):
                f = rec.split("\t")
                i = index.get(f[0])
                if i is None:
                    continue
                try:
                    beta[i], se[i] = float(f[6]), float(f[7])
                    ns.append(float(f[17]))
                except ValueError:
                    pass
    return beta, se, float(np.median(ns)) if ns else np.nan


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(mvp_dir, out_dir, n_windows, window_bp, seed, n_workers):
    t0 = time.time()
    traits = sorted(f for f in glob.glob(f"{mvp_dir}/*.sumstats.hg19.txt.gz") if "munged" not in f)
    if not traits:
        print(f"ERROR: No trait sumstats files found in {mvp_dir}")
        sys.exit(1)

    windows = sample_windows(n_windows, window_bp, seed)
    ids, chrom, pos, aaf = anchor_variants(traits[0], windows, window_bp)
    index = {v: i for i, v in enumerate(ids)}
    print(f"NOTE: {len(windows)} windows, {len(ids)} anchor variants, "
          f"{len(traits)} traits ({time.time() - t0:.0f}s)", flush=True)

    beta = np.empty((len(ids), len(traits)), dtype=np.float32)
    se = np.empty_like(beta)
    n = np.empty(len(traits))
    work = [(path, windows, window_bp, index) for path in traits]
    with get_context("fork").Pool(n_workers) as pool:
        for j, (b, s, med_n) in enumerate(pool.imap(pull_trait, work, chunksize=2)):
            beta[:, j], se[:, j], n[j] = b, s, med_n
            if (j + 1) % 100 == 0:
                print(f"NOTE: {j + 1}/{len(traits)} traits ({time.time() - t0:.0f}s)", flush=True)

    names = np.array([os.path.basename(p).split(".")[2] for p in traits])
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "mvp_trait_correlation_matrix.npz")
    np.savez(out_path, beta=beta, se=se, n=n, trait=names, id=np.array(ids),
             chrom=chrom, pos=pos, aaf=aaf)
    print(f"NOTE: done: {beta.shape}, missing {np.isnan(beta).mean():.4f} ({time.time() - t0:.0f}s)")
    print(f"NOTE: Written to {out_path}")


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, _ = getopt.getopt(argv[1:], "d:o:w:b:s:n:h")
    except getopt.GetoptError as e:
        print(f"ERROR: {e}")
        help(1)

    options_dict = dict(opts)
    if '-h' in options_dict:
        help()

    mvp_dir   = options_dict.get('-d', DEFAULT_MVP_DIR)
    out_dir   = options_dict.get('-o', DEFAULT_OUT_DIR)
    n_windows = int(options_dict.get('-w', DEFAULT_N_WINDOW))
    window_bp = int(options_dict.get('-b', DEFAULT_WINDOW_BP))
    seed      = int(options_dict.get('-s', DEFAULT_SEED))
    n_workers = int(options_dict.get('-n', DEFAULT_N_WORKERS))

    if not os.path.exists(mvp_dir):
        print(f"ERROR: MVP sumstats directory not found: {mvp_dir}")
        sys.exit(1)

    print("Acceptable Inputs Given")
    driver(mvp_dir, out_dir, n_windows, window_bp, seed, n_workers)


if __name__ == "__main__":
    main(sys.argv)
