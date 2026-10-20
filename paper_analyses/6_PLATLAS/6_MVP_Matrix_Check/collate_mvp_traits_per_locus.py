#!/usr/bin/env python3

'''

collate_mvp_traits_per_locus.py

The purpose of this script is to count, for each sampled locus and
population, how many of the MVP traits betanfer actually used as independent
(non-significant) evidence: a trait is excluded at a locus if it has any
record with P < 5e-8 inside the locus region.

Input:  sampled_loci_750kb.csv (CHR, START, END; from sample_loci.py).
        MVP_hg19 {AFR,EUR} sumstats, tabix-indexed (read-only).
Output: <out_dir>/mvp_traits_per_locus.tsv
          Columns: locus, population, n_total, n_gws, n_traits
          n_total  - number of trait sumstats files for the population
          n_gws    - traits with at least one P < 5e-8 record in the locus
          n_traits - n_total - n_gws (traits betanfer kept at this locus)

Usage:
  python collate_mvp_traits_per_locus.py [-l <locus_csv>] [-d <mvp_sumstats_dir>]
                                         [-a <populations>] [-o <out_dir>] [-n <n_workers>]

***Run in the betainfer conda environment on Aurora (needs pysam, which is
   not installed in the finemap / finemap_old environments)***

'''

import sys
import getopt
import glob
import os
from multiprocessing import Pool

import pandas as pd
import pysam

# ── Defaults ──────────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_RITM_DIR = ("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/"
             "mvp_results/RITM0322545")
DEFAULT_LOCUS_CSV   = os.path.join(_SCRIPT_DIR, "sampled_loci_750kb.csv")
DEFAULT_MVP_DIR     = "/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19"
DEFAULT_POPULATIONS = "AFR,EUR"
DEFAULT_OUT_DIR     = _RITM_DIR
DEFAULT_N_WORKERS   = 48
GWS_THRESHOLD       = 5e-8
COL_P               = 8   # 0-based column index of the P-value in MVP sumstats rows


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -l  Sampled locus CSV with columns CHR,START,END
        (default: sampled_loci_750kb.csv next to this script)
    -d  MVP sumstats base directory containing {population}/*.sumstats.hg19.txt.gz
        (default: /flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19)
    -a  Populations, comma-separated (default: AFR,EUR)
    -o  Output directory for mvp_traits_per_locus.tsv
        (default: RITM0322545/)
    -n  Number of parallel worker processes (default: 48)
    -h  Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
##########  PER-TRAIT WORKER  #################################################
###############################################################################

def trait_has_gws_hit(args):
    """
    Parameters: (path, regions) where regions is a list of (chrom, start, end).
    Returns: list[bool], one flag per region, True if that trait has any
    P < GWS_THRESHOLD record inside the region.
    """
    path, regions = args
    hits = []
    with pysam.TabixFile(path) as tf:
        for chrom, start, end in regions:
            flag = False
            for rec in tf.fetch(chrom, start - 1, end):
                try:
                    if float(rec.split("\t", COL_P + 1)[COL_P]) < GWS_THRESHOLD:
                        flag = True
                        break
                except ValueError:
                    pass
            hits.append(flag)
    return hits


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(locus_csv, mvp_dir, populations, out_dir, n_workers):
    os.makedirs(out_dir, exist_ok=True)

    loci = pd.read_csv(locus_csv)
    labels = [f"chr{int(r.CHR)}_{int(r.START)}_{int(r.END)}" for _, r in loci.iterrows()]
    regions = [(str(int(r.CHR)), int(r.START), int(r.END)) for _, r in loci.iterrows()]
    print(f"NOTE: {len(labels)} loci loaded from {locus_csv}")

    rows = []
    for population in populations:
        traits = sorted(f for f in glob.glob(f"{mvp_dir}/{population}/*.sumstats.hg19.txt.gz")
                        if "munged" not in f)
        with Pool(n_workers) as pool:
            hits = pool.map(trait_has_gws_hit, [(t, regions) for t in traits], chunksize=2)
        n_gws = pd.DataFrame(hits, columns=labels).sum()
        rows += [dict(locus=label, population=population, n_total=len(traits),
                      n_gws=int(n_gws[label]), n_traits=len(traits) - int(n_gws[label]))
                 for label in labels]
        print(f"NOTE: {population}: {len(traits)} traits done", flush=True)

    out_path = os.path.join(out_dir, "mvp_traits_per_locus.tsv")
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
    print(f"NOTE: Written to {out_path}")


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, _ = getopt.getopt(argv[1:], "l:d:a:o:n:h")
    except getopt.GetoptError as e:
        print(f"ERROR: {e}")
        help(1)

    options_dict = dict(opts)
    if '-h' in options_dict:
        help()

    locus_csv   = options_dict.get('-l', DEFAULT_LOCUS_CSV)
    mvp_dir     = options_dict.get('-d', DEFAULT_MVP_DIR)
    populations = [p.strip().upper() for p in options_dict.get('-a', DEFAULT_POPULATIONS).split(',')]
    out_dir     = options_dict.get('-o', DEFAULT_OUT_DIR)
    n_workers   = int(options_dict.get('-n', DEFAULT_N_WORKERS))

    for path, label in [(locus_csv, "locus CSV"), (mvp_dir, "MVP sumstats directory")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    print("Acceptable Inputs Given")
    driver(locus_csv, mvp_dir, populations, out_dir, n_workers)


if __name__ == "__main__":
    main(sys.argv)
