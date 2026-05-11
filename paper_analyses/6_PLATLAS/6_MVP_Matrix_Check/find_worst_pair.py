'''

find_worst_pair.py

Stream through the mvp_matrix_comparison.tsv and identify the variant pair
with the largest absolute delta_r where both variants are common (MAF > threshold)
in the relevant population.

MAF is looked up from maf_lookup.tsv (built by build_maf_lookup.py from MVP hg19
sumstats).  For EUR population rows, EUR_maf is used; for AFR rows, AFR_maf.

Variant pairs with missing MAF for either variant are silently skipped.

Usage:
  python find_worst_pair.py [-c <comparison_tsv>] [-m <maf_lookup_tsv>]
                             [-t <maf_threshold>] [-o <output_tsv>]

Output columns (TSV):
  locus, population, variant_1, variant_2, betanfer_r, mvp_r, delta_r,
  maf_v1, maf_v2

***Run in Conda finemap_old environment on Aurora***

'''

import sys
import getopt
import os
import pandas as pd
import numpy as np

# ── Defaults ─────────────────────────────────────────────────────────────────
_RITM_DIR = ("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/"
             "mvp_results/RITM0322545")
DEFAULT_COMPARISON = os.path.join(_RITM_DIR, "mvp_matrix_comparison.tsv")
DEFAULT_MAF_LOOKUP = os.path.join(_RITM_DIR, "maf_lookup.tsv")
DEFAULT_THRESHOLD  = 0.05
CHUNK_SIZE         = 100_000


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -c  Path to mvp_matrix_comparison.tsv
        (default: RITM0322545/mvp_matrix_comparison.tsv)
    -m  Path to maf_lookup.tsv
        (default: RITM0322545/maf_lookup.tsv)
    -t  MAF threshold; both variants must exceed this value
        (default: 0.05)
    -o  Output TSV path (result row written here in addition to stdout)
    -h  Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


def main():
    comparison = DEFAULT_COMPARISON
    maf_lookup = DEFAULT_MAF_LOOKUP
    threshold  = DEFAULT_THRESHOLD
    output     = None

    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hc:m:t:o:")
    except getopt.GetoptError as e:
        print(f"ERROR: {e}")
        help(1)

    for opt, val in opts:
        if opt == "-h":
            help()
        elif opt == "-c":
            comparison = val
        elif opt == "-m":
            maf_lookup = val
        elif opt == "-t":
            threshold = float(val)
        elif opt == "-o":
            output = val

    # ── Load MAF lookup ───────────────────────────────────────────────────────
    print(f"NOTE: Loading MAF lookup from {maf_lookup}...", flush=True)
    maf_df = pd.read_csv(maf_lookup, sep="\t")
    # Build population-specific dicts for vectorised .map() calls
    maf_eur = maf_df.set_index("rsid")["EUR_maf"].to_dict()
    maf_afr = maf_df.set_index("rsid")["AFR_maf"].to_dict()
    print(f"NOTE: Loaded MAF for {len(maf_df):,} rsIDs.", flush=True)

    # ── Stream comparison TSV ─────────────────────────────────────────────────
    print(f"NOTE: Streaming {comparison}...", flush=True)
    best_row      = None
    best_abs_delta = -1.0
    chunk_num      = 0

    for chunk in pd.read_csv(comparison, sep="\t", chunksize=CHUNK_SIZE):
        chunk_num += 1
        if chunk_num % 100 == 0:
            print(f"NOTE: Processed {chunk_num * CHUNK_SIZE:,} rows...",
                  flush=True)

        eur_mask = chunk["population"] == "EUR"
        afr_mask = chunk["population"] == "AFR"

        # EUR sub-chunk
        if eur_mask.any():
            sub = chunk[eur_mask].copy()
            sub["maf_v1"] = sub["variant_1"].map(maf_eur)
            sub["maf_v2"] = sub["variant_2"].map(maf_eur)
            sub = sub.dropna(subset=["maf_v1", "maf_v2"])
            sub = sub[(sub["maf_v1"] > threshold) & (sub["maf_v2"] > threshold)]
            if len(sub):
                sub["abs_delta"] = sub["delta_r"].abs()
                cand = sub.loc[sub["abs_delta"].idxmax()]
                if cand["abs_delta"] > best_abs_delta:
                    best_abs_delta = cand["abs_delta"]
                    best_row = cand.copy()

        # AFR sub-chunk
        if afr_mask.any():
            sub = chunk[afr_mask].copy()
            sub["maf_v1"] = sub["variant_1"].map(maf_afr)
            sub["maf_v2"] = sub["variant_2"].map(maf_afr)
            sub = sub.dropna(subset=["maf_v1", "maf_v2"])
            sub = sub[(sub["maf_v1"] > threshold) & (sub["maf_v2"] > threshold)]
            if len(sub):
                sub["abs_delta"] = sub["delta_r"].abs()
                cand = sub.loc[sub["abs_delta"].idxmax()]
                if cand["abs_delta"] > best_abs_delta:
                    best_abs_delta = cand["abs_delta"]
                    best_row = cand.copy()

    if best_row is None:
        print("ERROR: No variant pairs passed the MAF filter.", flush=True)
        sys.exit(1)

    # ── Output ─────────────────────────────────────────────────────────────────
    cols = ["locus", "population", "variant_1", "variant_2",
            "betanfer_r", "mvp_r", "delta_r", "maf_v1", "maf_v2"]
    header = "\t".join(cols)
    row_vals = "\t".join([
        str(best_row["locus"]),
        str(best_row["population"]),
        str(best_row["variant_1"]),
        str(best_row["variant_2"]),
        str(best_row["betanfer_r"]),
        str(best_row["mvp_r"]),
        str(best_row["delta_r"]),
        f"{best_row['maf_v1']:.6f}",
        f"{best_row['maf_v2']:.6f}",
    ])

    print(f"\nNOTE: Worst pair (MAF > {threshold}, |delta_r| = {best_abs_delta:.6f}):")
    print(header)
    print(row_vals)

    if output:
        out_dir = os.path.dirname(os.path.abspath(output))
        os.makedirs(out_dir, exist_ok=True)
        with open(output, "w") as fh:
            fh.write(header + "\n")
            fh.write(row_vals + "\n")
        print(f"NOTE: Result written to {output}")


if __name__ == "__main__":
    main()
