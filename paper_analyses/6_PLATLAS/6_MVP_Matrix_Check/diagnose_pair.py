'''

diagnose_pair.py

Given a specific variant pair identified by find_worst_pair.py, extract their
betas from all sumstats files used for the Betanfer LD matrix and report
trait-coverage statistics.  This tests the hypothesis that at least one variant
was tested for relatively few traits, producing a noisy beta-correlation estimate.

The same MVP hg19 sumstats file lists used to generate the LD matrices are used
here:
  EUR population -> MVP_EUR.MVP.EUR.txt  (hg38-lifted, RSID+CHR:POS IDs)
  AFR population -> MVP_AFR.MVP.AFR.txt

Reported statistics:
  N_total_traits      - total sumstats files in the file list
  N_v1_nonmissing     - files where variant_1 was found with a non-NaN beta
  N_v2_nonmissing     - files where variant_2 was found with a non-NaN beta
  N_joint_nonmissing  - files where BOTH were found (what Betanfer correlated over)
  pearson_r           - Pearson r recomputed from joint betas (should match betanfer_r)

Usage:
  python diagnose_pair.py -1 <rsid1> -2 <rsid2> -l <locus> -p <population>
                           [-f <file_list>] [-n <num_processes>] [-o <output_tsv>]

  <locus> format: chr{CHR}_{START}_{END}  (e.g. chr1_26750001_27500000)

Optional output TSV (columns): trait, v1_beta, v2_beta, both_nonmissing

***Run in Conda finemap_old environment on Aurora***

'''

import sys
import getopt
import os
import math
import pandas as pd
import numpy as np

# ── Make betanfer importable from repo root ────────────────────────────────────
_REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from betanfer.utils import autodetect_variant_beta_cols, extract_betas_from_range

# ── Defaults ─────────────────────────────────────────────────────────────────
_LD_BASE     = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices"
_FILE_LISTS  = {
    "EUR": os.path.join(_LD_BASE, "file_lists", "MVP_EUR.MVP.EUR.txt"),
    "AFR": os.path.join(_LD_BASE, "file_lists", "MVP_AFR.MVP.AFR.txt"),
}
DEFAULT_NUM_PROCESSES = 12


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -1  rsID of variant_1 (required)
    -2  rsID of variant_2 (required)
    -l  Locus string in format chr{CHR}_{START}_{END} (required)
    -p  Population: EUR or AFR (required)
    -f  Path to sumstats file list
        (default: auto-resolved from population)
    -n  Number of parallel worker threads for beta extraction
        (default: 12)
    -o  Output TSV path for per-trait beta table (optional)
    -h  Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


def parse_locus(locus_str):
    """Parse 'chr1_26750001_27500000' -> ('1', 26750001, 27500000)."""
    parts = locus_str.split("_")
    chrom = parts[0].replace("chr", "")
    start = int(parts[1])
    end   = int(parts[2])
    return chrom, start, end


def pearson_r(x, y):
    """Pearson r from two equal-length 1-D arrays (no NaNs expected)."""
    n = len(x)
    if n < 2:
        return float("nan")
    xm = x - x.mean()
    ym = y - y.mean()
    denom = math.sqrt((xm**2).sum() * (ym**2).sum())
    if denom == 0:
        return float("nan")
    return float((xm * ym).sum() / denom)


def main():
    variant1      = None
    variant2      = None
    locus         = None
    population    = None
    file_list_path = None
    num_processes = DEFAULT_NUM_PROCESSES
    output        = None

    try:
        opts, _ = getopt.getopt(sys.argv[1:], "h1:2:l:p:f:n:o:")
    except getopt.GetoptError as e:
        print(f"ERROR: {e}")
        help(1)

    for opt, val in opts:
        if opt == "-h":
            help()
        elif opt == "-1":
            variant1 = val
        elif opt == "-2":
            variant2 = val
        elif opt == "-l":
            locus = val
        elif opt == "-p":
            population = val
        elif opt == "-f":
            file_list_path = val
        elif opt == "-n":
            num_processes = int(val)
        elif opt == "-o":
            output = val

    for flag, name in [("-1", variant1), ("-2", variant2),
                       ("-l", locus), ("-p", population)]:
        if name is None:
            print(f"ERROR: {flag} is required.")
            help(1)

    if population not in _FILE_LISTS and file_list_path is None:
        print(f"ERROR: Population must be EUR or AFR, or supply -f.")
        help(1)

    if file_list_path is None:
        file_list_path = _FILE_LISTS[population]

    # ── Load file list ─────────────────────────────────────────────────────────
    with open(file_list_path) as fh:
        file_list = [line.strip() for line in fh if line.strip()]
    n_total = len(file_list)
    print(f"NOTE: {n_total} files in {file_list_path}", flush=True)

    # ── Auto-detect column names from first file ───────────────────────────────
    first_file = file_list[0]
    variant_col, beta_col, _ = autodetect_variant_beta_cols(first_file)
    if variant_col is None or beta_col is None:
        print(f"ERROR: Could not auto-detect variant/beta columns from {first_file}")
        sys.exit(1)
    print(f"NOTE: variant_col='{variant_col}'  beta_col='{beta_col}'", flush=True)

    # ── Parse locus ────────────────────────────────────────────────────────────
    chrom, start, end = parse_locus(locus)
    print(f"NOTE: Locus {locus} -> chrom={chrom} start={start} end={end}", flush=True)

    # ── Extract betas for the full locus ─────────────────────────────────────
    print(f"NOTE: Extracting betas across {n_total} files (this may take a few minutes)...",
          flush=True)
    beta_df = extract_betas_from_range(
        file_list=file_list,
        variant_col=variant_col,
        beta_col=beta_col,
        chrom=chrom,
        start=start,
        end=end,
        num_processes=num_processes,
    )

    if beta_df.empty:
        print(f"ERROR: No variants extracted for locus {locus}.")
        sys.exit(1)

    print(f"NOTE: Beta matrix shape: {beta_df.shape} (variants x traits)", flush=True)

    # ── Slice to the two target variants ──────────────────────────────────────
    missing = []
    for v in [variant1, variant2]:
        if v not in beta_df.index:
            missing.append(v)
    if missing:
        print(f"ERROR: Variant(s) not found in beta matrix: {missing}")
        print(f"       Available variants (first 10): {list(beta_df.index[:10])}")
        sys.exit(1)

    b1 = beta_df.loc[variant1]  # Series indexed by file path
    b2 = beta_df.loc[variant2]

    # ── Coverage statistics ────────────────────────────────────────────────────
    v1_nonmissing   = b1.notna().sum()
    v2_nonmissing   = b2.notna().sum()
    joint_mask      = b1.notna() & b2.notna()
    n_joint         = joint_mask.sum()

    print(f"\n{'='*60}")
    print(f"DIAGNOSTIC REPORT: {variant1} vs {variant2}")
    print(f"Locus: {locus}  |  Population: {population}")
    print(f"{'='*60}")
    print(f"  N_total_traits      = {n_total}")
    print(f"  N_v1_nonmissing     = {v1_nonmissing}  ({100*v1_nonmissing/n_total:.1f}%)")
    print(f"  N_v2_nonmissing     = {v2_nonmissing}  ({100*v2_nonmissing/n_total:.1f}%)")
    print(f"  N_joint_nonmissing  = {n_joint}  ({100*n_joint/n_total:.1f}%)")

    if n_joint < 2:
        print(f"  pearson_r           = NaN (N_joint < 2)")
    else:
        j1  = b1[joint_mask].values.astype(float)
        j2  = b2[joint_mask].values.astype(float)
        r   = pearson_r(j1, j2)
        print(f"  pearson_r           = {r:.6f}")
        print(f"  (Compare to betanfer_r from comparison TSV to verify)")

    print(f"{'='*60}\n")

    # ── Optional per-trait output ──────────────────────────────────────────────
    if output:
        rows = []
        for fpath in file_list:
            trait = os.path.basename(fpath).replace(".sumstats.txt.gz", "").replace(".sumstats.hg19.txt.gz", "")
            v1b = float(b1[fpath]) if pd.notna(b1[fpath]) else float("nan")
            v2b = float(b2[fpath]) if pd.notna(b2[fpath]) else float("nan")
            both = not (math.isnan(v1b) or math.isnan(v2b))
            rows.append((trait, v1b, v2b, both))

        out_df = pd.DataFrame(rows, columns=["trait", "v1_beta", "v2_beta", "both_nonmissing"])
        os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
        out_df.to_csv(output, sep="\t", index=False)
        print(f"NOTE: Per-trait beta table written to {output}")


if __name__ == "__main__":
    main()
