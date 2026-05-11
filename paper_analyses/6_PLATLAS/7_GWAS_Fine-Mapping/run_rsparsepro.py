'''
run_rsparsepro.py

Per-locus RSparsePro fine-mapping wrapper for the PLATLAS GWAS meta-analysis pipeline.

Extracts the locus region from the bgzip+tabix-indexed meta-analysis summary statistics,
intersects variants with the merged population LD matrix, and runs RSparsePro.

Meta-analysis sumstats format (bgzip+tabix, GWAMA output):
  Columns: rs_number, reference_allele, other_allele, eaf, beta, se, beta_95L, beta_95U,
           z, p-value, -log10_p-value, q_statistic, q_p-value, i2, n_studies, n_samples,
           effects, CHR, POS
  Variant IDs in rs_number have the form CHR:POS:REF:ALT.

LD matrix format:
  .merged.merged.vcor1      tab-delimited square float matrix (plink2 --r-unphased)
  .merged.merged.vcor1.vars one variant ID per line, matching matrix row/col order

Output columns:
  VARIANT_ID, Z, PIP, z_estimated, cs

Usage:
  python run_rsparsepro.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import os
import getopt
import subprocess
import time
import numpy as np
import pandas as pd
from io import StringIO

TABIX_BIN         = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/tabix"
FINEMAP_PYTHON    = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/python"
RSPARSEPRO_SCRIPT = "/flare/GeomicVar/mconery/tools/RSparsePro_LD/src/rsparsepro_ld.py"


def help(exit_num=0):
    print("""Usage: python run_rsparsepro.py [OPTIONS]
Required:
  --sumstats    Path to bgzip+tabix meta-analysis sumstats (.out.gz)
  --ld-file     Path to merged LD matrix (.merged.merged.vcor1)
  --ld-vars     Path to merged LD variant list (.merged.merged.vcor1.vars)
  --chr         Chromosome number
  --start       Locus start position (bp)
  --end         Locus end position (bp)
  --output-file Full path for final output .txt file
Optional:
  --K           Max causal signals for RSparsePro (default: 10)
  -h, --help    Show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  LD UTILITIES  #####################################################
###############################################################################

def read_vars_file(vars_path):
    with open(vars_path, 'r') as fh:
        return [line.strip() for line in fh if line.strip()]


def load_ld_matrix(vcor1_path):
    mat = np.loadtxt(vcor1_path, dtype=np.float64)
    if mat.ndim == 1:
        mat = mat.reshape(1, 1)
    return mat


###############################################################################
##########  SUMSTATS EXTRACTION  ##############################################
###############################################################################

def tabix_extract(gz_path, chrom, start, end):
    """Extract locus region from bgzip+tabix meta-analysis sumstats.

    Returns a DataFrame with columns as named in the file header.
    Returns None if tabix fails or no variants are found.
    """
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        [TABIX_BIN, "-h", gz_path, region],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"  [TABIX] ERROR extracting {gz_path}: {result.stderr.strip()}")
        return None
    lines = result.stdout.splitlines()
    if not lines:
        return None
    header_line = lines[0].lstrip('#')
    data_lines  = lines[1:]
    if not data_lines:
        return None
    content = header_line + "\n" + "\n".join(data_lines)
    df = pd.read_csv(StringIO(content), sep='\t', low_memory=False)
    return df


###############################################################################
##########  FAILURE OUTPUT  ###################################################
###############################################################################

def write_failure_output(output_file, timing_file):
    """Write an empty output file with correct headers and timing=NA."""
    cols = ['VARIANT_ID', 'Z', 'PIP', 'z_estimated', 'cs']
    pd.DataFrame(columns=cols).to_csv(output_file, sep='\t', index=False)
    os.makedirs(os.path.dirname(timing_file), exist_ok=True)
    with open(timing_file, 'w') as fh:
        fh.write("rsparsepro_time_seconds\nNA\n")
    print(f"  Failure output written to {output_file}")


###############################################################################
##########  MAIN  #############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "sumstats=", "ld-file=", "ld-vars=",
             "chr=", "start=", "end=", "output-file=", "K="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    sumstats_file = None
    ld_file       = None
    ld_vars_file  = None
    chrom         = None
    start         = None
    end           = None
    output_file   = None
    K             = 10

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--sumstats":
            sumstats_file = a
        elif o == "--ld-file":
            ld_file = a
        elif o == "--ld-vars":
            ld_vars_file = a
        elif o == "--chr":
            chrom = a
        elif o == "--start":
            start = int(a)
        elif o == "--end":
            end = int(a)
        elif o == "--output-file":
            output_file = a
        elif o == "--K":
            K = int(a)
        else:
            assert False, f"Unhandled option: {o}"

    required = {
        "--sumstats":    sumstats_file,
        "--ld-file":     ld_file,
        "--ld-vars":     ld_vars_file,
        "--chr":         chrom,
        "--start":       start,
        "--end":         end,
        "--output-file": output_file,
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        print(f"ERROR: Missing required arguments: {', '.join(missing)}")
        sys.exit(2)

    # ── Sentinel check ─────────────────────────────────────────────────────────
    if os.path.exists(output_file):
        print(f"NOTE: Output already exists, skipping: {output_file}")
        sys.exit(0)

    # ── Timing and path setup ──────────────────────────────────────────────────
    t_start = time.time()
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    timing_dir  = os.path.join(os.path.dirname(output_file), "timing")
    timing_file = os.path.join(timing_dir, os.path.basename(output_file))
    os.makedirs(timing_dir, exist_ok=True)
    output_prefix = os.path.splitext(output_file)[0]

    # ── Tabix-extract meta-analysis sumstats ───────────────────────────────────
    print(f"NOTE: Extracting chr{chrom}:{start}-{end} from {os.path.basename(sumstats_file)}")
    df = tabix_extract(sumstats_file, chrom, start, end)
    if df is None or df.empty:
        print("NOTE: No variants extracted from sumstats — writing empty output")
        write_failure_output(output_file, timing_file)
        sys.exit(0)

    # ── Read LD variant list ───────────────────────────────────────────────────
    ld_vars = read_vars_file(ld_vars_file)
    if not ld_vars:
        print("NOTE: LD vars file is empty — writing empty output")
        write_failure_output(output_file, timing_file)
        sys.exit(0)

    # ── Intersect sumstats with LD matrix variants ─────────────────────────────
    # Both rs_number and LD vars use CHR:POS:REF:ALT format; direct matching applies.
    ld_vars_set = set(ld_vars)
    df = df[df['rs_number'].astype(str).isin(ld_vars_set)].copy()
    df = df.drop_duplicates(subset='rs_number', keep='first')

    if df.empty:
        print("NOTE: No variants overlap LD matrix — writing empty output")
        write_failure_output(output_file, timing_file)
        sys.exit(0)

    df['z'] = pd.to_numeric(df['z'], errors='coerce')
    df = df.dropna(subset=['z'])

    if df.empty:
        print("NOTE: No valid z-scores after filtering — writing empty output")
        write_failure_output(output_file, timing_file)
        sys.exit(0)

    # Preserve LD vars ordering for the intersection
    sumstats_set = set(df['rs_number'].astype(str))
    ordered_vars = [v for v in ld_vars if v in sumstats_set]
    z_map        = dict(zip(df['rs_number'].astype(str), df['z']))

    print(f"NOTE: {len(ordered_vars)} variants in intersection with LD matrix")

    # ── Write temp z-score file ────────────────────────────────────────────────
    tmp_z  = f"{output_prefix}.tmp.z.txt"
    tmp_ld = f"{output_prefix}.tmp.ld.txt"

    z_df = pd.DataFrame({
        'RSID': ordered_vars,
        'Z':    [z_map[v] for v in ordered_vars],
    })
    z_df.to_csv(tmp_z, sep='\t', index=False)

    # ── Load and subset LD matrix ──────────────────────────────────────────────
    print(f"NOTE: Loading LD matrix from {os.path.basename(ld_file)}")
    try:
        ld_mat = load_ld_matrix(ld_file)
    except Exception as ex:
        print(f"ERROR: LD matrix load failed: {ex}")
        if os.path.exists(tmp_z):
            os.remove(tmp_z)
        write_failure_output(output_file, timing_file)
        sys.exit(0)

    idx_map     = {v: i for i, v in enumerate(ld_vars)}
    row_indices = [idx_map[v] for v in ordered_vars]
    ld_sub      = ld_mat[np.ix_(row_indices, row_indices)]

    np.savetxt(tmp_ld, ld_sub, fmt='%.6f', delimiter=' ')

    # ── Run RSparsePro ─────────────────────────────────────────────────────────
    tmp_save = f"{output_prefix}.tmp"
    print(f"NOTE: Running RSparsePro (K={K}, {len(ordered_vars)} variants)...")
    t_rs = time.time()
    proc = subprocess.run(
        [FINEMAP_PYTHON, RSPARSEPRO_SCRIPT,
         "--z", tmp_z, "--ld", tmp_ld, "--save", tmp_save, "--K", str(K)],
        capture_output=True, text=True
    )
    elapsed_rs = time.time() - t_rs

    for f in [tmp_z, tmp_ld]:
        if os.path.exists(f):
            os.remove(f)

    if proc.returncode != 0:
        print(f"ERROR: RSparsePro failed:\n{proc.stderr.strip()}")
        write_failure_output(output_file, timing_file)
        sys.exit(0)

    print(f"NOTE: RSparsePro completed in {elapsed_rs:.1f}s")

    # ── Read RSparsePro output and write final result ──────────────────────────
    tmp_result  = f"{tmp_save}.rsparsepro.txt"
    tmp_log     = f"{tmp_save}.rsparsepro.log"
    try:
        res_df = pd.read_csv(tmp_result, sep='\t')
        res_df = res_df.rename(columns={'RSID': 'VARIANT_ID'})
        res_df.to_csv(output_file, sep='\t', index=False, na_rep='NA')
        print(f"  Output written: {output_file} ({len(res_df)} variants)")
    except Exception as ex:
        print(f"ERROR: Could not read RSparsePro output {tmp_result}: {ex}")
        write_failure_output(output_file, timing_file)
    finally:
        for f in [tmp_result, tmp_log]:
            if os.path.exists(f):
                os.remove(f)

    # ── Write timing file ──────────────────────────────────────────────────────
    total_elapsed = time.time() - t_start
    with open(timing_file, 'w') as fh:
        fh.write("rsparsepro_time_seconds\n")
        fh.write(f"{total_elapsed:.2f}\n")

    print(f"NOTE: run_rsparsepro.py complete ({total_elapsed:.1f}s total)")


if __name__ == "__main__":
    main()
