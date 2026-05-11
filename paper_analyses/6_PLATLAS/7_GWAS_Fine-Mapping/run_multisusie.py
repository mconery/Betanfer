'''
run_multisusie.py

Per-locus MultiSuSiE multi-population fine-mapping wrapper for the PLATLAS
GWAS meta-analysis pipeline.

For each contributing cohort, this script:
  1. Extracts the locus region from the bgzip+tabix-indexed summary statistics.
  2. Reads the cohort's LD matrix (.vcor1) and variant list (.vcor1.vars).
  3. Filters sumstats to variants present in the LD matrix.

After per-cohort extraction, the script finds the strict variant intersection
across all cohorts' sumstats AND LD matrices, then calls MultiSuSiE.multisusie_rss
with the aligned beta / SE / LD inputs.

Sumstats format (bgzip+tabix, PLATLAS convention):
  Columns: #ID CHR POS REF ALT BETA SE P LOG10P ... AAF ... N ...
  Variant IDs in #ID have the form CHR:POS:REF:ALT.

var_y is set to 1.0 (appropriate for standardised continuous traits).
Use --estimate-residual-variance to let MultiSuSiE estimate it from the data
instead (recommended for unstandardised traits).

Output columns:
  VARIANT_ID, Z_{ancestry_label} (one per cohort), PIP, converged, cs

Usage:
  python run_multisusie.py [OPTIONS]

***Run in MultiSuSiE conda environment***
'''

import sys
import os
import getopt
import subprocess
import time
import re
import numpy as np
import pandas as pd

# ── Tool paths ────────────────────────────────────────────────────────────────
TABIX_BIN = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/tabix"


def help(exit_num=0):
    print("""Usage: python run_multisusie.py [OPTIONS]
Required:
  --gwas-files               Comma-separated bgzip sumstats paths (one per cohort)
  --ld-files                 Comma-separated .vcor1 LD file paths (one per cohort)
  --n-gwas                   Comma-separated sample sizes (one per cohort)
  --ancestry-labels          Comma-separated ancestry label strings (one per cohort)
  --output-file              Full path for final output .txt file
  --chr                      Chromosome number
  --start                    Locus start position (bp)
  --end                      Locus end position (bp)
Optional:
  --rho                      Off-diagonal effect-size correlation (default: 0.75)
  --L                        Max number of causal signals (default: 10)
  --coverage                 Credible set coverage (default: 0.95)
  --min-purity               Min credible set purity (default: 0.5)
  --max-iter                 Max iterations (default: 100)
  --tol                      Convergence tolerance (default: 1e-3)
  --estimate-residual-variance  Estimate residual variance from data (flag)
  -h, --help                 Show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  LD MATRIX UTILITIES  ##############################################
###############################################################################

def read_vars_file(vars_path):
    """Return ordered list of variant IDs from a .vcor1.vars file."""
    with open(vars_path, 'r') as fh:
        return [line.strip() for line in fh if line.strip()]


def load_ld_matrix(vcor1_path, n_vars):
    """Load a .vcor1 tab-delimited LD matrix as a float64 numpy array."""
    mat = np.loadtxt(vcor1_path, dtype=np.float64)
    if mat.ndim == 1:
        mat = mat.reshape(1, 1)
    return mat


###############################################################################
##########  SUMSTATS EXTRACTION  ##############################################
###############################################################################

def tabix_extract(gz_path, chrom, start, end):
    """Extract locus region from bgzip+tabix sumstats.

    Returns a DataFrame with the region data, with #ID renamed to SNP.
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
    from io import StringIO
    content = header_line + "\n" + "\n".join(data_lines)
    df = pd.read_csv(StringIO(content), sep='\t', low_memory=False)
    if 'ID' in df.columns:
        df = df.rename(columns={'ID': 'SNP'})
    return df


###############################################################################
##########  PER-COHORT PREPARATION  ###########################################
###############################################################################

def prepare_cohort(gz_path, ld_prefix, label, chrom, start, end):
    """Extract and prepare all inputs for one cohort.

    Returns a dict with:
      'vars':      ordered list of variant IDs from .vcor1.vars
      'sumstats':  DataFrame filtered to vars-present rows (no duplicates)
      'ld_mat':    numpy LD matrix (full size matching vars list)
      'z_scores':  dict SNP -> Z-score for output column
    or None on failure.
    """
    vars_path = ld_prefix + ".vars"
    if not os.path.exists(vars_path):
        print(f"  ERROR: .vars file not found: {vars_path}")
        return None
    if not os.path.exists(ld_prefix):
        print(f"  ERROR: LD matrix not found: {ld_prefix}")
        return None

    vars_list = read_vars_file(vars_path)
    if not vars_list:
        print(f"  ERROR: .vars file is empty: {vars_path}")
        return None

    # Tabix-extract the locus region
    df = tabix_extract(gz_path, chrom, start, end)
    if df is None or df.empty:
        print(f"  ERROR: no sumstats extracted for cohort {label} at "
              f"chr{chrom}:{start}-{end}")
        return None

    # Filter to variants present in LD matrix; remove duplicates
    vars_set = set(vars_list)
    df = df[df['SNP'].astype(str).isin(vars_set)].copy()
    df = df.drop_duplicates(subset='SNP', keep='first')
    if df.empty:
        print(f"  WARNING: no sumstats variants overlap LD matrix for {label}")
        return None

    # Validate BETA and SE columns
    try:
        df['BETA'] = pd.to_numeric(df['BETA'], errors='coerce')
        df['SE']   = pd.to_numeric(df['SE'],   errors='coerce')
        df = df.dropna(subset=['BETA', 'SE'])
        df = df[df['SE'] > 0]
    except Exception as ex:
        print(f"  ERROR: BETA/SE validation failed for {label}: {ex}")
        return None

    if df.empty:
        print(f"  WARNING: no valid BETA/SE variants for {label}")
        return None

    # Z-scores for output column
    z_scores = dict(zip(df['SNP'].astype(str),
                        df['BETA'] / df['SE']))

    # Read N: use the max N observed across variants in the locus
    try:
        n_val = int(df['N'].astype(float).max())
    except (KeyError, ValueError):
        n_val = 0

    # Load LD matrix
    try:
        ld_mat = load_ld_matrix(ld_prefix, len(vars_list))
    except Exception as ex:
        print(f"  ERROR: LD matrix load failed for {ld_prefix}: {ex}")
        return None

    if ld_mat.shape[0] != len(vars_list):
        print(f"  ERROR: LD matrix size {ld_mat.shape[0]} != "
              f"vars count {len(vars_list)} for {label}")
        return None

    return {
        'vars':     vars_list,
        'vars_set': vars_set,
        'sumstats': df,
        'ld_mat':   ld_mat,
        'z_scores': z_scores,
        'n':        n_val,
    }


###############################################################################
##########  ALIGN COHORTS TO SHARED VARIANT SET  ##############################
###############################################################################

def align_cohorts(cohort_data_list):
    """Find the strict variant intersection and align all cohort data to it.

    Returns (intersection_vars, aligned_list) where aligned_list has each
    cohort's Bhat, Shat, and R subsetted and reordered to the intersection.
    Returns (None, None) if intersection is empty.
    """
    # Intersection: variants present in ALL cohorts' sumstats AND LD matrices
    sets = [set(d['sumstats']['SNP'].astype(str)) & d['vars_set']
            for d in cohort_data_list]
    intersection = sets[0]
    for s in sets[1:]:
        intersection = intersection & s

    if not intersection:
        return None, None

    # Use the order from the first cohort's vars list for consistency
    first_vars = cohort_data_list[0]['vars']
    ordered_vars = [v for v in first_vars if v in intersection]

    if not ordered_vars:
        return None, None

    aligned = []
    for d in cohort_data_list:
        # Reorder sumstats to match ordered_vars using dict lookups (safe for duplicates)
        beta_map = dict(zip(d['sumstats']['SNP'].astype(str), d['sumstats']['BETA']))
        se_map   = dict(zip(d['sumstats']['SNP'].astype(str), d['sumstats']['SE']))
        bhat = np.array([float(beta_map[v]) for v in ordered_vars], dtype=np.float64)
        shat = np.array([float(se_map[v])   for v in ordered_vars], dtype=np.float64)

        # Subset and reorder LD matrix
        full_vars  = d['vars']
        idx_map    = {v: i for i, v in enumerate(full_vars)}
        row_indices = [idx_map[v] for v in ordered_vars]
        r_sub = d['ld_mat'][np.ix_(row_indices, row_indices)]

        # Validate
        if np.any(np.isnan(bhat)) or np.any(np.isnan(shat)):
            bhat = np.nan_to_num(bhat, nan=0.0)
            shat = np.where(shat == 0, 1.0, shat)
        if np.any(np.isnan(r_sub)):
            r_sub = np.nan_to_num(r_sub, nan=0.0)
            np.fill_diagonal(r_sub, 1.0)

        aligned.append({'Bhat': bhat, 'Shat': shat, 'R': r_sub, 'n': d['n']})

    return ordered_vars, aligned


###############################################################################
##########  FAILURE OUTPUT  ###################################################
###############################################################################

def write_failure_output(output_file, ancestry_labels, timing_file):
    """Write an empty output file with correct headers and timing=NA."""
    cols = (['VARIANT_ID']
            + [f'Z_{lbl}' for lbl in ancestry_labels]
            + ['PIP', 'converged', 'cs'])
    pd.DataFrame(columns=cols).to_csv(output_file, sep='\t', index=False)
    os.makedirs(os.path.dirname(timing_file), exist_ok=True)
    with open(timing_file, 'w') as fh:
        fh.write("multisusie_time_seconds\nNA\n")
    print(f"  Failure output written to {output_file}")


###############################################################################
##########  MAIN  #############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "gwas-files=", "ld-files=", "n-gwas=", "ancestry-labels=",
             "output-file=", "chr=", "start=", "end=",
             "rho=", "L=", "coverage=", "min-purity=",
             "max-iter=", "tol=", "estimate-residual-variance"]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    gwas_files_str  = None
    ld_files_str    = None
    n_gwas_str      = None
    ancestry_str    = None
    output_file     = None
    chrom           = None
    start           = None
    end             = None
    rho             = 0.75
    L               = 10
    coverage        = 0.95
    min_purity      = 0.5
    max_iter        = 100
    tol             = 1e-3
    estimate_resid  = False

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--gwas-files":
            gwas_files_str = a
        elif o == "--ld-files":
            ld_files_str = a
        elif o == "--n-gwas":
            n_gwas_str = a
        elif o == "--ancestry-labels":
            ancestry_str = a
        elif o == "--output-file":
            output_file = a
        elif o == "--chr":
            chrom = a
        elif o == "--start":
            start = int(a)
        elif o == "--end":
            end = int(a)
        elif o == "--rho":
            rho = float(a)
        elif o == "--L":
            L = int(a)
        elif o == "--coverage":
            coverage = float(a)
        elif o == "--min-purity":
            min_purity = float(a)
        elif o == "--max-iter":
            max_iter = int(a)
        elif o == "--tol":
            tol = float(a)
        elif o == "--estimate-residual-variance":
            estimate_resid = True
        else:
            assert False, f"Unhandled option: {o}"

    required = {
        "--gwas-files":      gwas_files_str,
        "--ld-files":        ld_files_str,
        "--n-gwas":          n_gwas_str,
        "--ancestry-labels": ancestry_str,
        "--output-file":     output_file,
        "--chr":             chrom,
        "--start":           start,
        "--end":             end,
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        print(f"ERROR: Missing required arguments: {', '.join(missing)}")
        sys.exit(2)

    gwas_files      = [f.strip() for f in gwas_files_str.split(',')]
    ld_prefixes     = [f.strip() for f in ld_files_str.split(',')]
    n_gwas_list     = [int(n.strip()) for n in n_gwas_str.split(',')]
    ancestry_labels = [lbl.strip() for lbl in ancestry_str.split(',')]

    n_pops = len(gwas_files)
    if not (len(ld_prefixes) == len(n_gwas_list) == len(ancestry_labels) == n_pops):
        print("ERROR: --gwas-files, --ld-files, --n-gwas, --ancestry-labels "
              "must all have the same number of entries")
        sys.exit(2)

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    timing_dir  = os.path.join(os.path.dirname(output_file), "timing")
    timing_file = os.path.join(timing_dir, os.path.basename(output_file))
    os.makedirs(timing_dir, exist_ok=True)

    t_start = time.time()

    # ── Per-cohort preparation ─────────────────────────────────────────────────
    print(f"NOTE: Preparing {n_pops} cohort(s) for chr{chrom}:{start}-{end}")
    cohort_data = []
    for gz, ld_prefix, label in zip(gwas_files, ld_prefixes, ancestry_labels):
        print(f"  Cohort {label}: {os.path.basename(gz)}")
        d = prepare_cohort(gz, ld_prefix, label, chrom, start, end)
        if d is None:
            write_failure_output(output_file, ancestry_labels, timing_file)
            sys.exit(0)
        cohort_data.append(d)

    # ── Align all cohorts to variant intersection ──────────────────────────────
    ordered_vars, aligned = align_cohorts(cohort_data)
    if ordered_vars is None:
        print("NOTE: No variants in intersection — writing empty output")
        write_failure_output(output_file, ancestry_labels, timing_file)
        sys.exit(0)

    n_vars = len(ordered_vars)
    print(f"NOTE: {n_vars} variants in intersection across {n_pops} cohort(s)")

    # ── Build MultiSuSiE inputs ────────────────────────────────────────────────
    Bhat_list = [a['Bhat'] for a in aligned]
    Shat_list = [a['Shat'] for a in aligned]
    R_list    = [a['R']    for a in aligned]
    n_list    = [a['n']    for a in aligned]

    # Effect-size correlation matrix (identity with off-diagonal = rho)
    rho_mat = np.full((n_pops, n_pops), rho)
    np.fill_diagonal(rho_mat, 1.0)

    # ── Run MultiSuSiE ─────────────────────────────────────────────────────────
    try:
        import MultiSuSiE
    except ImportError:
        print("ERROR: MultiSuSiE package not found in current Python environment")
        write_failure_output(output_file, ancestry_labels, timing_file)
        sys.exit(1)

    print("NOTE: Running MultiSuSiE.multisusie_rss ...")
    t_ms = time.time()
    try:
        result = MultiSuSiE.multisusie_rss(
            R_list=R_list,
            population_sizes=n_list,
            b_list=Bhat_list,
            s_list=Shat_list,
            varY_list=[1.0] * n_pops,
            rho=rho_mat,
            L=L,
            coverage=coverage,
            min_abs_corr=min_purity,
            max_iter=max_iter,
            tol=tol,
            estimate_residual_variance=estimate_resid,
        )
    except Exception as ex:
        elapsed = time.time() - t_start
        print(f"ERROR: MultiSuSiE failed: {ex}")
        write_failure_output(output_file, ancestry_labels, timing_file)
        sys.exit(0)

    elapsed_ms = time.time() - t_ms
    print(f"NOTE: MultiSuSiE completed in {elapsed_ms:.1f}s")

    # ── Parse result and write output ──────────────────────────────────────────
    try:
        # PIP: sum of alpha across causal components
        if hasattr(result, 'pip'):
            pip_arr = np.array(result.pip)
        elif hasattr(result, 'alpha'):
            pip_arr = np.array(result.alpha).sum(axis=0)
        else:
            pip_arr = np.zeros(n_vars)

        converged = bool(result.converged) if hasattr(result, 'converged') else False

        # Credible set membership
        # result.sets is a list: [cs_index_arrays, purity, claimed_coverage, include_mask]
        # include_mask[i] is True if CS i passes the purity filter.
        cs_assignments = np.zeros(n_vars, dtype=int)
        if hasattr(result, 'sets') and result.sets is not None and len(result.sets) >= 4:
            cs_list = result.sets[0]       # list of index arrays, one per L component
            include_mask = result.sets[3]  # bool array; True = CS passes purity
            cs_num = 0
            for i, cs_indices in enumerate(cs_list):
                if include_mask[i] and len(cs_indices) > 0:
                    cs_num += 1
                    for idx in cs_indices:
                        if 0 <= idx < n_vars:
                            cs_assignments[idx] = cs_num

    except Exception as ex:
        print(f"WARNING: result parsing encountered an error: {ex}; "
              "PIP/cs set to 0")
        pip_arr        = np.zeros(n_vars)
        converged      = False
        cs_assignments = np.zeros(n_vars, dtype=int)

    out = pd.DataFrame({'VARIANT_ID': ordered_vars})

    # Z-scores from original sumstats
    for label, d in zip(ancestry_labels, cohort_data):
        out[f'Z_{label}'] = out['VARIANT_ID'].map(d['z_scores']).fillna('NA')

    out['PIP']       = pip_arr
    out['converged'] = converged
    out['cs']        = cs_assignments

    out.to_csv(output_file, sep='\t', index=False, na_rep='NA')
    print(f"  Output written: {output_file} ({len(out)} variants)")

    # ── Write timing file ──────────────────────────────────────────────────────
    total_elapsed = time.time() - t_start
    with open(timing_file, 'w') as fh:
        fh.write("multisusie_time_seconds\n")
        fh.write(f"{total_elapsed:.2f}\n")

    print(f"NOTE: run_multisusie.py complete ({total_elapsed:.1f}s total)")


if __name__ == "__main__":
    main()
