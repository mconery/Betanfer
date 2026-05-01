'''
run_susiex.py

Per-locus SuSiEx multi-population fine-mapping wrapper for the PLATLAS GWAS
meta-analysis pipeline.

For each contributing cohort, this script:
  1. Extracts the locus region from the bgzip+tabix-indexed summary statistics.
  2. Filters to variants present in the cohort's LD matrix (.vcor1.vars).
  3. Writes a SuSiEx-format summary statistics temp file.
  4. Synthesises a PLINK frequency file (_frq.frq) from the AAF column in the
     sumstats — no reference panel or PLINK binary is needed for this step.
  5. Synthesises a PLINK BIM file (_ref.bim) from the .vcor1.vars variant list
     and the allele/position information extracted from variant IDs.
  6. Converts the .vcor1 LD matrix to float32 binary (.vcor1.ld.bin).

After per-cohort preparation, the script finds the variant intersection across
all cohorts, applies a MAF filter, and calls the SuSiEx binary.  SuSiEx output
files (.snp, .cs, .summary) are parsed into a standard tab-delimited output,
then all temporary files are removed.

Sumstats format (bgzip+tabix, PLATLAS convention):
  Columns: #ID CHR POS REF ALT BETA SE P LOG10P ... AAF ... N ...
  Variant IDs in #ID have the form CHR:POS:REF:ALT.

Output columns:
  VARIANT_ID, Z_{ancestry_label} (one per cohort), PIP, converged, cs,
  LogBF_CS{n}_{ancestry_label} (one per credible set per cohort)

Usage:
  python run_susiex.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import os
import getopt
import subprocess
import tempfile
import shutil
import time
import re
import numpy as np
import pandas as pd

# ── Tool paths ────────────────────────────────────────────────────────────────
TABIX_BIN   = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/tabix"
SUSIEX_BIN  = "/flare/GeomicVar/mconery/tools/SuSiEx/bin/SuSiEx"
PLINK_BIN   = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/plink"

# Column indices (1-based) in the temp sumstats file we write for SuSiEx
# File has columns: SNP CHR POS REF ALT BETA SE P
_COL_SNP  = 1
_COL_CHR  = 2
_COL_POS  = 3
_COL_A2   = 4   # REF = other allele
_COL_A1   = 5   # ALT = effect allele
_COL_BETA = 6
_COL_SE   = 7
_COL_P    = 8


def help(exit_num=0):
    print("""Usage: python run_susiex.py [OPTIONS]
Required:
  --gwas-files        Comma-separated bgzip sumstats paths (one per cohort)
  --ld-files          Comma-separated .vcor1 LD file paths (one per cohort)
  --n-gwas            Comma-separated sample sizes (one per cohort)
  --ancestry-labels   Comma-separated ancestry label strings (one per cohort)
  --output-file       Full path for final output .txt file
  --chr               Chromosome number
  --start             Locus start position (bp)
  --end               Locus end position (bp)
Optional:
  --maf               MAF threshold (default: 0.005)
  --n-sig             Max number of causal signals (default: 10)
  --level             Credible set coverage (default: 0.95)
  --min-purity        Min credible set purity (default: 0.5)
  --pval-thresh       P-value threshold (default: 1e-5)
  --max-iter          Max iterations (default: 100)
  --tol               Convergence tolerance (default: 1e-4)
  -h, --help          Show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  LD MATRIX UTILITIES  ##############################################
###############################################################################

def read_vars_file(vars_path):
    """Return ordered list of variant IDs from a .vcor1.vars file."""
    with open(vars_path, 'r') as fh:
        return [line.strip() for line in fh if line.strip()]


def convert_vcor1_to_bin(vcor1_path, bin_path):
    """Convert a tab-delimited .vcor1 matrix to float32 binary .ld.bin."""
    mat = np.loadtxt(vcor1_path, dtype=np.float64)
    if mat.ndim == 1:
        mat = mat.reshape(1, 1)
    mat_f32 = mat.astype(np.float32)
    with open(bin_path, 'wb') as fh:
        fh.write(mat_f32.tobytes('C'))


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
    # The header line may start with '#'; strip it for pandas
    header_line = lines[0].lstrip('#')
    data_lines  = lines[1:]
    if not data_lines:
        return None
    from io import StringIO
    content = header_line + "\n" + "\n".join(data_lines)
    df = pd.read_csv(StringIO(content), sep='\t', low_memory=False)
    # Rename the ID column (was '#ID' before stripping '#')
    if 'ID' in df.columns:
        df = df.rename(columns={'ID': 'SNP'})
    return df


###############################################################################
##########  SYNTHETIC REFERENCE FILE CREATION  ################################
###############################################################################

def parse_variant_id(vid):
    """Parse CHR:POS:REF:ALT from a variant ID string.

    Returns (chrom, pos, ref, alt) or None on failure.
    """
    m = re.match(r'^(\w+):(\d+):([^:]+):([^:]+)$', vid)
    if m:
        return m.group(1), int(m.group(2)), m.group(3), m.group(4)
    return None


def write_frq_file(frq_path, sumstats_df, vars_list):
    """Write a synthetic PLINK-format .frq file.

    Uses AAF from sumstats (effect allele = ALT frequency).
    Variants in vars_list but absent from sumstats are assigned MAF=0.5
    (conservative; SuSiEx will filter these if below --maf threshold).
    """
    # Build lookup from sumstats: SNP -> (CHR, A1, A2, AAF, N)
    lookup = {}
    for _, row in sumstats_df.iterrows():
        snp = str(row['SNP'])
        try:
            aaf = float(row['AAF'])
            n   = float(row['N'])
        except (KeyError, ValueError, TypeError):
            aaf, n = 0.5, 0
        lookup[snp] = (row['CHR'], row['ALT'], row['REF'], aaf, n)

    rows = []
    for vid in vars_list:
        if vid in lookup:
            chrom_val, a1, a2, aaf, n = lookup[vid]
            maf     = min(aaf, 1.0 - aaf)
            nchrobs = int(2 * n) if n > 0 else 2
        else:
            parsed = parse_variant_id(vid)
            if parsed:
                chrom_val, _, a2, a1 = parsed  # REF=A2, ALT=A1
            else:
                chrom_val, a1, a2 = "0", "A", "G"
            maf, nchrobs = 0.5, 2
        rows.append(f"{chrom_val}\t{vid}\t{a1}\t{a2}\t{maf:.6f}\t{nchrobs}")

    with open(frq_path, 'w') as fh:
        fh.write("CHR\tSNP\tA1\tA2\tMAF\tNCHROBS\n")
        fh.write("\n".join(rows) + "\n")


def write_bim_file(bim_path, vars_list, sumstats_df):
    """Write a synthetic PLINK-format .bim file (no header).

    Columns: CHR SNP CM BP A1(ALT) A2(REF)
    Allele and position info from sumstats where available; parsed from
    variant ID otherwise.
    """
    lookup = {}
    for _, row in sumstats_df.iterrows():
        snp = str(row['SNP'])
        try:
            lookup[snp] = (row['CHR'], int(row['POS']), row['ALT'], row['REF'])
        except (KeyError, ValueError, TypeError):
            pass

    rows = []
    for vid in vars_list:
        if vid in lookup:
            chrom_val, pos, a1, a2 = lookup[vid]
        else:
            parsed = parse_variant_id(vid)
            if parsed:
                chrom_val, pos, a2, a1 = parsed  # REF=A2, ALT=A1
            else:
                chrom_val, pos, a1, a2 = "0", 0, "A", "G"
        rows.append(f"{chrom_val}\t{vid}\t0\t{pos}\t{a1}\t{a2}")

    with open(bim_path, 'w') as fh:
        fh.write("\n".join(rows) + "\n")


###############################################################################
##########  PER-COHORT PREPARATION  ###########################################
###############################################################################

def prepare_cohort(gz_path, ld_prefix, n_gwas, label, chrom, start, end):
    """Extract and prepare all inputs for one cohort.

    Returns a dict with:
      'vars':       ordered list of variant IDs from .vcor1.vars
      'sumstats':   DataFrame filtered to vars-present rows
      'ld_bin':     path to .ld.bin (created here if not existing)
      'frq_file':   path to synthetic _frq.frq
      'bim_file':   path to synthetic _ref.bim
      'sst_file':   path to SuSiEx-format sumstats temp file
      'z_scores':   dict of SNP -> Z-score (for output column)
    or None on failure.
    """
    vars_path = ld_prefix + ".vars"
    bin_path  = ld_prefix + ".ld.bin"
    frq_path  = ld_prefix + "_frq.frq"
    bim_path  = ld_prefix + "_ref.bim"

    if not os.path.exists(vars_path):
        print(f"  ERROR: .vars file not found: {vars_path}")
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

    # Filter sumstats to variants present in the LD matrix; remove duplicates
    vars_set = set(vars_list)
    df = df[df['SNP'].astype(str).isin(vars_set)].copy()
    df = df.drop_duplicates(subset='SNP', keep='first')
    if df.empty:
        print(f"  WARNING: no sumstats variants overlap LD matrix for {label}")
        return None

    # Compute Z-scores for output column
    try:
        z_scores = dict(zip(df['SNP'].astype(str),
                            df['BETA'].astype(float) / df['SE'].astype(float)))
    except Exception:
        z_scores = {}

    # Write SuSiEx-format sumstats temp file: SNP CHR POS REF ALT BETA SE P
    sst_path = ld_prefix + "_susiex_sst.txt"
    sst_df = df[['SNP', 'CHR', 'POS', 'REF', 'ALT', 'BETA', 'SE', 'P']].copy()
    sst_df.to_csv(sst_path, sep='\t', index=False, na_rep='NA')

    # Synthetic reference files (needed by SuSiEx binary)
    write_frq_file(frq_path, df, vars_list)
    write_bim_file(bim_path, vars_list, df)

    # Convert LD matrix to binary (skip if already done by a prior run)
    if not os.path.exists(bin_path):
        try:
            convert_vcor1_to_bin(ld_prefix, bin_path)
        except Exception as ex:
            print(f"  ERROR: LD conversion failed for {ld_prefix}: {ex}")
            return None

    return {
        'vars':     vars_list,
        'vars_set': set(vars_list),
        'sumstats': df,
        'ld_bin':   bin_path,
        'frq_file': frq_path,
        'bim_file': bim_path,
        'sst_file': sst_path,
        'z_scores': z_scores,
    }


###############################################################################
##########  INTERSECT AND MAF FILTER  #########################################
###############################################################################

def intersect_and_filter(cohort_data_list, maf_thresh):
    """Find variant intersection across all cohorts and apply MAF filter.

    Returns the set of variants that pass, or an empty set if none remain.
    """
    # Intersection across all cohorts' sumstats-present variants
    present_sets = [set(d['sumstats']['SNP'].astype(str)) for d in cohort_data_list]
    intersection = present_sets[0]
    for s in present_sets[1:]:
        intersection = intersection & s

    if not intersection:
        return set()

    # Build per-cohort MAF lookup dicts (SNP -> MAF) for vectorised filtering
    maf_lookups = []
    for d in cohort_data_list:
        df = d['sumstats']
        try:
            aaf_series = pd.to_numeric(df['AAF'], errors='coerce')
            maf_series = aaf_series.apply(lambda a: min(a, 1.0 - a) if pd.notna(a) else 0.5)
            maf_lookups.append(dict(zip(df['SNP'].astype(str), maf_series)))
        except KeyError:
            maf_lookups.append({})

    # Keep variant only if ALL cohorts have MAF >= threshold
    passing = set()
    for vid in intersection:
        keep = all(
            lookup.get(vid, 0.5) >= maf_thresh
            for lookup in maf_lookups
        )
        if keep:
            passing.add(vid)

    return passing


###############################################################################
##########  OUTPUT PARSING  ###################################################
###############################################################################

def parse_susiex_output(snp_file, ancestry_labels, cohort_data_list, output_file):
    """Parse SuSiEx .snp output and write standard-format output file.

    Output columns:
      VARIANT_ID, Z_{label} (one per cohort), PIP, converged, cs,
      LogBF_CS{n}_{label} (one per CS per cohort)
    """
    try:
        snp_df = pd.read_csv(snp_file, sep='\t')
    except Exception as ex:
        print(f"  ERROR: could not read SuSiEx .snp file {snp_file}: {ex}")
        return False

    pip_cols   = [c for c in snp_df.columns if c.startswith('PIP(CS')]
    logbf_cols = [c for c in snp_df.columns if c.startswith('LogBF(CS')]

    # Overall PIP = sum of per-CS PIPs
    snp_df['PIP'] = snp_df[pip_cols].sum(axis=1) if pip_cols else 0.0

    out = pd.DataFrame()
    out['VARIANT_ID'] = snp_df['SNP']

    # Z-scores from original sumstats (BETA/SE)
    for i, (label, d) in enumerate(zip(ancestry_labels, cohort_data_list)):
        z_col = f"Z_{label}"
        out[z_col] = out['VARIANT_ID'].map(d['z_scores']).fillna('NA')

    out['PIP']       = snp_df['PIP']
    out['converged'] = 'NA'

    # Credible set membership: use the first CS with PIP > 0.01
    cs_assignments = []
    for _, row in snp_df.iterrows():
        cs_num = 0
        for pc in pip_cols:
            if pd.notna(row[pc]) and row[pc] > 0.01:
                m = re.match(r'PIP\(CS(\d+)\)', pc)
                if m:
                    cs_num = int(m.group(1))
                    break
        cs_assignments.append(cs_num)
    out['cs'] = cs_assignments

    # LogBF columns renamed by population
    for lc in logbf_cols:
        m = re.match(r'LogBF\(CS(\d+),Pop(\d+)\)', lc)
        if m:
            cs_num  = m.group(1)
            pop_idx = int(m.group(2)) - 1
            label   = ancestry_labels[pop_idx] if pop_idx < len(ancestry_labels) else f"Pop{pop_idx+1}"
            out[f"LogBF_CS{cs_num}_{label}"] = snp_df[lc]

    out.to_csv(output_file, sep='\t', index=False, na_rep='NA')
    print(f"  Output written: {output_file} ({len(out)} variants)")
    return True


###############################################################################
##########  FAILURE OUTPUT  ###################################################
###############################################################################

def write_failure_output(output_file, ancestry_labels, timing_file, elapsed):
    """Write an empty output file with correct headers and timing=NA."""
    cols = (['VARIANT_ID']
            + [f'Z_{lbl}' for lbl in ancestry_labels]
            + ['PIP', 'converged', 'cs'])
    pd.DataFrame(columns=cols).to_csv(output_file, sep='\t', index=False)
    os.makedirs(os.path.dirname(timing_file), exist_ok=True)
    with open(timing_file, 'w') as fh:
        fh.write("susiex_time_seconds\nNA\n")
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
             "maf=", "n-sig=", "level=", "min-purity=",
             "pval-thresh=", "max-iter=", "tol="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    gwas_files_str    = None
    ld_files_str      = None
    n_gwas_str        = None
    ancestry_str      = None
    output_file       = None
    chrom             = None
    start             = None
    end               = None
    maf               = 0.005
    n_sig             = 10
    level             = 0.95
    min_purity        = 0.5
    pval_thresh       = 1e-5
    max_iter          = 100
    tol               = 1e-4

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
        elif o == "--maf":
            maf = float(a)
        elif o == "--n-sig":
            n_sig = int(a)
        elif o == "--level":
            level = float(a)
        elif o == "--min-purity":
            min_purity = float(a)
        elif o == "--pval-thresh":
            pval_thresh = float(a)
        elif o == "--max-iter":
            max_iter = int(a)
        elif o == "--tol":
            tol = float(a)
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
    for gz, ld_prefix, n_val, label in zip(
        gwas_files, ld_prefixes, n_gwas_list, ancestry_labels
    ):
        print(f"  Cohort {label}: {os.path.basename(gz)}")
        d = prepare_cohort(gz, ld_prefix, n_val, label, chrom, start, end)
        if d is None:
            write_failure_output(output_file, ancestry_labels, timing_file,
                                 time.time() - t_start)
            sys.exit(0)
        cohort_data.append(d)

    # ── Variant intersection and MAF filter ───────────────────────────────────
    passing_variants = intersect_and_filter(cohort_data, maf)
    if not passing_variants:
        print("NOTE: No variants remain after intersection/MAF filter — writing empty output")
        write_failure_output(output_file, ancestry_labels, timing_file,
                             time.time() - t_start)
        # Cleanup temp files
        for d in cohort_data:
            for f in [d['sst_file'], d['frq_file'], d['bim_file'], d['ld_bin']]:
                if os.path.exists(f):
                    os.remove(f)
        sys.exit(0)

    print(f"NOTE: {len(passing_variants)} variants pass intersection + MAF filter")

    # ── Filter temp sumstats files to passing variants only ───────────────────
    sst_files_filtered = []
    for d, ld_prefix in zip(cohort_data, ld_prefixes):
        sst_df   = pd.read_csv(d['sst_file'], sep='\t')
        sst_filt = sst_df[sst_df['SNP'].astype(str).isin(passing_variants)]
        filt_path = ld_prefix + "_susiex_sst_filt.txt"
        sst_filt.to_csv(filt_path, sep='\t', index=False, na_rep='NA')
        sst_files_filtered.append(filt_path)

    # ── Call SuSiEx binary ─────────────────────────────────────────────────────
    output_base = os.path.splitext(output_file)[0]
    out_dir  = os.path.dirname(output_base)
    out_name = os.path.basename(output_base)
    bp_range = f"{start},{end}"

    # Column numbers in the temp sumstats files (1-based)
    n_pops_str = ",".join(str(n) for n in n_gwas_list)

    cmd = [
        SUSIEX_BIN,
        "--sst_file",   ",".join(sst_files_filtered),
        "--n_gwas",     ",".join(str(n) for n in n_gwas_list),
        "--ld_file",    ",".join(ld_prefixes),
        "--out_dir",    out_dir if out_dir else ".",
        "--out_name",   out_name,
        "--chr",        str(chrom),
        "--bp",         bp_range,
        "--chr_col",    ",".join([str(_COL_CHR)]  * n_pops),
        "--snp_col",    ",".join([str(_COL_SNP)]  * n_pops),
        "--bp_col",     ",".join([str(_COL_POS)]  * n_pops),
        "--a1_col",     ",".join([str(_COL_A1)]   * n_pops),
        "--a2_col",     ",".join([str(_COL_A2)]   * n_pops),
        "--eff_col",    ",".join([str(_COL_BETA)] * n_pops),
        "--se_col",     ",".join([str(_COL_SE)]   * n_pops),
        "--pval_col",   ",".join([str(_COL_P)]    * n_pops),
        "--maf",        str(maf),
        "--n_sig",      str(n_sig),
        "--level",      str(level),
        "--min_purity", str(min_purity),
        "--pval_thresh", str(pval_thresh),
        "--max_iter",   str(max_iter),
        "--tol",        str(tol),
        "--threads",    "1",
        "--plink",      PLINK_BIN,
    ]

    print(f"NOTE: Running SuSiEx:\n  {' '.join(cmd)}")
    t_susiex = time.time()
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        elapsed = time.time() - t_susiex
        print(f"NOTE: SuSiEx completed in {elapsed:.1f}s")
        if result.stdout:
            print(result.stdout)
    except subprocess.CalledProcessError as ex:
        elapsed = time.time() - t_start
        print(f"ERROR: SuSiEx failed (return code {ex.returncode})")
        if ex.stdout:
            print(ex.stdout)
        if ex.stderr:
            print(ex.stderr)
        write_failure_output(output_file, ancestry_labels, timing_file, elapsed)
        # Cleanup
        for d in cohort_data:
            for f in [d['sst_file'], d['frq_file'], d['bim_file'], d['ld_bin']]:
                if os.path.exists(f):
                    os.remove(f)
        for f in sst_files_filtered:
            if os.path.exists(f):
                os.remove(f)
        sys.exit(0)

    # ── Parse SuSiEx output and write standard format ─────────────────────────
    snp_file = os.path.join(out_dir, f"{out_name}.snp")
    if not os.path.exists(snp_file):
        print(f"ERROR: SuSiEx .snp file not found: {snp_file}")
        write_failure_output(output_file, ancestry_labels, timing_file,
                             time.time() - t_start)
    else:
        success = parse_susiex_output(snp_file, ancestry_labels, cohort_data, output_file)
        if not success:
            write_failure_output(output_file, ancestry_labels, timing_file,
                                 time.time() - t_start)

    # ── Write timing file ──────────────────────────────────────────────────────
    total_elapsed = time.time() - t_start
    with open(timing_file, 'w') as fh:
        fh.write("susiex_time_seconds\n")
        fh.write(f"{total_elapsed:.2f}\n")

    # ── Cleanup temp files ─────────────────────────────────────────────────────
    for f in [
        os.path.join(out_dir, f"{out_name}.snp"),
        os.path.join(out_dir, f"{out_name}.cs"),
        os.path.join(out_dir, f"{out_name}.summary"),
    ]:
        if os.path.exists(f):
            os.remove(f)

    for d, ld_prefix in zip(cohort_data, ld_prefixes):
        for f in [d['sst_file'], d['frq_file'], d['bim_file'], d['ld_bin']]:
            if os.path.exists(f):
                os.remove(f)

    for f in sst_files_filtered:
        if os.path.exists(f):
            os.remove(f)

    print(f"NOTE: run_susiex.py complete ({total_elapsed:.1f}s total)")


if __name__ == "__main__":
    main()
