#!/usr/bin/env python3

"""
SuSiEx Multi-Population Fine-Mapping Script
Updated to auto-detect chromosome and BP range from *.vars files
and generate PLINK-format frequency files
"""

import argparse
import os
import sys
import subprocess
import pandas as pd
import numpy as np
import re

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Multi-population fine-mapping using SuSiEx (auto chr/bp from LD matrix .vars)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--gwas_files', '-g', required=True,
                       help='Comma-delimited list of GWAS summary statistics files')
    parser.add_argument('--ld_files', '-l', required=True,
                       help='Comma-delimited list of LD matrix file prefixes (each must have a corresponding .vars file)')
    parser.add_argument('--output_file', '-o', required=True,
                       help='Full path for final output file')
    parser.add_argument('--n_gwas', '-n', required=True,
                       help='Comma-delimited list of GWAS sample sizes, same order as GWAS files')
    parser.add_argument('--maf', type=float, default=0.005)
    parser.add_argument('--n_sig', type=int, default=10)
    parser.add_argument('--level', type=int, default=95)
    parser.add_argument('--min_purity', type=float, default=0.5)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--keep_ambig', action='store_true')
    parser.add_argument('--mult_step', action='store_true')
    parser.add_argument('--pval_thresh', type=float, default=1e-5)
    parser.add_argument('--max_iter', type=int, default=100)
    parser.add_argument('--tol', type=float, default=1e-4)
    parser.add_argument('--susiex_bin', required=False,
        default='/lustre/orion/bif154/proj-shared/mconery/tools/SuSiEx/bin/SuSiEx',
        help='Path to SuSiEx executable')
    return parser.parse_args()

def read_vars_file(vars_file):
    """
    Reads variant IDs from .vars file. Returns chromosome, min_bp, max_bp.
    """
    variant_ids = []
    with open(vars_file, 'r') as f:
        for line in f:
            v = line.strip()
            if v:
                variant_ids.append(v)
    chrom, positions = None, []
    for vid in variant_ids:
        m = re.match(r'^(\w+):(\d+):[^:]+:[^:]+$', vid)
        if m:
            c, pos = m.group(1), int(m.group(2))
            chrom = c if chrom is None else chrom
            positions.append(pos)
        else:
            print(f'WARNING: Variant ID does not match CHR:POS:REF:ALT: {vid}')
    if not positions or chrom is None:
        print(f'ERROR: Could not parse chromosome and positions from {vars_file}')
        sys.exit(1)
    return chrom, min(positions), max(positions), variant_ids

def detect_columns(gwas_file, verbose=True):
    try:
        df = pd.read_csv(gwas_file, nrows=0, sep=' ')
        headers = df.columns.tolist()
        if verbose:
            print(f"  Headers in {os.path.basename(gwas_file)}: {headers}")

        column_patterns = {
            'chr': ['CHR', 'chr', 'chromosome', '#CHROM', 'chrom'],
            'snp': ['SNP', 'snp', 'rsid', 'RSID', 'ID', 'id', 'variant_id', 'VARIANT_ID', 'MarkerName', 'rs_number', 'MarkerID'],
            'bp': ['BP', 'bp', 'pos', 'POS', 'position', 'POSITION', 'base_pair_location'],
            'a1': ['A1', 'a1', 'effect_allele', 'EFFECT_ALLELE', 'ALT', 'alt', 'ea', 'Allele1'],
            'a2': ['A2', 'a2', 'other_allele', 'OTHER_ALLELE', 'REF', 'ref', 'nea', 'Allele2'],
            'eff': ['BETA', 'beta', 'Beta', 'OR', 'or', 'effect', 'EFFECT', 'b'],
            'se': ['SE', 'se', 'standard_error', 'STANDARD_ERROR', 'stderr', 'STDERR'],
            'pval': ['P', 'p', 'pval', 'PVAL', 'p_value', 'P_VALUE', 'pvalue', 'PVALUE', 'p.value'],
            # Add allele frequency and missing rate patterns
            'af': ['AF_Allele2', 'af', 'AF', 'allele_freq', 'ALLELE_FREQ', 'freq', 'FREQ', 'EAF', 'eaf', 'A1_FREQ', 'FRQ_A1'],
            'missing': ['MissingRate', 'missing_rate', 'MISSING_RATE', 'missing', 'MISSING', 'miss', 'MISS', 'callrate', 'CALLRATE']
        }

        detected_cols = {}
        for col_type, patterns in column_patterns.items():
            col_num = None
            for i, header in enumerate(headers):
                if header in patterns:
                    col_num = i + 1
                    break
            if col_num is None:
                defaults = {'chr': 1, 'snp': 2, 'bp': 3, 'a1': 4, 'a2': 5, 'eff': 6, 'se': 7, 'pval': 8, 'af': None, 'missing': None}
                col_num = defaults.get(col_type, None)
                if col_type in ['af', 'missing'] and col_num is None:
                    print(f'    WARNING: Could not detect {col_type} column')
                elif col_type not in ['af', 'missing']:
                    print(f'    WARNING: Could not detect {col_type} column, using default {col_num}')
            detected_cols[col_type] = col_num

        if verbose:
            print(f'  Detected columns: {detected_cols}')
        return detected_cols
    except Exception as e:
        print(f"ERROR: Failed to read GWAS file {gwas_file}: {e}")
        sys.exit(1)

def generate_frequency_file(gwas_file, vars_file, ld_prefix, sample_size, detected_cols):
    """
    Generate PLINK-format frequency file from GWAS summary statistics.
    """
    print(f"  Generating frequency file for {os.path.basename(gwas_file)}")

    # Read variant IDs from .vars file
    with open(vars_file, 'r') as f:
        variant_ids = [line.strip() for line in f if line.strip()]

    # Read GWAS data
    try:
        gwas_df = pd.read_csv(gwas_file, sep=' ')
    except:
        try:
            gwas_df = pd.read_csv(gwas_file, sep='\t')
        except:
            print(f"ERROR: Could not read GWAS file {gwas_file}")
            sys.exit(1)

    # Get column names from detected columns (convert back to 0-indexed)
    snp_col = list(gwas_df.columns)[detected_cols['snp'] - 1]
    a1_col = list(gwas_df.columns)[detected_cols['a1'] - 1]
    a2_col = list(gwas_df.columns)[detected_cols['a2'] - 1]

    # Get AF and missing rate columns if available
    af_col = None
    missing_col = None

    if detected_cols['af'] is not None:
        af_col = list(gwas_df.columns)[detected_cols['af'] - 1]

    if detected_cols['missing'] is not None:
        missing_col = list(gwas_df.columns)[detected_cols['missing'] - 1]

    # Filter GWAS data to only variants in the LD matrix
    gwas_filtered = gwas_df[gwas_df[snp_col].isin(variant_ids)].copy()

    if len(gwas_filtered) == 0:
        print(f"ERROR: No overlapping variants found between GWAS and LD matrix")
        sys.exit(1)

    # Create frequency dataframe in PLINK .frq format
    freq_df = pd.DataFrame()
    freq_df['CHR'] = gwas_filtered.iloc[:, detected_cols['chr'] - 1]
    freq_df['SNP'] = gwas_filtered[snp_col]
    freq_df['A1'] = gwas_filtered[a1_col]
    freq_df['A2'] = gwas_filtered[a2_col]

    # Calculate MAF from AF_Allele2 if available
    if af_col is not None and af_col in gwas_filtered.columns:
        af_values = gwas_filtered[af_col].astype(float)
        # Convert to MAF (minor allele frequency)
        freq_df['MAF'] = np.where(af_values > 0.5, 1 - af_values, af_values)
    else:
        print(f"    WARNING: No allele frequency column found, setting MAF to 0.5")
        freq_df['MAF'] = 0.5

    # Calculate NCHROBS from missing rate if available
    if missing_col is not None and missing_col in gwas_filtered.columns:
        missing_values = gwas_filtered[missing_col].astype(float)

        # Check if missing values are rates (< 1) or counts
        if missing_values.max() <= 1.0:
            # Missing rates: NCHROBS = 2 * sample_size * (1 - missing_rate)
            freq_df['NCHROBS'] = np.round(2 * sample_size * (1 - missing_values)).astype(int)
        else:
            # Missing counts: NCHROBS = 2 * sample_size - missing_count
            freq_df['NCHROBS'] = np.round(2 * sample_size - missing_values).astype(int)
    else:
        print(f"    WARNING: No missing rate column found, setting NCHROBS to {2 * sample_size}")
        freq_df['NCHROBS'] = 2 * sample_size

    # Save frequency file
    frq_file = f"{ld_prefix}_frq.frq"
    freq_df.to_csv(frq_file, sep='\t', index=False)

    print(f"    Frequency file saved: {frq_file} ({len(freq_df)} variants)")
    return frq_file

def validate_inputs(args):
    print("=== Validating inputs ===")
    if not os.path.isfile(args.susiex_bin):
        print(f"ERROR: SuSiEx executable not found: {args.susiex_bin}")
        sys.exit(1)

    gwas_files = [f.strip() for f in args.gwas_files.split(',')]
    ld_files = [f.strip() for f in args.ld_files.split(',')]
    n_gwas_list = [int(n.strip()) for n in args.n_gwas.split(',')]
    n_pops = len(gwas_files)

    if len(ld_files) != n_pops or len(n_gwas_list) != n_pops:
        print("ERROR: Number of populations must match across all inputs")
        sys.exit(1)

    for i, gwas_file in enumerate(gwas_files):
        if not os.path.isfile(gwas_file):
            print(f"ERROR: GWAS file not found: {gwas_file}")
            sys.exit(1)
        print(f"Population {i+1} GWAS file: {gwas_file} ✓")

    for i, ld_prefix in enumerate(ld_files):
        vars_file = f"{ld_prefix}.vars"
        if not os.path.isfile(f"{ld_prefix}") or not os.path.isfile(vars_file):
            print(f"ERROR: LD matrix ({ld_prefix}) or vars file ({vars_file}) not found")
            sys.exit(1)
        print(f"Population {i+1} LD matrix: {ld_prefix}, .vars ✓")

    output_dir = os.path.dirname(args.output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    print(f"Validated input files for {n_pops} populations.")
    return gwas_files, ld_files, n_gwas_list, n_pops

def detect_all_columns(gwas_files):
    print("\n=== Detecting column numbers ===")
    all_cols = {'chr': [], 'snp': [], 'bp': [], 'a1': [], 'a2': [], 'eff': [], 'se': [], 'pval': []}
    all_detected_cols = []

    for i, gwas_file in enumerate(gwas_files):
        print(f"Population {i+1}:")
        detected_cols = detect_columns(gwas_file)
        all_detected_cols.append(detected_cols)

        for col_type in all_cols.keys():
            all_cols[col_type].append(str(detected_cols[col_type]))

    col_strings = {k: ','.join(v) for k,v in all_cols.items()}
    return col_strings, all_detected_cols

def get_chr_and_bp_range(ld_files):
    chrs, min_pos, max_pos = [], [], []
    all_variant_ids = []

    for ld_prefix in ld_files:
        vars_file = f"{ld_prefix}.vars"
        c, mn, mx, variant_ids = read_vars_file(vars_file)
        chrs.append(c)
        min_pos.append(mn)
        max_pos.append(mx)
        all_variant_ids.append(variant_ids)

    # For SuSiEx, require all chromosomes are the same
    if len(set(chrs)) > 1:
        print(f"ERROR: Chromosomes differ across populations: {chrs}")
        sys.exit(1)

    return chrs[0], str(min(min_pos)), str(max(max_pos)), all_variant_ids

def generate_all_frequency_files(gwas_files, ld_files, n_gwas_list, all_detected_cols):
    """
    Generate frequency files for all populations.
    """
    print("\n=== Generating frequency files ===")
    frq_files = []

    for i, (gwas_file, ld_prefix, sample_size, detected_cols) in enumerate(zip(gwas_files, ld_files, n_gwas_list, all_detected_cols)):
        print(f"Population {i+1}:")
        vars_file = f"{ld_prefix}.vars"
        frq_file = generate_frequency_file(gwas_file, vars_file, ld_prefix, sample_size, detected_cols)
        frq_files.append(frq_file)

    return frq_files

def run_susiex(args, gwas_files, ld_files, n_gwas_list, col_strings, chrom, bp_start, bp_end):
    print("\n=== Running SuSiEx fine-mapping ===")
    output_dir = os.path.dirname(args.output_file)
    output_name = os.path.basename(args.output_file).replace('.susiex','').replace('.txt','')
    bp_range = f"{bp_start},{bp_end}"

    cmd = [
        args.susiex_bin,
        '--sst_file', args.gwas_files,
        '--n_gwas', args.n_gwas,
        '--ld_file', args.ld_files,
        '--out_dir', output_dir if output_dir else '.',
        '--out_name', output_name,
        '--chr', chrom,
        '--bp', bp_range,
        '--chr_col', col_strings['chr'],
        '--snp_col', col_strings['snp'],
        '--bp_col', col_strings['bp'],
        '--a1_col', col_strings['a1'],
        '--a2_col', col_strings['a2'],
        '--eff_col', col_strings['eff'],
        '--se_col', col_strings['se'],
        '--pval_col', col_strings['pval'],
        '--maf', str(args.maf),
        '--n_sig', str(args.n_sig),
        '--level', str(args.level),
        '--min_purity', str(args.min_purity),
        '--pval_thresh', str(args.pval_thresh),
        '--max_iter', str(args.max_iter),
        '--tol', str(args.tol),
        '--threads', str(args.threads)
    ]

    if args.keep_ambig:
        cmd.append('--keep-ambig')
    if args.mult_step:
        cmd.append('--mult-step')

    print(' '.join(cmd))
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("SuSiEx completed successfully")
        if result.stdout:
            print("STDOUT:\n", result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: SuSiEx failed with return code {e.returncode}")
        if e.stdout:
            print("STDOUT:\n", e.stdout)
        if e.stderr:
            print("STDERR:\n", e.stderr)
        sys.exit(1)

def cleanup_frequency_files(frq_files):
    """
    Delete frequency files after SuSiEx completion.
    """
    print("\n=== Cleaning up frequency files ===")
    for frq_file in frq_files:
        if os.path.exists(frq_file):
            os.remove(frq_file)
            print(f"Deleted: {frq_file}")

def main():
    print("SuSiEx Multi-Population Fine-Mapping: Auto chr/bp from .vars")
    args = parse_arguments()

    # Validate inputs
    gwas_files, ld_files, n_gwas_list, n_pops = validate_inputs(args)

    # Detect columns
    col_strings, all_detected_cols = detect_all_columns(gwas_files)

    # Get chromosome and BP range
    chrom, bp_start, bp_end, all_variant_ids = get_chr_and_bp_range(ld_files)
    print(f'Chromosome: {chrom}, BP range: {bp_start}-{bp_end}')

    # Generate frequency files
    frq_files = generate_all_frequency_files(gwas_files, ld_files, n_gwas_list, all_detected_cols)

    try:
        # Run SuSiEx
        run_susiex(args, gwas_files, ld_files, n_gwas_list, col_strings, chrom, bp_start, bp_end)
    finally:
        # Clean up frequency files
        cleanup_frequency_files(frq_files)

if __name__ == "__main__":
    main()
