#!/usr/bin/env python3
"""
SuSiEx Multi-Population Fine-Mapping Script
Updated to auto-detect chromosome and BP range from *.vars files
and generate PLINK-format frequency files with A1 frequency

EXTENDED: Now also includes automatic conversion of text LD matrices to PLINK bin4 format for SuSiEx.

Usage:
 - Provide text-based LD matrices as LD prefix (no .ld.bin file)
 - Script will convert .txt LD matrix to .ld.bin format automatically
 - PLINK-format frequency files are auto-generated as before
 - After running SuSiEx, bin4 files can optionally be deleted (code provided, commented out)
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
        description='Multi-population fine-mapping using SuSiEx (auto chr/bp from LD matrix .vars, LD text-to-bin4 conversion added)',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--gwas_files', '-g', required=True,
                        help='Comma-delimited list of GWAS summary statistics files')
    parser.add_argument('--ld_files', '-l', required=True,
                        help='Comma-delimited list of LD matrix file prefixes (each must have a corresponding .vars file, .ld.bin or .txt allowed)')
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
    parser.add_argument('--susiex_bin', required=False, default='/lustre/orion/bif154/proj-shared/mconery/tools/SuSiEx/bin/SuSiEx',
                        help='Path to SuSiEx executable')
    return parser.parse_args()

def convert_text_ld_to_bin4(text_file, output_prefix):
    """
    Convert a text LD matrix to PLINK bin4 format (single-precision float32 row-major).
    Output is output_prefix.ld.bin. Returns True if successful, False otherwise.
    """
    print(f"[LD CONVERT] Reading LD matrix from: {text_file}")
    if not os.path.exists(text_file):
        print(f"[LD CONVERT] Error: Input file {text_file} not found!")
        return False
    try:
        ld_matrix = np.loadtxt(text_file, dtype=np.float64)
        print(f"[LD CONVERT] Successfully loaded matrix of size: {ld_matrix.shape}")
    except Exception as e:
        print(f"[LD CONVERT] Error reading file: {e}")
        return False
    if ld_matrix.shape[0] != ld_matrix.shape[1]:
        print(f"[LD CONVERT] Error: Matrix is not square! Shape: {ld_matrix.shape}")
        return False
    # Basic checks
    diag_values = np.diag(ld_matrix)
    if not np.allclose(diag_values, 1.0, atol=0.01):
        print("[LD CONVERT] Warning: Diagonal values are not close to 1.0")
    if ld_matrix.min() < -1.1 or ld_matrix.max() > 1.1:
        print("[LD CONVERT] Warning: LD values outside expected range [-1, 1]")
    ld_matrix_f32 = ld_matrix.astype(np.float32)
    output_bin = f"{output_prefix}.ld.bin"
    print(f"[LD CONVERT] Writing binary file: {output_bin}")
    try:
        with open(output_bin, 'wb') as f:
            f.write(ld_matrix_f32.tobytes('C'))
        print("[LD CONVERT] Successfully converted to bin4 format!")
        return True
    except Exception as e:
        print(f"[LD CONVERT] Error writing binary file: {e}")
        return False

def delete_bin4_files(ld_bin_files):
    """Delete bin4 LD matrix files."""
    print("\n=== Cleaning up bin4 LD files ===")
    for binf in ld_bin_files:
        if os.path.exists(binf):
            #os.remove(binf)  # Uncomment to actually delete files
            print(f"Would delete: {binf}")
        else:
            print(f"Not found (skip): {binf}")

def read_vars_file(vars_file):
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
        df = pd.read_csv(gwas_file, nrows=0, sep='\t')
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
            'missing': ['MissingRate', 'missing_rate', 'MISSING_RATE', 'missing', 'MISSING', 'miss', 'MISS', 'callrate', 'CALLRATE']
        }

        # Enhanced allele frequency detection - look for allele-specific frequencies
        af_patterns = {
            'af_a1': ['AF_Allele1', 'A1_FREQ', 'FRQ_A1', 'AF_A1', 'freq_a1', 'FREQ_A1', 'EAF', 'eaf'],
            'af_a2': ['AF_Allele2', 'A2_FREQ', 'FRQ_A2', 'AF_A2', 'freq_a2', 'FREQ_A2'],
            'af_generic': ['af', 'AF', 'allele_freq', 'ALLELE_FREQ', 'freq', 'FREQ', 'frequency', 'FREQUENCY']
        }

        detected_cols = {}

        # Detect standard columns
        for col_type, patterns in column_patterns.items():
            col_num = None
            for i, header in enumerate(headers):
                if header in patterns:
                    col_num = i + 1
                    break
            if col_num is None:
                defaults = {'chr': 1, 'snp': 2, 'bp': 3, 'a1': 4, 'a2': 5, 'eff': 6, 'se': 7, 'pval': 8, 'missing': None}
                col_num = defaults.get(col_type, None)
                if col_type == 'missing' and col_num is None:
                    print(f'    WARNING: Could not detect {col_type} column')
                elif col_type != 'missing':
                    print(f'    WARNING: Could not detect {col_type} column, using default {col_num}')
            detected_cols[col_type] = col_num

        # Detect allele frequency columns with priority for allele-specific ones
        af_col = None
        af_type = None

        # First, look for A1-specific frequency
        for i, header in enumerate(headers):
            if header in af_patterns['af_a1']:
                af_col = i + 1
                af_type = 'a1'
                break

        # If no A1 frequency found, look for A2-specific frequency
        if af_col is None:
            for i, header in enumerate(headers):
                if header in af_patterns['af_a2']:
                    af_col = i + 1
                    af_type = 'a2'
                    break

        # Finally, look for generic frequency columns
        if af_col is None:
            for i, header in enumerate(headers):
                if header in af_patterns['af_generic']:
                    af_col = i + 1
                    af_type = 'generic'
                    print(f'    WARNING: Found generic frequency column {header}, assuming it refers to A1')
                    break

        detected_cols['af'] = af_col
        detected_cols['af_type'] = af_type

        if af_col is None:
            print('    WARNING: Could not detect any allele frequency column')
        else:
            print(f'    Found allele frequency column: {headers[af_col-1]} (type: {af_type})')

        if verbose:
            print(f'  Detected columns: {detected_cols}')
        return detected_cols
    except Exception as e:
        print(f"ERROR: Failed to read GWAS file {gwas_file}: {e}")
        sys.exit(1)

def generate_frequency_file(gwas_file, vars_file, ld_prefix, sample_size, detected_cols):
    """
    Generate PLINK-format frequency file from GWAS summary statistics.
    The MAF column will contain the frequency of the A1 allele (not minor allele frequency).
    """
    print(f"  Generating frequency file for {os.path.basename(gwas_file)}")

    # Read variant IDs from .vars file
    with open(vars_file, 'r') as f:
        variant_ids = [line.strip() for line in f if line.strip()]

    # Read GWAS data
    try:
        gwas_df = pd.read_csv(gwas_file, sep='\t')
    except:
        try:
            gwas_df = pd.read_csv(gwas_file, sep=' ')
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
        print("ERROR: No overlapping variants found between GWAS and LD matrix")
        sys.exit(1)

    # Create frequency dataframe in PLINK .frq format
    freq_df = pd.DataFrame()
    freq_df['CHR'] = gwas_filtered.iloc[:, detected_cols['chr'] - 1]
    freq_df['SNP'] = gwas_filtered[snp_col]
    freq_df['A1'] = gwas_filtered[a1_col]
    freq_df['A2'] = gwas_filtered[a2_col]

    # Calculate A1 frequency (not MAF)
    if af_col is not None and af_col in gwas_filtered.columns:
        af_values = gwas_filtered[af_col].astype(float)
        af_type = detected_cols['af_type']

        if af_type == 'a1' or af_type == 'generic':
            # Frequency is for A1, use as-is
            freq_df['MAF'] = af_values
            print(f"    Using A1 frequency from column {af_col}")
        elif af_type == 'a2':
            # Frequency is for A2, convert to A1 frequency
            freq_df['MAF'] = 1 - af_values
            print(f"    Converting A2 frequency to A1 frequency from column {af_col}")
    else:
        print("    WARNING: No allele frequency column found, setting A1 frequency to 0.5")
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
    
    #Create bin4 LD files if needed
    bin4_created = []
    for ld_prefix in ld_files:
        ld_bin = f"{ld_prefix}.ld.bin"
        ld_txt = f"{ld_prefix}"
        if not os.path.exists(ld_bin):
            if os.path.exists(ld_txt):
                print(f"Converting {ld_txt} to {ld_bin} ...")
                if convert_text_ld_to_bin4(ld_txt, ld_prefix):
                    bin4_created.append(ld_bin)
                else:
                    print(f"Error converting {ld_txt} to bin4 format. Exiting.")
                    sys.exit(1)
            else:
                print(f"Missing both {ld_bin} and {ld_txt} for LD matrix prefix {ld_prefix}")
                sys.exit(1)
        else:
            print(f"Found existing bin4 LD matrix: {ld_bin}")
            bin4_created.append(ld_bin)
    
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
        #cleanup_frequency_files(frq_files)
        #delete_bin4_files()
        pass

if __name__ == "__main__":
    main()
