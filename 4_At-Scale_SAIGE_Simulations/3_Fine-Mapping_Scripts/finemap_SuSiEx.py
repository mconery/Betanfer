#!/usr/bin/env python3
"""
SuSiEx Multi-Population Fine-Mapping Script
Updated to auto-detect chromosome and BP range from *.vars files
"""
import argparse
import os
import sys
import subprocess
import pandas as pd
import re
from pathlib import Path

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
    return chrom, min(positions), max(positions)

def detect_columns(gwas_file, verbose=True):
    try:
        df = pd.read_csv(gwas_file, nrows=0, sep='	')
        headers = df.columns.tolist()
        if verbose:
            print(f"  Headers in {os.path.basename(gwas_file)}: {headers}")
        column_patterns = {
            'chr': ['CHR', 'chr', 'chromosome', '#CHROM', 'chrom'],
            'snp': ['SNP', 'snp', 'rsid', 'RSID', 'ID', 'id', 'variant_id', 'VARIANT_ID', 'MarkerName'],
            'bp': ['BP', 'bp', 'pos', 'POS', 'position', 'POSITION', 'base_pair_location'],
            'a1': ['A1', 'a1', 'effect_allele', 'EFFECT_ALLELE', 'ALT', 'alt', 'ea'],
            'a2': ['A2', 'a2', 'other_allele', 'OTHER_ALLELE', 'REF', 'ref', 'nea'],
            'eff': ['BETA', 'beta', 'Beta', 'OR', 'or', 'effect', 'EFFECT', 'b'],
            'se': ['SE', 'se', 'standard_error', 'STANDARD_ERROR', 'stderr', 'STDERR'],
            'pval': ['P', 'p', 'pval', 'PVAL', 'p_value', 'P_VALUE', 'pvalue', 'PVALUE']
        }
        detected_cols = {}
        for col_type, patterns in column_patterns.items():
            col_num = None
            for i, header in enumerate(headers):
                if header in patterns:
                    col_num = i + 1
                    break
            if col_num is None:
                defaults = {'chr': 1, 'snp': 2, 'bp': 3, 'a1': 4, 'a2': 5, 'eff': 6, 'se': 7, 'pval': 8}
                col_num = defaults.get(col_type, 1)
                print(f'    WARNING: Could not detect {col_type} column, using default {col_num}')
            detected_cols[col_type] = col_num
        if verbose:
            print(f'  Detected columns: {detected_cols}')
        return detected_cols
    except Exception as e:
        print(f"ERROR: Failed to read GWAS file {gwas_file}: {e}")
        sys.exit(1)

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
        print(f"ERROR: Number of populations must match across all inputs")
        sys.exit(1)
    for i, gwas_file in enumerate(gwas_files):
        if not os.path.isfile(gwas_file):
            print(f"ERROR: GWAS file not found: {gwas_file}")
            sys.exit(1)
        print(f"Population {i+1} GWAS file: {gwas_file} ✓")
    for i, ld_prefix in enumerate(ld_files):
        vars_file = f"{ld_prefix}.vars"
        if not os.path.isfile(f"{ld_prefix}.npy") or not os.path.isfile(vars_file):
            print(f"ERROR: LD matrix ({ld_prefix}.npy) or vars file ({vars_file}) not found")
            sys.exit(1)
        print(f"Population {i+1} LD matrix: {ld_prefix}.npy, .vars ✓")
    output_dir = os.path.dirname(args.output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    print(f"Validated input files for {n_pops} populations.")
    return gwas_files, ld_files, n_gwas_list, n_pops

def detect_all_columns(gwas_files):
    print("=== Detecting column numbers ===")
    all_cols = {'chr': [], 'snp': [], 'bp': [], 'a1': [], 'a2': [], 'eff': [], 'se': [], 'pval': []}
    for i, gwas_file in enumerate(gwas_files):
        print(f"Population {i+1}:")
        detected_cols = detect_columns(gwas_file)
        for col_type in all_cols.keys():
            all_cols[col_type].append(str(detected_cols[col_type]))
    col_strings = {k: ','.join(v) for k,v in all_cols.items()}
    return col_strings

def get_chr_and_bp_range(ld_files):
    chrs, min_pos, max_pos = [], [], []
    for ld_prefix in ld_files:
        vars_file = f"{ld_prefix}.vars"
        c, mn, mx = read_vars_file(vars_file)
        chrs.append(c)
        min_pos.append(mn)
        max_pos.append(mx)
    # For SuSiEx, require all chromosomes are the same
    if len(set(chrs)) > 1:
        print(f"ERROR: Chromosomes differ across populations: {chrs}")
        sys.exit(1)
    return chrs[0], str(min(min_pos)), str(max(max_pos))

def run_susiex(args, gwas_files, ld_files, n_gwas_list, col_strings, chrom, bp_start, bp_end):
    print("=== Running SuSiEx fine-mapping ===")
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
            print("STDOUT:", result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: SuSiEx failed with return code {e.returncode}")
        if e.stdout:
            print("STDOUT:", e.stdout)
        if e.stderr:
            print("STDERR:", e.stderr)
        sys.exit(1)

def main():
    print("SuSiEx Multi-Population Fine-Mapping: Auto chr/bp from .vars")
    args = parse_arguments()
    gwas_files, ld_files, n_gwas_list, n_pops = validate_inputs(args)
    col_strings = detect_all_columns(gwas_files)
    chrom, bp_start, bp_end = get_chr_and_bp_range(ld_files)
    print(f'Chromosome: {chrom}, BP range: {bp_start}-{bp_end}')
    run_susiex(args, gwas_files, ld_files, n_gwas_list, col_strings, chrom, bp_start, bp_end)
if __name__ == "__main__":
    main()
