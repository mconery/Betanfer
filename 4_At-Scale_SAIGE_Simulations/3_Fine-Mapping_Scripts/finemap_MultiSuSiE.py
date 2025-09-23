#!/usr/bin/env python3

"""
MultiSuSiE Multi-Population Fine-Mapping Script

Python script designed to take input files and parameters and use them to run MultiSuSiE
Updated: Subsets the summary statistics for each population to only match its own LD variant list (from the .vars file), NOT the intersection across all populations.
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import re
from typing import List, Dict, Tuple, Any

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Multi-population fine-mapping using MultiSuSiE (each population subset only to its LD matrix variants)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--gwas_files', '-g', required=True, help='Comma-delimited list of GWAS summary statistics files')
    parser.add_argument('--ld_files', '-l', required=True, help='Comma-delimited list of LD matrix file prefixes (each must have a corresponding .vars file)')
    parser.add_argument('--output_file', '-o', required=True, help='Full path for output txt file')
    parser.add_argument('--n_gwas', '-n', required=True, help='Comma-delimited list of GWAS sample sizes, same order as GWAS files')
    parser.add_argument('--ancestry_labels', '-a', required=True, help='Comma-delimited list of ancestry labels')
    parser.add_argument('--maf', type=float, default=0.005, help='Minor allele frequency threshold (default: 0.005)')
    parser.add_argument('--L', type=int, default=10, help='Maximum number of causal effects (default: 10)')
    parser.add_argument('--coverage', type=float, default=0.95, help='Coverage for credible sets (default: 0.95)')
    parser.add_argument('--min_purity', type=float, default=0.5, help='Minimum purity for credible sets')
    parser.add_argument('--max_iter', type=int, default=100, help='Maximum number of iterations')
    parser.add_argument('--tol', type=float, default=1e-3, help='Convergence tolerance')
    parser.add_argument('--estimate_residual_variance', action='store_true', help='Estimate residual variance')
    parser.add_argument('--standardize', action='store_true', default=True, help='Standardize inputs')
    parser.add_argument('--prior_weights', default=None, help='Prior weights for variants')
    parser.add_argument('--var_y', default=None, help='Comma-delimited phenotype variances (optional)')
    return parser.parse_args()

def read_vars_file(vars_file: str) -> Tuple[str, int, int, List[str]]:
    if not os.path.exists(vars_file):
        raise FileNotFoundError(f"ERROR: .vars file not found: {vars_file}")
    with open(vars_file) as f:
        variant_ids = [line.strip() for line in f if line.strip()]
    chrom, positions = None, []
    for vid in variant_ids:
        m = re.match(r'^(\w+):(\d+):[^:]+:[^:]+$', vid)
        if m:
            c, pos = m.group(1), int(m.group(2))
            if chrom is None: chrom = c
            positions.append(pos)
    if not positions or chrom is None:
        raise ValueError(f'Could not parse chromosome and positions from {vars_file}')
    return chrom, min(positions), max(positions), variant_ids

def detect_columns(gwas_file: str) -> Dict[str, str]:
    try:
        df = pd.read_csv(gwas_file, nrows=0, sep='\t')
        headers = df.columns.tolist()
    except:
        df = pd.read_csv(gwas_file, nrows=0, sep=' ')
        headers = df.columns.tolist()
    column_patterns = {
        'chr': ['CHR', 'chr', 'chromosome', '#CHROM', 'chrom'],
        'snp': ['SNP', 'snp', 'rsid', 'RSID', 'ID', 'id', 'variant_id', 'VARIANT_ID', 'MarkerName'],
        'bp': ['BP', 'bp', 'pos', 'POS', 'position', 'POSITION', 'base_pair_location'],
        'a1': ['A1', 'a1', 'effect_allele', 'ALT', 'Allele1'],
        'a2': ['A2', 'a2', 'other_allele', 'REF', 'Allele2'],
        'beta': ['BETA', 'beta', 'Beta', 'OR', 'effect'],
        'se': ['SE', 'se', 'standard_error', 'stderr', 'STDERR'],
        'pval': ['P', 'p', 'pval', 'PVAL', 'p_value'],
        'n': ['N', 'n', 'sample_size', 'SAMPLE_SIZE', 'nobs', 'NOBS']
    }
    detected_cols = {}
    for col_type, patterns in column_patterns.items():
        col_name = next((h for h in headers if h in patterns), None)
        detected_cols[col_type] = col_name
    return detected_cols

def load_gwas_data_to_ld_variants(gwas_file: str, detected_cols: Dict[str, str], ld_variant_ids: List[str], n_gwas: int) -> pd.DataFrame:
    try:
        gwas_df = pd.read_csv(gwas_file, sep='\t')
    except:
        gwas_df = pd.read_csv(gwas_file, sep=' ')
    snp_col = detected_cols['snp']
    if snp_col is None or snp_col not in gwas_df.columns:
        raise ValueError(f'SNP column not found in GWAS file: {gwas_file}')
    gwas_sub = gwas_df[gwas_df[snp_col].isin(ld_variant_ids)]
    if gwas_sub.empty:
        raise ValueError(f'No overlapping variants found in GWAS and LD set')
    gwas_sub = gwas_sub.set_index(snp_col).reindex(ld_variant_ids).reset_index()
    result_df = pd.DataFrame({'SNP': gwas_sub[snp_col]})
    if detected_cols['a1']: result_df['A1'] = gwas_sub[detected_cols['a1']]
    if detected_cols['a2']: result_df['A2'] = gwas_sub[detected_cols['a2']]
    if detected_cols['beta']: result_df['Beta'] = gwas_sub[detected_cols['beta']]
    if detected_cols['se']: result_df['Se'] = gwas_sub[detected_cols['se']]
    if detected_cols['pval']: result_df['P'] = gwas_sub[detected_cols['pval']]
    if detected_cols['beta'] and detected_cols['se']:
        result_df['Z'] = result_df['Beta'] / result_df['Se']
    if detected_cols['n']:
        result_df['N'] = gwas_sub[detected_cols['n']]
    else:
        result_df['N'] = n_gwas
    return result_df

def load_ld_matrix(ld_file: str, vars_file: str, variant_ids: List[str]) -> np.ndarray:
    if not os.path.exists(ld_file):
        raise FileNotFoundError(f'LD matrix file not found: {ld_file}')
    try:
        ld_matrix = pd.read_csv(ld_file, header=None, sep='\t').values
    except:
        ld_matrix = pd.read_csv(ld_file, header=None, sep=' ').values
    _, _, _, all_variant_ids = read_vars_file(vars_file)
    if ld_matrix.shape[0] != len(all_variant_ids):
        raise ValueError('LD matrix size does not match its .vars list')
    # Prune to just variant_ids (which is all_variant_ids here)
    return ld_matrix

def validate_inputs(args):
    gwas_files = [f.strip() for f in args.gwas_files.split(',')]
    ld_files = [f.strip() for f in args.ld_files.split(',')]
    n_gwas_list = [int(n.strip()) for n in args.n_gwas.split(',')]
    ancestry_labels = [a.strip() for a in args.ancestry_labels.split(',')]
    if not (len(gwas_files) == len(ld_files) == len(n_gwas_list) == len(ancestry_labels)):
        raise ValueError('Mismatch between number of populations across input lists.')
    return gwas_files, ld_files, n_gwas_list, ancestry_labels

def run_multisusie(gwas_files, ld_files, n_gwas_list, ancestry_labels, args):
    try:
        import MultiSuSiE
    except ImportError:
        raise ImportError("Please install MultiSuSiE: pip install git+https://github.com/jordanero/MultiSuSiE.git")
    b_list, s_list, R_list, summ_stat_list, varY_list = [], [], [], [], []
    for i in range(len(gwas_files)):
        vars_file = ld_files[i] + ".vars"
        _, _, _, ld_variant_ids = read_vars_file(vars_file)
        detected_cols = detect_columns(gwas_files[i])
        gwas_df = load_gwas_data_to_ld_variants(gwas_files[i], detected_cols, ld_variant_ids, n_gwas_list[i])
        ld_matrix = load_ld_matrix(ld_files[i], vars_file, ld_variant_ids)
        summ_stat_list.append(gwas_df)
        b_list.append(gwas_df['Beta'].values)
        s_list.append(gwas_df['Se'].values)
        R_list.append(ld_matrix)
        varY_list.append(1.0 if args.var_y is None else float(args.var_y.split(',')[i]))
    result = MultiSuSiE.multisusie_rss(
        b_list=b_list,s_list=s_list,R_list=R_list,varY_list=varY_list,
        population_sizes=n_gwas_list,L=args.L,max_iter=args.max_iter,
        tol=args.tol,estimate_residual_variance=args.estimate_residual_variance,
        standardize=args.standardize,prior_weights=args.prior_weights,coverage=args.coverage)
    return result, summ_stat_list

def save_results(multisusie_result, summ_stat_list, ancestry_labels, output_file):
    variant_ids = summ_stat_list[0]['SNP'].tolist()
    df = pd.DataFrame({'VARIANT_ID': variant_ids})
    for i, pop in enumerate(ancestry_labels):
        df[f'Z_{pop}'] = summ_stat_list[i]['Z'].values
    df['PIP'] = getattr(multisusie_result, 'pip', np.zeros(len(df)))
    df['converged'] = getattr(multisusie_result, 'converged', 'NA')
    df['cs'] = 0
    if hasattr(multisusie_result, 'sets') and hasattr(multisusie_result.sets, 'cs'):
        for i, cslist in enumerate(multisusie_result.sets.cs, 1):
            for indx in cslist:
                if indx < len(df):
                    df.loc[indx, 'cs'] = i
    for i, pop in enumerate(ancestry_labels):
        df[f'PIP_{pop}'] = 0  # Placeholder, MultiSuSiE-per-pop PIPs if available
    if len(ancestry_labels) >= 2:
        df['PIP_' + '_'.join(ancestry_labels)] = 0  # Placeholder
    df.to_csv(output_file, sep='\t', index=False, na_rep='NA')

def main():
    args = parse_arguments()
    gwas_files, ld_files, n_gwas_list, ancestry_labels = validate_inputs(args)
    multisusie_result, summ_stat_list = run_multisusie(gwas_files, ld_files, n_gwas_list, ancestry_labels, args)
    save_results(multisusie_result, summ_stat_list, ancestry_labels, args.output_file)
    print(f"Results written to {args.output_file}")

if __name__ == '__main__':
    main()
