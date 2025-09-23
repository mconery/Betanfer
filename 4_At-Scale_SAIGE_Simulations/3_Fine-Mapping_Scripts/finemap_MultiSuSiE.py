#!/usr/bin/env python3

"""
MultiSuSiE Multi-Population Fine-Mapping Script

Python script designed to take input files and parameters and use them to run MultiSuSiE
Updated: Uses phenotype files to extract sample sizes and phenotype variances automatically
Added: Effect-size correlation parameter (rho) for off-diagonal values in rho matrix
Fixed: Ensures all populations have the same set of variants (intersection approach)
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
        description='Multi-population fine-mapping using MultiSuSiE (intersection of variants across all populations)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--gwas_files', '-g', required=True, help='Comma-delimited list of GWAS summary statistics files')
    parser.add_argument('--ld_files', '-l', required=True, help='Comma-delimited list of LD matrix file prefixes (each must have a corresponding .vars file)')
    parser.add_argument('--pheno_files', '-p', required=True, help='Comma-delimited list of phenotype files with PHENO column')
    parser.add_argument('--output_file', '-o', required=True, help='Full path for output txt file')
    parser.add_argument('--ancestry_labels', '-a', required=True, help='Comma-delimited list of ancestry labels')
    
    # MultiSuSiE parameters
    parser.add_argument('--maf', type=float, default=0.005, help='Minor allele frequency threshold (default: 0.005)')
    parser.add_argument('--L', type=int, default=10, help='Maximum number of causal effects (default: 10)')
    parser.add_argument('--coverage', type=float, default=0.95, help='Coverage for credible sets (default: 0.95)')
    parser.add_argument('--min_purity', type=float, default=0.5, help='Minimum purity for credible sets')
    parser.add_argument('--max_iter', type=int, default=100, help='Maximum number of iterations')
    parser.add_argument('--tol', type=float, default=1e-3, help='Convergence tolerance')
    parser.add_argument('--rho', type=float, default=0.75, help='Effect-size correlation between populations (off-diagonal values in rho matrix, default: 0.75)')
    
    # Boolean flags
    parser.add_argument('--estimate_residual_variance', action='store_true', help='Estimate residual variance')
    parser.add_argument('--standardize', action='store_true', default=True, help='Standardize inputs')
    
    # Optional parameters
    parser.add_argument('--prior_weights', default=None, help='Prior weights for variants')
    
    return parser.parse_args()

def read_vars_file(vars_file: str) -> Tuple[str, int, int, List[str]]:
    """Read variant IDs from .vars file and extract chromosome/position info."""
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
    """Detect column names in GWAS file."""
    try:
        df = pd.read_csv(gwas_file, nrows=0, sep='\t')
        headers = df.columns.tolist()
    except:
        df = pd.read_csv(gwas_file, nrows=0, sep=' ')
        headers = df.columns.tolist()
    column_patterns = {
        'chr': ['CHR', 'chr', 'chromosome', '#CHROM', 'chrom'],
        'snp': ['SNP', 'snp', 'rsid', 'RSID', 'ID', 'id', 'variant_id', 'VARIANT_ID', 'MarkerName', 'MarkerID'],
        'bp': ['BP', 'bp', 'pos', 'POS', 'position', 'POSITION', 'base_pair_location'],
        'a1': ['A1', 'a1', 'effect_allele', 'ALT', 'Allele1'],
        'a2': ['A2', 'a2', 'other_allele', 'REF', 'Allele2'],
        'beta': ['BETA', 'beta', 'Beta', 'OR', 'effect'],
        'se': ['SE', 'se', 'standard_error', 'stderr', 'STDERR'],
        'pval': ['P', 'p', 'pval', 'PVAL', 'p_value', "p.value"],
        'n': ['N', 'n', 'sample_size', 'SAMPLE_SIZE', 'nobs', 'NOBS']
    }
    detected_cols = {}
    for col_type, patterns in column_patterns.items():
        col_name = next((h for h in headers if h in patterns), None)
        detected_cols[col_type] = col_name
    return detected_cols

def load_phenotype_data(pheno_file: str, ancestry_label: str) -> Tuple[int, float]:
    """
    Load phenotype data and extract sample size and phenotype variance.
    
    Returns:
        n_samples: Number of samples (rows in phenotype file)
        var_y: Variance of the PHENO column
    """
    print(f"Loading phenotype data for {ancestry_label} from {os.path.basename(pheno_file)}")
    if not os.path.exists(pheno_file):
        raise FileNotFoundError(f"ERROR: Phenotype file not found: {pheno_file}")
    
    # Try different separators
    try:
        pheno_df = pd.read_csv(pheno_file, sep='\t')
    except:
        try:
            pheno_df = pd.read_csv(pheno_file, sep=' ')
        except Exception as e:
            raise ValueError(f"ERROR: Could not read phenotype file {pheno_file}: {e}")
    
    # Check for PHENO column
    if 'PHENO' not in pheno_df.columns:
        # Try other common phenotype column names
        pheno_cols = [col for col in pheno_df.columns if col.upper() in ['PHENO', 'PHENOTYPE', 'Y', 'VALUE']]
        if not pheno_cols:
            raise ValueError(f"ERROR: No PHENO column found in {pheno_file}. Available columns: {list(pheno_df.columns)}")
        pheno_col = pheno_cols[0]
        print(f"  Using column '{pheno_col}' as phenotype")
    else:
        pheno_col = 'PHENO'
    
    # Extract phenotype values (remove missing values)
    pheno_values = pheno_df[pheno_col].dropna()
    
    if len(pheno_values) == 0:
        raise ValueError(f"ERROR: No valid phenotype values found in {pheno_file}")
    
    n_samples = len(pheno_values)
    var_y = np.var(pheno_values, ddof=1)  # Sample variance
    
    print(f"  Sample size: {n_samples}")
    print(f"  Phenotype variance: {var_y:.6f}")
    print(f"  Phenotype mean: {np.mean(pheno_values):.6f}")
    print(f"  Phenotype std: {np.std(pheno_values, ddof=1):.6f}")
    
    return n_samples, var_y

def find_common_variants(ld_files: List[str], gwas_files: List[str], ancestry_labels: List[str]) -> List[str]:
    """
    Find intersection of variants across ALL populations (both LD matrices and GWAS files).
    This ensures all populations have exactly the same variants.
    """
    print("\n=== Finding common variants across all populations ===")
    
    # Get variants from all LD matrices
    all_ld_variants = {}
    for i, ld_file in enumerate(ld_files):
        vars_file = f"{ld_file}.vars"
        _, _, _, variant_ids = read_vars_file(vars_file)
        all_ld_variants[ancestry_labels[i]] = set(variant_ids)
        print(f"Population {ancestry_labels[i]} LD matrix: {len(variant_ids)} variants")
    
    # Get variants from all GWAS files
    all_gwas_variants = {}
    for i, gwas_file in enumerate(gwas_files):
        # Detect columns
        detected_cols = detect_columns(gwas_file)
        
        # Read GWAS file
        try:
            gwas_df = pd.read_csv(gwas_file, sep='\t')
        except:
            gwas_df = pd.read_csv(gwas_file, sep=' ')
        
        snp_col = detected_cols['snp']
        if snp_col is None or snp_col not in gwas_df.columns:
            raise ValueError(f"ERROR: Could not find SNP column in GWAS file: {gwas_file}")
        
        gwas_variants = set(gwas_df[snp_col].dropna().tolist())
        all_gwas_variants[ancestry_labels[i]] = gwas_variants
        print(f"Population {ancestry_labels[i]} GWAS: {len(gwas_variants)} variants")
    
    # Find intersection of all LD variants
    common_ld_variants = all_ld_variants[ancestry_labels[0]]
    for i in range(1, len(ancestry_labels)):
        common_ld_variants = common_ld_variants.intersection(all_ld_variants[ancestry_labels[i]])
    
    print(f"Common LD variants across all populations: {len(common_ld_variants)}")
    
    # Find intersection of all GWAS variants
    common_gwas_variants = all_gwas_variants[ancestry_labels[0]]
    for i in range(1, len(ancestry_labels)):
        common_gwas_variants = common_gwas_variants.intersection(all_gwas_variants[ancestry_labels[i]])
    
    print(f"Common GWAS variants across all populations: {len(common_gwas_variants)}")
    
    # Find final intersection of LD and GWAS variants
    final_common_variants = list(common_ld_variants.intersection(common_gwas_variants))
    print(f"Final common variants (LD ∩ GWAS): {len(final_common_variants)}")
    
    if len(final_common_variants) == 0:
        raise ValueError("ERROR: No variants are common across all LD matrices and GWAS files")
    
    # Sort variants to ensure consistent ordering
    final_common_variants.sort()
    
    return final_common_variants

def load_gwas_data_to_common_variants(gwas_file: str, detected_cols: Dict[str, str], common_variants: List[str], n_gwas: int) -> pd.DataFrame:
    """Load GWAS data filtered to match common variants across all populations."""
    try:
        gwas_df = pd.read_csv(gwas_file, sep='\t')
    except:
        gwas_df = pd.read_csv(gwas_file, sep=' ')
    snp_col = detected_cols['snp']
    if snp_col is None or snp_col not in gwas_df.columns:
        raise ValueError(f'SNP column not found in GWAS file: {gwas_file}')
    
    # Filter to common variants
    gwas_sub = gwas_df[gwas_df[snp_col].isin(common_variants)]
    if gwas_sub.empty:
        raise ValueError(f'No overlapping variants found in GWAS and common variant set')
    
    # Order by common_variants to ensure consistent ordering across populations
    gwas_sub = gwas_sub.set_index(snp_col).reindex(common_variants).reset_index()
    gwas_sub = gwas_sub.dropna(subset=[snp_col])  # Remove variants not in GWAS
    
    # Create standardized output dataframe
    result_df = pd.DataFrame({'SNP': gwas_sub[snp_col]})
    
    # Add other columns
    if detected_cols['a1']: result_df['A1'] = gwas_sub[detected_cols['a1']]
    if detected_cols['a2']: result_df['A2'] = gwas_sub[detected_cols['a2']]
    if detected_cols['beta']: result_df['Beta'] = gwas_sub[detected_cols['beta']]
    if detected_cols['se']: result_df['Se'] = gwas_sub[detected_cols['se']]
    if detected_cols['pval']: result_df['P'] = gwas_sub[detected_cols['pval']]
    
    # Calculate Z-scores
    if detected_cols['beta'] and detected_cols['se']:
        result_df['Z'] = result_df['Beta'] / result_df['Se']
    
    # Add sample size
    if detected_cols['n']:
        result_df['N'] = gwas_sub[detected_cols['n']]
    else:
        result_df['N'] = n_gwas
    
    return result_df

def load_ld_matrix_to_common_variants(ld_file: str, vars_file: str, common_variants: List[str]) -> np.ndarray:
    """Load LD matrix pruned to match common variants across all populations."""
    if not os.path.exists(ld_file):
        raise FileNotFoundError(f'LD matrix file not found: {ld_file}')
    try:
        ld_matrix = pd.read_csv(ld_file, header=None, sep='\t').values
    except:
        ld_matrix = pd.read_csv(ld_file, header=None, sep=' ').values
    
    _, _, _, all_variant_ids = read_vars_file(vars_file)
    if ld_matrix.shape[0] != len(all_variant_ids):
        raise ValueError('LD matrix size does not match its .vars list')
    
    # Find indices of common variants in this population's LD matrix
    common_indices = []
    for cv in common_variants:
        try:
            idx = all_variant_ids.index(cv)
            common_indices.append(idx)
        except ValueError:
            # This should not happen if find_common_variants worked correctly
            raise ValueError(f'Common variant {cv} not found in LD matrix .vars file')
    
    # Subset the LD matrix to common variants in the correct order
    ld_matrix_pruned = ld_matrix[np.ix_(common_indices, common_indices)]
    
    return ld_matrix_pruned

def create_rho_matrix(n_populations: int, rho_value: float) -> np.ndarray:
    """
    Create the effect-size correlation matrix (rho matrix) for MultiSuSiE.
    
    Args:
        n_populations: Number of populations
        rho_value: Correlation value for off-diagonal elements
        
    Returns:
        rho_matrix: n_populations x n_populations correlation matrix
    """
    rho_matrix = np.eye(n_populations)  # Start with identity matrix (1s on diagonal)
    
    # Fill off-diagonal elements with rho_value
    for i in range(n_populations):
        for j in range(n_populations):
            if i != j:
                rho_matrix[i, j] = rho_value
    
    print(f"Created rho matrix with off-diagonal correlation {rho_value}:")
    print(rho_matrix)
    
    return rho_matrix

def validate_inputs(args):
    """Validate input files and arguments."""
    gwas_files = [f.strip() for f in args.gwas_files.split(',')]
    ld_files = [f.strip() for f in args.ld_files.split(',')]
    pheno_files = [f.strip() for f in args.pheno_files.split(',')]
    ancestry_labels = [a.strip() for a in args.ancestry_labels.split(',')]
    
    if not (len(gwas_files) == len(ld_files) == len(pheno_files) == len(ancestry_labels)):
        raise ValueError('ERROR: Number of populations must match across GWAS files, LD files, phenotype files, and ancestry labels')
    
    if len(gwas_files) < 2:
        raise ValueError('ERROR: At least 2 populations are required for multi-ancestry analysis')
    
    # Validate rho parameter
    if not -1 <= args.rho <= 1:
        raise ValueError(f'ERROR: Effect-size correlation (rho) must be between -1 and 1, got {args.rho}')
    
    # Check file existence
    for i, (gwas_file, ld_file, pheno_file) in enumerate(zip(gwas_files, ld_files, pheno_files)):
        ancestry = ancestry_labels[i]
        
        if not os.path.exists(gwas_file):
            raise FileNotFoundError(f"ERROR: GWAS file not found for {ancestry}: {gwas_file}")
        
        vars_file = f"{ld_file}.vars"
        if not os.path.exists(ld_file) or not os.path.exists(vars_file):
            raise FileNotFoundError(f"ERROR: LD matrix ({ld_file}) or .vars file ({vars_file}) not found for {ancestry}")
        
        if not os.path.exists(pheno_file):
            raise FileNotFoundError(f"ERROR: Phenotype file not found for {ancestry}: {pheno_file}")
        
        print(f"Population {i+1} ({ancestry}) files validated ✓")
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    return gwas_files, ld_files, pheno_files, ancestry_labels

def run_multisusie(gwas_files, ld_files, pheno_files, ancestry_labels, args):
    """Run MultiSuSiE fine-mapping analysis."""
    print("\n=== Running MultiSuSiE fine-mapping ===")
    
    try:
        import MultiSuSiE
    except ImportError:
        raise ImportError("ERROR: Please install MultiSuSiE: pip install git+https://github.com/jordanero/MultiSuSiE.git")
    
    # CRITICAL: Find common variants across ALL populations first
    common_variants = find_common_variants(ld_files, gwas_files, ancestry_labels)
    
    b_list, s_list, R_list, summ_stat_list, varY_list, n_gwas_list = [], [], [], [], [], []
    
    for i in range(len(gwas_files)):
        ancestry = ancestry_labels[i]
        print(f"\nProcessing population {i+1} ({ancestry}):")
        
        # Load phenotype data to get sample size and variance
        n_samples, var_y = load_phenotype_data(pheno_files[i], ancestry)
        n_gwas_list.append(n_samples)
        varY_list.append(var_y)
        
        # Load and process GWAS data (filtered to common variants)
        detected_cols = detect_columns(gwas_files[i])
        gwas_df = load_gwas_data_to_common_variants(gwas_files[i], detected_cols, common_variants, n_samples)
        print(f"  GWAS data filtered to {len(gwas_df)} common variants")
        
        # Load LD matrix (pruned to common variants)
        vars_file = ld_files[i] + ".vars"
        ld_matrix = load_ld_matrix_to_common_variants(ld_files[i], vars_file, common_variants)
        print(f"  LD matrix pruned to: {ld_matrix.shape}")
        
        # Verify all populations have same number of variants
        if len(common_variants) != ld_matrix.shape[0] or len(common_variants) != len(gwas_df):
            raise ValueError(f"ERROR: Variant count mismatch for {ancestry}: common={len(common_variants)}, LD={ld_matrix.shape[0]}, GWAS={len(gwas_df)}")
        
        # Store data for MultiSuSiE
        summ_stat_list.append(gwas_df)
        b_list.append(gwas_df['Beta'].values)
        s_list.append(gwas_df['Se'].values)
        R_list.append(ld_matrix)
    
    # Create rho matrix for effect-size correlation between populations
    rho_matrix = create_rho_matrix(len(ancestry_labels), args.rho)
    
    print(f"\n=== Running MultiSuSiE analysis ===")
    print(f"Populations: {ancestry_labels}")
    print(f"Common variants: {len(common_variants)}")
    print(f"Sample sizes: {n_gwas_list}")
    print(f"Phenotype variances: {[f'{v:.6f}' for v in varY_list]}")
    print(f"Effect-size correlation (rho): {args.rho}")
    print(f"LD matrix shapes: {[R.shape for R in R_list]}")
    
    # Verify all LD matrices have same shape (should be guaranteed now)
    if not all(R.shape == R_list[0].shape for R in R_list):
        raise ValueError("ERROR: LD matrices have different shapes after filtering")
    
    # Parse prior weights if provided
    prior_weights = None
    if args.prior_weights is not None:
        if os.path.exists(args.prior_weights):
            prior_weights = np.loadtxt(args.prior_weights)
        else:
            prior_weights = np.array([float(w.strip()) for w in args.prior_weights.split(',')])
        print(f"Using prior weights: {prior_weights.shape if hasattr(prior_weights, 'shape') else len(prior_weights)} elements")
    
    # Run MultiSuSiE
    result = MultiSuSiE.multisusie_rss(
        b_list=b_list,
        s_list=s_list,
        R_list=R_list,
        varY_list=varY_list,
        population_sizes=n_gwas_list,
        rho=rho_matrix,  # Pass the rho matrix to MultiSuSiE
        L=args.L,
        max_iter=args.max_iter,
        tol=args.tol,
        estimate_residual_variance=args.estimate_residual_variance,
        standardize=args.standardize,
        prior_weights=prior_weights,
        coverage=args.coverage,
        min_abs_corr=args.min_purity
    )
    
    print("MultiSuSiE analysis completed successfully!")
    return result, summ_stat_list

def save_results(multisusie_result, summ_stat_list, ancestry_labels, output_file):
    """Save results in format compatible with other fine-mapping scripts."""
    print("\n=== Saving results ===")
    
    variant_ids = summ_stat_list[0]['SNP'].tolist()
    df = pd.DataFrame({'VARIANT_ID': variant_ids})
    
    # Add Z-score columns for each population
    for i, pop in enumerate(ancestry_labels):
        df[f'Z_{pop}'] = summ_stat_list[i]['Z'].values
    
    # Add overall PIP column
    df['PIP'] = getattr(multisusie_result, 'pip', np.zeros(len(df)))
    
    # Add convergence status
    df['converged'] = getattr(multisusie_result, 'converged', 'NA')
    
    # Add credible set assignments
    df['cs'] = 0
    if hasattr(multisusie_result, 'sets'):
        for i in range(len(multisusie_result.sets)):
            cslist=multisusie_result.sets[0][i]
            for indx in cslist:
                if indx < len(df):
                    df.loc[indx, 'cs'] = i + 1
    
    # Write output
    df.to_csv(output_file, sep='\t', index=False, na_rep='NA')
    
    print(f"Results saved to: {output_file}")
    print(f"Output contains {len(df)} variants and {len(df.columns)} columns")

def main():
    """Main function."""
    print("=== MultiSuSiE Multi-Population Fine-Mapping ===")
    
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Validate inputs
        gwas_files, ld_files, pheno_files, ancestry_labels = validate_inputs(args)
        
        # Run MultiSuSiE analysis
        multisusie_result, summ_stat_list = run_multisusie(gwas_files, ld_files, pheno_files, ancestry_labels, args)
        
        # Save results
        save_results(multisusie_result, summ_stat_list, ancestry_labels, args.output_file)
        
        print("\n=== Analysis completed successfully! ===")
        
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()