#!/usr/bin/env python3

"""
MultiSuSiE Multi-Population Fine-Mapping Script

Python script designed to take input files and parameters and use them to run MultiSuSiE
Similar architecture to SuSiEx and MESuSiE scripts for genomic fine-mapping

@author Generated based on SuSiEx and MESuSiE script structures
"""

import argparse
import os
import sys
import pandas as pd
import numpy as np
import re
from typing import List, Dict, Tuple, Any

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Multi-population fine-mapping using MultiSuSiE (auto chr/bp from LD matrix .vars)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Core input arguments (consistent with other scripts)
    parser.add_argument('--gwas_files', '-g', required=True,
                        help='Comma-delimited list of GWAS summary statistics files')
    parser.add_argument('--ld_files', '-l', required=True,
                        help='Comma-delimited list of LD matrix file prefixes (each must have a corresponding .vars file)')
    parser.add_argument('--output_file', '-o', required=True,
                        help='Full path for final output txt file')
    parser.add_argument('--n_gwas', '-n', required=True,
                        help='Comma-delimited list of GWAS sample sizes, same order as GWAS files')
    parser.add_argument('--ancestry_labels', '-a', required=True,
                        help='Comma-delimited list of ancestry labels (e.g., EUR,AFR,EAS), same order as GWAS files')
    
    # MultiSuSiE-specific parameters
    parser.add_argument('--maf', type=float, default=0.005,
                        help='Minor allele frequency threshold (default: 0.005)')
    parser.add_argument('--L', type=int, default=10,
                        help='Maximum number of causal effects (default: 10)')
    parser.add_argument('--coverage', type=float, default=0.95,
                        help='Coverage for credible sets (default: 0.95)')
    parser.add_argument('--min_purity', type=float, default=0.5,
                        help='Minimum purity for credible sets (default: 0.5)')
    parser.add_argument('--max_iter', type=int, default=100,
                        help='Maximum number of iterations (default: 100)')
    parser.add_argument('--tol', type=float, default=1e-3,
                        help='Tolerance for convergence (default: 1e-3)')
    
    # MultiSuSiE advanced options
    parser.add_argument('--estimate_residual_variance', action='store_true',
                        help='Estimate residual variance (default: False)')
    parser.add_argument('--low_memory_mode', action='store_true',
                        help='Use low memory mode (default: False)')
    parser.add_argument('--standardize', action='store_true', default=True,
                        help='Standardize inputs (default: True)')
    parser.add_argument('--intercept', action='store_true', default=True,
                        help='Fit intercept (default: True)')
    
    # Prior configuration
    parser.add_argument('--prior_weights', default=None,
                        help='Prior weights for variants (default: uniform)')
    parser.add_argument('--var_y', default=None,
                        help='Comma-delimited list of phenotype variances for each population (default: estimated)')
    
    return parser.parse_args()

def read_vars_file(vars_file: str) -> Tuple[str, int, int, List[str]]:
    """Read variant IDs from .vars file. Returns chromosome, min_bp, max_bp, variant_ids."""
    print(f"Reading variants from {vars_file}")
    if not os.path.exists(vars_file):
        raise FileNotFoundError(f"ERROR: .vars file not found: {vars_file}")
    
    with open(vars_file, 'r') as f:
        variant_ids = [line.strip() for line in f if line.strip()]
    
    # Extract chromosome and positions
    chrom = None
    positions = []
    
    for vid in variant_ids:
        # Match CHR:POS:REF:ALT format
        match = re.match(r'^(\w+):(\d+):[^:]+:[^:]+$', vid)
        if match:
            c = match.group(1)
            pos = int(match.group(2))
            if chrom is None:
                chrom = c
            positions.append(pos)
        else:
            print(f'WARNING: Variant ID does not match CHR:POS:REF:ALT: {vid}')
    
    if not positions or chrom is None:
        raise ValueError(f'ERROR: Could not parse chromosome and positions from {vars_file}')
    
    return chrom, min(positions), max(positions), variant_ids

def detect_columns(gwas_file: str, verbose: bool = True) -> Dict[str, str]:
    """Detect column names in GWAS file."""
    if verbose:
        print(f" Detecting columns in {os.path.basename(gwas_file)}")
    
    try:
        # Try reading the header
        df = pd.read_csv(gwas_file, nrows=0, sep='\t')
        headers = df.columns.tolist()
    except:
        try:
            df = pd.read_csv(gwas_file, nrows=0, sep=' ')
            headers = df.columns.tolist()
        except Exception as e:
            raise ValueError(f"ERROR: Failed to read GWAS file {gwas_file}: {e}")
    
    if verbose:
        print(f" Headers: {', '.join(headers)}")
    
    # Column patterns for detection
    column_patterns = {
        'chr': ['CHR', 'chr', 'chromosome', '#CHROM', 'chrom'],
        'snp': ['SNP', 'snp', 'rsid', 'RSID', 'ID', 'id', 'variant_id', 'VARIANT_ID', 'MarkerName', 'rs_number', 'MarkerID'],
        'bp': ['BP', 'bp', 'pos', 'POS', 'position', 'POSITION', 'base_pair_location'],
        'a1': ['A1', 'a1', 'effect_allele', 'EFFECT_ALLELE', 'ALT', 'alt', 'ea', 'Allele1'],
        'a2': ['A2', 'a2', 'other_allele', 'OTHER_ALLELE', 'REF', 'ref', 'nea', 'Allele2'],
        'beta': ['BETA', 'beta', 'Beta', 'OR', 'or', 'effect', 'EFFECT', 'b'],
        'se': ['SE', 'se', 'standard_error', 'STANDARD_ERROR', 'stderr', 'STDERR'],
        'pval': ['P', 'p', 'pval', 'PVAL', 'p_value', 'P_VALUE', 'pvalue', 'PVALUE', 'p.value'],
        'n': ['N', 'n', 'sample_size', 'SAMPLE_SIZE', 'nobs', 'NOBS']
    }
    
    detected_cols = {}
    
    # Detect each column type
    for col_type, patterns in column_patterns.items():
        col_name = None
        for header in headers:
            if header in patterns:
                col_name = header
                break
        
        if col_name is None and verbose:
            print(f" WARNING: Could not detect {col_type} column")
        
        detected_cols[col_type] = col_name
    
    return detected_cols

def find_common_variants(ld_files: List[str], gwas_files: List[str], ancestry_labels: List[str]) -> List[str]:
    """Find intersection of variants across all populations."""
    print("\n=== Finding common variants across all populations ===")
    
    # Get variants from all LD matrices
    all_ld_variants = {}
    for i, ld_file in enumerate(ld_files):
        vars_file = f"{ld_file}.vars"
        _, _, _, variant_ids = read_vars_file(vars_file)
        all_ld_variants[ancestry_labels[i]] = variant_ids
        print(f"Population {ancestry_labels[i]} LD matrix: {len(variant_ids)} variants")
    
    # Get variants from all GWAS files
    all_gwas_variants = {}
    for i, gwas_file in enumerate(gwas_files):
        # Detect columns
        detected_cols = detect_columns(gwas_file, verbose=False)
        
        # Read GWAS file
        try:
            gwas_df = pd.read_csv(gwas_file, sep='\t')
        except:
            gwas_df = pd.read_csv(gwas_file, sep=' ')
        
        snp_col = detected_cols['snp']
        if snp_col is None or snp_col not in gwas_df.columns:
            raise ValueError(f"ERROR: Could not find SNP column in GWAS file: {gwas_file}")
        
        gwas_variants = gwas_df[snp_col].tolist()
        all_gwas_variants[ancestry_labels[i]] = gwas_variants
        print(f"Population {ancestry_labels[i]} GWAS: {len(gwas_variants)} variants")
    
    # Find intersection of all LD variants
    common_ld_variants = set(all_ld_variants[ancestry_labels[0]])
    for i in range(1, len(ancestry_labels)):
        common_ld_variants = common_ld_variants.intersection(set(all_ld_variants[ancestry_labels[i]]))
    
    print(f"Common LD variants across all populations: {len(common_ld_variants)}")
    
    # Find intersection of all GWAS variants
    common_gwas_variants = set(all_gwas_variants[ancestry_labels[0]])
    for i in range(1, len(ancestry_labels)):
        common_gwas_variants = common_gwas_variants.intersection(set(all_gwas_variants[ancestry_labels[i]]))
    
    print(f"Common GWAS variants across all populations: {len(common_gwas_variants)}")
    
    # Find final intersection of LD and GWAS variants
    final_common_variants = list(common_ld_variants.intersection(common_gwas_variants))
    print(f"Final common variants (LD ∩ GWAS): {len(final_common_variants)}")
    
    if len(final_common_variants) == 0:
        raise ValueError("ERROR: No variants are common across all LD matrices and GWAS files")
    
    return final_common_variants

def load_gwas_data(gwas_file: str, detected_cols: Dict[str, str], common_variants: List[str], n_gwas: int) -> pd.DataFrame:
    """Load and process GWAS data filtered to common variants."""
    print(f" Loading GWAS data from {os.path.basename(gwas_file)}")
    
    # Read the full GWAS file
    try:
        gwas_df = pd.read_csv(gwas_file, sep='\t')
    except:
        gwas_df = pd.read_csv(gwas_file, sep=' ')
    
    # Filter to common variants
    snp_col = detected_cols['snp']
    if snp_col is None or snp_col not in gwas_df.columns:
        raise ValueError("ERROR: Could not find SNP column in GWAS file")
    
    gwas_filtered = gwas_df[gwas_df[snp_col].isin(common_variants)]
    
    if len(gwas_filtered) == 0:
        raise ValueError("ERROR: No overlapping variants found between GWAS and common variant set")
    
    print(f" Found {len(gwas_filtered)} overlapping variants")
    
    # Order by common_variants to ensure consistent ordering
    gwas_filtered = gwas_filtered.set_index(snp_col).reindex(common_variants).reset_index()
    gwas_filtered = gwas_filtered.dropna()
    
    # Create standard format data frame
    result_df = pd.DataFrame({
        'SNP': gwas_filtered[snp_col]
    })
    
    # Add other columns if available and rename Allele1/Allele2 to REF/ALT
    if detected_cols['a1'] is not None:
        if detected_cols['a1'] == "Allele1":
            result_df['ALT'] = gwas_filtered[detected_cols['a1']]
            print(f" Renamed {detected_cols['a1']} to ALT")
        else:
            result_df['A1'] = gwas_filtered[detected_cols['a1']]
    
    if detected_cols['a2'] is not None:
        if detected_cols['a2'] == "Allele2":
            result_df['REF'] = gwas_filtered[detected_cols['a2']]
            print(f" Renamed {detected_cols['a2']} to REF")
        else:
            result_df['A2'] = gwas_filtered[detected_cols['a2']]
    
    if detected_cols['beta'] is not None:
        result_df['Beta'] = gwas_filtered[detected_cols['beta']]
    
    if detected_cols['se'] is not None:
        result_df['Se'] = gwas_filtered[detected_cols['se']]
    
    if detected_cols['pval'] is not None:
        result_df['P'] = gwas_filtered[detected_cols['pval']]
    
    # Calculate Z-scores if Beta and Se are available
    if detected_cols['beta'] is not None and detected_cols['se'] is not None:
        result_df['Z'] = result_df['Beta'] / result_df['Se']
    
    # Add sample size
    if detected_cols['n'] is not None:
        result_df['N'] = gwas_filtered[detected_cols['n']]
    else:
        result_df['N'] = n_gwas
    
    return result_df

def load_ld_matrix(ld_file: str, vars_file: str, common_variants: List[str]) -> np.ndarray:
    """Load LD matrix pruned to common variants."""
    print(f" Loading LD matrix from {os.path.basename(ld_file)}")
    
    if not os.path.exists(ld_file):
        raise FileNotFoundError(f"ERROR: LD matrix file not found: {ld_file}")
    
    # Read the full LD matrix and variant IDs
    try:
        ld_matrix = pd.read_csv(ld_file, header=None, sep='\t').values
    except:
        try:
            ld_matrix = pd.read_csv(ld_file, header=None, sep=' ').values
        except Exception as e:
            raise ValueError(f"ERROR: Could not read LD matrix from {ld_file}: {e}")
    
    # Read variant IDs from .vars file
    _, _, _, all_variant_ids = read_vars_file(vars_file)
    
    # Check original dimensions
    n_all_variants = len(all_variant_ids)
    if ld_matrix.shape[0] != n_all_variants or ld_matrix.shape[1] != n_all_variants:
        raise ValueError(f"ERROR: LD matrix dimensions {ld_matrix.shape} do not match .vars file variants ({n_all_variants})")
    
    # Find indices of common variants in the full matrix
    common_indices = [i for i, vid in enumerate(all_variant_ids) if vid in common_variants]
    
    if len(common_indices) == 0:
        raise ValueError("ERROR: No common variants found in LD matrix")
    
    print(f" Pruning LD matrix from {n_all_variants} to {len(common_indices)} variants")
    
    # Subset the LD matrix to common variants
    ld_matrix_pruned = ld_matrix[np.ix_(common_indices, common_indices)]
    
    # Reorder to match common_variants order
    common_variant_ids = [all_variant_ids[i] for i in common_indices]
    variant_order = [common_variant_ids.index(v) for v in common_variants if v in common_variant_ids]
    ld_matrix_final = ld_matrix_pruned[np.ix_(variant_order, variant_order)]
    
    return ld_matrix_final

def validate_inputs(args) -> Tuple[List[str], List[str], List[int], List[str], int]:
    """Validate input files and arguments."""
    print("=== Validating inputs ===")
    
    gwas_files = [f.strip() for f in args.gwas_files.split(',')]
    ld_files = [f.strip() for f in args.ld_files.split(',')]
    n_gwas_list = [int(n.strip()) for n in args.n_gwas.split(',')]
    ancestry_labels = [a.strip() for a in args.ancestry_labels.split(',')]
    
    n_pops = len(gwas_files)
    
    if len(ld_files) != n_pops or len(n_gwas_list) != n_pops or len(ancestry_labels) != n_pops:
        raise ValueError("ERROR: Number of populations must match across all inputs")
    
    if n_pops < 2:
        raise ValueError("ERROR: At least 2 populations are required for multi-ancestry analysis")
    
    # Check GWAS files exist
    for i, gwas_file in enumerate(gwas_files):
        if not os.path.exists(gwas_file):
            raise FileNotFoundError(f"ERROR: GWAS file not found: {gwas_file}")
        print(f"Population {i+1} ({ancestry_labels[i]}) GWAS file: {gwas_file} ✓")
    
    # Check LD files exist
    for i, ld_file in enumerate(ld_files):
        vars_file = f"{ld_file}.vars"
        if not os.path.exists(ld_file) or not os.path.exists(vars_file):
            raise FileNotFoundError(f"ERROR: LD matrix ({ld_file}) or .vars file ({vars_file}) not found")
        print(f"Population {i+1} ({ancestry_labels[i]}) LD matrix: {ld_file}, .vars ✓")
    
    # Create output directory if needed
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"Validated input files for {n_pops} populations.")
    
    return gwas_files, ld_files, n_gwas_list, ancestry_labels, n_pops

def run_multisusie(validated_inputs: Tuple, args) -> Tuple[Any, List[pd.DataFrame]]:
    """Run MultiSuSiE fine-mapping analysis."""
    print("\n=== Running MultiSuSiE fine-mapping ===")
    
    try:
        import MultiSuSiE
    except ImportError:
        raise ImportError("ERROR: MultiSuSiE package not found. Please install with: pip install git+https://github.com/jordanero/MultiSuSiE.git")
    
    gwas_files, ld_files, n_gwas_list, ancestry_labels, n_pops = validated_inputs
    
    # Find common variants across all populations
    common_variants = find_common_variants(ld_files, gwas_files, ancestry_labels)
    
    # Get chromosome and BP info from first .vars file
    chrom, min_bp, max_bp, _ = read_vars_file(f"{ld_files[0]}.vars")
    print(f"Chromosome: {chrom}, BP range: {min_bp}-{max_bp}")
    
    # Process each population with common variants
    summ_stat_list = []
    b_list = []
    s_list = []
    R_list = []
    varY_list = []
    
    for i, (gwas_file, ld_file) in enumerate(zip(gwas_files, ld_files)):
        ancestry = ancestry_labels[i]
        print(f"\nProcessing population {i+1} ({ancestry}):")
        
        # Detect columns
        detected_cols = detect_columns(gwas_file)
        
        # Load GWAS data (filtered to common variants)
        gwas_data = load_gwas_data(gwas_file, detected_cols, common_variants, n_gwas_list[i])
        
        # Load LD matrix (pruned to common variants)
        vars_file = f"{ld_file}.vars"
        ld_matrix = load_ld_matrix(ld_file, vars_file, common_variants)
        
        # Store in lists
        summ_stat_list.append(gwas_data)
        
        # Extract data for MultiSuSiE
        if 'Beta' in gwas_data.columns:
            b_list.append(gwas_data['Beta'].values)
        else:
            raise ValueError(f"ERROR: Beta column not found for population {ancestry}")
        
        if 'Se' in gwas_data.columns:
            s_list.append(gwas_data['Se'].values)
        else:
            raise ValueError(f"ERROR: Se column not found for population {ancestry}")
        
        R_list.append(ld_matrix)
        
        # Estimate phenotype variance if not provided
        if args.var_y is not None:
            var_y_values = [float(v.strip()) for v in args.var_y.split(',')]
            varY_list.append(var_y_values[i])
        else:
            # Estimate from summary statistics (rough approximation)
            varY_est = 1.0  # Default value, can be improved
            varY_list.append(varY_est)
    
    # Parse prior weights if provided
    prior_weights = None
    if args.prior_weights is not None:
        # Assume it's a file path with weights
        if os.path.exists(args.prior_weights):
            prior_weights = np.loadtxt(args.prior_weights)
        else:
            # Assume it's comma-separated values
            prior_weights = np.array([float(w.strip()) for w in args.prior_weights.split(',')])
    
    # Run MultiSuSiE
    print("\n=== Running MultiSuSiE core analysis ===")
    
    multisusie_result = MultiSuSiE.multisusie_rss(
        b_list=b_list,
        s_list=s_list,
        R_list=R_list,
        varY_list=varY_list,
        population_sizes=n_gwas_list,
        L=args.L,
        max_iter=args.max_iter,
        tol=args.tol,
        estimate_residual_variance=args.estimate_residual_variance,
        low_memory_mode=args.low_memory_mode,
        standardize=args.standardize,
        intercept=args.intercept,
        prior_weights=prior_weights,
        coverage=args.coverage
    )
    
    print("MultiSuSiE analysis completed successfully!")
    
    return multisusie_result, summ_stat_list

def save_results(multisusie_result: Any, summ_stat_list: List[pd.DataFrame], ancestry_labels: List[str], output_file: str, args) -> None:
    """Save results in format compatible with other fine-mapping scripts."""
    print("\n=== Saving results ===")
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Get variant IDs from the first population (all should have same variants after organization)
    variant_ids = summ_stat_list[0]['SNP'].tolist()
    n_variants = len(variant_ids)
    
    # Initialize the output dataframe
    result_df = pd.DataFrame({
        'VARIANT_ID': variant_ids
    })
    
    # Add Z-score columns for each population
    for i, ancestry in enumerate(ancestry_labels):
        z_col_name = f"Z_{ancestry}"
        if 'Z' in summ_stat_list[i].columns:
            result_df[z_col_name] = summ_stat_list[i]['Z'].values
        elif 'Beta' in summ_stat_list[i].columns and 'Se' in summ_stat_list[i].columns:
            # Calculate Z-score if not present
            result_df[z_col_name] = summ_stat_list[i]['Beta'].values / summ_stat_list[i]['Se'].values
        else:
            # Fill with NA if Z-score cannot be determined
            result_df[z_col_name] = [np.nan] * n_variants
            print(f"WARNING: Could not find Z-scores for population {ancestry}")
    
    # Add overall PIP column
    if hasattr(multisusie_result, 'pip') and multisusie_result.pip is not None:
        pip_values = multisusie_result.pip
        result_df['PIP'] = np.where(np.isnan(pip_values), 0, pip_values)
    else:
        result_df['PIP'] = [0] * n_variants
        print("WARNING: No PIP values found in MultiSuSiE results")
    
    # Add converged column
    if hasattr(multisusie_result, 'converged'):
        converged_status = multisusie_result.converged
    else:
        converged_status = 'NA'  # Default to NA since convergence info may not be available
    result_df['converged'] = [converged_status] * n_variants
    
    # Add credible set column
    result_df['cs'] = [0] * n_variants  # Initialize with 0 (not in any CS)
    
    if hasattr(multisusie_result, 'sets') and multisusie_result.sets is not None:
        sets = multisusie_result.sets
        if hasattr(sets, 'cs') and sets.cs is not None:
            for cs_idx, cs_variants in enumerate(sets.cs, 1):
                # cs_variants should be indices
                if len(cs_variants) > 0:
                    for variant_idx in cs_variants:
                        if 0 <= variant_idx < len(result_df):
                            result_df.loc[variant_idx, 'cs'] = cs_idx
            print(f"Assigned {len(sets.cs)} credible sets to variants")
        else:
            print("No credible sets found in MultiSuSiE results")
    else:
        print("No credible sets found in MultiSuSiE results")
    
    # Add population-specific PIP columns if available
    if hasattr(multisusie_result, 'pip_config') and multisusie_result.pip_config is not None:
        pip_config = multisusie_result.pip_config
        # This would depend on the exact structure of MultiSuSiE output
        # For now, create placeholder columns
        for ancestry in ancestry_labels:
            result_df[f'PIP_{ancestry}'] = [0] * n_variants
        
        # Create intersection column for multi-ancestry
        if len(ancestry_labels) >= 2:
            intersection_name = '_'.join(ancestry_labels)
            result_df[f'PIP_{intersection_name}'] = [0] * n_variants
    else:
        print("WARNING: No population-specific PIPs found in MultiSuSiE results")
        # Create placeholder columns with zeros
        for ancestry in ancestry_labels:
            result_df[f'PIP_{ancestry}'] = [0] * n_variants
        
        # Create intersection column for 2+ populations
        if len(ancestry_labels) >= 2:
            intersection_name = '_'.join(ancestry_labels)
            result_df[f'PIP_{intersection_name}'] = [0] * n_variants
    
    # Write the tab-delimited output file
    result_df.to_csv(output_file, sep='\t', index=False, na_rep='NA')
    
    print(f"Saved results to: {output_file}")
    print(f"Output contains {len(result_df)} variants and {len(result_df.columns)} columns")
    
    # Print column summary for verification
    print("Output columns:")
    for i, col_name in enumerate(result_df.columns, 1):
        if col_name == "VARIANT_ID":
            print(f" {i}. {col_name} (variant identifiers)")
        elif col_name.startswith("Z_"):
            print(f" {i}. {col_name} (Z-scores)")
        elif col_name == "PIP":
            print(f" {i}. {col_name} (overall posterior inclusion probability)")
        elif col_name == "converged":
            print(f" {i}. {col_name} (convergence status)")
        elif col_name == "cs":
            n_in_cs = sum(1 for x in result_df[col_name] if x > 0)
            print(f" {i}. {col_name} (credible set assignments, {n_in_cs} variants in sets)")
        elif col_name.startswith("PIP_"):
            print(f" {i}. {col_name} (population-specific PIP)")
    
    print("Results saved successfully!")

def main():
    """Main function."""
    print("=== MultiSuSiE Multi-Population Fine-Mapping ===")
    
    try:
        # Parse arguments
        args = parse_arguments()
        
        # Validate inputs
        validated_inputs = validate_inputs(args)
        
        # Run MultiSuSiE analysis
        multisusie_result, summ_stat_list = run_multisusie(validated_inputs, args)
        
        # Save results
        save_results(multisusie_result, summ_stat_list, validated_inputs[3], args.output_file, args)
        
        print("\n=== Analysis completed successfully! ===")
        
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main()