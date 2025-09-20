'''
plot_susiex_results.multi-pop.py

This script plots the results of a SuSiEx multi-population fine-mapping given a results file
and a list of the causal variants used to simulate the GWAS.

Creates facet plots with:
- Separate facets for each population's p-values (from Z-scores)
- Single facet for PIPs
- Highlighting of causal variants and credible sets

***This needs to be run in an environment with matplotlib, pandas, numpy, and scipy ***
'''

import matplotlib
matplotlib.use('Agg')  # Set backend before importing pyplot
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import norm
import argparse
import re

def identify_populations(df):
    """Identify population names from Z-score columns (Z_AFR, Z_EUR, etc.)"""
    z_cols = [col for col in df.columns if col.startswith('Z_')]
    populations = [col.replace('Z_', '') for col in z_cols]
    
    if not populations:
        # Fallback: look for single Z column (single population)
        if 'Z' in df.columns:
            populations = ['Pop1']  # Default name for single population
            df['Z_Pop1'] = df['Z']  # Rename for consistency
        else:
            raise ValueError("No Z-score columns found in the data")
    
    return populations

def calculate_pvalues_from_zscore(z_scores):
    """Calculate -log10 p-values from Z-scores"""
    # Handle missing/NA values
    valid_mask = pd.notna(z_scores) & (z_scores != 'NA')
    z_numeric = pd.to_numeric(z_scores, errors='coerce')
    
    # Calculate p-values for valid Z-scores
    log_sf = norm.logsf(np.abs(z_numeric))  # log(1 - CDF(|Z|))
    log_p = log_sf + np.log(2)  # Account for two-tailed test
    neg_log10_p = -log_p / np.log(10)  # Convert to -log10(p)
    
    # Set invalid values to NaN
    neg_log10_p[~valid_mask] = np.nan
    
    return neg_log10_p

def plot_susiex_results(results_file, causal_file, output_file):
    """Plot SuSiEx multi-population fine-mapping results"""
    
    # Load data
    print("Loading results file...")
    df = pd.read_csv(results_file, sep="\t")
    
    if causal_file:
        print("Loading causal variants...")
        causal_snps = pd.read_csv(causal_file, header=None, names=['SNP'])
        causal_set = set(causal_snps['SNP'])
    else:
        causal_set = set()
        print("No causal variants file provided")
    
    # Extract genomic position from VARIANT_ID column
    print("Processing genomic positions...")
    df['pos'] = df['VARIANT_ID'].apply(
        lambda x: int(x.split(':')[1]) if pd.notna(x) and ':' in str(x) else np.nan
    )
    
    # Get chromosome and position range for labeling
    first_variant = df['VARIANT_ID'].dropna().iloc[0]
    chrom = first_variant.split(':')[0]
    min_pos = df['pos'].min()
    max_pos = df['pos'].max()
    
    # Identify populations
    populations = identify_populations(df)
    n_pops = len(populations)
    print(f"Found {n_pops} populations: {populations}")
    
    # Calculate p-values for each population
    for pop in populations:
        z_col = f'Z_{pop}'
        if z_col in df.columns:
            df[f'neg_log10_p_{pop}'] = calculate_pvalues_from_zscore(df[z_col])
        else:
            print(f"Warning: Z-score column {z_col} not found")
    
    # Create color mapping for credible sets (exclude cs=0)
    unique_cs = df.loc[df['cs'] != 0, 'cs'].dropna().unique()
    if len(unique_cs) > 0:
        color_map = plt.cm.get_cmap('tab10', len(unique_cs))
        cs_colors = {cs: color_map(i) for i, cs in enumerate(unique_cs)}
    else:
        cs_colors = {}
    
    # Assign colors (treat cs=0 as non-credible set)
    df['color'] = df['cs'].apply(
        lambda x: cs_colors[x] if (pd.notna(x) and x != 0 and x in cs_colors) else 'lightgrey'
    )
    
    # Create subplot layout: n_pops rows for p-values + 1 row for PIPs
    n_rows = n_pops + 1
    fig, axes = plt.subplots(n_rows, 1, figsize=(14, 3 * n_rows), sharex=True)
    
    # Handle case where there's only one subplot
    if n_rows == 1:
        axes = [axes]
    
    # Plot p-values for each population
    for i, pop in enumerate(populations):
        ax = axes[i]
        neg_log10_p_col = f'neg_log10_p_{pop}'
        
        if neg_log10_p_col in df.columns:
            # Plot association signals
            valid_data = df.dropna(subset=[neg_log10_p_col])
            ax.scatter(valid_data['pos'], valid_data[neg_log10_p_col], 
                      c=valid_data['color'], edgecolor='none', s=40, alpha=0.7)
            
            ax.set_ylabel(f'-log10(p) {pop}')
            ax.set_title(f'Association Signals - {pop} Population')
            
            # Highlight causal variants if provided
            if causal_set:
                causal_df = df[df['VARIANT_ID'].isin(causal_set) & df[neg_log10_p_col].notna()]
                if len(causal_df) > 0:
                    ax.scatter(causal_df['pos'], causal_df[neg_log10_p_col], 
                             facecolors='none', edgecolors='black', s=100, linewidth=2)
        else:
            ax.text(0.5, 0.5, f'No data for {pop}', transform=ax.transAxes, 
                   ha='center', va='center', fontsize=12)
            ax.set_ylabel(f'-log10(p) {pop}')
    
    # Plot PIPs in the bottom panel
    pip_ax = axes[-1]
    
    # Convert PIP to numeric, handle 'NA' values
    df['PIP_numeric'] = pd.to_numeric(df['PIP'], errors='coerce')
    valid_pip_data = df.dropna(subset=['PIP_numeric'])
    
    pip_ax.scatter(valid_pip_data['pos'], valid_pip_data['PIP_numeric'], 
                  c=valid_pip_data['color'], edgecolor='none', s=40, alpha=0.7)
    pip_ax.set_ylabel('PIP')
    pip_ax.set_xlabel(f'Genomic Position (chr{chrom}:{min_pos:,}-{max_pos:,})')
    pip_ax.set_title('Posterior Inclusion Probabilities (PIPs)')
    
    # Highlight causal variants in PIP plot
    if causal_set:
        causal_df = df[df['VARIANT_ID'].isin(causal_set) & df['PIP_numeric'].notna()]
        if len(causal_df) > 0:
            pip_ax.scatter(causal_df['pos'], causal_df['PIP_numeric'], 
                         facecolors='none', edgecolors='black', s=100, linewidth=2)
    
    # Create legend
    legend_elements = []
    
    # Add credible set colors to legend
    for cs in sorted(unique_cs):
        if cs in cs_colors:
            legend_elements.append(
                plt.Line2D([0], [0], marker='o', color='w', label=f'CS {int(cs)}',
                          markerfacecolor=cs_colors[cs], markersize=8)
            )
    
    # Add non-credible set and causal SNP to legend
    legend_elements.extend([
        plt.Line2D([0], [0], marker='o', color='w', label='Non-CS',
                   markerfacecolor='lightgrey', markersize=8),
    ])
    
    if causal_set:
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='black', label='Causal SNP',
                      markerfacecolor='none', markersize=10, markeredgewidth=2)
        )
    
    # Place legend on the top subplot
    if legend_elements:
        axes[0].legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add overall title
    fig.suptitle(f'SuSiEx Multi-Population Fine-Mapping Results\\nchr{chrom}:{min_pos:,}-{max_pos:,}', 
                 fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)  # Make room for suptitle
    
    # Save plot
    print(f"Saving plot to {output_file}")
    plt.savefig(output_file, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Plot completed successfully!")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot SuSiEx multi-population fine-mapping results with population-specific p-values and PIPs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--results_file',
        type=str, 
        required=True,
        help='Path to the SuSiEx results file (tab-separated, .txt format from create_susie_format_output)'
    )
    
    parser.add_argument(
        '--causal_file',
        type=str,
        required=False,
        help='Path to the causal SNPs file (one SNP per line, optional)'
    )
    
    parser.add_argument(
        '--output_file',
        type=str,
        required=True,
        help='Output path for JPEG plot (e.g., susiex_results_plot.jpg)'
    )
    
    args = parser.parse_args()
    
    plot_susiex_results(args.results_file, args.causal_file, args.output_file)