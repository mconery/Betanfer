'''
plot_rsparsepro_results.py

This script plots the results of an RSparsePro fine-mapping given a results file
and a list of the causal variants used to simulate the GWAS. 

***This needs to be run in the rsparsepro conda environment ***
'''

import matplotlib
matplotlib.use('Agg')  # Set backend before importing pyplot
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import norm
import argparse

def plot_finemapping_results(results_file, causal_file, output_file):
    # Load data
    df = pd.read_csv(results_file)
    causal_snps = pd.read_csv(causal_file, header=None, names=['SNP'])
    causal_set = set(causal_snps['SNP'])

    # Calculate p-values from Z-scores
    df['neg_log10_p'] = -np.log10(2 * (1 - norm.cdf(np.abs(df['Z']))))
    
    # Create color mapping for credible sets
    unique_cs = df['cs'].dropna().unique()
    color_map = plt.cm.get_cmap('tab10', len(unique_cs))
    cs_colors = {cs: color_map(i) for i, cs in enumerate(unique_cs)}
    df['color'] = df['cs'].apply(lambda x: cs_colors[x] if pd.notna(x) else 'lightgrey')

    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    # Top plot: -log10 p-values
    ax1.scatter(df['pos'], df['neg_log10_p'], c=df['color'], edgecolor='none', s=40)
    ax1.set_ylabel('-log10(p-value)')
    ax1.set_title('Association Signals and Credible Sets')
    
    # Bottom plot: PIP values
    ax2.scatter(df['pos'], df['PIP'], c=df['color'], edgecolor='none', s=40)
    ax2.set_ylabel('PIP')
    ax2.set_xlabel('Genomic Position')
    
    # Highlight causal variants
    causal_df = df[df['SNP'].isin(causal_set)]
    for ax in [ax1, ax2]:
        y = causal_df['neg_log10_p'] if ax == ax1 else causal_df['PIP']
        ax.scatter(causal_df['pos'], y, facecolors='none', 
                   edgecolors='black', s=100, linewidth=1.5)
    
    # Create legend
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label=f'CS {cs}',
                  markerfacecolor=cs_colors[cs], markersize=8) for cs in unique_cs
    ]
    legend_elements += [
        plt.Line2D([0], [0], marker='o', color='w', label='Non-CS',
                  markerfacecolor='lightgrey', markersize=8),
        plt.Line2D([0], [0], marker='o', color='black', label='Causal SNP',
                  markerfacecolor='none', markersize=10, markeredgewidth=1.5)
    ]
    ax1.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(output_file, format='jpeg', dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plot RSparsePro fine-mapping results with credible sets and causal SNPs.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'results_file',
        type=str,
        help='Path to the RSparsePro results CSV file'
    )
    parser.add_argument(
        'causal_file',
        type=str,
        help='Path to the causal SNPs CSV file (one SNP per line)'
    )
    parser.add_argument(
        'output_file',
        type=str,
        help='Output path for JPEG plot (e.g., results_plot.jpg)'
    )
    args = parser.parse_args()

    plot_finemapping_results(args.results_file, args.causal_file, args.output_file)
