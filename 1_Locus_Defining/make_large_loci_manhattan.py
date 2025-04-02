# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 17:08:08 2025

@author: mitch
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

def generate_manhattan(summary_stats_path, genomic_region, 
                      chr_col="CHR", pos_col="POS", p_col="P",
                      sig_threshold=5e-8, point_size=8, 
                      figsize=(14, 6), save_path=None):
    """
    Generate a Manhattan plot for a specific genomic region
    
    Parameters:
    summary_stats_path (str): Path to summary statistics file (CSV)
    genomic_region (str): Genomic region in "chr:start-end" format
    chr_col (str): Chromosome column name (default: "CHR")
    pos_col (str): Position column name (default: "POS")
    p_col (str): P-value column name (default: "P")
    sig_threshold (float): Genome-wide significance threshold (default: 5e-8)
    point_size (int): Marker size for points (default: 8)
    figsize (tuple): Figure dimensions (default: (14, 6))
    save_path (str): Optional path to save plot (default: None)
    
    Returns:
    matplotlib.figure.Figure: Generated figure object
    """
    # Read and prepare data
    df = pd.read_csv(summary_stats_path)
    df = df.rename(columns={chr_col: "CHR", pos_col: "POS", p_col: "P"})
    
    # Parse genomic region
    chrom, region = genomic_region.split(":")
    start, end = map(int, region.split("-"))
    
    # Filter data for specified region
    region_df = df[
        (df["CHR"].astype(str) == chrom.lstrip("chr")) & 
        (df["POS"].between(start, end))
    ].copy()
    
    # Calculate transformed p-values
    region_df["-log10P"] = -np.log10(region_df["P"])
    
    # Create plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Color chromosomes alternately
    chroms = region_df["CHR"].unique()
    colors_cycle = ["#4E79A7", "#F28E2B"]  # Alternating colors
    
    for i, chrom in enumerate(chroms):
        chrom_df = region_df[region_df["CHR"] == chrom]
        ax.scatter(
            x=chrom_df["POS"],
            y=chrom_df["-log10P"],
            color=colors_cycle[i % len(colors_cycle)],
            s=point_size,
            alpha=0.7
        )
    
    # Add significance threshold line
    ax.axhline(-np.log10(sig_threshold), color="red", linestyle="--")
    
    # Format plot
    ax.set_xlabel(f"Chromosome {chrom} Position", fontsize=12)
    ax.set_ylabel("-log10(p-value)", fontsize=12)
    ax.set_title(f"Manhattan Plot for {genomic_region}", fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    
    # Save if requested
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        
    return fig
