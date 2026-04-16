'''

plot_matrices.py

For each of the 100 sampled PLATLAS loci and 2 populations (200 pairs total),
this script reads the shared-variant subset matrices produced by
subset_matrices.py and generates a side-by-side heatplot saved as a PNG.

Each plot shows the betanfer LD matrix and the raw MVP LD matrix for the same
locus and population side-by-side, using a diverging colour scale (RdBu_r,
range -1 to 1) with a shared colourbar.

Output:
  {plot_dir}/{POP}/chr{CHR}_{START}_{END}.png

Usage:
  python plot_matrices.py -i <subset_dir> -o <plot_dir> -l <locus_csv>
                          [-a <populations>] [-n <n_workers>]

***Run on Frontier (Python/3.9)***

'''

import matplotlib
matplotlib.use('Agg')   # non-interactive backend — must precede pyplot import

import sys
import getopt
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from multiprocessing import Pool

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_POPULATIONS = "AFR,EUR"
DEFAULT_N_WORKERS   = 20


def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
-i => Subset matrices base directory (output of subset_matrices.py) REQUIRED
      Expected structure: {subset_dir}/{POP}/chr*_{START}_{END}.betanfer.vcor1
-o => Plot output directory REQUIRED
-l => Sampled locus CSV with columns CHR,START,END REQUIRED
-a => Populations comma-separated (default: AFR,EUR) OPTIONAL
-n => Number of parallel worker processes (default: 20) OPTIONAL
-h => Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:l:a:n:h")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    if '-h' in options_dict:
        help(0)

    # Required arguments
    try:
        subset_dir = options_dict['-i']
        plot_dir   = options_dict['-o']
        locus_csv  = options_dict['-l']
    except KeyError:
        print("ERROR: One or more required arguments are missing.")
        help()

    # Optional arguments
    populations = [p.strip().upper()
                   for p in options_dict.get('-a', DEFAULT_POPULATIONS).split(',')]
    n_workers   = int(options_dict.get('-n', DEFAULT_N_WORKERS))

    for path, label in [
        (subset_dir, "subset matrices directory"),
        (locus_csv,  "locus CSV"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    print("Acceptable Inputs Given")

    driver(subset_dir, plot_dir, locus_csv, populations, n_workers)


###############################################################################
################################## WORKER ####################################
###############################################################################

def process_locus_pop(work_item):
    """
    Load the betanfer and MVP subset matrices for one (locus, population) pair
    and save a side-by-side heatplot as a PNG.

    Returns (success, status_message).
    """
    locus_label = work_item['locus_label']
    pop         = work_item['pop']
    b_vcor1     = work_item['b_vcor1']
    m_vcor1     = work_item['m_vcor1']
    plot_path   = work_item['plot_path']

    tag = f"{locus_label}.{pop}"

    # Check input files exist
    for path, label in [(b_vcor1, "betanfer subset vcor1"),
                        (m_vcor1, "MVP subset vcor1")]:
        if not os.path.exists(path):
            return (False, f"WARNING: Missing {label} for {tag} — skipping")

    # Load matrices
    sub_b = np.loadtxt(b_vcor1, delimiter='\t', dtype=np.float64)
    sub_m = np.loadtxt(m_vcor1, delimiter='\t', dtype=np.float64)

    # Ensure 2-D and non-trivial
    if sub_b.ndim < 2 or sub_b.shape[0] < 2:
        return (False, f"NOTE: Matrix too small to plot for {tag} — skipping")

    n_vars = sub_b.shape[0]

    # ── Build side-by-side figure ──────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    cmap = 'RdBu_r'
    vmin, vmax = -1, 1

    im = axes[0].imshow(sub_b, cmap=cmap, vmin=vmin, vmax=vmax,
                        aspect='auto', interpolation='nearest')
    axes[0].set_title('Betanfer', fontsize=12)
    axes[0].set_xticks([])
    axes[0].set_yticks([])

    axes[1].imshow(sub_m, cmap=cmap, vmin=vmin, vmax=vmax,
                   aspect='auto', interpolation='nearest')
    axes[1].set_title('Raw MVP (plink)', fontsize=12)
    axes[1].set_xticks([])
    axes[1].set_yticks([])

    # Shared colourbar on the right
    fig.colorbar(im, ax=axes, label='r', fraction=0.046, pad=0.04)

    fig.suptitle(f"{locus_label}  |  {pop}  ({n_vars} shared variants)",
                 fontsize=13)

    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

    return (True, f"NOTE: Saved plot for {tag} ({n_vars} variants)")


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(subset_dir, plot_dir, locus_csv, populations, n_workers):

    # ── Read loci ──────────────────────────────────────────────────────────
    print(f"NOTE: Reading locus file: {locus_csv}")
    loci = pd.read_csv(locus_csv)
    print(f"NOTE: {len(loci)} loci to process")

    # ── Create output population subdirectories ────────────────────────────
    for pop in populations:
        os.makedirs(os.path.join(plot_dir, pop), exist_ok=True)

    # ── Build work item list (one per locus × population) ─────────────────
    work_items = []
    for _, locus in loci.iterrows():
        chr_num = int(locus['CHR'])
        start   = int(locus['START'])
        end     = int(locus['END'])
        locus_label = f"chr{chr_num}_{start}_{end}"

        for pop in populations:
            pop_sub  = os.path.join(subset_dir, pop)
            pop_plot = os.path.join(plot_dir,   pop)
            work_items.append({
                'locus_label': locus_label,
                'pop':         pop,
                'b_vcor1':     os.path.join(pop_sub,  f"{locus_label}.betanfer.vcor1"),
                'm_vcor1':     os.path.join(pop_sub,  f"{locus_label}.mvp.vcor1"),
                'plot_path':   os.path.join(pop_plot, f"{locus_label}.png"),
            })

    print(f"NOTE: {len(work_items)} plots to generate using {n_workers} workers")

    # ── Run parallel pool ──────────────────────────────────────────────────
    saved   = 0
    skipped = 0

    with Pool(processes=n_workers) as pool:
        for success, msg in pool.imap_unordered(process_locus_pop, work_items):
            print(msg)
            if success:
                saved += 1
            else:
                skipped += 1

    print(f"\nSUMMARY: {saved}/{len(work_items)} plots saved to {plot_dir}")
    if skipped > 0:
        print(f"WARNING: {skipped} plot(s) skipped — check warnings above")
    else:
        print("SUCCESS: All plots generated without skipping")


###############################################################################
################################ CALL MAIN ###################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv)
