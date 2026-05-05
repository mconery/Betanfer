#!/usr/bin/env python3
'''
plot_mvp_loci.py

Randomly selects N loci per MVP population (AFR, EUR) that appear in both the
original analysis (Data_S1.xlsx) and the new fine-mapping results, then
produces a 3-panel locus plot for each selected locus:
  Panel 1: Original SuSiE credible sets (from Data_S1.xlsx)
  Panel 2: New SuSiE credible sets
  Panel 3: New RSparsePro credible sets

Y-axis: -log10(GWAS p-value) from MVP summary statistics (via tabix).
X-axis: Genomic position (hg19).
CS variants colored by tab10 palette; non-CS variants in lightgrey.
'''

import sys
import getopt
import os
import re
import random
import subprocess
import io
import multiprocessing

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D

EXCEL_DEFAULT    = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/Data_S1.xlsx'
RESULTS_DEFAULT  = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/MVP'
SUMSTATS_DEFAULT = '/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19'
TABIX_DEFAULT    = '/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/tabix'

TAB10 = cm.get_cmap('tab10').colors


def usage():
    print('Usage: plot_mvp_loci.py --output-dir <dir> [options]', file=sys.stderr)
    print('  --excel-file  <path>   default: Data_S1.xlsx', file=sys.stderr)
    print('  --results-dir <path>   default: fine-mapping_results/MVP', file=sys.stderr)
    print('  --sumstats-dir <path>  default: sumstats/MVP_hg19', file=sys.stderr)
    print('  --n-loci <int>         loci per population (default: 100)', file=sys.stderr)
    print('  --seed <int>           random seed (default: 42)', file=sys.stderr)
    print('  --tabix <path>         tabix binary', file=sys.stderr)
    print('  --workers <int>        parallel workers (default: 16)', file=sys.stderr)
    sys.exit(1)


def cs_color(cs_id):
    if int(cs_id) == 0:
        return 'lightgrey'
    return TAB10[(int(cs_id) - 1) % 10]


def parse_stem(stem):
    """'chr10_114250001_115250000' -> ('10', 114250001, 115250000)"""
    m = re.match(r'chr(\w+)_(\d+)_(\d+)', stem)
    if not m:
        raise ValueError(f'Cannot parse locus stem: {stem}')
    return m.group(1), int(m.group(2)), int(m.group(3))


def tabix_extract(tabix_bin, sumstats_path, chrom, start, end):
    """Return DataFrame with RSID, POS, LOG10P for the region."""
    region = f'{chrom}:{start}-{end}'
    result = subprocess.run(
        [tabix_bin, '-h', sumstats_path, region],
        capture_output=True, text=True
    )
    if result.returncode != 0 or not result.stdout.strip():
        return pd.DataFrame(columns=['RSID', 'POS', 'LOG10P'])
    df = pd.read_csv(io.StringIO(result.stdout), sep='\t', low_memory=False)
    df.columns = [c.lstrip('#') for c in df.columns]
    df = df[['RSID', 'POS', 'LOG10P']].copy()
    df['POS']    = pd.to_numeric(df['POS'],    errors='coerce')
    df['LOG10P'] = pd.to_numeric(df['LOG10P'], errors='coerce')
    return df.dropna(subset=['POS', 'LOG10P']).reset_index(drop=True)


def load_finemap_cs(path, id_col):
    """Return {rsid: cs_int} from a SuSiE or RSparsePro output file."""
    try:
        df = pd.read_csv(path, sep='\t', usecols=[id_col, 'cs'])
        return dict(zip(df[id_col].astype(str), df['cs'].fillna(0).astype(int)))
    except Exception as e:
        print(f'Warning: could not read {path}: {e}', file=sys.stderr)
        return {}


def plot_locus(args):
    (trait, pop, stem, orig_cs,
     results_dir, sumstats_dir, output_dir, tabix_bin) = args

    out_path = os.path.join(output_dir, pop, f'{trait}.{stem}.locus_plot.jpg')
    if os.path.exists(out_path):
        return

    try:
        chrom, start, end = parse_stem(stem)
    except ValueError as e:
        print(f'Warning: {e}', file=sys.stderr)
        return

    sumstats_path = os.path.join(
        sumstats_dir, pop,
        f'MVP_R4.1000G_AGR.{trait}.{pop}.GIA.dbGaP.sumstats.hg19.txt.gz'
    )
    susie_path = os.path.join(results_dir, 'susie', pop,
                              f'{trait}.{stem}.susie.txt')
    rsp_path   = os.path.join(results_dir, 'rsparsepro', pop,
                              f'{trait}.{stem}.rsparsepro.txt')

    tabix_df = tabix_extract(tabix_bin, sumstats_path, chrom, start, end)
    if tabix_df.empty:
        print(f'Warning: no tabix data for {trait} {pop} {stem}', file=sys.stderr)
        return

    susie_cs = load_finemap_cs(susie_path, 'VARIANT_ID')
    rsp_cs   = load_finemap_cs(rsp_path,   'RSID')

    panels = [
        ('Original SuSiE', orig_cs),
        ('New SuSiE',       susie_cs),
        ('New RSparsePro',  rsp_cs),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    fig.suptitle(
        f'{trait}  |  {pop}  |  chr{chrom}:{start:,}–{end:,}',
        fontsize=11, y=0.998
    )

    for ax, (label, cs_dict) in zip(axes, panels):
        df = tabix_df.copy()
        df['cs'] = df['RSID'].map(cs_dict).fillna(0).astype(int)
        colors = df['cs'].map(cs_color).tolist()

        ax.scatter(df['POS'], df['LOG10P'],
                   c=colors, s=8, alpha=0.7, linewidths=0, rasterized=True)
        ax.set_ylabel(r'$-\log_{10}(P)$', fontsize=9)
        ax.set_title(label, loc='left', fontsize=9, pad=3)
        ax.spines[['top', 'right']].set_visible(False)

        cs_ids = sorted(df.loc[df['cs'] > 0, 'cs'].unique())
        if cs_ids:
            handles = [
                Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=cs_color(c), markersize=6,
                       label=f'CS {c}')
                for c in cs_ids
            ]
            ax.legend(handles=handles, fontsize=7, frameon=False,
                      loc='upper right', ncol=min(len(cs_ids), 5))

    axes[-1].set_xlabel('Genomic Position (hg19)', fontsize=9)
    plt.tight_layout(rect=[0, 0, 1, 0.997])
    fig.savefig(out_path, dpi=150, format='jpeg', bbox_inches='tight')
    plt.close(fig)


def main():
    try:
        opts, _ = getopt.getopt(sys.argv[1:], '', [
            'excel-file=', 'results-dir=', 'sumstats-dir=', 'output-dir=',
            'n-loci=', 'seed=', 'tabix=', 'workers=',
        ])
    except getopt.GetoptError as e:
        print(f'Error: {e}', file=sys.stderr)
        usage()

    excel_file   = EXCEL_DEFAULT
    results_dir  = RESULTS_DEFAULT
    sumstats_dir = SUMSTATS_DEFAULT
    output_dir   = None
    n_loci       = 100
    seed         = 42
    tabix_bin    = TABIX_DEFAULT
    workers      = 16

    for opt, val in opts:
        if   opt == '--excel-file':   excel_file   = val
        elif opt == '--results-dir':  results_dir  = val
        elif opt == '--sumstats-dir': sumstats_dir = val
        elif opt == '--output-dir':   output_dir   = val
        elif opt == '--n-loci':       n_loci       = int(val)
        elif opt == '--seed':         seed         = int(val)
        elif opt == '--tabix':        tabix_bin    = val
        elif opt == '--workers':      workers      = int(val)

    if output_dir is None:
        print('Error: --output-dir is required', file=sys.stderr)
        usage()

    os.makedirs(os.path.join(output_dir, 'AFR'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'EUR'), exist_ok=True)

    # Load Excel and build per-locus original CS map
    print('Loading Excel...', flush=True)
    df = pd.read_excel(excel_file, header=1)
    df = df[df['Population'].isin(['AFR', 'EUR'])].copy()
    df['locus_stem'] = df['Locus'].astype(str).str.replace('.', '_', regex=False)
    df['Merged Signal'] = pd.to_numeric(
        df['Merged Signal'], errors='coerce').fillna(0).astype(int)

    orig_cs_map = {}
    for (trait, pop, stem), grp in df.groupby(
            ['Trait', 'Population', 'locus_stem']):
        orig_cs_map[(trait, pop, stem)] = dict(
            zip(grp['RSID'].astype(str), grp['Merged Signal'])
        )
    print(f'  {len(orig_cs_map)} unique (Trait, Pop, locus) entries', flush=True)

    # Build per-population candidate list (both new result files must exist)
    candidates = {'AFR': [], 'EUR': []}
    for (trait, pop, stem) in orig_cs_map:
        susie = os.path.join(results_dir, 'susie', pop,
                             f'{trait}.{stem}.susie.txt')
        rsp   = os.path.join(results_dir, 'rsparsepro', pop,
                             f'{trait}.{stem}.rsparsepro.txt')
        if os.path.exists(susie) and os.path.exists(rsp):
            candidates[pop].append((trait, pop, stem))

    for pop in ('AFR', 'EUR'):
        print(f'  {pop} candidates: {len(candidates[pop])}', flush=True)

    # Random sample
    random.seed(seed)
    selected = []
    for pop in ('AFR', 'EUR'):
        k = min(n_loci, len(candidates[pop]))
        chosen = random.sample(candidates[pop], k)
        selected.extend(chosen)
        print(f'  Selected {k} {pop} loci', flush=True)

    pool_args = [
        (trait, pop, stem,
         orig_cs_map[(trait, pop, stem)],
         results_dir, sumstats_dir, output_dir, tabix_bin)
        for trait, pop, stem in selected
    ]

    print(f'Generating {len(pool_args)} plots with {workers} workers...', flush=True)
    with multiprocessing.Pool(workers) as pool:
        pool.map(plot_locus, pool_args)

    print('Done.', flush=True)


if __name__ == '__main__':
    main()
