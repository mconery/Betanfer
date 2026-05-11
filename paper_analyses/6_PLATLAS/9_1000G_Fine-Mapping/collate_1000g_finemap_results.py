#!/usr/bin/env python3
'''
collate_1000g_finemap_results.py

Compares original MVP fine-mapping results (Data_S1.xlsx) against new SuSiE
and RSparsePro results produced with 1000G LD matrices.

Produces two output files:
  {output_prefix}.signal_recovery.tsv  -- one row per original signal
  {output_prefix}.pip_comparison.tsv   -- one row per original variant

Only AFR and EUR populations are compared (AMR and EAS were not re-run).

Variant matching uses the MVP ID column (CHR:POS:REF:ALT) from the Excel file,
with an allele-flip fallback (CHR:POS:ALT:REF) to handle cases where the 1000G
LD matrix uses the opposite allele orientation.
'''

import sys
import getopt
import os
import pandas as pd

EXCEL_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/Data_S1.xlsx'
RESULTS_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/MVP_1000G'


def usage():
    print('Usage: collate_1000g_finemap_results.py --output-prefix <prefix> '
          '[--excel-file <path>] [--results-dir <path>]', file=sys.stderr)
    sys.exit(1)


def in_cs(cs_val):
    """Return True if cs_val indicates credible set membership (nonzero, non-NA)."""
    try:
        return pd.notna(cs_val) and int(cs_val) != 0
    except (ValueError, TypeError):
        return False


def flip_alleles(variant_id):
    """Swap REF and ALT in a CHR:POS:REF:ALT ID string."""
    parts = variant_id.split(':')
    if len(parts) == 4:
        return f'{parts[0]}:{parts[1]}:{parts[3]}:{parts[2]}'
    return None


def lookup_variant(data, vid):
    """Look up a variant in a result dict, trying allele-flip if not found directly."""
    if vid in data:
        return data[vid]
    flipped = flip_alleles(vid)
    if flipped and flipped in data:
        return data[flipped]
    return None


def load_finemap_file(path):
    """Load a SuSiE/RSparsePro output file. Returns {variant_id: (pip, cs)} or None."""
    if not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep='\t')
        id_col = 'VARIANT_ID' if 'VARIANT_ID' in df.columns else 'RSID'
        return dict(zip(df[id_col].astype(str), zip(df['PIP'], df['cs'])))
    except Exception as e:
        print(f'Warning: could not read {path}: {e}', file=sys.stderr)
        return None


def main():
    try:
        opts, _ = getopt.getopt(sys.argv[1:], '',
                                ['excel-file=', 'results-dir=', 'output-prefix='])
    except getopt.GetoptError as e:
        print(f'Error: {e}', file=sys.stderr)
        usage()

    excel_file = EXCEL_DEFAULT
    results_dir = RESULTS_DEFAULT
    output_prefix = None

    for opt, val in opts:
        if opt == '--excel-file':
            excel_file = val
        elif opt == '--results-dir':
            results_dir = val
        elif opt == '--output-prefix':
            output_prefix = val

    if output_prefix is None:
        print('Error: --output-prefix is required', file=sys.stderr)
        usage()

    # Load Excel and filter to AFR/EUR only
    print(f'Loading {excel_file}...', flush=True)
    df = pd.read_excel(excel_file, header=1)
    df = df[df['Population'].isin(['AFR', 'EUR'])].copy()
    # Locus in Excel uses dots (chr10.114250001.115250000); filenames use underscores
    df['locus_stem'] = df['Locus'].astype(str).str.replace('.', '_', regex=False)
    print(f'  {len(df)} AFR/EUR rows loaded', flush=True)

    # Load result files for all unique (Trait, Population, locus) combinations
    print('Loading new fine-mapping result files...', flush=True)
    unique_loci = df[['Trait', 'Population', 'locus_stem']].drop_duplicates()
    results_cache = {}
    for _, row in unique_loci.iterrows():
        trait, pop, stem = row['Trait'], row['Population'], row['locus_stem']
        for method in ('susie', 'rsparsepro'):
            fname = f'{trait}.{stem}.{method}.txt'
            fpath = os.path.join(results_dir, method, pop, fname)
            results_cache[(trait, pop, stem, method)] = load_finemap_file(fpath)

    n_total = len(unique_loci) * 2
    n_found = sum(1 for v in results_cache.values() if v is not None)
    print(f'  Found {n_found}/{n_total} expected result files', flush=True)

    # Part 1: Signal Recovery — one row per original (Trait, Population, Locus, Signal)
    print('Computing signal recovery...', flush=True)
    signal_rows = []
    for (trait, pop, locus, signal), grp in df.groupby(
            ['Trait', 'Population', 'Locus', 'Merged Signal'], sort=False):
        stem = grp['locus_stem'].iloc[0]
        mvp_ids = set(grp['MVP ID'].astype(str))
        mafs = pd.to_numeric(grp['EAF Population'], errors='coerce').apply(
            lambda e: min(e, 1 - e) if pd.notna(e) else float('nan'))
        susie_data = results_cache.get((trait, pop, stem, 'susie'))
        rsp_data = results_cache.get((trait, pop, stem, 'rsparsepro'))
        any_susie = (susie_data is not None and
                     any(in_cs(lookup_variant(susie_data, m)[1])
                         for m in mvp_ids
                         if lookup_variant(susie_data, m) is not None))
        any_rsp = (rsp_data is not None and
                   any(in_cs(lookup_variant(rsp_data, m)[1])
                       for m in mvp_ids
                       if lookup_variant(rsp_data, m) is not None))
        signal_rows.append({
            'Trait': trait,
            'Population': pop,
            'Locus': locus,
            'Merged_Signal': signal,
            'n_original_variants': len(mvp_ids),
            'mean_maf': mafs.mean(),
            'median_maf': mafs.median(),
            'susie_file_exists': susie_data is not None,
            'rsparsepro_file_exists': rsp_data is not None,
            'any_variant_in_susie_cs': any_susie,
            'any_variant_in_rsparsepro_cs': any_rsp,
        })

    signal_out = f'{output_prefix}.signal_recovery.tsv'
    pd.DataFrame(signal_rows).to_csv(signal_out, sep='\t', index=False)
    print(f'  Wrote {len(signal_rows)} signals to {signal_out}', flush=True)

    # Part 2: PIP Comparison — one row per original variant
    print('Computing PIP comparison...', flush=True)
    pip_rows = []
    for _, row in df.iterrows():
        trait, pop, stem = row['Trait'], row['Population'], row['locus_stem']
        mvp_id = str(row['MVP ID'])
        susie_data = results_cache.get((trait, pop, stem, 'susie'))
        rsp_data = results_cache.get((trait, pop, stem, 'rsparsepro'))

        s_result = lookup_variant(susie_data, mvp_id) if susie_data is not None else None
        susie_pip, susie_cs = s_result if s_result is not None else (float('nan'), float('nan'))

        r_result = lookup_variant(rsp_data, mvp_id) if rsp_data is not None else None
        rsp_pip, rsp_cs = r_result if r_result is not None else (float('nan'), float('nan'))

        eaf = pd.to_numeric(row['EAF Population'], errors='coerce')
        maf = min(eaf, 1 - eaf) if pd.notna(eaf) else float('nan')
        pip_rows.append({
            'Trait': trait,
            'Population': pop,
            'Locus': row['Locus'],
            'Merged_Signal': row['Merged Signal'],
            'RSID': row['RSID'],
            'MVP_ID': mvp_id,
            'maf': maf,
            'original_pip': row['Overall PIP'],
            'susie_pip': susie_pip,
            'susie_cs': susie_cs,
            'rsparsepro_pip': rsp_pip,
            'rsparsepro_cs': rsp_cs,
        })

    pip_out = f'{output_prefix}.pip_comparison.tsv'
    pd.DataFrame(pip_rows).to_csv(pip_out, sep='\t', index=False, na_rep='NA')
    print(f'  Wrote {len(pip_rows)} variants to {pip_out}', flush=True)


if __name__ == '__main__':
    main()
