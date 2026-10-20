#!/usr/bin/env python3
'''
build_cohort_sample_sizes.py

Extracts the per-cohort sample size (N, N_CASE, N_CONTROL) for every fine-mapped
trait x ancestry group, by reading the first data row of each cohort-level GWAS
summary statistics file linked to the meta-analysis.

Input:
  --finemapped-variants  gwas_finemapped_variants.tsv from step 11 (credible-set
                          collation); only its trait and pop columns are used.
  --linker                meta_sumstats.linker.betanfer.tsv: one row per
                          "<trait>.<pop>" key, mapping to the cohort sumstats file
                          relative path "<cohort>/<cohort_ancestry>/<file>".
  --sumstats-dir          directory holding the files named in --linker; each file
                          is bgzip-compressed with a '#'-prefixed header that
                          includes N and N_CASE.
  --pipeline-sample-sizes optional cross-check: the meta-analysis pipeline's own
                          per-cohort N table (TRAIT, COHORT, COHORT_POP, N). Only
                          used to print how many rows agree; never asserted.

Output (--output, default cohort_sample_sizes.tsv):
  trait, pop, cohort, cohort_ancestry, n, n_case, n_control, file
  One row per trait x ancestry group x contributing cohort.

Usage:
  python build_cohort_sample_sizes.py --finemapped-variants <path> --output <path>

***Run in finemap conda environment (adds openpyxl: finemap_old)***
'''

import sys
import getopt
import gzip
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

FINEMAPPED_VARIANTS_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/credible_sets/gwas_finemapped_variants.tsv'
LINKER_DEFAULT = '/flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv'
SUMSTATS_DIR_DEFAULT = '/flare/GeomicVar/multi-biobank-meta/sumstats'
PIPELINE_SAMPLE_SIZES_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices/sample_sizes.tsv'
OUTPUT_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/supplementary_tables/cohort_sample_sizes.tsv'
THREADS_DEFAULT = 16


def help(exit_num=0):
    print("""Usage: python build_cohort_sample_sizes.py [OPTIONS]
Options:
  --finemapped-variants    gwas_finemapped_variants.tsv (default: application_results/credible_sets/...)
  --linker                 meta_sumstats.linker.betanfer.tsv (default: """ + LINKER_DEFAULT + """)
  --sumstats-dir            directory of cohort-level summary statistics (default: """ + SUMSTATS_DIR_DEFAULT + """)
  --pipeline-sample-sizes   pipeline N table used only for a cross-check NOTE (default: """ + PIPELINE_SAMPLE_SIZES_DEFAULT + """)
  --output                  output TSV path (default: """ + OUTPUT_DEFAULT + """)
  --threads                 parallel file reads (default: 16)
  -h, --help                show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  COHORT FILE READING  #############################################
###############################################################################
def first_data_row(sumstats_dir, rel_path):
    """Read N and N_CASE from the first data row of one cohort sumstats file.

    Parameters: sumstats_dir (str), rel_path (str, relative to sumstats_dir).
    Returns: (rel_path, n, n_case).
    """
    with gzip.open(f'{sumstats_dir}/{rel_path}', 'rt') as fh:
        header = fh.readline().rstrip('\n').lstrip('#').split('\t')
        row = dict(zip(header, fh.readline().rstrip('\n').split('\t')))
    return rel_path, int(float(row['N'])), int(float(row['N_CASE']))


def load_cohort_links(linker_path, keys):
    """Read the meta-analysis linker and keep only the given (trait, pop) keys.

    Parameters: linker_path (str), keys (set of (trait, pop)).
    Returns: DataFrame with trait, pop, cohort, cohort_ancestry, file.
    """
    link = pd.read_csv(linker_path, sep='\t', header=None, names=['key', 'file'])
    link[['trait', 'pop']] = link.key.str.rsplit('.', n=1, expand=True)
    link = link[[k in keys for k in zip(link.trait, link['pop'])]].copy()
    parts = link.file.str.split('/', expand=True)
    link['cohort'], link['cohort_ancestry'] = parts[0], parts[1]
    assert set(zip(link.trait, link['pop'])) == keys, 'every fine-mapped trait x pop needs linker rows'
    return link


def check_against_pipeline(out, pipeline_sample_sizes_path):
    """Print how many (trait, cohort, cohort_ancestry) rows agree with the
    pipeline's own N table. Diagnostic only; never raises.
    """
    ref = pd.read_csv(pipeline_sample_sizes_path, sep='\t').rename(
        columns={'TRAIT': 'trait', 'COHORT': 'cohort', 'COHORT_POP': 'cohort_ancestry', 'N': 'n_pipeline'})
    ref['cohort'] = ref.cohort.str.upper()
    chk = out.drop_duplicates(['trait', 'cohort', 'cohort_ancestry']).merge(
        ref, on=['trait', 'cohort', 'cohort_ancestry'], how='left')
    mismatches = (chk.n_pipeline.notna() & (chk.n != chk.n_pipeline)).sum()
    print(f'NOTE: pipeline N available for {chk.n_pipeline.notna().sum()}/{len(chk)} rows; '
          f'{mismatches} mismatches')


def main():
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'h', ['help', 'finemapped-variants=', 'linker=', 'sumstats-dir=',
                                                     'pipeline-sample-sizes=', 'output=', 'threads='])
    except getopt.GetoptError as err:
        print(f'ERROR: {err}')
        help(2)

    finemapped_variants = FINEMAPPED_VARIANTS_DEFAULT
    linker_path = LINKER_DEFAULT
    sumstats_dir = SUMSTATS_DIR_DEFAULT
    pipeline_sample_sizes = PIPELINE_SAMPLE_SIZES_DEFAULT
    output = OUTPUT_DEFAULT
    threads = THREADS_DEFAULT
    for opt, val in opts:
        if opt in ('-h', '--help'):
            help()
        elif opt == '--finemapped-variants':
            finemapped_variants = val
        elif opt == '--linker':
            linker_path = val
        elif opt == '--sumstats-dir':
            sumstats_dir = val
        elif opt == '--pipeline-sample-sizes':
            pipeline_sample_sizes = val
        elif opt == '--output':
            output = val
        elif opt == '--threads':
            threads = int(val)

    fm = pd.read_csv(finemapped_variants, sep='\t', usecols=['trait', 'pop'])
    keys = set(zip(fm.trait, fm['pop']))
    link = load_cohort_links(linker_path, keys)

    print(f'NOTE: reading {link.file.nunique()} cohort sumstats files with {threads} threads', flush=True)
    with ThreadPoolExecutor(threads) as pool:
        rows = pd.DataFrame(
            list(pool.map(lambda f: first_data_row(sumstats_dir, f), sorted(link.file.unique()))),
            columns=['file', 'n', 'n_case'])
    out = link.merge(rows, on='file', validate='many_to_one')
    out['n_control'] = out.n - out.n_case
    assert (out.n_case > 0).all() and (out.n_control > 0).all(), 'every cohort needs cases and controls'

    check_against_pipeline(out, pipeline_sample_sizes)

    cols = ['trait', 'pop', 'cohort', 'cohort_ancestry', 'n', 'n_case', 'n_control', 'file']
    out[cols].sort_values(cols[:4]).to_csv(output, sep='\t', index=False)
    print(f'NOTE: {len(out)} trait x pop x cohort rows, {out.file.nunique()} files -> {output}', flush=True)


if __name__ == '__main__':
    main()
