#!/usr/bin/env python3
'''
build_supplementary_data.py

Assembles Supplementary Data 1 (the complete PLATLAS fine-mapping results) from
the credible-set collation (step 11) and gene nomination (step 12) outputs.

Input:
  --credible-sets-dir   directory with gwas_finemapped_variants.tsv,
                         gwas_cs_members.tsv and cs_match_status_native.tsv.gz
                         (step 11; "pop"/"locus"/"tool" columns, SuSiEx credible
                         sets already following the SuSiEx 95% rule).
  --gene-nomination-dir  directory with gene_nominations_comprehensive.tsv
                          (step 12).
  --cohort-sample-sizes  per-cohort N/N_CASE/N_CONTROL table
                          (build_cohort_sample_sizes.py).
  --phecode-catalog, --panukb-manifest, --phecode-manual
                          phenotype description/category lookups (see each
                          file's header for its source).

Output (--output-dir, one package directory plus its tar):
  <package-name>/phenotypes.tsv.gz             one row per trait x ancestry group
  <package-name>/cohort_sample_sizes.tsv.gz    one row per trait x ancestry x cohort
  <package-name>/regions.tsv.gz                one row per trait x ancestry x region
  <package-name>/credible_sets.tsv.gz          one row per credible set (CS)
  <package-name>/credible_set_members.tsv.gz   one row per CS member
  <package-name>/variant_pips.tsv.gz           one row per reported variant
  <package-name>/gene_nominations.tsv.gz       one row per gene nomination
  <package-name>/README.txt, md5sums.txt       filled from --readme-template
  <package-name>.tar                           the package directory, tarred

Column definitions are in --readme-template (data_package_readme_template.txt),
which becomes the package's own README.txt.

Usage:
  python build_supplementary_data.py --output-dir <path>

***Run in finemap conda environment (adds openpyxl: finemap_old)***
'''

import sys
import getopt
import hashlib
import json
import shutil
from pathlib import Path
import tarfile

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent

CREDIBLE_SETS_DIR_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/credible_sets'
GENE_NOMINATION_DIR_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination'
COHORT_SAMPLE_SIZES_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/supplementary_tables/cohort_sample_sizes.tsv'
PHECODE_CATALOG_DEFAULT = str(HERE / 'phecode_catalog.tsv')
PANUKB_MANIFEST_DEFAULT = str(HERE / 'panukb_phecode_manifest.tsv')
PHECODE_MANUAL_DEFAULT = str(HERE / 'phecode_manual.tsv')
README_TEMPLATE_DEFAULT = str(HERE / 'data_package_readme_template.txt')
OUTPUT_DIR_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/supplementary_tables'
PACKAGE_NAME_DEFAULT = 'betainfer_platlas_finemapping'

TOOLS = {'rsparsepro': 'RSparsePro', 'susiex': 'SuSiEx', 'multisusie': 'MultiSuSiE'}
EVIDENCE = {'nonsynonymous_coding_variant': 'Coding', 'abc_interaction': 'ABC link', 'abc_promoter_overlap': 'ABC promoter'}
POPS = ['EUR', 'AFR', 'ALL']
GZIP_DETERMINISTIC = {'method': 'gzip', 'mtime': 0}


def help(exit_num=0):
    print("""Usage: python build_supplementary_data.py [OPTIONS]
Options:
  --credible-sets-dir     step 11 output directory (default: """ + CREDIBLE_SETS_DIR_DEFAULT + """)
  --gene-nomination-dir   step 12 output directory (default: """ + GENE_NOMINATION_DIR_DEFAULT + """)
  --cohort-sample-sizes   per-cohort N table (default: """ + COHORT_SAMPLE_SIZES_DEFAULT + """)
  --phecode-catalog       phecode 1.2 PheWAS Catalog map (default: next to this script)
  --panukb-manifest       Pan-UK Biobank phenotype manifest (default: next to this script)
  --phecode-manual        manually curated phecodes (default: next to this script)
  --readme-template       package README.txt template (default: next to this script)
  --output-dir            output directory (default: """ + OUTPUT_DIR_DEFAULT + """)
  --package-name          package directory / tar basename (default: """ + PACKAGE_NAME_DEFAULT + """)
  -h, --help              show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  PHENOTYPE METADATA  ###############################################
###############################################################################
def phecode_of(trait):
    """PLATLAS trait identifier (Phe_250_2) -> phecode string (250.2)."""
    parts = trait.removeprefix('Phe_').split('_')
    return parts[0] + ('.' + ''.join(parts[1:]) if len(parts) > 1 else '')


def phenotype_catalog(traits, phecode_catalog_path, panukb_manifest_path, phecode_manual_path):
    """Phecode -> description/category, preferring the PheWAS Catalog, then the
    Pan-UK Biobank manifest, then the manually curated table.

    Returns: DataFrame with trait, phecode, description, category, non_catalog.
    """
    cat = pd.read_csv(phecode_catalog_path, sep='\t', dtype=str).rename(columns={'group': 'category'})
    pan = pd.read_csv(panukb_manifest_path, sep='\t', dtype=str)
    manual = pd.read_csv(phecode_manual_path, sep='\t', dtype=str)
    out = pd.DataFrame({'trait': sorted(traits)})
    out['phecode'] = out.trait.map(phecode_of)
    out['description'], out['category'] = None, None
    for src in (cat, pan, manual):
        idx = src.drop_duplicates('phecode').set_index('phecode')
        for col in ('description', 'category'):
            out[col] = out[col].fillna(out.phecode.map(idx[col]))
    assert out.description.notna().all() and out.category.notna().all()
    out['non_catalog'] = ~out.phecode.isin(cat.phecode)
    return out


def order_key(df):
    return df.assign(_p=df.phecode.astype(float), _a=df.ancestry.map(POPS.index))


###############################################################################
##########  CREDIBLE-SET TABLE  ################################################
###############################################################################
def credible_sets(members, variants, jaccard_path):
    """Credible-set members and one row per CS (lead variant, size, ABF sharing).

    CS members keep their own method's PIP; the lead variant is the member with
    the highest PIP. ABF sharing comes from the PIP-weighted Jaccard comparison
    (step 11), joined on trait/ancestry/region/method/cs.

    Parameters: members, variants (DataFrames); jaccard_path (str).
    Returns: (members without cs_size, one row per CS).
    """
    m = members.merge(variants[['trait', 'region', 'ancestry', 'variant_id', *[f'{t}_pip' for t in TOOLS]]],
                      on=['trait', 'region', 'ancestry', 'variant_id'], how='left', validate='many_to_one')
    m['pip'] = np.select([m.method == t for t in TOOLS], [m[f'{t}_pip'] for t in TOOLS], np.nan)
    assert m.pip.notna().all(), 'every CS member needs its method PIP'
    m = m[['trait', 'ancestry', 'region', 'method', 'cs', 'variant_id', 'pip', 'cs_size']]
    keys = ['trait', 'ancestry', 'region', 'method', 'cs']
    lead = m.sort_values(keys + ['pip', 'variant_id'], ascending=[True] * 5 + [False, True]).drop_duplicates(keys)
    cs = lead.rename(columns={'variant_id': 'lead_variant', 'pip': 'lead_pip'})
    size = m.groupby(keys).size().rename('n').reset_index()
    cs = cs.merge(size, on=keys, validate='one_to_one')
    assert (cs.n == cs.cs_size).all()
    cs = cs.drop(columns='n').rename(columns={'cs_size': 'n_variants'})

    jac = pd.read_csv(jaccard_path, sep='\t').rename(columns={'pop': 'ancestry', 'tool': 'method', 'locus': 'region',
                                                               'best_pip_weighted_jaccard': 'best_jaccard_abf', 'shared': 'shared_with_abf'})
    cs = cs.merge(jac, on=keys, how='left', validate='one_to_one', indicator=True)
    assert cs._merge.value_counts()['both'] == len(jac), 'every compared CS must map to a collated CS'
    cs['shared_with_abf'] = cs.shared_with_abf.map({True: 'yes', False: 'no'})
    cols = keys + ['n_variants', 'lead_variant', 'lead_pip', 'best_jaccard_abf', 'shared_with_abf']
    return m.drop(columns='cs_size'), cs[cols]


###############################################################################
##########  REGION-LEVEL TABLE  ################################################
###############################################################################
def regions(variants, cs, noms):
    """One row per trait x ancestry x region: variant/CS counts, ABF sharing,
    the top variant by median PIP, and nominated genes.
    """
    keys = ['trait', 'ancestry', 'region']
    top = variants.sort_values(keys + ['median_pip', 'variant_id'], ascending=[True] * 3 + [False, True]).drop_duplicates(keys)
    out = variants.groupby(keys).size().rename('n_variants').reset_index()
    out = out.merge(top[keys + ['variant_id', 'median_pip']].rename(columns={'variant_id': 'top_variant', 'median_pip': 'top_median_pip'}), on=keys)
    n_cs = cs.pivot_table(index=keys, columns='method', values='cs', aggfunc='size').add_prefix('n_cs_')
    compared = cs[cs.best_jaccard_abf.notna()]
    n_shared = compared.assign(s=compared.shared_with_abf == 'yes').pivot_table(index=keys, columns='method', values='s', aggfunc='sum').add_prefix('n_cs_shared_abf_')
    out = out.merge(n_cs, on=keys, how='left').merge(n_shared, on=keys, how='left')
    for t in TOOLS:
        out[f'n_cs_{t}'] = out[f'n_cs_{t}'].fillna(0).astype(int)
        out[f'n_cs_shared_abf_{t}'] = out[f'n_cs_shared_abf_{t}'].astype('Int64')
    out = out.merge(cs.groupby(keys).n_variants.min().rename('smallest_cs_size').reset_index(), on=keys, how='left')
    out['smallest_cs_size'] = out.smallest_cs_size.astype('Int64')
    genes = noms.groupby(keys).gene.agg(lambda g: ';'.join(sorted(set(g)))).rename('nominated_genes').reset_index()
    out = out.merge(genes, on=keys, how='left')
    cols = [f'n_cs_{t}' for t in TOOLS] + [f'n_cs_shared_abf_{t}' for t in TOOLS]
    return out[keys + ['n_variants'] + cols + ['smallest_cs_size', 'top_variant', 'top_median_pip', 'nominated_genes']]


###############################################################################
##########  GENE NOMINATION TABLE  #############################################
###############################################################################
def nominations(variants, gene_nomination_path):
    """Gene nominations with their region and per-method PIP, duplicates dropped.

    Returns: (DataFrame, number of exact-duplicate rows removed).
    """
    raw = pd.read_csv(gene_nomination_path, sep='\t').rename(columns={'pop': 'ancestry'})
    noms = raw.drop_duplicates()
    noms = noms.merge(variants[['trait', 'ancestry', 'variant_id', 'region', *[f'{t}_pip' for t in TOOLS]]],
                      on=['trait', 'ancestry', 'variant_id'], how='left', validate='many_to_one')
    assert noms.region.notna().all()
    noms['evidence'] = noms.implication_type.map(EVIDENCE)
    assert noms.evidence.notna().all()
    noms = noms.rename(columns={'abc_max_score': 'abc_score', 'abc_n_celltypes': 'abc_n_cell_types'})
    noms['abc_n_cell_types'] = noms.abc_n_cell_types.astype('Int64')
    cols = ['gene', 'trait', 'ancestry', 'region', 'variant_id', 'evidence', 'consequence', 'abc_score', 'abc_n_cell_types',
            'median_pip', *[f'{t}_pip' for t in TOOLS]]
    return noms[cols], len(raw) - len(noms)


###############################################################################
##########  PHENOTYPE TABLE  ###################################################
###############################################################################
def phenotypes(catalog, variants, cs, noms, cohort_sample_sizes_path):
    """One row per trait x ancestry group: cohorts, cases/controls, region/CS/
    gene-nomination counts.

    Returns: (phenotype table, per-cohort sample-size table without helper columns).
    """
    sizes = pd.read_csv(cohort_sample_sizes_path, sep='\t').rename(columns={'pop': 'ancestry'})
    sizes['cohort'] = sizes.cohort.map({'MVP': 'MVP', 'UKB': 'UKB', 'FINNGEN': 'FinnGen'})
    sizes['label'] = sizes.cohort + '-' + sizes.cohort_ancestry
    sizes = sizes.sort_values(['trait', 'ancestry', 'label'])
    sizes['detail'] = sizes.label + ' ' + sizes.n_case.map('{:,}'.format) + '/' + sizes.n_control.map('{:,}'.format)
    keys = ['trait', 'ancestry']
    out = sizes.groupby(keys).agg(cohorts=('label', ';'.join), n_cases=('n_case', 'sum'), n_controls=('n_control', 'sum'),
                                  cases_controls_by_cohort=('detail', '; '.join)).reset_index()
    out = catalog.drop(columns='non_catalog').merge(out, on='trait')
    out = out.merge(variants.groupby(keys).region.nunique().rename('n_regions').reset_index(), on=keys, validate='one_to_one')
    n_cs = cs.pivot_table(index=keys, columns='method', values='cs', aggfunc='size').add_prefix('n_cs_').reset_index()
    out = out.merge(n_cs, on=keys, how='left').merge(noms.groupby(keys).gene.nunique().rename('n_nominated_genes').reset_index(), on=keys, how='left')
    for c in [f'n_cs_{t}' for t in TOOLS] + ['n_nominated_genes']:
        out[c] = out[c].fillna(0).astype(int)
    out = order_key(out).sort_values(['_p', 'trait', '_a'])
    cols = ['trait', 'phecode', 'description', 'category', 'ancestry', 'cohorts', 'n_cases', 'n_controls',
            'cases_controls_by_cohort', 'n_regions', *[f'n_cs_{t}' for t in TOOLS], 'n_nominated_genes']
    out = out[cols]
    sizes = sizes.drop(columns=['label', 'detail', 'file'])
    return out, sizes


def attach_catalog(df, catalog):
    return df.merge(catalog[['trait', 'phecode', 'description']], on='trait', how='left', validate='many_to_one')


###############################################################################
##########  PACKAGE OUTPUT  ####################################################
###############################################################################
def write_table(df, package_dir, name, sort):
    df = df if sort is None else df.sort_values(sort)
    df.to_csv(package_dir / name, sep='\t', index=False, float_format='%.6g', compression=GZIP_DETERMINISTIC)
    return len(df)


def md5(path):
    return hashlib.md5(path.read_bytes()).hexdigest()


def main():
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'h', ['help', 'credible-sets-dir=', 'gene-nomination-dir=',
                                                     'cohort-sample-sizes=', 'phecode-catalog=', 'panukb-manifest=',
                                                     'phecode-manual=', 'readme-template=', 'output-dir=',
                                                     'package-name='])
    except getopt.GetoptError as err:
        print(f'ERROR: {err}')
        help(2)

    credible_sets_dir = CREDIBLE_SETS_DIR_DEFAULT
    gene_nomination_dir = GENE_NOMINATION_DIR_DEFAULT
    cohort_sample_sizes_path = COHORT_SAMPLE_SIZES_DEFAULT
    phecode_catalog_path = PHECODE_CATALOG_DEFAULT
    panukb_manifest_path = PANUKB_MANIFEST_DEFAULT
    phecode_manual_path = PHECODE_MANUAL_DEFAULT
    readme_template_path = README_TEMPLATE_DEFAULT
    output_dir = OUTPUT_DIR_DEFAULT
    package_name = PACKAGE_NAME_DEFAULT
    for opt, val in opts:
        if opt in ('-h', '--help'):
            help()
        elif opt == '--credible-sets-dir':
            credible_sets_dir = val
        elif opt == '--gene-nomination-dir':
            gene_nomination_dir = val
        elif opt == '--cohort-sample-sizes':
            cohort_sample_sizes_path = val
        elif opt == '--phecode-catalog':
            phecode_catalog_path = val
        elif opt == '--panukb-manifest':
            panukb_manifest_path = val
        elif opt == '--phecode-manual':
            phecode_manual_path = val
        elif opt == '--readme-template':
            readme_template_path = val
        elif opt == '--output-dir':
            output_dir = val
        elif opt == '--package-name':
            package_name = val

    output_dir = Path(output_dir)
    package_dir = output_dir / package_name
    tar_path = output_dir / f'{package_name}.tar'
    if package_dir.exists():
        shutil.rmtree(package_dir)
    package_dir.mkdir(parents=True)

    credible_sets_dir = Path(credible_sets_dir)
    variants = pd.read_csv(credible_sets_dir / 'gwas_finemapped_variants.tsv', sep='\t').rename(columns={'pop': 'ancestry', 'locus': 'region'})
    members = pd.read_csv(credible_sets_dir / 'gwas_cs_members.tsv', sep='\t').rename(columns={'pop': 'ancestry', 'locus': 'region', 'tool': 'method'})
    gene_nomination_path = Path(gene_nomination_dir) / 'gene_nominations_comprehensive.tsv'
    jaccard_path = credible_sets_dir / 'cs_match_status_native.tsv.gz'

    catalog = phenotype_catalog(variants.trait.unique(), phecode_catalog_path, panukb_manifest_path, phecode_manual_path)
    members, cs = credible_sets(members, variants, jaccard_path)
    noms, n_dup = nominations(variants, gene_nomination_path)
    region_tab = regions(variants, cs, noms)
    pheno, sizes = phenotypes(catalog, variants, cs, noms, cohort_sample_sizes_path)

    # ---- display names (shared by the data package and the workbook) ----
    for df in (members, cs):
        df['method'] = df.method.map(TOOLS)
    region_tab = attach_catalog(region_tab, catalog)
    noms = attach_catalog(noms, catalog)

    rows = {}
    for df, name, sort in [
        (pheno, 'phenotypes.tsv.gz', None),
        (sizes, 'cohort_sample_sizes.tsv.gz', ['trait', 'ancestry', 'cohort', 'cohort_ancestry']),
        (region_tab[['trait', 'phecode', 'description'] + [c for c in region_tab if c not in ('trait', 'phecode', 'description')]], 'regions.tsv.gz', ['trait', 'ancestry', 'region']),
        (cs, 'credible_sets.tsv.gz', ['trait', 'ancestry', 'region', 'method', 'cs']),
        (members, 'credible_set_members.tsv.gz', ['trait', 'ancestry', 'region', 'method', 'cs', 'variant_id']),
        (variants[['trait', 'ancestry', 'region', 'variant_id', *[f'{t}_pip' for t in TOOLS], 'median_pip', 'in_cs_any']], 'variant_pips.tsv.gz', ['trait', 'ancestry', 'region', 'variant_id']),
        (noms[['gene', 'trait', 'phecode', 'description'] + [c for c in noms if c not in ('gene', 'trait', 'phecode', 'description')]], 'gene_nominations.tsv.gz', ['gene', 'trait', 'ancestry', 'variant_id', 'evidence']),
    ]:
        rows[name] = write_table(df, package_dir, name, sort)

    # ---- README and checksums ----
    readme = Path(readme_template_path).read_text().format(**{k.split('.')[0]: f'{v:,}' for k, v in rows.items()})
    (package_dir / 'README.txt').write_text(readme)
    files = sorted(p for p in package_dir.iterdir())
    (package_dir / 'md5sums.txt').write_text(''.join(f'{md5(p)}  {p.name}\n' for p in files))
    with tarfile.open(tar_path, 'w') as tar:
        for p in sorted(package_dir.iterdir()):
            tar.add(p, arcname=f'{package_name}/{p.name}')

    # ---- consistency checks ----
    region_counts = variants.groupby('ancestry')[['trait', 'region']].apply(lambda d: len(d.drop_duplicates())).to_dict()
    checks = {
        'credible_sets': len(cs), 'credible_sets_by_method': cs.method.value_counts().to_dict(),
        'members': len(members), 'variants': len(variants),
        'regions': len(region_tab), 'regions_by_ancestry': region_counts,
        'regions_with_cs': int((region_tab[[f'n_cs_{t}' for t in TOOLS]].sum(axis=1) > 0).sum()),
        'cs_compared_with_abf': int(cs.best_jaccard_abf.notna().sum()),
        'cs_shared_with_abf': int((cs.shared_with_abf == 'yes').sum()),
        'phenotype_ancestry_rows': len(pheno), 'traits': int(pheno.trait.nunique()),
        'traits_described_outside_phewas_catalog': int(catalog.non_catalog.sum()),
        'gene_nominations': len(noms), 'exact_duplicate_nominations_removed': int(n_dup),
        'nominated_genes': int(noms.gene.nunique()),
        'susiex_pip_above_1': int((variants.susiex_pip > 1).sum()),
    }
    assert checks['regions'] == len(variants[['trait', 'ancestry', 'region']].drop_duplicates())
    assert pheno.n_regions.sum() == checks['regions']
    for t, label in TOOLS.items():
        assert region_tab[f'n_cs_{t}'].sum() == pheno[f'n_cs_{t}'].sum() == checks['credible_sets_by_method'][label]
    assert noms.groupby(['trait', 'ancestry']).gene.nunique().sum() == pheno.n_nominated_genes.sum()
    checks['files'] = {p.name: {'bytes': p.stat().st_size, 'md5': md5(p)} for p in files}
    checks['tar'] = {'path': tar_path.name, 'bytes': tar_path.stat().st_size}
    summary_path = output_dir / f'{package_name}_summary.json'
    summary_path.write_text(json.dumps(checks, indent=2) + '\n')
    print(json.dumps({k: v for k, v in checks.items() if k != 'files'}, indent=2))


if __name__ == '__main__':
    main()
