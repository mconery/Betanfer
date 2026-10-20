#!/usr/bin/env python3
'''
build_supplementary_workbook.py

Builds Supplementary Tables 1-4 (Phenotypes/cohorts, Methods, Regions, Gene
nominations) as a single styled Excel workbook from the Supplementary Data 1
package (build_supplementary_data.py) and the Supplementary Table 2 methods
TSV (s2_methods.tsv). Title/header/body fonts, the title fill and the header
rule are the explicit colors and fonts of the published reference
supplementary workbook this table format follows (see the STYLE section).

Input:
  --package-dir       Supplementary Data 1 package directory (phenotypes.tsv.gz,
                       regions.tsv.gz, gene_nominations.tsv.gz).
  --summary-json       the *_summary.json written by build_supplementary_data.py
                        (counts quoted in the table legends).
  --s2-methods         Supplementary Table 2 content (Software/Input/Settings/
                        Output per analysis step; static file next to this
                        script).

Output (--output-dir):
  SupplementaryTables_S1-S4.xlsx  Legend + S1 + S2 + S3 + S4 sheets.
  legend.md                       the same legend text as the Legend sheet.

The script re-reads the workbook it just wrote and asserts every cell matches
what was written, before reporting success; it exits with an assertion error
if not.

Usage:
  python build_supplementary_workbook.py --package-dir <path> --summary-json <path> \\
      --output-dir <path>

***Run in finemap conda environment (adds openpyxl: finemap_old)***
'''

import sys
import getopt
import json
from pathlib import Path
import textwrap

import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Border, Font, NamedStyle, PatternFill, Side
from openpyxl.utils import get_column_letter

HERE = Path(__file__).resolve().parent

PACKAGE_DIR_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/supplementary_tables/betainfer_platlas_finemapping'
SUMMARY_JSON_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/supplementary_tables/betainfer_platlas_finemapping_summary.json'
S2_METHODS_DEFAULT = str(HERE / 's2_methods.tsv')
OUTPUT_DIR_DEFAULT = '/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/supplementary_tables'

POP_ORDER = {'EUR': 0, 'AFR': 1, 'ALL': 2}
CHROM = {str(i): i for i in range(1, 24)}
THIN = Side(style='thin', color='000000')
FORMATS = {'text': ('General', 'left', False), 'wrap': ('General', 'left', True),
           'int': ('#,##0', 'right', False), 'dec': ('0.0000', 'right', False)}
METHODS = [('rsparsepro', 'RSparsePro'), ('susiex', 'SuSiEx'), ('multisusie', 'MultiSuSiE')]


def help(exit_num=0):
    print("""Usage: python build_supplementary_workbook.py [OPTIONS]
Options:
  --package-dir         Supplementary Data 1 package directory (default: """ + PACKAGE_DIR_DEFAULT + """)
  --summary-json        build_supplementary_data.py's summary JSON (default: """ + SUMMARY_JSON_DEFAULT + """)
  --s2-methods           Supplementary Table 2 content TSV (default: next to this script)
  --output-dir           output directory (default: """ + OUTPUT_DIR_DEFAULT + """)
  -h, --help              show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  COLUMN SPECS: (header, source column, width, format)  ############
###############################################################################
S1 = [('Trait ID', 'trait', 11, 'text'), ('Phecode', 'phecode', 9, 'text'), ('Phenotype', 'description', 36, 'wrap'),
      ('Category', 'category', 22, 'wrap'), ('Ancestry group', 'ancestry', 10, 'text'), ('Cohorts', 'cohorts', 27, 'wrap'),
      ('Cases', 'n_cases', 11, 'int'), ('Controls', 'n_controls', 12, 'int'),
      ('Cases/controls by cohort', 'cases_controls_by_cohort', 58, 'wrap'), ('Fine-mapped regions', 'n_regions', 11, 'int'),
      *[(f'{name} CS', f'n_cs_{m}', 11, 'int') for m, name in METHODS], ('Nominated genes', 'n_nominated_genes', 11, 'int')]
S2_WIDTHS = [20, 30, 45, 80, 40]
S3 = [('Trait ID', 'trait', 11, 'text'), ('Phenotype', 'description', 40, 'text'), ('Ancestry group', 'ancestry', 10, 'text'),
      ('Region (GRCh38)', 'region', 25, 'text'), ('Reported variants', 'n_variants', 10, 'int'),
      *[(f'{name} CS', f'n_cs_{m}', 11, 'int') for m, name in METHODS],
      *[(f'{name} CS shared with ABF', f'n_cs_shared_abf_{m}', 13, 'int') for m, name in METHODS],
      ('Smallest CS size', 'smallest_cs_size', 10, 'int'), ('Top variant (GRCh38)', 'top_variant', 22, 'text'),
      ('Top variant median PIP', 'top_median_pip', 11, 'dec'), ('Nominated genes', 'nominated_genes', 40, 'text')]
S4 = [('Gene', 'gene', 12, 'text'), ('Trait ID', 'trait', 11, 'text'), ('Phenotype', 'description', 40, 'text'),
      ('Ancestry group', 'ancestry', 10, 'text'), ('Region (GRCh38)', 'region', 25, 'text'),
      ('Variant (GRCh38)', 'variant_id', 22, 'text'), ('Evidence', 'evidence', 13, 'text'),
      ('VEP consequence', 'consequence', 22, 'text'), ('ABC score', 'abc_score', 9, 'dec'),
      ('ABC cell types', 'abc_n_cell_types', 9, 'int'), ('Median PIP', 'median_pip', 9, 'dec'),
      *[(f'{name} PIP', f'{m}_pip', 11, 'dec') for m, name in METHODS]]


###############################################################################
##########  DATA: read the package's gz TSVs, display-only transforms  #######
###############################################################################
def read(package_dir, name):
    """Read one gzipped TSV from the Supplementary Data 1 package."""
    return pd.read_csv(package_dir / f'{name}.tsv.gz', sep='\t', dtype={'phecode': str}, keep_default_na=False,
                       na_values=[''], float_precision='round_trip')


def ordered(df, position, first=None):
    """Sort rows by phecode, ancestry, then chromosome:position parsed from
    the GRCh38 position-like column (variant_id or region).
    """
    parts = position.str.removeprefix('chr').str.split(r'[_:]', regex=True)
    keys = pd.DataFrame({'phe': df.phecode.astype(float), 'pop': df.ancestry.map(POP_ORDER),
                         'chr': parts.str[0].map(CHROM), 'pos': parts.str[1].astype(int)}, index=df.index)
    if first:
        keys.insert(0, 'first', df[first].str.upper())
    assert keys.notna().all().all()
    return df.loc[keys.sort_values(list(keys.columns), kind='stable').index].reset_index(drop=True)


def load(package_dir, s2_methods_path):
    """Load and order the four table contents: S1 phenotypes, S2 methods
    (static), S3 regions, S4 gene nominations.
    """
    s1 = read(package_dir, 'phenotypes').assign(cohorts=lambda d: d.cohorts.str.replace(';', '; '))
    s2 = pd.read_csv(s2_methods_path, sep='\t', dtype=str)
    s3 = read(package_dir, 'regions').assign(nominated_genes=lambda d: d.nominated_genes.str.replace(';', '; '))
    s3 = ordered(s3, s3.region)
    s4 = read(package_dir, 'gene_nominations')
    s4 = ordered(s4, s4.variant_id, first='gene')
    return s1, s2, s3, s4


###############################################################################
##########  LEGEND: the same text goes into the Legend sheet and legend.md  ##
###############################################################################
def legends(s1, s3, s4, meta):
    """Build the (label, title, description) entries for all five tables."""
    shared_note = ('blank where no comparison was made (SuSiEx and MultiSuSiE in EUR and AFR, and phenotypes outside '
                   'the comparison set; Supplementary Table 2)')
    return [
        ('Supplementary Table 1', 'Phenotypes, cohorts and fine-mapping yield',
         f'One row per phenotype and ancestry group ({len(s1):,} rows; {s1.trait.nunique():,} phenotypes). Phenotypes are '
         f'phecodes; descriptions and categories follow the phecode 1.2 map of the PheWAS Catalog, except for '
         f'{meta["traits_described_outside_phewas_catalog"]} phecodes absent from that map, which are described from the '
         'Pan-UK Biobank phenotype manifest (19) or the MVP phenotype list (1). Cohorts lists the cohort–ancestry inputs to '
         'the meta-analysis. Cases and controls are summed across cohorts and also given per cohort (cases/controls). '
         'Fine-mapped regions is the number of regions fine-mapped (Supplementary Table 3). Method columns give the number '
         'of 95% credible sets (CS). Nominated genes is the number of distinct genes nominated (Supplementary Table 4).'),
        ('Supplementary Table 2', 'Analysis steps, software and settings',
         'Software, versions, inputs, settings and outputs of each analysis step, from meta-analysis to gene nomination '
         'and downstream analyses. Commit identifiers are given for software installed from source repositories.'),
        ('Supplementary Table 3', 'Region-level fine-mapping results',
         f'One row per phenotype, ancestry group and fine-mapped region ({len(s3):,} rows; {meta["regions_with_cs"]:,} '
         'with at least one CS). Region identifiers (chr<chromosome>_<start>_<end>, GRCh38; chromosome X is coded as 23) match those in Supplementary '
         'Data 1. Reported variants is the number of variants that belong to a CS from any method or have PIP ≥ 0.01 from '
         'at least one method. Method CS columns give the number of 95% CS. CS shared with ABF gives the number of those '
         'CS with PIP-weighted Jaccard similarity > 0.10 to at least one PLATLAS approximate Bayes factor (ABF) CS for the '
         f'same phenotype and ancestry group; it is {shared_note}. Smallest CS size is the number of variants in the '
         'smallest CS from any method (blank if none). Top variant is the variant with the highest median PIP across '
         'methods. Nominated genes lists genes nominated in the region (Supplementary Table 4).'),
        ('Supplementary Table 4', 'Gene nominations',
         f'One row per nominated gene, phenotype, ancestry group, variant and evidence type ({len(s4):,} rows; '
         f'{s4.gene.nunique():,} genes). Variants with median PIP ≥ 0.01 across the available methods were eligible. '
         'Evidence: Coding, a missense, frameshift, stop-gained, stop-lost, start-lost, splice-acceptor or splice-donor '
         'consequence in a protein-coding gene (VEP consequence); ABC link, the variant lies in an activity-by-contact '
         '(ABC) element linked to the gene (ABC score > 0.1); ABC promoter, the variant lies in the promoter element of '
         'the gene. Coding evidence takes precedence for the same variant and gene. ABC score is the maximum across cell '
         'types and ABC cell types the number of cell types with a link. Method PIPs are blank where the method did not '
         'model the variant; SuSiEx PIPs are sums of single-effect probabilities and can exceed 1.'),
        ('Supplementary Data 1', 'Complete fine-mapping results (betainfer_platlas_finemapping)',
         f'Tab-separated files with a README describing all columns: phenotypes, per-cohort sample sizes, regions, CS '
         f'({meta["credible_sets"]:,}), CS members ({meta["members"]:,}), variant PIPs ({meta["variants"]:,}) and gene '
         'nominations. CS-level fields include the lead variant and the best PIP-weighted Jaccard similarity to ABF CS.'),
    ]


ABBREVIATIONS = ('Abbreviations: ABC, activity-by-contact; ABF, approximate Bayes factor; AFR, African; ALL, all ancestry '
                 'groups combined (EUR, AFR and UK Biobank South Asian); CS, credible set; EUR, European; MVP, Million '
                 'Veteran Program; PIP, posterior inclusion probability; SAS, South Asian; UKB, UK Biobank; VEP, Variant '
                 'Effect Predictor.')


###############################################################################
##########  STYLE: explicit title/header/body fonts, fill and border  #######
###############################################################################
# Colors are resolved once from the published reference supplementary workbook
# this table format follows, so no copy of that file is needed at run time:
# its theme is the default Office theme, where every font color is dk1 at
# tint 0 (black, #000000), and the title band's fill is accent3 (#A5A5A5)
# lightened by tint 0.6 toward white, which resolves to #DBDBDB. Borders use
# the automatic color, i.e. black. Everything else (alignment, number_format)
# the reference supplied was either 'General'/unset or is overridden right
# after by this script, so only these four values carry the reference's
# appearance forward.
TITLE_FONT = Font(name='Calibri', size=11, bold=True, color='000000')
TITLE_FILL = PatternFill(fill_type='solid', fgColor='DBDBDB')
HEADER_FONT = Font(name='Calibri', size=11, bold=True, color='000000')
HEADER_BORDER = Border(bottom=THIN)
BODY_FONT = Font(name='Calibri', size=11, bold=False, color='000000')
LEGEND_HEADER_ALIGN = [Alignment(horizontal='right'), Alignment(horizontal='center'), Alignment(horizontal='center')]


def lines(text, width):
    return len(textwrap.wrap(str(text), max(8, int(width) - 3))) or 1


def add_styles(wb):
    """Register one named style per column format, using the fixed body font."""
    for key, (fmt, horizontal, wrap) in FORMATS.items():
        wb.add_named_style(NamedStyle(name=f'body_{key}', font=BODY_FONT, border=Border(), number_format=fmt,
                                      alignment=Alignment(horizontal=horizontal, vertical='top', wrap_text=wrap)))


def make_sheet(wb, name, title, widths):
    """Create a sheet with the title-row style and print/freeze settings."""
    ws = wb.create_sheet(name)
    for i, width in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = width
        cell = ws.cell(1, i)
        cell.font = TITLE_FONT
        cell.fill = TITLE_FILL
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(widths))
    ws.cell(1, 1, title).alignment = Alignment(horizontal='left', vertical='center', wrap_text=True)
    ws.row_dimensions[1].height = 22
    ws.row_dimensions[2].height = 15
    ws.freeze_panes = 'A4'
    ws.sheet_properties.pageSetUpPr.fitToPage = True
    ws.page_setup.orientation = 'landscape'
    ws.page_setup.paperSize = ws.PAPERSIZE_A3
    ws.page_setup.fitToWidth = 1
    ws.page_setup.fitToHeight = 0
    ws.print_options.gridLines = False
    ws.print_title_rows = '1:3'
    return ws


def table_sheet(wb, name, title, df, spec):
    """Write one data table sheet; returns the written rows for a read-back check."""
    ws = make_sheet(wb, name, title, [w for _, _, w, _ in spec])
    for i, (label, _, _, _) in enumerate(spec, 1):
        c = ws.cell(3, i, label)
        c.font = HEADER_FONT
        c.border = HEADER_BORDER
        c.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
    ws.row_dimensions[3].height = 15 * max(lines(label, w) for label, _, w, _ in spec) + 8
    cols = [col for _, col, _, _ in spec]
    names = [f'body_{kind}' for _, _, _, kind in spec]
    wrapped = [(i, w) for i, (_, _, w, kind) in enumerate(spec) if kind == 'wrap']
    values = df[cols].astype(object).where(df[cols].notna(), None)
    for r, row in enumerate(values.itertuples(index=False, name=None), 4):
        for c, (value, named) in enumerate(zip(row, names), 1):
            ws.cell(r, c, value).style = named
        if wrapped:
            ws.row_dimensions[r].height = max(18, 15 * max(lines(row[i] or '', w) for i, w in wrapped) + 5)
    for c in range(1, len(spec) + 1):
        ws.cell(3 + len(df), c).border = Border(bottom=THIN)
    ws.print_area = f'A1:{get_column_letter(len(spec))}{3 + len(df)}'
    return [tuple(r) for r in values.itertuples(index=False, name=None)]


def legend_sheet(wb, entries):
    """Write the Legend sheet (Table/Title/Description) plus the abbreviations footer."""
    widths = [22, 44, 120]
    ws = make_sheet(wb, 'Legend', 'Supplementary Tables 1–4 and Supplementary Data 1', widths)
    for i, label in enumerate(['Table', 'Title', 'Description'], 1):
        c = ws.cell(3, i, label)
        c.font = HEADER_FONT
        c.border = HEADER_BORDER
        c.alignment = LEGEND_HEADER_ALIGN[i - 1]
    for r, values in enumerate(entries, 4):
        for c, value in enumerate(values, 1):
            cell = ws.cell(r, c, value)
            cell.font = BODY_FONT
            cell.alignment = Alignment(horizontal='left', vertical='top', wrap_text=True)
        ws.row_dimensions[r].height = 15 * max(lines(v, w) for v, w in zip(values, widths)) + 6
    row = 5 + len(entries)
    ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=3)
    c = ws.cell(row, 1, ABBREVIATIONS)
    c.font = BODY_FONT
    c.alignment = Alignment(horizontal='left', vertical='top', wrap_text=True)
    ws.row_dimensions[row].height = 15 * lines(ABBREVIATIONS, sum(widths) * .85) + 6
    ws.freeze_panes = None
    ws.print_area = f'A1:C{row}'


###############################################################################
##########  MAIN: build -> save -> read back and check -> legend.md  #########
###############################################################################
def main():
    try:
        opts, _ = getopt.getopt(sys.argv[1:], 'h', ['help', 'package-dir=', 'summary-json=', 's2-methods=',
                                                     'output-dir='])
    except getopt.GetoptError as err:
        print(f'ERROR: {err}')
        help(2)

    package_dir = PACKAGE_DIR_DEFAULT
    summary_json_path = SUMMARY_JSON_DEFAULT
    s2_methods_path = S2_METHODS_DEFAULT
    output_dir = OUTPUT_DIR_DEFAULT
    for opt, val in opts:
        if opt in ('-h', '--help'):
            help()
        elif opt == '--package-dir':
            package_dir = val
        elif opt == '--summary-json':
            summary_json_path = val
        elif opt == '--s2-methods':
            s2_methods_path = val
        elif opt == '--output-dir':
            output_dir = val

    package_dir = Path(package_dir)
    output_dir = Path(output_dir)

    meta = json.loads(Path(summary_json_path).read_text())
    s1, s2, s3, s4 = load(package_dir, s2_methods_path)
    entries = legends(s1, s3, s4, meta)
    wb = Workbook()
    wb.remove(wb.active)
    add_styles(wb)
    legend_sheet(wb, entries)
    s2_spec = [(col, col, w, 'wrap') for col, w in zip(s2.columns, S2_WIDTHS)]
    expected = {}
    for sheet, (label, title, _), df, spec in zip(['S1', 'S2', 'S3', 'S4'], entries, [s1, s2, s3, s4], [S1, s2_spec, S3, S4]):
        expected[sheet] = ([h for h, _, _, _ in spec], table_sheet(wb, sheet, f'{label}: {title}', df, spec))
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / 'SupplementaryTables_S1-S4.xlsx'
    wb.save(path)

    # ---- read back every cell and compare ----
    check = load_workbook(path, read_only=True)
    for sheet, (headers, rows) in expected.items():
        got = list(check[sheet].iter_rows(min_row=3, values_only=True))
        assert list(got[0]) == headers, sheet
        assert [tuple(r) for r in got[1:]] == rows, sheet
    check.close()

    md = ['# Supplementary Tables 1–4 and Supplementary Data 1', '']
    for label, title, text in entries:
        md += [f'**{label}. {title}.** {text}', '']
    (output_dir / 'legend.md').write_text('\n'.join(md + [ABBREVIATIONS, '']))

    cells = sum(len(r) * len(h) for h, r in expected.values())
    print(f'NOTE: wrote {path} ({path.stat().st_size:,} bytes); {cells:,} cells across '
          f'{len(expected)} table sheets plus Legend read back identical', flush=True)


if __name__ == '__main__':
    main()
