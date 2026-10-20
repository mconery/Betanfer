#!/usr/bin/env python3
'''
compute_domain_enrichment.py

GO/KEGG pathway enrichment of nominated genes within each phecode domain
(step 12, gene nomination). A nominated gene's domain is the phecode
"group" of the trait it was nominated for; enrichment is tested per domain,
per GO-BP/KEGG term, against two backgrounds -- GENCODE protein-coding genes,
and the full set of nominated genes (a sensitivity check for background
choice). Benjamini-Hochberg FDR is applied twice: once across every test run
for a background (q_global), and once within each (domain, library) family
(q_within_domain_library).

Figure-specific term selection (choosing a non-redundant display set of
terms for the Supplementary Figure 5 enrichment panel) is not done here --
see build_figureS5_tables.py in 14_Application_Figures, which reads this
script's domain_enrichment_results.tsv and domain_enrichment_term_genes.tsv.

Input:
  --nominations  gene_nominations_comprehensive.tsv (nominate_genes_comprehensive.py),
                 columns trait, gene (one row per trait x gene nomination)
  --gtf          GENCODE annotation GTF (used only for its protein_coding gene set)
  --phecodes     phecode metadata table, columns phecode, group, description
                 (default: phecode_catalog.tsv in 13_Supplementary_Tables)
  --gmt-cache    Directory with <library>.gmt gene-set files for LIBRARIES;
                 any missing library is downloaded from Enrichr into it

Output (written to --out-dir):
  base_trait_domain_map.tsv        base_trait, domain, member_domains -- the
                                    phecode domain(s) each base trait maps to
  domain_enrichment_results.tsv    background, domain, library, term,
                                    n_universe, n_query, n_term, n_overlap,
                                    expected, fold, p, overlap_genes, q_global,
                                    q_within_domain_library -- one row per
                                    (background, domain, library, term) test
  domain_enrichment_term_genes.tsv library, term, n_genes, genes -- each
                                    term's gene set intersected with the
                                    nominated-genes background (redundancy
                                    input for downstream term selection)
  domain_sizes.tsv                 domain, n_traits, n_genes_raw, n_genes, tested

Usage:
  python compute_domain_enrichment.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import csv
import gzip
import getopt
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

HERE = Path(__file__).resolve().parent
NOMINATIONS_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination/gene_nominations_comprehensive.tsv"
OUT_DIR_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination"
# GENCODE is an external resource only available in our tree, same convention as
# nominate_genes_comprehensive.py's --gtf.
GTF_DEFAULT = "/lus/flare/projects/GeomicVar/jianhuawang/betainfer/data/raw/annotation/gencode.v19.annotation.gtf.gz"
PHECODES_DEFAULT = str(HERE.parent / "13_Supplementary_Tables" / "phecode_catalog.tsv")
GMT_CACHE_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination/gmt"
ENRICHR_URL = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={}"

LIBRARIES = ("KEGG_2021_Human", "GO_Biological_Process_2023")
MIN_QUERY = 20              # minimum query genes (in-background) for a domain to be tested
MIN_TERM_SIZE, MAX_TERM_SIZE = 10, 500   # term size range (in-background) kept from the GMT

# Raw phecode "group" category -> short domain slug. Unlisted categories (e.g.
# "endocrine/metabolic", "dermatologic") pass through unchanged as their own slug.
DOMAIN_NAMES = {
    "circulatory system": "circulatory", "hematopoietic": "blood",
    "infectious diseases": "infectious", "injuries & poisonings": "injury",
    "mental disorders": "mental", "pregnancy complications": "pregnancy",
    "congenital anomalies": "congenital", "unclassified": "other",
}


def help(exit_num=0):
    print("""Usage: python compute_domain_enrichment.py [OPTIONS]
Optional:
  --nominations  Comprehensive gene-nomination table (default: {n})
  --gtf          GENCODE annotation GTF (default: {g})
  --phecodes     Phecode metadata table (default: {p})
  --gmt-cache    Directory with KEGG/GO-BP .gmt files; missing ones are
                 downloaded from Enrichr (default: {c})
  --out-dir      Output directory (default: {o})
  -h, --help     Show this help message""".format(
        n=NOMINATIONS_DEFAULT, g=GTF_DEFAULT, p=PHECODES_DEFAULT,
        c=GMT_CACHE_DEFAULT, o=OUT_DIR_DEFAULT))
    sys.exit(exit_num)


###############################################################################
##########  REFERENCE LOADING  ################################################
###############################################################################

def protein_coding_genes(gtf_path):
    """Return the set of gene symbols GENCODE marks gene_type protein_coding."""
    keep = set()
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if "\tgene\t" in line and 'gene_type "protein_coding"' in line:
                keep.add(line.split('gene_name "')[1].split('"')[0])
    return keep


def base_trait(trait):
    """Roll a phecode subtype ("Phe_250_2") up to its parent ("Phe_250")."""
    parts = trait.split("_")
    return "_".join(parts[:2]) if len(parts) == 3 and parts[0] == "Phe" and parts[2].isdigit() else trait


def phecode_group(trait, metadata):
    """Raw phecode "group" category for a trait, or "unclassified" if untracked."""
    code = trait.removeprefix("Phe_").replace("_", ".")
    return metadata.loc[code, "group"] if code in metadata.index else "unclassified"


###############################################################################
##########  DOMAIN QUERIES AND BASE-TRAIT MAP  #################################
###############################################################################

def load_queries(nominations_path, metadata, background):
    """Per-domain nominated gene sets. Returns (genes_by_domain, nominated genes
    restricted to `background`, a per-domain size/eligibility table)."""
    genes, traits, all_genes = defaultdict(set), defaultdict(set), set()
    with open(nominations_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene, trait = row["gene"], row["trait"]
            all_genes.add(gene)
            domain = phecode_group(trait, metadata)
            genes[domain].add(gene)
            traits[domain].add(trait)
    sizes = []
    for domain in sorted(genes):
        n_genes = len(genes[domain] & background)
        sizes.append(dict(domain=domain, n_traits=len(traits[domain]),
                           n_genes_raw=len(genes[domain]), n_genes=n_genes, tested=n_genes >= MIN_QUERY))
    return genes, all_genes & background, pd.DataFrame(sizes)


def domain_map(nominations_path, metadata):
    """base_trait -> the phecode domain(s) its (sub)traits carry, short-slugged.
    A base trait with exactly one known domain across its subtypes takes that
    domain; zero or more than one known domain both map to "other"."""
    subtraits = defaultdict(set)
    with open(nominations_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            subtraits[base_trait(row["trait"])].add(row["trait"])

    rows = []
    for base in sorted(subtraits):
        known = set()
        for trait in subtraits[base]:
            raw = phecode_group(trait, metadata)
            if raw != "unclassified":
                known.add(DOMAIN_NAMES.get(raw, raw))
        domain = next(iter(known)) if len(known) == 1 else "other"
        rows.append(dict(base_trait=base, domain=domain, member_domains=";".join(sorted(known))))
    return rows


###############################################################################
##########  GO/KEGG TERMS AND HYPERGEOMETRIC TEST  #############################
###############################################################################

def fetch_gmt(gmt_cache):
    """Download any LIBRARIES gene-set file missing from `gmt_cache` from Enrichr."""
    Path(gmt_cache).mkdir(parents=True, exist_ok=True)
    for library in LIBRARIES:
        path = Path(gmt_cache) / f"{library}.gmt"
        if not path.exists():
            print(f"NOTE: downloading {library} from Enrichr -> {path}")
            with urllib.request.urlopen(ENRICHR_URL.format(library), timeout=300) as response:
                path.write_bytes(response.read())


def load_terms(gmt_cache, background):
    """KEGG/GO-BP gene sets, intersected with `background` and kept only if
    their in-background size falls in [MIN_TERM_SIZE, MAX_TERM_SIZE]."""
    terms = []
    for library in LIBRARIES:
        with (Path(gmt_cache) / f"{library}.gmt").open() as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                members = {g.split(",")[0].strip() for g in fields[2:] if g.strip()} & background
                if MIN_TERM_SIZE <= len(members) <= MAX_TERM_SIZE:
                    terms.append((library, fields[0], members))
    return terms


def term_genes_table(terms, nominated):
    """Each term's gene set intersected with the nominated-genes background --
    the redundancy input for the downstream term-selection step."""
    rows = []
    for library, term, members in terms:
        genes = sorted(members & nominated)
        rows.append(dict(library=library, term=term, n_genes=len(genes), genes=";".join(genes)))
    return pd.DataFrame(rows)


def bh(p):
    """Benjamini-Hochberg q values for the complete supplied family, preserving row order."""
    p = np.asarray(p, dtype=float)
    order = np.argsort(p, kind="stable")
    ranked = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    q = np.empty(len(p))
    q[order] = np.minimum(ranked, 1)
    return q


def calculate(queries, backgrounds, terms):
    """Hypergeometric enrichment of each (domain, term) query against each
    named background; q_global is BH-corrected within a background across
    every domain/library/term test, q_within_domain_library within each
    (background, domain, library) family."""
    rows = []
    for background_name, background in backgrounds.items():
        for domain, genes in sorted(queries.items()):
            query = genes & background
            if len(query) < MIN_QUERY:
                continue
            for library, term, members in terms:
                pathway = members & background
                overlap = query & pathway
                expected = len(query) * len(pathway) / len(background)
                rows.append(dict(
                    background=background_name, domain=domain, library=library, term=term,
                    n_universe=len(background), n_query=len(query), n_term=len(pathway),
                    n_overlap=len(overlap), expected=expected,
                    fold=len(overlap) / expected if expected else np.nan,
                    p=float(hypergeom.sf(len(overlap) - 1, len(background), len(pathway), len(query))),
                    overlap_genes=";".join(sorted(overlap))))
    results = pd.DataFrame(rows)
    for _, idx in results.groupby("background").groups.items():
        results.loc[idx, "q_global"] = bh(results.loc[idx, "p"])
    for _, idx in results.groupby(["background", "domain", "library"]).groups.items():
        results.loc[idx, "q_within_domain_library"] = bh(results.loc[idx, "p"])
    return results


###############################################################################
##########  MAIN  ##############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "nominations=", "gtf=", "phecodes=", "gmt-cache=", "out-dir="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    nominations_path = NOMINATIONS_DEFAULT
    gtf_path = GTF_DEFAULT
    phecodes_path = PHECODES_DEFAULT
    gmt_cache = GMT_CACHE_DEFAULT
    out_dir = Path(OUT_DIR_DEFAULT)

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--nominations":
            nominations_path = a
        elif o == "--gtf":
            gtf_path = a
        elif o == "--phecodes":
            phecodes_path = a
        elif o == "--gmt-cache":
            gmt_cache = a
        elif o == "--out-dir":
            out_dir = Path(a)
        else:
            assert False, f"Unhandled option: {o}"

    for path in [nominations_path, gtf_path, phecodes_path]:
        if not Path(path).exists():
            print(f"ERROR: input not found: {path}")
            sys.exit(1)
    fetch_gmt(gmt_cache)
    out_dir.mkdir(parents=True, exist_ok=True)

    pc_genes = protein_coding_genes(gtf_path)
    metadata = pd.read_csv(phecodes_path, sep="\t", dtype=str).set_index("phecode")
    print(f"NOTE: GENCODE protein-coding genes: {len(pc_genes)}")

    queries, nominated, sizes = load_queries(nominations_path, metadata, pc_genes)
    print(f"NOTE: nominated genes (any domain): {len(nominated)}; "
          f"domains tested (n_genes>={MIN_QUERY}): {int(sizes.tested.sum())}/{len(sizes)}")
    sizes.to_csv(out_dir / "domain_sizes.tsv", sep="\t", index=False)

    terms = load_terms(gmt_cache, pc_genes)
    print(f"NOTE: KEGG/GO-BP terms kept (size {MIN_TERM_SIZE}-{MAX_TERM_SIZE} in background): {len(terms)}")

    results = calculate(queries, {"gencode_protein_coding": pc_genes, "nominated_genes": nominated}, terms)
    results.to_csv(out_dir / "domain_enrichment_results.tsv", sep="\t", index=False)
    print(f"NOTE: enrichment tests: {len(results)} -> {out_dir / 'domain_enrichment_results.tsv'}")

    term_genes_table(terms, nominated).to_csv(out_dir / "domain_enrichment_term_genes.tsv", sep="\t", index=False)

    rows = domain_map(nominations_path, metadata)
    pd.DataFrame(rows).to_csv(out_dir / "base_trait_domain_map.tsv", sep="\t", index=False)
    print(f"NOTE: base traits mapped to a domain: {len(rows)} -> {out_dir / 'base_trait_domain_map.tsv'}")


if __name__ == "__main__":
    main()
