#!/usr/bin/env python3
'''
build_figureS5_tables.py

Builds the source-data tables consumed by plot_supplementary_figure5.py:
the credible-set matching sensitivity curve (panel a), the exact VEP
annotation composition (panel b, reused unchanged from Figure 5's own
annotation table -- the display collapse to figure_style.VEP_GROUPS happens
at draw time in plot_supplementary_figure5.py, not here), the
gene-nomination evidence summary (panel c), and the GO/KEGG pathway
enrichment term selection and dot-plot cells (panel d).

Panel d selects, from every (domain, library, term) enrichment test computed
by 12_Gene_Nomination/compute_domain_enrichment.py, a non-redundant set of
up to MAX_TERMS significant terms with at least one representative per
phecode domain (round-robin by rank, skipping any term whose gene-set
overlap with an already-selected term exceeds REDUNDANCY_OC), then expands
that selection back across every domain to build the dot-plot matrix.

Input:
  CS_MATCH_STATUS            same credible-set collation table as
                              build_figure5_tables.py (pop, tool, trait,
                              locus, cs, best_pip_weighted_jaccard, shared)
  ANNOTATION_COMPOSITION     annotation_composition.tsv, as written by
                              build_figure5_tables.py (group, category,
                              count, total, percent)
  GENE_EVIDENCE_COMPOSITION  gene_evidence_composition.tsv, as written by
                              build_figure5_tables.py (evidence, genes)
  DOMAIN_ENRICHMENT_RESULTS  domain_enrichment_results.tsv, as written by
                              12_Gene_Nomination/compute_domain_enrichment.py
  DOMAIN_ENRICHMENT_TERM_GENES  domain_enrichment_term_genes.tsv, as written
                              by the same script (each term's gene set,
                              already intersected with the nominated genes)

Output (all tab-separated, written under OUT_DIR):
  a_matching_sensitivity.tsv  pop, tool, threshold, new_cs, new_matched
  b_annotation_exact.tsv      group, category, count, total, percent (copied through)
  c_gene_evidence.tsv         evidence, label, genes, percent
  d_enrichment_selected_terms.tsv  the selected non-redundant terms, with
                              selection_rank / selection_domain
  d_enrichment_all_cells.tsv  every (domain, selected term) cell -- the
                              selected terms' significance expanded across
                              every tested domain
  d_enrichment_visible_cells.tsv  d_enrichment_all_cells.tsv restricted to
                              significant (q_global < 0.05) cells -- what
                              the dot plot actually draws
  d_enrichment_display_order.tsv  the selected terms in row-display order
                              (selection_domain, q_global, term), with a
                              display_row index

Usage:
  python build_figureS5_tables.py [--cs-match-status PATH]
      [--annotation-composition PATH] [--gene-evidence-composition PATH]
      [--domain-enrichment-results PATH] [--domain-enrichment-term-genes PATH]
      [--out-dir DIR]

***Run in finemap conda environment***
'''
import getopt
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import figure_style as S

###############################################################################
##########  DEFAULT PATHS  ####################################################
###############################################################################

ROOT = Path("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping")
APPLICATION_RESULTS = ROOT / "application_results"
CS_MATCH_STATUS = APPLICATION_RESULTS / "credible_sets" / "cs_match_status_native.tsv.gz"
FIG5_SOURCE_DATA = APPLICATION_RESULTS / "figures" / "source_data" / "fig5"
ANNOTATION_COMPOSITION = FIG5_SOURCE_DATA / "annotation_composition.tsv"
GENE_EVIDENCE_COMPOSITION = FIG5_SOURCE_DATA / "gene_evidence_composition.tsv"
GENE_NOMINATION = APPLICATION_RESULTS / "gene_nomination"
DOMAIN_ENRICHMENT_RESULTS = GENE_NOMINATION / "domain_enrichment_results.tsv"
DOMAIN_ENRICHMENT_TERM_GENES = GENE_NOMINATION / "domain_enrichment_term_genes.tsv"
OUT_DIR = APPLICATION_RESULTS / "figures" / "source_data" / "figs5"

EVIDENCE_LABELS = {"coding": "Coding only", "abc": "ABC only", "both": "Both"}

# Term-selection constants (GO/KEGG panel d), ported unchanged.
MAX_TERMS = 12
REDUNDANCY_OC = 0.7
LIBRARY_LABELS = {"KEGG_2021_Human": "KEGG", "GO_Biological_Process_2023": "GO BP"}
# Raw phecode "group" -> short display name for a handful of long/awkward labels;
# every other raw group falls back to str.capitalize().
DOMAIN_DISPLAY_NAMES = {
    "circulatory system": "Circulatory", "hematopoietic": "Blood",
    "infectious diseases": "Infectious", "injuries & poisonings": "Injury/poisoning",
    "mental disorders": "Mental", "endocrine/metabolic": "Endocrine/metabolic",
}

###############################################################################
##########  PANEL A: matching sensitivity across PIP-weighted-Jaccard thresholds  ####
###############################################################################

def compute_matching_sensitivity(status):
    rows = []
    for pop, tool in S.PAIRS:
        d = status[status["pop"].eq(pop) & status.tool.eq(tool)]
        new_cs = len(d)
        for threshold in S.MATCH_THRESHOLDS:
            new_matched = int((d.best_pip_weighted_jaccard >= threshold).sum())
            rows.append(dict(pop=pop, tool=tool, threshold=threshold, new_cs=new_cs, new_matched=new_matched))
    return pd.DataFrame(rows)


###############################################################################
##########  PANEL C: gene-nomination evidence summary  ########################
###############################################################################

def compute_gene_evidence(evidence_composition):
    d = evidence_composition.set_index("evidence").loc[["coding", "abc", "both"]].reset_index()
    d["percent"] = 100 * d.genes / d.genes.sum()
    d["label"] = d.evidence.map(EVIDENCE_LABELS)
    return d[["evidence", "label", "genes", "percent"]]


###############################################################################
##########  PANEL D: GO/KEGG term selection and dot-plot matrix  ##############
###############################################################################

def overlap_coefficient(a, b):
    return len(a & b) / min(len(a), len(b)) if a and b else 0.0


def select(candidates, members):
    """Round-robin by phecode domain: each pass offers every domain its best
    remaining candidate (by q_global, p, domain, library, term order), skips
    any candidate too similar (overlap coefficient >= REDUNDANCY_OC) to an
    already-selected term's gene set, and stops at MAX_TERMS or an empty pass."""
    candidates = candidates.sort_values(["q_global", "p", "domain", "library", "term"])
    domains = candidates.drop_duplicates("domain").domain.tolist()
    queues = {domain: list(frame.to_dict("records")) for domain, frame in candidates.groupby("domain")}
    selected, excluded = [], []
    while len(selected) < MAX_TERMS:
        added = False
        for domain in domains:
            while queues[domain]:
                row = queues[domain].pop(0)
                key = (row["library"], row["term"])
                redundant = [prior for prior in selected if overlap_coefficient(
                    members[key], members[(prior["library"], prior["term"])]) >= REDUNDANCY_OC]
                if redundant:
                    prior = redundant[0]
                    excluded.append({"domain": domain, "library": key[0], "term": key[1],
                                      "represented_by_library": prior["library"],
                                      "represented_by_term": prior["term"],
                                      "overlap_coefficient": overlap_coefficient(
                                          members[key], members[(prior["library"], prior["term"])])})
                    continue
                row["selection_rank"] = len(selected) + 1
                row["selection_domain"] = domain
                selected.append(row)
                added = True
                break
            if len(selected) == MAX_TERMS:
                break
        if not added:
            break
    return pd.DataFrame(selected), pd.DataFrame(excluded), domains


def load_term_members(term_genes):
    """(library, term) -> gene set, from domain_enrichment_term_genes.tsv's
    semicolon-joined `genes` column (already intersected with nominated genes)."""
    members = {}
    for row in term_genes.itertuples():
        genes = set(row.genes.split(";")) if isinstance(row.genes, str) and row.genes else set()
        members[(row.library, row.term)] = genes
    return members


def compute_enrichment_selection(results, term_genes):
    """Select non-redundant significant terms (one per domain per round) and
    expand them into the dot-plot's full cell matrix, visible-cell subset,
    and display row order. Returns (selected, matrix, display_order)."""
    primary = results[results.background.eq("nominated_genes") & ~results.domain.eq("unclassified")].copy()
    members = load_term_members(term_genes)
    candidates = primary[primary.q_global < 0.05]
    selected, excluded, domains = select(candidates, members)
    if selected.empty:
        selected = pd.DataFrame(columns=list(primary.columns) + ["selection_rank", "selection_domain"])
    selected["term_label"] = selected.term.map(lambda value: re.sub(r"\s*\(GO:\d+\)$", "", value))
    selected["library_label"] = selected.library.map(LIBRARY_LABELS)
    keys = selected[["library", "term", "term_label", "library_label", "selection_rank", "selection_domain"]]
    matrix = primary.merge(keys, on=["library", "term"], validate="many_to_one")
    matrix["domain_label"] = matrix.domain.map(lambda value: DOMAIN_DISPLAY_NAMES.get(value, value.capitalize()))
    matrix["significant"], matrix["display"] = matrix.q_global < 0.05, matrix.q_global < 0.05
    matrix = matrix.sort_values(["domain", "selection_rank"])
    display_order = selected.sort_values(["selection_domain", "q_global", "term"]).reset_index(drop=True)
    display_order["display_row"] = np.arange(len(display_order))
    print(f"NOTE: GO/KEGG enrichment: {len(selected)} selected terms "
          f"({selected.library.value_counts().to_dict() if len(selected) else {}}), "
          f"{int(matrix.display.sum())} visible cells across {matrix.loc[matrix.display, 'domain'].nunique()} domains")
    return selected, matrix, display_order


###############################################################################
##########  CLI + MAIN  ########################################################
###############################################################################

def usage():
    print(__doc__)


def main(argv):
    cs_match_status = CS_MATCH_STATUS
    annotation_composition, gene_evidence_composition = ANNOTATION_COMPOSITION, GENE_EVIDENCE_COMPOSITION
    domain_enrichment_results, domain_enrichment_term_genes = DOMAIN_ENRICHMENT_RESULTS, DOMAIN_ENRICHMENT_TERM_GENES
    out_dir = OUT_DIR
    try:
        opts, _ = getopt.getopt(argv, "h", ["help", "cs-match-status=", "annotation-composition=",
                                             "gene-evidence-composition=", "domain-enrichment-results=",
                                             "domain-enrichment-term-genes=", "out-dir="])
    except getopt.GetoptError as e:
        print(f"ERROR: {e}"); usage(); sys.exit(2)
    for opt, val in opts:
        if opt in ("-h", "--help"):
            usage(); sys.exit(0)
        elif opt == "--cs-match-status":
            cs_match_status = Path(val)
        elif opt == "--annotation-composition":
            annotation_composition = Path(val)
        elif opt == "--gene-evidence-composition":
            gene_evidence_composition = Path(val)
        elif opt == "--domain-enrichment-results":
            domain_enrichment_results = Path(val)
        elif opt == "--domain-enrichment-term-genes":
            domain_enrichment_term_genes = Path(val)
        elif opt == "--out-dir":
            out_dir = Path(val)

    for path in [cs_match_status, annotation_composition, gene_evidence_composition]:
        if not path.exists():
            print(f"ERROR: input not found: {path}"); sys.exit(1)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"NOTE: reading credible-set matching status from {cs_match_status}")
    status = pd.read_csv(cs_match_status, sep="\t")
    compute_matching_sensitivity(status).to_csv(out_dir / "a_matching_sensitivity.tsv", sep="\t", index=False)

    print(f"NOTE: reusing annotation composition from {annotation_composition}")
    shutil.copy(annotation_composition, out_dir / "b_annotation_exact.tsv")

    print(f"NOTE: reusing gene evidence composition from {gene_evidence_composition}")
    evidence_composition = pd.read_csv(gene_evidence_composition, sep="\t")
    compute_gene_evidence(evidence_composition).to_csv(out_dir / "c_gene_evidence.tsv", sep="\t", index=False)

    if domain_enrichment_results.exists() and domain_enrichment_term_genes.exists():
        print(f"NOTE: reading domain enrichment results from {domain_enrichment_results}")
        results = pd.read_csv(domain_enrichment_results, sep="\t")
        term_genes = pd.read_csv(domain_enrichment_term_genes, sep="\t", keep_default_na=False)
        selected, matrix, display_order = compute_enrichment_selection(results, term_genes)
        selected.to_csv(out_dir / "d_enrichment_selected_terms.tsv", sep="\t", index=False)
        matrix.to_csv(out_dir / "d_enrichment_all_cells.tsv", sep="\t", index=False)
        matrix[matrix.display].to_csv(out_dir / "d_enrichment_visible_cells.tsv", sep="\t", index=False)
        display_order.to_csv(out_dir / "d_enrichment_display_order.tsv", sep="\t", index=False)
    else:
        print(f"NOTE: {domain_enrichment_results} or {domain_enrichment_term_genes} not found; "
              f"panel d (GO/KEGG enrichment) will be skipped by plot_supplementary_figure5.py")

    print(f"NOTE: wrote Supplementary Figure 5 source-data tables to {out_dir}")


if __name__ == "__main__":
    main(sys.argv[1:])
