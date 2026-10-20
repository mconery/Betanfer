#!/usr/bin/env python3
'''
build_figure5_tables.py

Builds the source-data tables consumed by plot_figure5.py: credible-set (CS)
matching against the ABF baseline, fine-mapping resolution by CS size, VEP
annotation and coding-variant composition by CS size, per-gene nomination
ranking / domain matrix / evidence breakdown, the pleiotropy distribution,
and the new "genes reached only through CS with no ABF counterpart" panel.

Every table is derived straight from the credible-set collation and gene
nomination outputs; no exploratory analysis logic is included.

Input:
  CS_MATCH_STATUS   credible-set collation: pop, tool, trait, locus, cs,
                    best_pip_weighted_jaccard, shared (shared = matched an
                    ABF credible set at the 0.10 PIP-weighted-Jaccard threshold)
  CS_MEMBERS        credible-set collation: trait, locus, pop, tool, cs,
                    cs_size, variant_id (one row per CS member variant, every
                    tool/pop; used both for per-tool CS resolution and for each
                    variant's minimum CS size across every tool/pop it is a
                    member of)
  GWAS_FINEMAPPED_VARIANTS  credible-set collation: variant_id (one row per
                    distinct finemapped variant across all tools/pops; the
                    "kept" set for the annotation-composition panel)
  GENE_NOMINATIONS  gene nomination: trait, pop, gene, variant_id,
                    implication_type, consequence, median_pip, abc_max_score,
                    abc_n_celltypes
  VEP_MOST_SEVERE   VEP --most_severe tab output (run_vep_annotation.sh in
                    12_Gene_Nomination): #Uploaded_variation, Consequence
                    ('##' header lines are skipped)
  DOMAIN_MAP        base_trait -> domain slug lookup, written by
                    compute_domain_enrichment.py (12_Gene_Nomination):
                    base_trait, domain, member_domains

Output (all tab-separated, written under OUT_DIR):
  figure5a_cs_matching.tsv     pop, tool, new_cs, new_matched, new_unmatched
  cs_resolution_four_bins.tsv  tool, size_bin, count, total, percent
  annotation_composition.tsv   group, category, count, total, percent -- one
                                row per distinct variant that is both VEP-
                                annotated and in the kept set, counted once in
                                "All" and once more under the size_bin of its
                                minimum CS size across every tool/pop CS it
                                belongs to (not restricted to one tool)
  coding_proportion.tsv        group, coding_variants, annotated_variants, percent
  gene_ranking.tsv             gene, phecode_groups, named_domains, evidence,
                                main_selected
  gene_domain_matrix.tsv       gene, domain, phecode_groups
  gene_evidence_composition.tsv  evidence, genes
  f_top_bar_counts.tsv         gene, evidence, phecode_groups -- top TOP_N_GENES
                                genes only (gene_ranking.tsv main_selected)
  pleiotropy_distribution.tsv  phecode_groups, genes
  d_genes_by_cs_class.tsv      view, Shared CS only, Shared and additional CS,
                                Additional CS only, total

Usage:
  python build_figure5_tables.py [--cs-match-status PATH] [--cs-members PATH]
      [--gwas-finemapped-variants PATH] [--gene-nominations PATH]
      [--vep-most-severe PATH] [--domain-map PATH] [--out-dir DIR]

***Run in finemap conda environment***
'''
import getopt
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
CS_MEMBERS = APPLICATION_RESULTS / "credible_sets" / "gwas_cs_members.tsv"
GWAS_FINEMAPPED_VARIANTS = APPLICATION_RESULTS / "credible_sets" / "gwas_finemapped_variants.tsv"
GENE_NOMINATIONS = APPLICATION_RESULTS / "gene_nomination" / "gene_nominations_comprehensive.tsv"
VEP_MOST_SEVERE = APPLICATION_RESULTS / "gene_nomination" / "vep_all_most_severe.tsv"
DOMAIN_MAP = APPLICATION_RESULTS / "gene_nomination" / "base_trait_domain_map.tsv"
OUT_DIR = APPLICATION_RESULTS / "figures" / "source_data" / "fig5"

EVIDENCE_MAP = {"nonsynonymous_coding_variant": "coding", "abc_interaction": "abc", "abc_promoter_overlap": "abc"}
CS_MEMBERS_CHUNKSIZE = 1_000_000
TOP_N_GENES = 12

###############################################################################
##########  CS MATCHING (panel a)  ############################################
###############################################################################

def compute_cs_matching(status):
    rows = []
    for pop, tool in S.PAIRS:
        d = status[status["pop"].eq(pop) & status.tool.eq(tool)]
        rows.append(dict(pop=pop, tool=tool, new_cs=len(d), new_matched=int(d.shared.sum()),
                          new_unmatched=int((~d.shared).sum())))
    return pd.DataFrame(rows)


###############################################################################
##########  CS MEMBERS: resolution (panel b) + annotation (panel c) 共用一次扫描  ####
###############################################################################

def read_cs_members_derived(path):
    """One pass over the (large) CS-member table: dedupe to one row per credible
    set for the resolution panel, and track each variant's minimum CS size across
    every tool/pop/trait/locus/cs it is a member of (not restricted to one tool)
    for the annotation-composition panel."""
    cs_cols = ["pop", "tool", "trait", "locus", "cs", "cs_size"]
    cs_rows = []
    min_cs_size = {}
    chunks = pd.read_csv(path, sep="\t", chunksize=CS_MEMBERS_CHUNKSIZE,
                          usecols=cs_cols + ["variant_id"])
    for chunk in chunks:
        cs_rows.append(chunk[cs_cols].drop_duplicates(subset=cs_cols[:-1]))
        for vid, size in chunk.groupby("variant_id").cs_size.min().items():
            if vid not in min_cs_size or size < min_cs_size[vid]:
                min_cs_size[vid] = size
    cs = pd.concat(cs_rows, ignore_index=True).drop_duplicates(subset=cs_cols[:-1])
    return cs, min_cs_size


def load_kept_variants(path):
    """Distinct finemapped variant ids across every tool/pop (the "kept" set
    that the annotation-composition panel is restricted to)."""
    return set(pd.read_csv(path, sep="\t", usecols=["variant_id"]).variant_id)


def compute_resolution(cs):
    cs = cs.assign(size_bin=cs.cs_size.map(S.resolution_bin))
    counts = cs.groupby(["tool", "size_bin"]).size().rename("count").reset_index()
    totals = cs.groupby("tool").size().rename("total")
    counts = counts.join(totals, on="tool")
    counts["percent"] = 100 * counts["count"] / counts.total
    order = {b: i for i, b in enumerate(S.RESOLUTION_BINS)}
    return counts.sort_values(["tool", "size_bin"], key=lambda s: s.map(order) if s.name == "size_bin" else s)


###############################################################################
##########  VEP ANNOTATION (panel c / S5 panel b)  ############################
###############################################################################

def read_vep_categories(path):
    with open(path) as fh:
        skip = 0
        for line in fh:
            if not line.startswith("##"):
                break
            skip += 1
    vep = pd.read_csv(path, sep="\t", skiprows=skip, usecols=["#Uploaded_variation", "Consequence"])
    vep.columns = ["variant_id", "consequence"]
    # Unrecognized consequence terms fall back to "Non-Gene Unknown" rather than
    # being dropped, matching the category map's own default.
    vep["category"] = vep.consequence.str.split(",").str[0].map(lambda c: S.SO2CAT.get(c, "Non-Gene Unknown"))
    return vep.drop_duplicates("variant_id").set_index("variant_id").category


def compute_annotation(kept, min_cs_size, vep_category):
    """One row per distinct variant that is both VEP-annotated and in the kept
    set, counted once in "All" and once more under the size_bin of its minimum
    CS size across every tool/pop CS it belongs to."""
    annotated_kept = vep_category.index.intersection(kept)
    categories = vep_category.loc[annotated_kept]
    bins = pd.Series({vid: S.size_bin(min_cs_size[vid]) for vid in annotated_kept if vid in min_cs_size})
    cats = [c for c, _ in S.MVP_CATS]
    rows = []
    for group, vids in [("All", annotated_kept)] + [(g, bins[bins.eq(g)].index) for g in ["11+", "6-10", "2-5", "1"]]:
        total = len(vids)
        counts = categories.loc[vids].value_counts().reindex(cats, fill_value=0)
        for cat, count in counts.items():
            rows.append(dict(group=group, category=cat, count=int(count), total=total,
                              percent=100 * count / total if total else 0.0))
    return pd.DataFrame(rows)


def compute_coding_proportion(annotation):
    rows = []
    for group, gd in annotation.groupby("group"):
        coding = int(gd[gd.category.isin(S.CODING_CATS)]["count"].sum())
        total = int(gd.total.iloc[0])
        rows.append(dict(group=group, coding_variants=coding, annotated_variants=total,
                          percent=100 * coding / total if total else 0.0))
    return pd.DataFrame(rows)


###############################################################################
##########  GENE NOMINATION: ranking, domain matrix, evidence (panels e, f)  ##
###############################################################################

def load_gene_links(gene_nominations_path, domain_map_path):
    nom = pd.read_csv(gene_nominations_path, sep="\t")
    nom["evidence"] = nom.implication_type.map(EVIDENCE_MAP)
    unknown = nom.evidence.isna().sum()
    if unknown:
        print(f"NOTE: {unknown} nomination rows have an unrecognized implication_type; dropped")
    nom = nom.dropna(subset=["evidence"])
    nom["group"] = nom.trait.map(S.base_trait)
    domain = pd.read_csv(domain_map_path, sep="\t", usecols=["base_trait", "domain"]).set_index("base_trait").domain
    nom["domain"] = nom.group.map(domain).fillna("unmapped")
    return nom


def gene_group_evidence(nom):
    """One row per (gene, phecode group) with the set of evidence kinds seen for it."""
    g = nom.groupby(["gene", "group"]).evidence.agg(lambda s: frozenset(s)).reset_index()
    g["kind"] = np.where(g.evidence.map(len).eq(2), "both", g.evidence.map(lambda s: next(iter(s))))
    domains = nom[["gene", "group", "domain"]].drop_duplicates(["gene", "group"])
    return g.merge(domains, on=["gene", "group"])


def compute_gene_tables(nom):
    links = gene_group_evidence(nom)
    per_gene = links.groupby("gene")
    phecode_groups = per_gene.size().rename("phecode_groups")
    named_domains = per_gene.domain.agg(lambda d: d[~d.isin(["other", "unmapped", "unclassified"])].nunique()).rename("named_domains")
    # A gene can have one group whose own kind is already "both" alongside other
    # groups that are purely "coding" or purely "abc", so the set of per-group
    # kinds can have 3 distinct members; any mix at all means "both" overall.
    evidence_union = per_gene.kind.agg(lambda s: frozenset(s))
    evidence = evidence_union.map(lambda s: "both" if len(s) > 1 else next(iter(s))).rename("evidence")

    genes = pd.concat([phecode_groups, named_domains, evidence], axis=1).reset_index()
    genes = genes.sort_values(["named_domains", "phecode_groups", "gene"], ascending=[False, False, True]).reset_index(drop=True)
    genes["main_selected"] = genes.index < TOP_N_GENES

    domain_matrix = links.groupby(["gene", "domain"]).size().rename("phecode_groups").reset_index()

    top_genes = genes.loc[genes.main_selected, "gene"]
    bar_counts = (links.groupby(["gene", "kind"]).size().rename("phecode_groups").reset_index()
                  .rename(columns={"kind": "evidence"}))
    bar_counts = (bar_counts.set_index(["gene", "evidence"])
                  .reindex(pd.MultiIndex.from_product([top_genes, ["coding", "abc", "both"]], names=["gene", "evidence"]), fill_value=0)
                  .reset_index())

    evidence_composition = evidence.value_counts().rename("genes").rename_axis("evidence").reset_index()
    pleiotropy = phecode_groups.value_counts().rename("genes").rename_axis("phecode_groups").reset_index().sort_values("phecode_groups")
    return genes, domain_matrix, bar_counts, evidence_composition, pleiotropy


###############################################################################
##########  NEW PANEL D: genes reached only through CS with no ABF match  #####
###############################################################################

def compute_class_counts(nom, cs_members_path, status):
    classes = ["Shared CS only", "Shared and additional CS", "Additional CS only"]
    rename = {"Shared CS only": "Shared CS only", "Shared and specific CS": "Shared and additional CS",
              "Specific CS only": "Additional CS only"}
    chunks = pd.read_csv(cs_members_path, sep="\t", chunksize=CS_MEMBERS_CHUNKSIZE,
                          usecols=["trait", "locus", "pop", "tool", "cs", "variant_id"])
    members = pd.concat(c[c.tool.eq("rsparsepro")].drop(columns="tool") for c in chunks)
    status_r = status[status.tool.eq("rsparsepro")][["pop", "trait", "locus", "cs", "shared"]]
    d = nom.merge(members, on=["trait", "pop", "variant_id"]).merge(status_r, on=["pop", "trait", "locus", "cs"])

    def classify(view):
        g = view.groupby("gene").shared.agg(any_shared="any", all_shared="all")
        raw = np.select([g.all_shared, g.any_shared], ["Shared CS only", "Shared and specific CS"], "Specific CS only")
        return pd.Series(raw, index=g.index).map(rename).value_counts().reindex(classes, fill_value=0)

    views = {"All": d, "Coding": d[d.evidence.eq("coding")], "ABC": d[d.evidence.eq("abc")]}
    out = pd.DataFrame({view: classify(v) for view, v in views.items()}).T.rename_axis("view")
    out["total"] = out.sum(axis=1)
    return out.loc[["All", "ABC", "Coding"]]


###############################################################################
##########  CLI + MAIN  ########################################################
###############################################################################

def usage():
    print(__doc__)


def main(argv):
    cs_match_status, cs_members = CS_MATCH_STATUS, CS_MEMBERS
    gwas_finemapped_variants = GWAS_FINEMAPPED_VARIANTS
    gene_nominations, vep_most_severe, domain_map = GENE_NOMINATIONS, VEP_MOST_SEVERE, DOMAIN_MAP
    out_dir = OUT_DIR
    try:
        opts, _ = getopt.getopt(argv, "h", ["help", "cs-match-status=", "cs-members=",
                                             "gwas-finemapped-variants=", "gene-nominations=",
                                             "vep-most-severe=", "domain-map=", "out-dir="])
    except getopt.GetoptError as e:
        print(f"ERROR: {e}"); usage(); sys.exit(2)
    for opt, val in opts:
        if opt in ("-h", "--help"):
            usage(); sys.exit(0)
        elif opt == "--cs-match-status":
            cs_match_status = Path(val)
        elif opt == "--cs-members":
            cs_members = Path(val)
        elif opt == "--gwas-finemapped-variants":
            gwas_finemapped_variants = Path(val)
        elif opt == "--gene-nominations":
            gene_nominations = Path(val)
        elif opt == "--vep-most-severe":
            vep_most_severe = Path(val)
        elif opt == "--domain-map":
            domain_map = Path(val)
        elif opt == "--out-dir":
            out_dir = Path(val)

    for path in [cs_match_status, cs_members, gwas_finemapped_variants, gene_nominations, vep_most_severe, domain_map]:
        if not path.exists():
            print(f"ERROR: input not found: {path}"); sys.exit(1)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"NOTE: reading credible-set matching status from {cs_match_status}")
    status = pd.read_csv(cs_match_status, sep="\t")
    compute_cs_matching(status).to_csv(out_dir / "figure5a_cs_matching.tsv", sep="\t", index=False)

    print(f"NOTE: reading credible-set members from {cs_members}")
    cs, min_cs_size = read_cs_members_derived(cs_members)
    compute_resolution(cs).to_csv(out_dir / "cs_resolution_four_bins.tsv", sep="\t", index=False)

    print(f"NOTE: reading kept variant set from {gwas_finemapped_variants}")
    kept = load_kept_variants(gwas_finemapped_variants)

    print(f"NOTE: reading VEP most-severe consequences from {vep_most_severe}")
    vep_category = read_vep_categories(vep_most_severe)
    annotation = compute_annotation(kept, min_cs_size, vep_category)
    annotation.to_csv(out_dir / "annotation_composition.tsv", sep="\t", index=False)
    compute_coding_proportion(annotation).to_csv(out_dir / "coding_proportion.tsv", sep="\t", index=False)

    print(f"NOTE: reading gene nominations from {gene_nominations}")
    nom = load_gene_links(gene_nominations, domain_map)
    genes, domain_matrix, bar_counts, evidence_composition, pleiotropy = compute_gene_tables(nom)
    genes[["gene", "phecode_groups", "named_domains", "evidence", "main_selected"]].to_csv(
        out_dir / "gene_ranking.tsv", sep="\t", index=False)
    domain_matrix.to_csv(out_dir / "gene_domain_matrix.tsv", sep="\t", index=False)
    bar_counts.to_csv(out_dir / "f_top_bar_counts.tsv", sep="\t", index=False)
    evidence_composition.to_csv(out_dir / "gene_evidence_composition.tsv", sep="\t", index=False)
    pleiotropy.to_csv(out_dir / "pleiotropy_distribution.tsv", sep="\t", index=False)

    print("NOTE: computing new panel d (genes reached only through unmatched credible sets)")
    class_counts = compute_class_counts(nom, cs_members, status)
    class_counts.to_csv(out_dir / "d_genes_by_cs_class.tsv", sep="\t")

    print(f"NOTE: wrote Figure 5 source-data tables to {out_dir}")


if __name__ == "__main__":
    main(sys.argv[1:])
