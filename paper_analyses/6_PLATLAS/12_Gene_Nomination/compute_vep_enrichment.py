#!/usr/bin/env python3
'''
compute_vep_enrichment.py

VEP consequence composition by fine-mapping resolution, and its Fisher-exact
coding enrichment (step 12 enrichment tables; feed the enrichment panel of
the application figures).

A variant's resolution group is the size of the smallest credible set it
belongs to, across every tool/pop/trait membership (All fine-mapped ->
11+ -> 6-10 -> 2-5 -> 1 variant). "Coding" = Splice/Start/Stop Gain-Loss +
Missense + Synonymous. Fisher's exact test compares the coding fraction
between adjacent resolution groups, plus the extreme All-vs-1 comparison.

Input:
  --vep      VEP --most_severe tab output (run_vep_annotation.sh), columns
             Uploaded_variation, Consequence
  --members  gwas_cs_members.tsv (step 11), columns variant_id, cs_size

Output:
  vep_enrichment.tsv     category x resolution-group percentage table
  vep_coding_fisher.tsv  columns: group1, group2, coding_pct1, coding_pct2, odds_ratio, p
                         (written next to --out-tab)

Usage:
  python compute_vep_enrichment.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import csv
import getopt
from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact

VEP_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination/vep_all_most_severe.tsv"
MEMBERS_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/credible_sets/gwas_cs_members.tsv"
OUT_TAB_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination/vep_enrichment.tsv"

# "coding" = Splice/Start/Stop Gain-Loss + Missense + Synonymous (Fisher coding enrichment)
CODING_CATS = ("Splice/Start/Stop Gain/Loss", "Missense", "Synonymous")

# Stacking order (bottom -> top)
CAT_NAMES = [
    "Regulatory Region", "TF Binding Site", "Non-Gene Unknown", "Intron",
    "Non-Coding Exon", "3/5 prime UTR", "Mature miRNA", "Synonymous",
    "Missense", "Splice/Start/Stop Gain/Loss",
]

# VEP SO term -> consequence category
SO2CAT = {}
for t in ("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost",
          "start_lost", "frameshift_variant", "transcript_ablation", "transcript_amplification"):
    SO2CAT[t] = "Splice/Start/Stop Gain/Loss"
for t in ("missense_variant", "protein_altering_variant", "inframe_insertion", "inframe_deletion"):
    SO2CAT[t] = "Missense"
for t in ("synonymous_variant", "stop_retained_variant", "start_retained_variant",
          "incomplete_terminal_codon_variant"):
    SO2CAT[t] = "Synonymous"
SO2CAT["mature_miRNA_variant"] = "Mature miRNA"
for t in ("5_prime_UTR_variant", "3_prime_UTR_variant"):
    SO2CAT[t] = "3/5 prime UTR"
for t in ("non_coding_transcript_exon_variant", "non_coding_transcript_variant", "NMD_transcript_variant"):
    SO2CAT[t] = "Non-Coding Exon"
for t in ("intron_variant", "splice_region_variant", "splice_donor_region_variant",
          "splice_polypyrimidine_tract_variant", "splice_donor_5th_base_variant"):
    SO2CAT[t] = "Intron"
for t in ("intergenic_variant", "upstream_gene_variant", "downstream_gene_variant"):
    SO2CAT[t] = "Non-Gene Unknown"
for t in ("TF_binding_site_variant", "TFBS_ablation", "TFBS_amplification"):
    SO2CAT[t] = "TF Binding Site"
for t in ("regulatory_region_variant", "regulatory_region_ablation", "regulatory_region_amplification"):
    SO2CAT[t] = "Regulatory Region"

GROUPS = ["All", "11+", "6-10", "2-5", "1"]


def help(exit_num=0):
    print("""Usage: python compute_vep_enrichment.py [OPTIONS]
Optional:
  --vep       VEP --most_severe tab output (default: {vep})
  --members   Credible-set members table, step 11 output (default: {m})
  --out-tab   Enrichment table output path; vep_coding_fisher.tsv is written
              next to it (default: {o})
  -h, --help  Show this help message""".format(vep=VEP_DEFAULT, m=MEMBERS_DEFAULT, o=OUT_TAB_DEFAULT))
    sys.exit(exit_num)


###############################################################################
##########  LOADING  ###########################################################
###############################################################################

def size_bin(sz):
    if sz == 1:
        return "1"
    if sz <= 5:
        return "2-5"
    if sz <= 10:
        return "6-10"
    return "11+"


def load_consequence(path):
    """VEP --most_severe output -> {variant_id: category}."""
    vc, idx = {}, None
    for line in open(path):
        if line.startswith("#"):
            if line.startswith("#Uploaded"):
                hdr = line[1:].rstrip("\n").split("\t")
                idx = {c: i for i, c in enumerate(hdr)}
            continue
        if idx is None:
            continue
        f = line.rstrip("\n").split("\t")
        csq = f[idx["Consequence"]].split(",")[0]
        vc[f[idx["Uploaded_variation"]]] = SO2CAT.get(csq, "Non-Gene Unknown")
    return vc


def load_min_cs_size(path):
    """Smallest credible-set size each variant belongs to, across all tool/pop/trait memberships."""
    mn = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            v, sz = row["variant_id"], int(row["cs_size"])
            if v not in mn or sz < mn[v]:
                mn[v] = sz
    return mn


###############################################################################
##########  FISHER ENRICHMENT  #################################################
###############################################################################

def coding_nc(counts, tot, g):
    cod = sum(counts[g][c] for c in CODING_CATS)
    return cod, tot[g] - cod


def fisher_coding(counts, tot, g1, g2):
    a, b = coding_nc(counts, tot, g1)
    c, d = coding_nc(counts, tot, g2)
    odds, p = fisher_exact([[a, b], [c, d]])   # OR > 1: g1 is more coding-enriched than g2
    return p, odds


###############################################################################
##########  MAIN  ##############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "h", ["help", "vep=", "members=", "out-tab="])
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    vep_path = VEP_DEFAULT
    members_path = MEMBERS_DEFAULT
    out_tab = OUT_TAB_DEFAULT

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--vep":
            vep_path = a
        elif o == "--members":
            members_path = a
        elif o == "--out-tab":
            out_tab = a
        else:
            assert False, f"Unhandled option: {o}"

    vc = load_consequence(vep_path)
    mn = load_min_cs_size(members_path)

    # ── Group counts: unit = unique variant; background = all annotated variants,
    #    resolution bucket = each variant's smallest CS size (buckets are mutually exclusive) ──
    counts = {g: defaultdict(int) for g in GROUPS}
    for v, cat in vc.items():
        counts["All"][cat] += 1
        sz = mn.get(v)
        if sz is not None:
            counts[size_bin(sz)][cat] += 1
    tot = {g: sum(counts[g].values()) for g in GROUPS}

    Path(out_tab).parent.mkdir(parents=True, exist_ok=True)
    with open(out_tab, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["category"] + [f"{g}(n={tot[g]})" for g in GROUPS])
        for cat in CAT_NAMES:
            w.writerow([cat] + [f"{100 * counts[g][cat] / tot[g]:.2f}%" if tot[g] else "0" for g in GROUPS])

    # ── Fisher exact: coding enrichment between adjacent groups + the extreme All-vs-1 ──
    pairs = [("All", "11+"), ("11+", "6-10"), ("6-10", "2-5"), ("2-5", "1"), ("All", "1")]
    fres = {p: fisher_coding(counts, tot, *p) for p in pairs}

    fisher_tab = str(Path(out_tab).with_name("vep_coding_fisher.tsv"))
    with open(fisher_tab, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["group1", "group2", "coding_pct1", "coding_pct2", "odds_ratio", "p"])
        for g1, g2 in pairs:
            c1, _ = coding_nc(counts, tot, g1)
            c2, _ = coding_nc(counts, tot, g2)
            p, odds = fres[(g1, g2)]
            w.writerow([g1, g2, f"{100 * c1 / tot[g1]:.4f}", f"{100 * c2 / tot[g2]:.4f}",
                        f"{odds:.4f}", f"{p:.4e}"])

    print("NOTE: group sizes:", {g: tot[g] for g in GROUPS})
    for g in GROUPS:
        cod, _ = coding_nc(counts, tot, g)
        print(f"NOTE:   {g:5s} coding={100 * cod / tot[g]:.2f}%  (n={tot[g]})")
    print("NOTE: Fisher (coding vs non-coding): " +
          "; ".join(f"{a}~{b}: p={v[0]:.1e} OR={v[1]:.2f}" for (a, b), v in fres.items()))
    print(f"NOTE: -> {out_tab} + {fisher_tab}")


if __name__ == "__main__":
    main()
