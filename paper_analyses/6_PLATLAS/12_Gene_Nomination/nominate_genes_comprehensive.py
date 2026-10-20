#!/usr/bin/env python3
'''
nominate_genes_comprehensive.py

Comprehensive gene nomination = coding (nominate_coding_genes.py) +
non-coding ABC enhancer-gene links (step 12).

Two evidence lines, matching MVP's gene-nomination approach:
  1. nonsynonymous_coding_variant -- from nominate_coding_genes.py (VEP missense/LoF).
  2. ABC enhancer-gene            -- a fine-mapped variant (median PIP >= 0.01) whose
     hg19 position falls inside a Nasser et al. ABC-Max enhancer is linked to that
     enhancer's TargetGene. implication_type is abc_promoter_overlap for a self-
     promoter interval, abc_interaction otherwise; keeps the max ABC score and the
     supporting cell-type count; requires max ABC score > 0.1 and a protein-coding
     target gene.
Coding evidence wins over ABC evidence for the same (trait, pop, gene, variant_id).

Enhancer-interval membership uses the half-open convention consistent with the ABC
file being BED-derived (BED intervals are [start, end), so a 1-based position p is
inside the interval iff start < p <= end). Interval overlap is computed once per
chromosome via two np.searchsorted calls against a sorted array of fine-mapped
variant positions (no bedtools).

Data source (hg19): ABC-Max predictions AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3
from Nasser et al., Nature 2021.

Input:
  --variants  gwas_finemapped_variants.tsv (step 11): trait, pop, variant_id, median_pip
  --coding    gene_nominations.tsv (nominate_coding_genes.py)
  --lift-vcf  GRCh37 VCF annotated by VEP (prep_vep_input.py, normalized by
              run_vep_annotation.sh), ID = hg38 variant_id, POS = hg19 --
              variants are hg38 but the ABC file is hg19, so overlap uses this
              build-matched position, not the hg38 position encoded in variant_id
  --abc       ABC-Max predictions file (hg19)
  --gtf       GENCODE annotation GTF (used only for its protein_coding gene set)

Output:
  gene_nominations_comprehensive.tsv  columns: trait, pop, gene, variant_id,
    implication_type, consequence, median_pip, abc_max_score, abc_n_celltypes

Usage:
  python nominate_genes_comprehensive.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import csv
import gzip
import getopt
from collections import defaultdict
from pathlib import Path

import numpy as np

VARIANTS_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/credible_sets/gwas_finemapped_variants.tsv"
CODING_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination/gene_nominations.tsv"
LIFT_VCF_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination/vep_input_all_hg19.normalized.vcf"
OUT_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination/gene_nominations_comprehensive.tsv"
# ABC-Max and GENCODE are external resources only available in our tree.
ABC_DEFAULT = "/lus/flare/projects/GeomicVar/jianhuawang/betainfer/data/raw/abc/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"
GTF_DEFAULT = "/lus/flare/projects/GeomicVar/jianhuawang/betainfer/data/raw/annotation/gencode.v19.annotation.gtf.gz"

PIP_MIN, ABC_MIN = 0.01, 0.1
# fine-mapping chromosome codes -> ABC "chrN" naming (X/Y/MT collapsed as in prep_vep_input.py)
CHR_MAP = {str(i): f"chr{i}" for i in range(1, 23)}
CHR_MAP.update({"23": "chrX", "24": "chrY", "25": "chrX", "26": "chrM", "X": "chrX", "Y": "chrY"})


def help(exit_num=0):
    print("""Usage: python nominate_genes_comprehensive.py [OPTIONS]
Optional:
  --variants   Fine-mapped variants table, step 11 output (default: {v})
  --coding     Coding gene nominations, nominate_coding_genes.py output (default: {c})
  --lift-vcf   GRCh37 VCF annotated by VEP, run_vep_annotation.sh output (default: {l})
  --abc        ABC-Max predictions file, hg19 (default: {a})
  --gtf        GENCODE annotation GTF (default: {g})
  --out        Output path (default: {o})
  -h, --help   Show this help message""".format(
        v=VARIANTS_DEFAULT, c=CODING_DEFAULT, l=LIFT_VCF_DEFAULT,
        a=ABC_DEFAULT, g=GTF_DEFAULT, o=OUT_DEFAULT))
    sys.exit(exit_num)


###############################################################################
##########  REFERENCE LOADING  ################################################
###############################################################################

def protein_coding_genes(gtf_path):
    """Return the set of gene symbols GENCODE marks gene_type protein_coding."""
    keep = set()
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line[0] == "#" or "\tgene\t" not in line:
                continue
            if 'gene_type "protein_coding"' in line:
                keep.add(line.split('gene_name "')[1].split('"')[0])
    return keep


def load_lift(lift_vcf_path):
    """variant_id -> hg19 position, from the lifted VCF (ID col = original hg38 variant_id)."""
    lift = {}
    with open(lift_vcf_path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            lift[f[2]] = int(f[1])
    return lift


def load_finemapped(variants_path, lift):
    """Return (per_variant, chr_index).
    per_variant: variant_id -> [(trait, pop, median_pip), ...] (median_pip >= PIP_MIN).
    chr_index: ABC chrom -> (sorted hg19 position ndarray, aligned variant_id ndarray).
    Overlap uses the hg19 position from `lift` (build-matched to the ABC file), not
    the hg38 position encoded in variant_id.
    """
    per_variant = defaultdict(list)
    pos_by_chr = defaultdict(list)
    with open(variants_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mp = row["median_pip"]
            if not mp or float(mp) < PIP_MIN:
                continue
            vid = row["variant_id"]
            p = lift.get(vid)
            if p is None:
                continue
            per_variant[vid].append((row["trait"], row["pop"], float(mp)))
            ac = CHR_MAP.get(vid.split(":")[0])
            if ac:
                pos_by_chr[ac].append((p, vid))
    chr_index = {}
    for ac, lst in pos_by_chr.items():
        lst = sorted(set(lst))
        chr_index[ac] = (np.array([p for p, _ in lst]), np.array([v for _, v in lst], dtype=object))
    return per_variant, chr_index


###############################################################################
##########  ABC OVERLAP  #######################################################
###############################################################################

def abc_links(abc_path, chr_index):
    """Stream the ABC file once; return {variant_id: {gene: [max_abc, celltypes, is_promoter]}}
    for links with ABC score > ABC_MIN.

    Interval membership is start < pos <= end (BED half-open convention): the left
    bound uses searchsorted 'right' at start, not 'left', so a variant sitting exactly
    at an enhancer's start boundary is correctly excluded.
    """
    hits = defaultdict(lambda: defaultdict(lambda: [0.0, set(), False]))
    with gzip.open(abc_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(header)}
        c_chr, c_s, c_e = ci["chr"], ci["start"], ci["end"]
        c_gene, c_cls, c_self = ci["TargetGene"], ci["class"], ci["isSelfPromoter"]
        c_abc, c_ct = ci["ABC.Score"], ci["CellType"]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                score = float(f[c_abc])
            except (ValueError, IndexError):
                continue
            if score <= ABC_MIN:
                continue
            idx = chr_index.get(f[c_chr])
            if idx is None:
                continue
            pos, vids = idx
            lo = np.searchsorted(pos, int(f[c_s]), "right")   # start < pos
            hi = np.searchsorted(pos, int(f[c_e]), "right")   # pos <= end
            if lo == hi:
                continue
            gene, cls, is_self, ct = f[c_gene], f[c_cls], f[c_self], f[c_ct]
            promoter = (cls == "promoter" and is_self == "True")
            for vid in vids[lo:hi]:
                rec = hits[vid][gene]
                rec[0] = max(rec[0], score)
                rec[1].add(ct)
                rec[2] = rec[2] or promoter
    return hits


###############################################################################
##########  MAIN  ##############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "variants=", "coding=", "lift-vcf=", "abc=", "gtf=", "out="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    variants_path = VARIANTS_DEFAULT
    coding_path = CODING_DEFAULT
    lift_vcf_path = LIFT_VCF_DEFAULT
    abc_path = ABC_DEFAULT
    gtf_path = GTF_DEFAULT
    out_path = OUT_DEFAULT

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--variants":
            variants_path = a
        elif o == "--coding":
            coding_path = a
        elif o == "--lift-vcf":
            lift_vcf_path = a
        elif o == "--abc":
            abc_path = a
        elif o == "--gtf":
            gtf_path = a
        elif o == "--out":
            out_path = a
        else:
            assert False, f"Unhandled option: {o}"

    pc_genes = protein_coding_genes(gtf_path)
    lift = load_lift(lift_vcf_path)
    per_variant, chr_index = load_finemapped(variants_path, lift)
    n_var = sum(len(v[1]) for v in chr_index.values())
    print(f"NOTE: fine-mapped unique variants (median PIP>={PIP_MIN}): {n_var}; protein_coding genes: {len(pc_genes)}")

    hits = abc_links(abc_path, chr_index)
    print(f"NOTE: variants with >=1 ABC link (score>{ABC_MIN}): {len(hits)}")

    # ── Coding nominations win over ABC for the same (trait, pop, gene, variant) ──
    coding_keys = set()
    coding_rows = []
    with open(coding_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            coding_keys.add((r["trait"], r["pop"], r["gene"], r["variant_id"]))
            coding_rows.append([r["trait"], r["pop"], r["gene"], r["variant_id"],
                                "nonsynonymous_coding_variant", r["consequence"], r["median_pip"], "", ""])

    # ── ABC rows, expanded to (trait, pop) via fine-mapped membership ──
    abc_rows = []
    for vid, genes in hits.items():
        for gene, (max_abc, cts, promoter) in genes.items():
            if gene not in pc_genes:
                continue
            impl = "abc_promoter_overlap" if promoter else "abc_interaction"
            for trait, pop, mp in per_variant.get(vid, []):
                if (trait, pop, gene, vid) in coding_keys:   # coding wins
                    continue
                abc_rows.append([trait, pop, gene, vid, impl, "", f"{mp:.6g}",
                                 f"{max_abc:.4f}", str(len(cts))])

    header = ["trait", "pop", "gene", "variant_id", "implication_type",
              "consequence", "median_pip", "abc_max_score", "abc_n_celltypes"]
    all_rows = coding_rows + abc_rows
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(all_rows)

    ng = len({(r[0], r[1], r[2]) for r in all_rows})
    print(f"NOTE: comprehensive nominations: {len(all_rows)} rows "
          f"({len(coding_rows)} coding + {len(abc_rows)} ABC); {ng} distinct (trait,pop,gene)")
    print(f"NOTE:  -> {out_path}")


if __name__ == "__main__":
    main()
