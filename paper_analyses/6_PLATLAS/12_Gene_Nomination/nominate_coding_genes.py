#!/usr/bin/env python3
'''
nominate_coding_genes.py

Nominates genes for fine-mapped variants via coding VEP consequence
(the coding evidence line of gene nomination, step 12).

A variant nominates every protein-coding gene it overlaps when:
  - its cross-method median PIP is >= --pip-thresh (default 0.01);
  - VEP's consequence for that gene is coding (missense / frameshift /
    stop or start gain-loss / splice acceptor-donor);
  - VEP's biotype for that gene is protein_coding.
Run with VEP's --per_gene: a variant that overlaps more than one gene gets
one row per coding gene (rows are not deduplicated by variant).

Input:
  --vep       VEP --per_gene tab output (run_vep_annotation.sh), columns
              Uploaded_variation, Location, Consequence, SYMBOL, Gene, BIOTYPE
  --variants  gwas_finemapped_variants.tsv (step 11), column median_pip

Output:
  gene_nominations.tsv  columns: trait, pop, gene, variant_id, consequence, median_pip
                         one row per (trait, pop, gene, variant_id) coding nomination

Usage:
  python nominate_coding_genes.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import csv
import getopt
from collections import defaultdict
from pathlib import Path

VEP_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination/vep_output.tsv"
VARIANTS_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/credible_sets/gwas_finemapped_variants.tsv"
OUT_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/gene_nomination/gene_nominations.tsv"

# Coding consequence set: missense + splice/start/stop gain-loss + frameshift LoF
CODING_CSQ = {
    "missense_variant",
    "frameshift_variant",
    "stop_gained", "stop_lost", "start_lost",
    "splice_acceptor_variant", "splice_donor_variant",
}


def help(exit_num=0):
    print("""Usage: python nominate_coding_genes.py [OPTIONS]
Optional:
  --vep          VEP --per_gene tab output (default: {vep})
  --variants     Fine-mapped variants table, step 11 output (default: {v})
  --out          Output path (default: {o})
  --pip-thresh   Minimum median PIP to nominate (default: 0.01)
  -h, --help     Show this help message""".format(vep=VEP_DEFAULT, v=VARIANTS_DEFAULT, o=OUT_DEFAULT))
    sys.exit(exit_num)


###############################################################################
##########  VEP PARSING  #######################################################
###############################################################################

def iter_coding(path):
    """Yield (variant_id, consequence, gene) for coding + protein_coding rows only.

    VEP --per_gene gives one row per overlapping gene for a variant; all
    matching rows are kept (not deduplicated by variant).
    """
    idx = None
    for line in open(path):
        if line.startswith("#"):
            if line.startswith("#Uploaded"):
                hdr = line[1:].rstrip("\n").split("\t")
                idx = {c: i for i, c in enumerate(hdr)}
            continue
        if idx is None:
            continue
        f = line.rstrip("\n").split("\t")
        csq = f[idx["Consequence"]].split(",")[0]   # most-severe consequence comes first
        sym = f[idx["SYMBOL"]] if "SYMBOL" in idx else "-"
        bt = f[idx["BIOTYPE"]] if "BIOTYPE" in idx else "-"
        if csq in CODING_CSQ and bt == "protein_coding" and sym not in ("-", ""):
            yield f[idx["Uploaded_variation"]], csq, sym


def load_variant_traits(path, pip_thresh):
    """gwas_finemapped_variants.tsv -> {variant_id: [(trait, pop, median_pip), ...]}."""
    vt = defaultdict(list)
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                mp = float(row["median_pip"])
            except (ValueError, KeyError):
                continue
            if mp >= pip_thresh:
                vt[row["variant_id"]].append((row["trait"], row["pop"], mp))
    return vt


###############################################################################
##########  MAIN  ##############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h", ["help", "vep=", "variants=", "out=", "pip-thresh="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    vep_path = VEP_DEFAULT
    variants_path = VARIANTS_DEFAULT
    out_path = OUT_DEFAULT
    pip_thresh = 0.01

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--vep":
            vep_path = a
        elif o == "--variants":
            variants_path = a
        elif o == "--out":
            out_path = a
        elif o == "--pip-thresh":
            pip_thresh = float(a)
        else:
            assert False, f"Unhandled option: {o}"

    vt = load_variant_traits(variants_path, pip_thresh)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    n_nom, coding_vars, genes = 0, set(), set()
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["trait", "pop", "gene", "variant_id", "consequence", "median_pip"])
        for vid, csq, sym in iter_coding(vep_path):
            if vid not in vt:
                continue
            coding_vars.add(vid)
            genes.add(sym)
            for trait, pop, mp in vt[vid]:
                w.writerow([trait, pop, sym, vid, csq, f"{mp:.6g}"])
                n_nom += 1

    print(f"NOTE: coding-nominating variants: {len(coding_vars)}; unique genes: {len(genes)}")
    print(f"NOTE: nominations (trait x gene x variant rows): {n_nom} -> {out_path}")


if __name__ == "__main__":
    main()
