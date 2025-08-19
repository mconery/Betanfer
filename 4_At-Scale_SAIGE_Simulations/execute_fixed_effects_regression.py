#!/usr/bin/env python3.9

import argparse
import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import re
from statsmodels.api import Logit, OLS, add_constant
import multiprocessing
from functools import partial

def get_all_variants_with_info(plink_prefix):
    """Load variant IDs and annotation info from BIM file."""
    bim_file = plink_prefix + ".bim"
    if not os.path.exists(bim_file):
        raise FileNotFoundError(f"BIM file not found: {bim_file}")
    bim_df = pd.read_csv(
        bim_file, delimiter="\t", header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"]
    )
    return bim_df

def extract_genotypes(plink_prefix, snps_to_extract, tmp_dir):
    """Extract genotypes for given SNPs using PLINK2 --export A."""
    snp_file = os.path.join(tmp_dir, "snp_subset.txt")
    pd.Series(snps_to_extract).to_csv(snp_file, index=False, header=False)
    out_pref = os.path.join(tmp_dir, "extract")
    cmd = [
        "plink2",
        "--bfile", plink_prefix,
        "--extract", snp_file,
        "--export", "A",
        "--out", out_pref
    ]
    subprocess.run(cmd, check=True)
    geno = pd.read_csv(out_pref + ".raw", sep='\s+')

    # Fix SNP column names: remove allele suffix from IDs
    def normalize_snp_col(col):
        if ":" in col:
            return re.sub(r"_[ACGT]+$", "", col)
        else:
            return col

    baseline_cols = geno.columns[:6]
    snp_cols = [normalize_snp_col(c) for c in geno.columns[6:]]
    geno.columns = list(baseline_cols) + snp_cols
    return geno

def single_snp_regression(args):
    """Run regression for a single SNP (for parallel mapping)."""
    snp, merged, covars, bim_info, trait_type = args
    try:
        X = merged[[snp] + covars]
        X = add_constant(X, has_constant='add')
        y = merged["PHENO"]
        annot = bim_info.set_index("SNP").loc[snp]
        chr_num = annot["CHR"]
        pos = annot["BP"]
        ref_allele = annot["A1"]
        alt_allele = annot["A2"]
        geno_series = merged[snp].dropna()
        aaf = geno_series.mean() / 2

        if trait_type == "binary":
            model = Logit(y, X, missing="drop")
        else:
            model = OLS(y, X, missing="drop")
        res = model.fit(disp=0)
        beta = res.params[snp]
        se = res.bse[snp]
        p = res.pvalues[snp]
    except Exception:
        # If regression fails (monomorphic, singular, etc), report None
        beta, se, p, aaf = None, None, None, None

    return {
        "SNP": snp,
        "CHR": chr_num,
        "BP": pos,
        "REF": ref_allele,
        "ALT": alt_allele,
        "AAF": aaf,
        "Beta": beta,
        "StdErr": se,
        "Pval": p
    }

def main():
    parser = argparse.ArgumentParser(description="Parallel GWAS regression for all variants.")
    parser.add_argument('--plink_dir', required=True, help='Directory containing PLINK binary files')
    parser.add_argument('--plink_pattern', required=True, help='PLINK filename pattern (e.g., chr1.plink.maf)')
    parser.add_argument('--pheno_file', required=True, help='Phenotype file with #FID, IID, PHENO, PC1-PC10')
    parser.add_argument('--trait_type', choices=["binary", "continuous"], required=True)
    parser.add_argument('--out_results', required=True, help='Output TSV file for results')
    args = parser.parse_args()

    plink_prefix = os.path.join(args.plink_dir, args.plink_pattern)
    bim_df = get_all_variants_with_info(plink_prefix)
    all_variants = bim_df["SNP"].tolist()

    pheno = pd.read_csv(args.pheno_file, sep='\s+')
    covars = ["PC{}".format(i) for i in range(1, 11)]
    pheno = pheno[["#FID", "IID", "PHENO"] + covars]

    # PLINK2 is reasonably fast; keep extracting all at once for one batch unless you have memory issues.
    with tempfile.TemporaryDirectory() as tmpdir:
        geno = extract_genotypes(plink_prefix, all_variants, tmpdir)

    pheno = pheno.rename(columns={"#FID": "FID"})
    merged = pd.merge(pheno, geno, on=["FID", "IID"])

    num_cores = multiprocessing.cpu_count()
    print(f"Using {num_cores} CPU cores for regression.")

    # Prepare arguments for each process
    arglist = [(snp, merged, covars, bim_df, args.trait_type) for snp in all_variants]

    with multiprocessing.Pool(num_cores) as pool:
        results = pool.map(single_snp_regression, arglist, chunksize=10)

    results_df = pd.DataFrame(results)
    results_df.to_csv(args.out_results, sep="\t", index=False)
    print(f"Results written to {args.out_results}")

if __name__ == "__main__":
    main()
