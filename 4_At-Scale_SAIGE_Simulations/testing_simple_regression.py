'''
testing_simple_regression.py

The purpose of this script is to run a quick regression to see if we are getting
close to the correct significance level for the causal variants.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

import argparse
import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
from statsmodels.api import Logit, OLS, add_constant

def load_variants(variants_file):
    return pd.read_csv(variants_file, delimiter="\t", header=None, names=["SNP", "Other"])

def load_phenos(pheno_file):
    pheno = pd.read_csv(pheno_file, delim_whitespace=True)
    covars = ["PC{}".format(i) for i in range(1, 11)]
    return pheno[["#FID", "IID", "PHENO"] + covars]

def extract_genotypes(plink_prefix, snps, tmp_dir):
    snp_file = os.path.join(tmp_dir, "snp_subset.txt")
    snps["SNP"].to_csv(snp_file, index=False, header=False)
    out_pref = os.path.join(tmp_dir, "extract")
    cmd = [
        "plink2",
        "--bfile", plink_prefix,
        "--extract", snp_file,
        "--export", "A",
        "--out", out_pref
    ]
    subprocess.run(cmd, check=True)
    geno = pd.read_csv(out_pref + ".raw", delim_whitespace=True)
    return geno

def run_regressions(merged, snp_list, trait_type):
    covars = ["PC{}".format(i) for i in range(1, 11)]
    results = []
    for snp in snp_list["SNP"]:
        X = merged[[snp] + covars]
        X = add_constant(X)
        y = merged["PHENO"]
        if trait_type == "binary":
            model = Logit(y, X, missing="drop")
        else:
            model = OLS(y, X, missing="drop")
        res = model.fit(disp=0)
        beta = res.params[snp]
        se = res.bse[snp]
        p = res.pvalues[snp]
        results.append({"SNP": snp, "Beta": beta, "StdErr": se, "Pval": p})
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink_prefix', required=True)
    parser.add_argument('--pheno_file', required=True)
    parser.add_argument('--variants_file', required=True)
    parser.add_argument('--trait_type', choices=["binary", "continuous"], required=True)
    args = parser.parse_args()

    snps = load_variants(args.variants_file)
    pheno = load_phenos(args.pheno_file)

    with tempfile.TemporaryDirectory() as tmpdir:
        geno = extract_genotypes(args.plink_prefix, snps, tmpdir)

    merged = pd.merge(pheno, geno, on=["#FID", "IID"])
    results = run_regressions(merged, snps, args.trait_type)
    
    print(results.to_string(index=False))

if __name__ == "__main__":
    main()
