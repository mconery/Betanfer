'''

testing_simple_regression.py

The purpose of this script is to run a quick regression to see if we are getting
close to the correct significance level for the causal variants.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python. ***

'''

import argparse
import pandas as pd
import numpy as np
import subprocess
import tempfile
import os
import re
from statsmodels.api import Logit, OLS, add_constant

def load_variants(variants_file):
    return pd.read_csv(variants_file, delimiter="\t", header=None, names=["SNP", "Other"])

def extract_chr_from_variant(snp_id):
    # Handles formats like "1:12345:A:G"
    match = re.match(r"^(\d+)[^:]*:", snp_id)
    if match:
        return match.group(1)
    else:
        raise ValueError("Could not parse chromosome from SNP ID: {}".format(snp_id))

def load_phenos(pheno_file):
    pheno = pd.read_csv(pheno_file, delim_whitespace=True)
    covars = ["PC{}".format(i) for i in range(1, 11)]
    return pheno[["#FID", "IID", "PHENO"] + covars]

def get_all_variants_from_plink(plink_prefix):
    """Read all variant IDs from the .bim file"""
    bim_file = plink_prefix + ".bim"
    if not os.path.exists(bim_file):
        raise FileNotFoundError(f"BIM file not found: {bim_file}")
    
    # Read the .bim file (columns: chr, rsid, genetic_distance, bp_pos, allele1, allele2)
    bim_df = pd.read_csv(bim_file, delimiter="\t", header=None, 
                         names=["CHR", "SNP", "CM", "BP", "A1", "A2"])
    return bim_df["SNP"].tolist()

def select_random_variants(all_variants, causal_variants, n_random):
    """Select n_random variants from all_variants that are not in causal_variants"""
    if n_random == 0:
        return []
    
    # Get variants that are not in the causal set
    causal_set = set(causal_variants)
    available_variants = [v for v in all_variants if v not in causal_set]
    
    if len(available_variants) < n_random:
        print(f"Warning: Only {len(available_variants)} non-causal variants available, "
              f"but {n_random} requested. Using all available.")
        n_random = len(available_variants)
    
    return np.random.choice(available_variants, size=n_random, replace=False).tolist()

def extract_genotypes(plink_prefix, all_snps_to_extract, tmp_dir):
    """Extract genotypes for all SNPs in the list"""
    snp_file = os.path.join(tmp_dir, "snp_subset.txt")
    pd.Series(all_snps_to_extract).to_csv(snp_file, index=False, header=False)
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

    #Fix genotype column names
    def normalize_snp_col(col):
        if ":" in col:
            return re.sub(r"_[ACGT]+$", "", col)
        else:
            return col

    baseline_cols = geno.columns[:6]
    snp_cols = [normalize_snp_col(c) for c in geno.columns[6:]]
    geno.columns = list(baseline_cols) + snp_cols

    return geno

def run_regressions(merged, causal_variants, random_variants, trait_type):
    """Run regressions on both causal and random variants, labeling each type"""
    covars = ["PC{}".format(i) for i in range(1, 11)]
    results = []
    
    # Process causal variants
    for snp in causal_variants:
        if snp in merged.columns:
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
            results.append({"SNP": snp, "Beta": beta, "StdErr": se, "Pval": p, "Type": "causal"})
    
    # Process random variants
    for snp in random_variants:
        if snp in merged.columns:
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
            results.append({"SNP": snp, "Beta": beta, "StdErr": se, "Pval": p, "Type": "non-causal"})
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink_dir', required=True)
    parser.add_argument('--plink_pattern', default='chr*.plink.maf')
    parser.add_argument('--pheno_file', required=True)
    parser.add_argument('--variants_file', required=True)
    parser.add_argument('--trait_type', choices=["binary", "continuous"], required=True)
    parser.add_argument('--n_random', type=int, default=0, 
                        help='Number of random variants to select (default: 0)')
    args = parser.parse_args()

    # Load causal variants
    causal_snps = load_variants(args.variants_file)
    causal_variant_list = causal_snps["SNP"].tolist()
    
    # Extract chromosome from the first variant
    first_snp = causal_snps.iloc[0]['SNP']
    chr_num = extract_chr_from_variant(first_snp)
    
    # Build plink file prefix
    plink_basename = args.plink_pattern.replace('*', chr_num)
    plink_prefix = os.path.join(args.plink_dir, plink_basename)

    # Get all variants from plink file and select random ones
    all_variants = get_all_variants_from_plink(plink_prefix)
    random_variants = select_random_variants(all_variants, causal_variant_list, args.n_random)
    
    # Combine all variants to extract
    all_variants_to_extract = causal_variant_list + random_variants
    
    # Load phenotypes
    pheno = pd.read_csv(args.pheno_file, delim_whitespace=True)
    covars = ["PC{}".format(i) for i in range(1, 11)]
    pheno = pheno[["#FID", "IID", "PHENO"] + covars]

    # Extract genotypes for all variants
    with tempfile.TemporaryDirectory() as tmpdir:
        geno = extract_genotypes(plink_prefix, all_variants_to_extract, tmpdir)

    # Merge data
    pheno.rename(columns={"#FID": "FID"}, inplace=True)
    merged = pd.merge(pheno, geno, on=["FID", "IID"])

    # Run regressions
    results = run_regressions(merged, causal_variant_list, random_variants, args.trait_type)
    
    # Print results
    print(results.to_string(index=False))

if __name__ == "__main__":
    main()
