'''
execute_fixed_effects_regression.py

The purpose of this script is to run a basic fixed-effects regression on a set 
of variants given an input phenotype file

'''

import argparse
import pandas as pd
import numpy as np
import multiprocessing
from statsmodels.api import Logit, OLS, add_constant
from pysnptools.snpreader import Bed
import sys

def worker(args):
    snp_idxs, bed, pheno_df, covar_list, trait_type, sid_list, chr_list, bp_list, a1_list, a2_list, progress_queue, total_snps = args
    results = []
    n_done = 0
    last_update = 0
    
    for snp_offset, snp_idx in enumerate(snp_idxs):
        snp_id = sid_list[snp_idx]
        chr_num = chr_list[snp_idx]
        bp_pos = bp_list[snp_idx]
        ref = a1_list[snp_idx]
        alt = a2_list[snp_idx]
        geno = bed[:, snp_idx].read().val.flatten()
        aaf = np.nanmean(geno) / 2
        
        df = pheno_df.copy()
        df['SNP'] = geno
        X = df[["SNP"] + covar_list]
        X = add_constant(X, has_constant='add')
        y = df["PHENO"]
        
        try:
            if trait_type == "binary":
                model = Logit(y, X, missing="drop")
            else:
                model = OLS(y, X, missing="drop")
            fit = model.fit(disp=0)
            beta = fit.params["SNP"]
            se = fit.bse["SNP"]
            p = fit.pvalues["SNP"]
        except Exception:
            beta = se = p = None
        
        results.append({"SNP": snp_id, "CHR": chr_num, "BP": bp_pos,
                        "REF": ref, "ALT": alt, "AAF": aaf,
                        "Beta": beta, "StdErr": se, "Pval": p})
        
        n_done += 1
        if progress_queue is not None:
            progress_queue.put(1)
    
    return results

def main():
    parser = argparse.ArgumentParser(description="GWAS fast regression, direct-from-bed, with progress.")
    parser.add_argument("--plink_prefix", required=True, help="Prefix of BED/BIM/FAM (no extension)")
    parser.add_argument("--pheno_file", required=True, help="Phenotype file with #FID IID PHENO PC1-PC10")
    parser.add_argument("--trait_type", choices=["binary", "continuous"], required=True)
    parser.add_argument("--out_results", required=True, help="Output file (TSV)")
    parser.add_argument("--n_cores", type=int, default=None, help="Cores to use (default: all)")
    args = parser.parse_args()

    # Load genotype info with Bed
    bed = Bed(args.plink_prefix, count_A1=True)       # Keep A1 dosage convention!
    sid_list = bed.sid                                # Variant IDs
    chr_list = bed.pos[:,0]                           # Chromosome numbers (as int or str)
    bp_list = bed.pos[:,2]                            # Physical positions
    n_snps = len(sid_list)
    
    # Read the .bim file that corresponds to .bed to get allele info
    bim = pd.read_csv(args.plink_prefix + '.bim', sep='\t', header=None)
    bim.columns = ['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2']
    a1_list = bim['A1'].tolist()              # A1 (effect/ref) allele
    a2_list = bim['A2'].tolist()              # A2 (other) allele
    
    # Load and merge phenotype+covariates
    pheno_df = pd.read_csv(args.pheno_file, sep=r'\s+')
    covar_list = [f"PC{i}" for i in range(1, 11)]
    pheno_df = pheno_df[["#FID", "IID", "PHENO"] + covar_list]
    pheno_df = pheno_df.rename(columns={"#FID": "FID"})
    fam_df = pd.read_csv(f"{args.plink_prefix}.fam", sep=r'\s+', header=None)
    fam_df.columns = ["FID","IID","PAT","MAT","SEX","PHENO_IGNORE"]
    pheno_df = pd.merge(pheno_df, fam_df[["FID","IID"]], on=["FID","IID"])
    pheno_df = pheno_df.sort_values(["FID", "IID"]).reset_index(drop=True)
    
    # Setup multiprocessing
    n_cores = args.n_cores or multiprocessing.cpu_count()
    chunk_size = int(np.ceil(n_snps / n_cores))
    chunks = [list(range(i*chunk_size, min((i+1)*chunk_size, n_snps))) for i in range(n_cores)]
    progress_queue = multiprocessing.Queue()

    pool_args = [
        (chunk, bed, pheno_df, covar_list, args.trait_type, sid_list, chr_list, bp_list, a1_list, a2_list, progress_queue, n_snps)
        for chunk in chunks
    ]

    # Launch processes
    pool = multiprocessing.Pool(n_cores)
    total_done = 0
    last_percent = -1

    # Start async
    results_async = [pool.apply_async(worker, (pa,)) for pa in pool_args]

    # Progress bar
    while total_done < n_snps:
        # Each sub-process puts a "1" into the queue for each SNP done.
        try:
            progress_queue.get(timeout=1)
            total_done += 1
            percent = int((total_done/n_snps)*100)
            if percent % 10 == 0 and percent != last_percent:
                print(f"{percent}% of variants analyzed ({total_done} / {n_snps})", flush=True)
                last_percent = percent
        except:
            if all([a.ready() for a in results_async]):
                break

    # Collect results
    pool.close()
    pool.join()
    all_results = []
    for res in results_async:
        all_results.extend(res.get())
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(args.out_results, sep="\t", index=False)
    print(f"Results written to {args.out_results}")

if __name__ == "__main__":
    main()