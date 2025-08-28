#!/usr/bin/env python3

import argparse
import os
import re
import pandas as pd
from glob import glob

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract causal variant results from SAIGE step2 output files."
    )
    
    parser.add_argument(
        "--output-dir", required=True,
        help="Directory containing population subdirectories with step2 output files."
    )
    
    parser.add_argument(
        "--causal-var-dir", required=True,
        help="Directory containing causal variant files."
    )
    
    parser.add_argument(
        "--out", required=True,
        help="Path for the combined output file (CSV)."
    )
    
    parser.add_argument(
        "--populations", nargs="+", default=None,
        help="List of populations to process (default: all directories in output-dir)."
    )
    
    parser.add_argument(
        "--legacy-filenames", default=False, action='store_true',
        help="Switch to expect output files WITHOUT heritability in filename."
    )
    
    return parser.parse_args()

def get_causal_variants_map(causal_var_dir):
    # Mapping: key = <pheno_name>.<num_causal>_causal.heritability_<heritability>, value = causal_var_file
    mapping = dict()
    causal_pattern = re.compile(r'(.+?\.\d+_causal)\.variant_set\.txt')
    
    for f in glob(os.path.join(causal_var_dir, "*.txt")):
        m = causal_pattern.match(os.path.basename(f))
        if m:
            mapping[m.group(1)] = f
    
    return mapping

def parse_output_filename(filename, legacy_filenames=False):
    """
    Matches:
      <Pheno_name>.<num_causal_variants>_causal.heritability_<heritability>.(binary_<prevalence>|continuous).<chr#>.step2.txt
    Returns a metadata dict with possible keys: pheno_name, causal_part, trait_type, chromosome, prevalence.
    """
    base = os.path.basename(filename)
    # Regex captures both binary_<prevalence> and continuous
    if legacy_filenames:
        # No heritability in filename
        pattern = re.compile(
            r"(.+?)\.(\d+_causal)\.(binary_(\d*\.?\d+)|continuous)\.(chr\d+)\.step2\.txt$"
        )
    else:
        pattern = re.compile(
            r"(.+?)\.(\d+_causal)\.(heritability_\d*\.?\d*)\.(binary_(\d*\.?\d+)|continuous)\.(chr\d+)\.step2\.txt$"
        )
    
    m = pattern.match(base)
    if not m:
        return None
    
    pheno_name, causal_part, heritability, trait_type_full, prevalence, chrom = m.groups()
    trait_type = "binary" if trait_type_full.startswith("binary") else "continuous"
    
    meta = {
        "pheno_name": pheno_name,
        "causal_part": causal_part,
        "heritability": heritability,
        "trait_type": trait_type,
        "chromosome": chrom,
    }
    
    if trait_type == "binary":
        meta["prevalence"] = float(prevalence) if prevalence else None
    
    return meta

def extract_heritability(causal_part):
    # causal_part is like '1_causal.heritability_5'
    m = re.search(r"heritability_(\d*\.?\d+)", causal_part)
    return float(m.group(1)) if m else None

def extract_num_causal_variants(causal_part):
    """
    Extract number of causal variants from causal_part
    causal_part examples: '1_causal', '5_causal.heritability_0.1'
    """
    m = re.match(r"(\d+)_causal", causal_part)
    return int(m.group(1)) if m else None

def extract_causal_results_for_file(out_file, variant_ids, population, meta, legacy_filenames=False):
    df = pd.read_csv(out_file, sep=None, engine="python", dtype=str)
    if not all(col in df.columns for col in ["MarkerID", "BETA", "p.value"]):
        return pd.DataFrame()
    
    filtered = df[df["MarkerID"].isin(variant_ids)].copy()
    if filtered.empty:
        return filtered
    
    filtered = filtered[["MarkerID", "BETA", "p.value"]]
    filtered["population"] = population
    filtered["pheno_name"] = meta["pheno_name"]
    
    if meta["trait_type"] == "binary" and "prevalence" in meta and meta["prevalence"] is not None:
        filtered["trait_type"] = f"binary_{meta['prevalence']}"
    else:
        filtered["trait_type"] = meta["trait_type"]
    
    if not legacy_filenames:
        filtered["heritability"] = extract_heritability(meta["causal_part"])
    
    # Add the number of causal variants column
    filtered["num_causal_variants"] = extract_num_causal_variants(meta["causal_part"])
    
    return filtered

def main():
    args = parse_args()
    output_dir = args.output_dir
    causal_var_dir = args.causal_var_dir
    out_file = args.out
    legacy_filenames = args.legacy_filenames
    
    if args.populations:
        populations = args.populations
    else:
        populations = [d for d in os.listdir(output_dir) if os.path.isdir(os.path.join(output_dir, d))]
    
    # All keys in this mapping include heritability
    causal_map = get_causal_variants_map(causal_var_dir)
    
    results = []
    
    for pop in populations:
        pop_dir = os.path.join(output_dir, pop)
        
        for file in glob(os.path.join(pop_dir, "*.step2.txt")):
            meta = parse_output_filename(file, legacy_filenames=legacy_filenames)
            if not meta:
                continue
            
            if legacy_filenames:
                # Find the first matching causal_map entry for this pheno_name and causal_part prefix (ignore heritability)
                # meta["causal_part"] will be e.g. "2_causal"
                candidates = [
                    k for k in causal_map.keys()
                    if k.startswith(f"{meta['pheno_name']}.{meta['causal_part']}")
                ]
                
                if not candidates:
                    continue
                
                # If exactly one, use it
                causal_key = candidates[0]
            else:
                causal_key = f"{meta['pheno_name']}.{meta['causal_part']}"
            
            causal_file = causal_map.get(causal_key)
            if not causal_file or not os.path.isfile(causal_file):
                continue
            
            cdf = pd.read_csv(causal_file, engine="python", header=None)
            variant_ids = set(cdf.iloc[:,0].astype(str))
            
            extracted = extract_causal_results_for_file(
                file, variant_ids, pop, meta, legacy_filenames=legacy_filenames
            )
            
            if not extracted.empty:
                results.append(extracted)
    
    if not results:
        print("No results found.")
        return
    
    combined = pd.concat(results, axis=0, ignore_index=True)
    combined.to_csv(out_file, index=False)
    print(f"Combined results written to {out_file}")

if __name__ == "__main__":
    main()