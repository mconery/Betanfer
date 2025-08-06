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
    return parser.parse_args()

def get_causal_variants_map(causal_var_dir):
    # Mapping: key = <Pheno_name>.<num_causal_variants>_causal.effect_size_<effect_size>, value = causal_var_file
    mapping = dict()
    causal_pattern = re.compile(r'(.+?\.\d+_causal\.effect_size_\d*\.?\d*)\.500000_space100000\.txt')
    for f in glob(os.path.join(causal_var_dir, "*.txt")):
        m = causal_pattern.match(os.path.basename(f))
        if m:
            mapping[m.group(1)] = f
    return mapping

def parse_output_filename(filename):
    """
    Matches:
      <Pheno_name>.<num_causal_variants>_causal.effect_size_<effect_size>.(binary_<prevalence>|continuous).<chr#>.step2.txt
    Returns a metadata dict with possible keys: pheno_name, causal_part, trait_type, chromosome, prevalence.
    """
    base = os.path.basename(filename)
    # Regex captures both binary_<prevalence> and continuous
    pattern = re.compile(
        r"(.+?)\.(\d+_causal\.effect_size_\d*\.?\d*)\.(binary_(\d*\.?\d+)|continuous)\.(chr\d+)\.step2\.txt$"
    )
    m = pattern.match(base)
    if not m:
        return None
    pheno_name, causal_part, trait_type_full, prevalence, chrom = m.groups()
    trait_type = "binary" if trait_type_full.startswith("binary") else "continuous"
    meta = {
        "pheno_name": pheno_name,
        "causal_part": causal_part,
        "trait_type": trait_type,
        "chromosome": chrom,
    }
    if trait_type == "binary":
        meta["prevalence"] = float(prevalence) if prevalence else None
    return meta

def extract_effect_size(causal_part):
    # causal_part is like '1_causal.effect_size_5'
    m = re.search(r"effect_size_(\d*\.?\d+)", causal_part)
    return float(m.group(1)) if m else None

def extract_causal_results_for_file(out_file, variant_ids, population, meta):
    df = pd.read_csv(out_file, sep=None, engine="python", dtype=str)
    if not all(col in df.columns for col in ["MarkerID", "BETA", "p.value"]):
        return pd.DataFrame()  # Return empty if columns missing
    
    # Filter rows with causal variant IDs
    filtered = df[df["MarkerID"].isin(variant_ids)].copy()
    if filtered.empty:
        return filtered

    # Select only the columns needed
    filtered = filtered[["MarkerID", "BETA", "p.value"]]

    # Add metadata columns
    filtered["population"] = population
    filtered["pheno_name"] = meta["pheno_name"]

    # For trait_type include prevalence if binary, e.g. "binary_0.25"
    if meta["trait_type"] == "binary" and "prevalence" in meta and meta["prevalence"] is not None:
        filtered["trait_type"] = f"binary_{meta['prevalence']}"
    else:
        filtered["trait_type"] = meta["trait_type"]

    # Extract effect_size from causal_part
    filtered["effect_size"] = extract_effect_size(meta["causal_part"])

    return filtered

def main():
    args = parse_args()
    output_dir = args.output_dir
    causal_var_dir = args.causal_var_dir
    out_file = args.out

    if args.populations:
        populations = args.populations
    else:
        # Use all directories in output_dir
        populations = [d for d in os.listdir(output_dir)
                         if os.path.isdir(os.path.join(output_dir, d))]

    causal_map = get_causal_variants_map(causal_var_dir)

    # Gather all output files and map to their experiment meta info
    results = []
    for pop in populations:
        pop_dir = os.path.join(output_dir, pop)
        for file in glob(os.path.join(pop_dir, "*.step2.txt")):
            meta = parse_output_filename(file)
            if not meta:
                continue
            causal_key = f"{meta['pheno_name']}.{meta['causal_part']}"
            causal_file = causal_map.get(causal_key)
            if not causal_file or not os.path.isfile(causal_file):
                continue
            # Load causal variant IDs for this set
            cdf = pd.read_csv(causal_file, sep=None, engine="python", header=None)
            variant_ids = set(cdf.iloc[:,0].astype(str))
            meta["trait_type"] = meta["trait_type"]  # e.g., 'binary'/'continuous'
            meta["chromosome"] = meta["chromosome"]
            meta["pheno_name"] = meta["pheno_name"]
            # Extract only relevant rows
            extracted = extract_causal_results_for_file(file, variant_ids, pop, meta)
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
