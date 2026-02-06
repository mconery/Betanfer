#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd
import numpy as np


def read_vars_file(vars_file):
    if not os.path.exists(vars_file):
        raise FileNotFoundError(f".vars file not found: {vars_file}")
    with open(vars_file) as f:
        variant_ids = [line.strip() for line in f if line.strip()]
    if not variant_ids:
        raise ValueError(f".vars file is empty: {vars_file}")
    return variant_ids


def combine_variant_sets(vars_files, mode="union"):
    """
    vars_files: list of paths
    mode: 'union' or 'intersection'
    returns: sorted list of variant IDs
    """
    if mode not in {"union", "intersection"}:
        raise ValueError(f"variant_mode must be 'union' or 'intersection', got {mode}")

    sets = []
    for vf in vars_files:
        ids = read_vars_file(vf)
        sets.append(set(ids))

    if mode == "union":
        combined = set().union(*sets)
    else:  # intersection
        combined = set.intersection(*sets) if sets else set()

    if not combined:
        raise ValueError(
            f"No variants found after {mode} over {len(vars_files)} variant files"
        )

    return sorted(combined)


def detect_snp_column(gwas_file):
    try:
        df = pd.read_csv(gwas_file, nrows=0, sep="\t")
    except Exception:
        df = pd.read_csv(gwas_file, nrows=0, sep=" ")
    headers = df.columns.tolist()
    candidates = [
        "SNP", "snp", "rsid", "RSID", "ID", "id",
        "variant_id", "VARIANT_ID", "MarkerName", "MarkerID"
    ]
    for c in candidates:
        if c in headers:
            return c
    raise ValueError(
        f"Could not detect SNP column in {gwas_file}; "
        f"found columns: {headers}"
    )


def detect_beta_se_columns(gwas_file):
    try:
        df = pd.read_csv(gwas_file, nrows=0, sep="\t")
    except Exception:
        try:
            df = pd.read_csv(gwas_file, nrows=0, sep=" ")
        except Exception:
            return None, None
    headers = df.columns.tolist()
    beta_candidates = ["BETA", "beta", "Beta", "effect_size", "EFFECT_SIZE"]
    se_candidates = ["SE", "se", "standard_error", "STDERR", "stderr"]
    beta = next((h for h in headers if h in beta_candidates), None)
    se = next((h for h in headers if h in se_candidates), None)
    return beta, se


def load_z_scores(gwas_file, snp_col, variant_ids):
    """
    Returns: list of Z-scores aligned to variant_ids (NaN if missing)
    """
    beta_col, se_col = detect_beta_se_columns(gwas_file)
    if beta_col is None or se_col is None:
        return [np.nan] * len(variant_ids)

    try:
        df = pd.read_csv(gwas_file, sep="\t",
                         usecols=[snp_col, beta_col, se_col])
    except Exception:
        df = pd.read_csv(gwas_file, sep=" ",
                         usecols=[snp_col, beta_col, se_col])

    df = df.drop_duplicates(subset=[snp_col]).set_index(snp_col)

    z = []
    for vid in variant_ids:
        if vid in df.index:
            beta = df.loc[vid, beta_col]
            se = df.loc[vid, se_col]
            if pd.notna(beta) and pd.notna(se) and se > 0:
                z.append(beta / se)
            else:
                z.append(np.nan)
        else:
            z.append(np.nan)
    return z


def write_dummy_results(
    output_file,
    variant_ids,
    z_dict,
    single_no_pop_name=False,
    create_timing=False,
    method_name="multisusie",
):
    """
    z_dict: mapping from label -> list of Zs (same length as variant_ids)
    single_no_pop_name: if True, write a single column 'Z' instead of 'Z_label'
    method_name: name for timing file header (e.g., "susie", "multisusie")
    """
    df = pd.DataFrame({"VARIANT_ID": variant_ids})

    if single_no_pop_name:
        # Only one entry in z_dict
        label, z_vals = next(iter(z_dict.items()))
        df["Z"] = z_vals
    else:
        for label, z_vals in z_dict.items():
            df[f"Z_{label}"] = z_vals

    df["PIP"] = np.nan
    df["converged"] = np.nan
    df["cs"] = 0

    out_dir = os.path.dirname(output_file)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    df.to_csv(output_file, sep="\t", index=False, na_rep="NA")

    if create_timing:
        timing_dir = os.path.join(out_dir, "timing")
        os.makedirs(timing_dir, exist_ok=True)
        timing_file = os.path.join(timing_dir, os.path.basename(output_file))
        header = f"{method_name}_time_seconds"
        with open(timing_file, "w") as f:
            f.write(f"{header}\nNA\n")


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Create dummy fine-mapping results for failed loci.\n"
            "Accepts a single GWAS file with one or more variant files.\n"
            "Combines variants by union (default) or intersection."
        )
    )
    p.add_argument(
        "--gwas_file", "-g", required=True,
        help="GWAS summary statistics file (single file only)"
    )
    p.add_argument(
        "--vars_files", "-v", required=True,
        help="Comma-delimited list of variant files (.vars-style, one ID per line)"
    )
    p.add_argument(
        "--variant_mode", "-m", choices=["union", "intersection"],
        default="union",
        help="How to combine variant sets from vars_files (default: union)"
    )
    p.add_argument(
        "--output_file", "-o", required=True,
        help="Output results file"
    )
    p.add_argument(
        "--populations", "-p", default="",
        help=(
            "Comma-delimited population labels; must match number of variant files "
            "if provided. If omitted and there is exactly one variant file, "
            "a single 'Z' column is written."
        )
    )
    p.add_argument(
        "--method", default="multisusie",
        help=(
            "Name of the fine-mapping method; used for timing file header "
            "(default: multisusie). Timing file will contain "
            "'<method>_time_seconds' as header."
        )
    )
    p.add_argument(
        "--calculate_z", action="store_true",
        help="Compute Z from Beta/SE if available"
    )
    p.add_argument(
        "--create_timing", action="store_true",
        help="Also write timing file with NA"
    )
    return p.parse_args()


def main():
    args = parse_args()

    # Single GWAS file
    gwas_file = args.gwas_file.strip()
    if not gwas_file:
        raise ValueError("GWAS file path cannot be empty")

    # Parse variant files
    vars_files = [v.strip() for v in args.vars_files.split(",") if v.strip()]
    if not vars_files:
        raise ValueError("No valid paths given in --vars_files")

    # Parse populations; may be empty string
    pop_raw = [p.strip() for p in args.populations.split(",") if p.strip()]
    single_no_pop = False

    if len(pop_raw) == 0:
        # Special case: single variant file, no populations => column named 'Z'
        if len(vars_files) == 1:
            single_no_pop = True
            labels = ["SINGLE"]  # internal label, will be written as 'Z'
        else:
            raise ValueError(
                "No populations provided, but multiple variant files given; "
                "please supply one population name per variant file."
            )
    else:
        labels = pop_raw

    # Validate counts match
    if len(labels) != len(vars_files):
        raise ValueError(
            "Number of populations must match number of variant files.\n"
            f"Got {len(labels)} populations and {len(vars_files)} variant files."
        )

    # Combine variants over all vars_files
    variant_ids = combine_variant_sets(vars_files, mode=args.variant_mode)

    # Build Z dict: use the single GWAS file for all populations
    z_dict = {}
    if args.calculate_z:
        snp_col = detect_snp_column(gwas_file)
        for label in labels:
            z_vals = load_z_scores(gwas_file, snp_col, variant_ids)
            z_dict[label] = z_vals
    else:
        # No Z calculation: fill NaNs
        for label in labels:
            z_dict[label] = [np.nan] * len(variant_ids)

    write_dummy_results(
        output_file=args.output_file,
        variant_ids=variant_ids,
        z_dict=z_dict,
        single_no_pop_name=single_no_pop,
        create_timing=args.create_timing,
        method_name=args.method,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
