#!/usr/bin/env python3
"""
Generate Multi-SuSiE multi-population fine-mapping commands
==========================================================

This script scans a directory tree that contains:

• GWAS summary-statistics files  
• phenotype files (with the **same population / locus naming scheme as the GWAS files**)  
• LD‐matrix files in `<matrix_type>/<POP>/` sub-directories

For every phenotype × matrix-type combination that is found in *all* requested
populations, the script writes one command line that will run the
Multi-SuSiE wrapper (`finemap_MultiSuSiE.py` **or** `finemap_MultiSuSiE.R`).

The output is a plain-text file that can be executed (e.g. with `bash` or
submitted as a SLURM array job).

Differences vs. the original SuSiEx version
-------------------------------------------
1. **Removed** all Plink-reference logic (`--geno-dir`, `--plink_conv`,
   `--n_gwas`, `--ref_panels`, etc.).
2. **Added** phenotype discovery and the `--pheno_files` and
   `--ancestry_labels` CLI arguments required by Multi-SuSiE.
3. **Wrapper script** can be either a `.py` or a `.R` file; the launcher
   (`python` vs `Rscript`) is chosen automatically.

Author: *adapted from generate_susiex_commands.py*  
Date  : 2025-09-25
"""

import argparse
import glob
import os
from pathlib import Path
from multiprocessing import Pool, cpu_count
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def detect_matrix_subdirs(ld_base_dir: str, matrix_filter=None) -> List[str]:
    """Return all matrix-type subdirs that contain at least one population."""
    if not os.path.exists(ld_base_dir):
        print(f"⚠️  LD base directory does not exist: {ld_base_dir}")
        return []
    
    common_pops = {"AFR", "EUR", "EAS", "SAS", "AMR"}
    matrix_types = []
    
    for entry in os.listdir(ld_base_dir):
        mdir = os.path.join(ld_base_dir, entry)
        if not os.path.isdir(mdir):
            continue
        if matrix_filter and entry not in matrix_filter:
            continue
        if any(os.path.isdir(os.path.join(mdir, p)) for p in common_pops):
            matrix_types.append(entry)
    
    return matrix_types

def normalize_heritability(heritability_str):
    """
    Normalize heritability string for comparison.
    Converts decimal representations to consistent format.
    
    Parameters:
    - heritability_str: Heritability string from filename (e.g., "heritability_0.0003")
    
    Returns:
    - float: Normalized heritability value
    """
    try:
        # Remove "heritability_" prefix if present
        if heritability_str.startswith("heritability_"):
            heritability_str = heritability_str.replace("heritability_", "")
        
        # Convert to float for consistent comparison
        return float(heritability_str)
    except (ValueError, TypeError):
        return None

def parse_gwas_filename(meta_filename, target_heritability=None):
    """
    Parse GWAS filename to extract components, handling decimal heritability values.
    Expected format: Pheno#.#_causal.heritability_#.#.continuous/binary_0.25.chr#.step2.txt
    
    Parameters:
    - meta_filename: The GWAS filename
    - target_heritability: Target heritability value to filter for (None = no filter)
    
    Returns:
    - dict: Dictionary with parsed components or None if parsing fails
    """
    # Remove the .step2.txt extension
    base_name = meta_filename.replace('.step2.txt', '')
    
    # Split by dots
    parts = base_name.split('.')
    if len(parts) < 5:  # Need at least: Pheno#, #_causal, heritability_#, #, trait_type, chr#
        return None
    
    try:
        # First two parts: Pheno# and #_causal
        pheno_causal = parts[0] + "." + parts[1]  # e.g., "Pheno1.2_causal"
        
        # Find heritability parts - they start with "heritability_" and continue until we hit trait type
        heritability_parts = []
        trait_type_index = None
        
        for i in range(2, len(parts)):
            part = parts[i]
            # Check if this is a trait type (continuous or binary_0.25)
            if part in ['continuous', 'binary_0'] or (part == '25' and i > 2 and parts[i-1] == 'binary_0'):
                trait_type_index = i
                break
            else:
                heritability_parts.append(part)
        
        if trait_type_index is None:
            return None
        
        # Reconstruct heritability (e.g., "heritability_0.01" from ["heritability_0", "01"])
        heritability = ".".join(heritability_parts)
        
        # Filter by heritability if specified
        if target_heritability is not None:
            parsed_heritability = normalize_heritability(heritability)
            if parsed_heritability is None:
                return None
            
            # Compare heritability values with small tolerance for floating point comparison
            if abs(parsed_heritability - target_heritability) > 1e-10:
                return None
        
        # Handle trait type - could be "continuous" or "binary_0.25"
        if parts[trait_type_index] == 'continuous':
            trait_type = 'continuous'
            chromosome_index = trait_type_index + 1
        elif parts[trait_type_index] == 'binary_0' and trait_type_index + 1 < len(parts) and parts[trait_type_index + 1] == '25':
            trait_type = 'binary_0.25'
            chromosome_index = trait_type_index + 2
        else:
            return None
        
        # Get chromosome
        if chromosome_index >= len(parts):
            return None
        chromosome = parts[chromosome_index]
        
        return {
            'pheno_causal': pheno_causal,
            'heritability': heritability,
            'trait_type': trait_type,
            'chromosome': chromosome
        }
    
    except (IndexError, ValueError):
        return None
    
def find_files(pop_dir: str, pattern: str) -> str | None:
    """Return a single file matching *pattern* inside *pop_dir*, or None."""
    matches = glob.glob(os.path.join(pop_dir, pattern))
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        print(f"⚠️  Multiple matches in {pop_dir}: {matches} – using first")
        return matches[0]
    return None

# ---------------------------------------------------------------------------
# Core processing
# ---------------------------------------------------------------------------

def process_one_pattern(args_tuple) -> Tuple[List[str], List[str]]:
    """
    Build Multi-SuSiE commands for a single phenotype/locus pattern.
    Returns (cmd_list, warnings)
    """
    gwas_example_file, gwas_base, pheno_base, ld_base, matrix_types, wrapper, out_base, populations, skip_existing, target_heritability = args_tuple
    
    warnings: List[str] = []
    cmds: List[str] = []
    meta = parse_gwas_filename(os.path.basename(gwas_example_file))
    if meta is None:
        warnings.append(f"Unparsable GWAS file name: {gwas_example_file}")
        return cmds, warnings
    pheno_causal = meta["pheno_causal"]
    trait_type   = meta["trait_type"]
    chrom        = meta["chromosome"]
    heritability = meta["heritability"]
    # ------------------------------------------------------------------ GWAS
    gwas_files: Dict[str, str] = {}
    for pop in populations:
        g_dir = os.path.join(gwas_base, pop)
        g_file = find_files(g_dir,
                            f"{pheno_causal}.{heritability}.{trait_type}.{chrom}.step2.txt")
        if g_file:
            gwas_files[pop] = g_file
        else:
            warnings.append(f"Missing GWAS for {pop}: {pheno_causal}, {chrom}")
    if len(gwas_files) != len(populations):
        return cmds, warnings
    # ------------------------------------------------------------- phenotype
    pheno_files: Dict[str, str] = {}
    for pop in populations:
        p_dir = os.path.join(pheno_base, pop)
        p_file = find_files(
            p_dir,
            f"{pheno_causal}.{heritability}.{trait_type}.txt"
        )
        if p_file:
            pheno_files[pop] = p_file
        else:
            warnings.append(f"Missing PHENO for {pop}: {pheno_causal}, {chrom}")
    if len(pheno_files) != len(populations):
        return cmds, warnings
    # -------------------------------------------------------------- matrices
    for m_type in matrix_types:
        ld_files: Dict[str, str] = {}
        for pop in populations:
            search_dir = os.path.join(ld_base, m_type, pop)
            ld_file = find_files(search_dir, f"{pheno_causal}.*.vcor1")
            if ld_file:
                ld_files[pop] = ld_file
            else:
                warnings.append(
                    f"Missing LD for {m_type}/{pop}: {pheno_causal}")
        if len(ld_files) != len(populations):
            continue
        # Output path
        out_dir = os.path.join(out_base, m_type)
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        out_file = os.path.join(
            out_dir, f"{pheno_causal}.{trait_type}.{m_type}.multisusie.txt"
        )
        if skip_existing and os.path.isfile(out_file):
            continue
        # ---------------------------------------------------- build command
        gwas_list   = ",".join(gwas_files[p]  for p in populations)
        ld_list     = ",".join(ld_files[p]    for p in populations)
        pheno_list  = ",".join(pheno_files[p] for p in populations)
        labels_list = ",".join(populations)
        launcher = "Rscript" if wrapper.endswith(".R") else "python"
        cmd = (
            f"{launcher} {wrapper} "
            f"--gwas_files {gwas_list} "
            f"--ld_files {ld_list} "
            f"--pheno_files {pheno_list} "
            f"--output_file {out_file} "
            f"--ancestry_labels {labels_list}"
        )
        cmds.append(cmd)
    return cmds, warnings

# ---------------------------------------------------------------------------
# High-level driver
# ---------------------------------------------------------------------------

def generate_multisusie_commands(
    gwas_dir: str,
    pheno_dir: str,
    ld_dir: str,
    wrapper_script: str,
    out_dir: str,
    cmd_file: str,
    populations: List[str],
    matrix_filter: List[str] | None,
    processes: int | None,
    skip_existing: bool,
    target_heritability=None
):
    # Print heritability filter information
    if target_heritability is not None:
        print(f"Heritability filter: {target_heritability}")
        print(f"Only processing traits with heritability = {target_heritability}")
    else:
        print(f"Heritability filter: DISABLED - processing all heritability values")
    # Discover matrix types
    matrix_types = detect_matrix_subdirs(ld_dir, matrix_filter)
    if not matrix_types:
        print("❌ No matrix types detected – check --ld-dir or --matrix-types.")
        return
    # Discover unique GWAS patterns using the *first* population folder
    first_pop_dir = os.path.join(gwas_dir, populations[0])
    example_files = glob.glob(os.path.join(first_pop_dir, "*_causal*.step2.txt"))
    example_files = [f for f in example_files 
                    if 'drop_causal' not in os.path.basename(f) 
                    and 'matched_to_' not in os.path.basename(f)]
    patterns: set[str] = set()
    gwas_example_list: List[str] = []
    for f in example_files:
        meta = parse_gwas_filename(os.path.basename(f), target_heritability=target_heritability)
        if meta:
            key = f"{meta['pheno_causal']}.{meta['trait_type']}.{meta['chromosome']}"
            if key not in patterns:
                patterns.add(key)
                gwas_example_list.append(f)
    # Parallel processing
    pool_size = processes or cpu_count()
    args_list = [
        (
            f,
            gwas_dir,
            pheno_dir,
            ld_dir,
            matrix_types,
            wrapper_script,
            out_dir,
            populations,
            skip_existing, 
            target_heritability
        )
        for f in gwas_example_list
    ]
    all_cmds: List[str] = []
    all_warnings: List[str] = []
    with Pool(pool_size) as pool:
        for cmd_list, warn_list in pool.map(process_one_pattern, args_list):
            all_cmds.extend(cmd_list)
            all_warnings.extend(warn_list)
    # Write commands
    Path(Path(cmd_file).parent).mkdir(parents=True, exist_ok=True)
    with open(cmd_file, "w") as fh:
        for c in all_cmds:
            fh.write(c + "\n")
    print(f"✅  Wrote {len(all_cmds)} commands → {cmd_file}")
    if all_warnings:
        print(f"⚠️  {len(all_warnings)} warnings – see stderr for details")
        for w in all_warnings:
            print("  ", w)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def cli() -> None:
    default_root   = "/path/to/project"        # ← edit these defaults
    default_gwas   = f"{default_root}/GWAS"
    default_pheno  = f"{default_root}/phenotypes"
    default_ld     = f"{default_root}/LD_matrices"
    default_out    = f"{default_root}/multisusie_results"
    default_cmds   = f"{default_out}/multisusie_commands.txt"
    default_wrap   = f"{default_root}/finemap_MultiSuSiE.py"
    p = argparse.ArgumentParser(
        description="Generate command file for Multi-SuSiE fine-mapping."
    )
    p.add_argument("--wrapper-script", "-w", default=default_wrap,
                   help="Path to finemap_MultiSuSiE.{py|R}")
    p.add_argument("--gwas-dir", "-g", default=default_gwas,
                   help="Base directory containing GWAS sub-dirs per population")
    p.add_argument("--pheno-dir", "-p", default=default_pheno,
                   help="Base directory containing phenotype sub-dirs per population")
    p.add_argument("--ld-dir", "-l", default=default_ld,
                   help="Base directory with LD matrices")
    p.add_argument("--output-dir", "-o", default=default_out,
                   help="Directory where Multi-SuSiE results will go")
    p.add_argument("--commands-file", "-c", default=default_cmds,
                   help="Destination *.txt file holding all commands")
    p.add_argument("--populations", "-P", default="AFR,EUR",
                   help="Comma-delimited population codes (default: AFR,EUR)")
    p.add_argument("--matrix-types", "-m", default=None,
                   help="Comma list of matrix types to keep (else: all)")
    p.add_argument("--processes", "-t", type=int, default=None,
                   help="Parallel processes (default: CPU count)")
    p.add_argument("--skip-existing", "-s", action="store_true",
                   help="Skip loci whose output already exists")
    p.add_argument('--heritability', type=float, default=0.0003,
                        help='Target heritability value to filter traits (default: 0.0003)')
    args = p.parse_args()
    pops = [s.strip().upper() for s in args.populations.split(",") if s.strip()]
    mfilter = (
        [m.strip() for m in args.matrix_types.split(",") if m.strip()]
        if args.matrix_types
        else None
    )
    generate_multisusie_commands(
        gwas_dir=args.gwas_dir,
        pheno_dir=args.pheno_dir,
        ld_dir=args.ld_dir,
        wrapper_script=args.wrapper_script,
        out_dir=args.output_dir,
        cmd_file=args.commands_file,
        populations=pops,
        matrix_filter=mfilter,
        processes=args.processes,
        skip_existing=args.skip_existing,
        target_heritability=args.heritability
    )


if __name__ == "__main__":
    cli()
