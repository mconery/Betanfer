#!/usr/bin/env python3

"""
Generate dummy results commands for failed MultiSuSiE fine-mapping runs.

This script:
1. Scans the expected output directory for MultiSuSiE results
2. Identifies which loci should have results but don't (failed runs)
3. Generates create_dummy_results.py commands to create dummy output files

The output is a shell script of create_dummy_results.py commands that can be
executed to fill in dummy results for all failed loci.

Author: Generated for error-handling workflow
Date: 2026-02-06
"""

import argparse
import glob
import os
from pathlib import Path
from typing import Dict, List, Tuple


def parse_gwas_filename(meta_filename):
    """
    Parse GWAS filename to extract pheno/causal/heritability/trait_type/chromosome.
    Expected format: Pheno#.#_causal.heritability_#.#.continuous/binary_0.25.chr#.step2.txt
    
    Returns:
        dict with keys: pheno_causal, heritability, trait_type, chromosome
    """
    base_name = meta_filename.replace('.step2.txt', '')
    parts = base_name.split('.')
    
    if len(parts) < 5:
        return None
    
    try:
        pheno_causal = parts[0] + "." + parts[1]
        heritability_parts = []
        trait_type_index = None
        
        for i in range(2, len(parts)):
            part = parts[i]
            if part in ['continuous', 'binary_0'] or (part == '25' and i > 2 and parts[i-1] == 'binary_0'):
                trait_type_index = i
                break
            else:
                heritability_parts.append(part)
        
        if trait_type_index is None:
            return None
        
        heritability = ".".join(heritability_parts)
        
        if parts[trait_type_index] == 'continuous':
            trait_type = 'continuous'
            chromosome_index = trait_type_index + 1
        elif parts[trait_type_index] == 'binary_0' and trait_type_index + 1 < len(parts) and parts[trait_type_index + 1] == '25':
            trait_type = 'binary_0.25'
            chromosome_index = trait_type_index + 2
        else:
            return None
        
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


def find_all_gwas_files(gwas_base: str, populations: List[str]) -> List[Tuple[str, Dict]]:
    """
    Find all GWAS files across all populations.
    
    Returns:
        List of (gwas_file_path, parsed_metadata) tuples
    """
    all_files = []
    
    for pop in populations:
        pop_dir = os.path.join(gwas_base, pop)
        if not os.path.isdir(pop_dir):
            continue
        
        # Find all .step2.txt files
        gwas_files = glob.glob(os.path.join(pop_dir, "*.step2.txt"))
        for gwas_file in gwas_files:
            meta = parse_gwas_filename(os.path.basename(gwas_file))
            if meta:
                all_files.append((gwas_file, meta))
    
    return all_files


def find_expected_output_files(gwas_files: List[Tuple[str, Dict]], 
                               ld_base: str,
                               output_base: str,
                               matrix_types: List[str],
                               method: str) -> List[Dict]:
    """
    Build list of expected output files based on GWAS files and available LD matrices.
    
    Returns:
        List of dicts with keys: pheno_causal, trait_type, matrix_type, output_file,
                                  gwas_files (dict), vars_files (dict), populations
    """
    expected = []
    
    # Group GWAS files by pheno/heritability/trait
    by_pattern = {}
    for gwas_file, meta in gwas_files:
        pop_dir = os.path.dirname(gwas_file)
        pop = os.path.basename(pop_dir)
        
        key = (meta['pheno_causal'], meta['heritability'], meta['trait_type'])
        if key not in by_pattern:
            by_pattern[key] = {'gwas': {}, 'vars': {}}
        
        by_pattern[key]['gwas'][pop] = gwas_file
        
        # Find corresponding .vars file
        ld_prefix = meta['pheno_causal']
        vars_file = os.path.join(ld_base, pop, 
                                f"{ld_prefix}.*.vars")
        vars_matches = glob.glob(vars_file)
        if vars_matches:
            by_pattern[key]['vars'][pop] = vars_matches[0]
    
    # For each pattern, check each matrix type
    for (pheno_causal, heritability, trait_type), files in by_pattern.items():
        for matrix_type in matrix_types:
            # Build output filename
            if method == 'multisusie':
                output_file = os.path.join(
                    output_base, matrix_type,
                    f"{pheno_causal}.{trait_type}.{matrix_type}.{method}.txt"
                )
            elif method in ['susiex', 'mesusie']:
                output_file = os.path.join(
                    output_base, matrix_type,
                    f"{pheno_causal}.{trait_type}.{matrix_type}.txt"
                )
            else:
                output_file = os.path.join(
                    output_base, matrix_type,
                    f"{pheno_causal}.{trait_type}.weighted.txt"
                )
            
            expected.append({
                'pheno_causal': pheno_causal,
                'trait_type': trait_type,
                'matrix_type': matrix_type,
                'output_file': output_file,
                'gwas_files': files['gwas'],
                'vars_files': files['vars'],
            })
    
    return expected


def find_missing_files(expected_outputs: List[Dict]) -> List[Dict]:
    """
    Filter to only expected outputs that don't exist.
    
    Returns:
        List of failed output specifications
    """
    missing = []
    for spec in expected_outputs:
        if not os.path.isfile(spec['output_file']):
            missing.append(spec)
    
    return missing


def generate_dummy_commands(missing_files: List[Dict],
                           populations: List[str],
                           method_name: str = "multisusie") -> List[str]:
    """
    Generate create_dummy_results.py commands for all missing files.
    
    Returns:
        List of command strings
    """
    commands = []
    
    #Get directory of current script
    script_path = Path(__file__).resolve()      # full path to this .py file
    script_dir  = script_path.parent            # directory containing this file
    
    for spec in missing_files:
        gwas_files_dict = spec['gwas_files']
        vars_files_dict = spec['vars_files']
        output_file = spec['output_file']
        
        # Check that we have files for all populations
        if set(gwas_files_dict.keys()) != set(populations):
            print(f"WARNING: Missing GWAS files for some populations in {spec['pheno_causal']}")
            continue
        
        if set(vars_files_dict.keys()) != set(populations):
            print(f"WARNING: Missing .vars files for some populations in {spec['pheno_causal']}")
            continue
        
        # Build command
        gwas_files = ",".join(gwas_files_dict[pop] for pop in populations)
        vars_files = ",".join(vars_files_dict[pop] for pop in populations)
        pop_labels = ",".join(populations)
        
        cmd = (
            f"python {script_dir}/create_dummy_results.py \\\n"
            f"  --gwas_file '{gwas_files}' \\\n"
            f"  --vars_files '{vars_files}' \\\n"
            f"  --output_file '{output_file}' \\\n"
            f"  --populations '{pop_labels}' \\\n"
            f"  --variant_mode union \\\n"
            f"  --method {method_name} \\\n"
            f"  --calculate_z \\\n"
            f"  --create_timing"
        )
        commands.append(cmd)
    
    return commands


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate create_dummy_results.py commands for failed MultiSuSiE runs.\n\n"
            "This script identifies MultiSuSiE output files that should exist but don't,\n"
            "then generates create_dummy_results.py commands to create dummy results."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    p.add_argument(
        "--gwas_dir", "-g", required=True,
        help="Base directory containing GWAS files organized by population subdirectories"
    )
    p.add_argument(
        "--ld_dir", "-l", required=True,
        help="Base directory containing LD matrices organized by matrix_type/population"
    )
    p.add_argument(
        "--output_dir", "-o", required=True,
        help="Base output directory where MultiSuSiE results should be stored"
    )
    p.add_argument(
        "--commands_file", "-c", required=True,
        help="Output file to write the shell commands"
    )
    p.add_argument(
        "--populations", "-p", default="AFR,EUR",
        help="Comma-delimited population labels (default: AFR,EUR)"
    )
    p.add_argument(
        "--matrix_types", "-m", default="mvp_external",
        help="Comma-delimited matrix types to check (default: mvp_external)"
    )
    p.add_argument(
        "--method", default="multisusie",
        help="Method name for timing file header (default: multisusie)"
    )
    
    return p.parse_args()


def main():
    args = parse_args()
    
    populations = [p.strip() for p in args.populations.split(",") if p.strip()]
    matrix_types = [m.strip() for m in args.matrix_types.split(",") if m.strip()]
    
    print(f"Scanning for failed MultiSuSiE runs...")
    print(f"  GWAS dir: {args.gwas_dir}")
    print(f"  LD dir: {args.ld_dir}")
    print(f"  Output dir: {args.output_dir}")
    print(f"  Populations: {populations}")
    print(f"  Matrix types: {matrix_types}")
    print()
    
    # Find all GWAS files
    print("Finding GWAS files...")
    gwas_files = find_all_gwas_files(args.gwas_dir, populations)
    print(f"  Found {len(gwas_files)} GWAS files")
    
    # Find expected output files
    print("Building list of expected output files...")
    expected = find_expected_output_files(gwas_files, args.ld_dir, 
                                         args.output_dir, matrix_types, args.method)
    print(f"  Expected {len(expected)} output files")
    
    # Find missing files
    print("Checking for missing output files...")
    missing = find_missing_files(expected)
    print(f"  Found {len(missing)} missing files")
    print()
    
    if len(missing) == 0:
        print("✓ All expected output files exist. Nothing to do.")
        return
    
    # Generate commands
    print(f"Generating create_dummy_results.py commands...")
    commands = generate_dummy_commands(missing, populations, args.method)
    print(f"  Generated {len(commands)} commands")
    print()
    
    # Write output file
    print(f"Writing commands to: {args.commands_file}")
    
    output_dir = os.path.dirname(args.commands_file)
    if output_dir and not os.path.exists(output_dir):
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    with open(args.commands_file, 'w') as f:
        f.write("#!/bin/bash\n\n")
        f.write("# Auto-generated commands to create dummy results for failed MultiSuSiE runs\n")
        f.write(f"# Generated: {args.populations}\n")
        f.write(f"# Method: {args.method}\n")
        f.write(f"# Total commands: {len(commands)}\n\n")
        
        for i, cmd in enumerate(commands, 1):
            f.write(f"# Command {i}/{len(commands)}\n")
            f.write(cmd)
            f.write("\n\n")
    
    print(f"✓ Commands file created: {args.commands_file}")
    print(f"✓ Execute with: bash {args.commands_file}")


if __name__ == "__main__":
    main()
