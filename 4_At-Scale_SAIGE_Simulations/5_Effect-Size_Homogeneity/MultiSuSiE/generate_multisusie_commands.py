#!/usr/bin/env python3

"""
Generate Multi-SuSiE multi-population fine-mapping commands
==========================================================
This script scans a directory tree that contains:
• GWAS summary-statistics files
• phenotype files (with the **same population / locus naming scheme as the GWAS files**)
• LD‐matrix files in `//` sub-directories
For every phenotype × matrix-type combination that is found in *all* requested
populations, the script writes one command line that will run the
Multi-SuSiE wrapper (`finemap_MultiSuSiE.py` **or** `finemap_MultiSuSiE.R`).
The output is a plain-text file that can be executed (e.g. with `bash` or
submitted as a SLURM array job).
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
        print(f"⚠️ LD base directory does not exist: {ld_base_dir}")
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

    if len(parts) < 5: # Need at least: Pheno#, #_causal, heritability_#, #, trait_type, chr#

        return None

    try:

        # First two parts: Pheno# and #_causal

        pheno_causal = parts[0] + "." + parts[1] # e.g., "Pheno1.2_causal"

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

        print(f"⚠️ Multiple matches in {pop_dir}: {matches} – using first")

        return matches[0]

    return None

# ---------------------------------------------------------------------------

# File discovery functions for different analysis types

# ---------------------------------------------------------------------------

def find_pheno_and_gwas_files_standard(gwas_base_dir, pheno_base_dir, populations, 
                                       pheno_causal, heritability, trait_type, chromosome):
    """
    Find matching GWAS and phenotype files for standard analysis (identical filenames across populations).
    
    Parameters:
    - gwas_base_dir: Base directory containing population subdirectories
    - pheno_base_dir: Base directory containing population subdirectories
    - populations: List of population names
    - pheno_causal: Phenotype pattern (e.g., "Pheno1.2_causal")
    - heritability: Heritability string
    - trait_type: Trait type (continuous or binary_0.25)
    - chromosome: Chromosome identifier
    
    Returns:
    - tuple: (gwas_files dict, pheno_files dict) or (empty, empty) if not all found
    """
    gwas_files = {}
    pheno_files = {}
    
    expected_gwas_filename = f"{pheno_causal}.{heritability}.{trait_type}.{chromosome}.step2.txt"
    expected_pheno_filename = f"{pheno_causal}.{heritability}.{trait_type}.txt"
    
    for pop in populations:
        pop_gwas_dir = os.path.join(gwas_base_dir, pop)
        pop_pheno_dir = os.path.join(pheno_base_dir, pop)
        
        # Find GWAS file
        gwas_path = os.path.join(pop_gwas_dir, expected_gwas_filename)
        if os.path.exists(gwas_path):
            gwas_files[pop] = gwas_path
        else:
            broad_pattern = f"{pheno_causal}.*.{trait_type}.{chromosome}.step2.txt"
            matching = glob.glob(os.path.join(pop_gwas_dir, broad_pattern))
            if matching:
                gwas_files[pop] = matching[0]
        
        # Find phenotype file
        pheno_path = os.path.join(pop_pheno_dir, expected_pheno_filename)
        if os.path.exists(pheno_path):
            pheno_files[pop] = pheno_path
        else:
            broad_pattern = f"{pheno_causal}.*.{trait_type}.txt"
            matching = glob.glob(os.path.join(pop_pheno_dir, broad_pattern))
            if matching:
                pheno_files[pop] = matching[0]
    
    return (gwas_files if len(gwas_files) == len(populations) else {}, 
            pheno_files if len(pheno_files) == len(populations) else {})

def find_pheno_and_gwas_files_dropped_causal(gwas_base_dir, pheno_base_dir, populations, central_pop,
                                             pheno_causal, heritability, trait_type, chromosome):
    """
    Find GWAS and phenotype files for dropped_causal analysis.
    Central population uses standard files. Other populations use drop_causal variant.
    
    Parameters:
    - gwas_base_dir: Base directory containing population subdirectories
    - pheno_base_dir: Base directory containing population subdirectories
    - populations: List of population names
    - central_pop: Reference population (uses standard files, others get drop_causal)
    - pheno_causal: Phenotype pattern
    - heritability: Heritability string
    - trait_type: Trait type
    - chromosome: Chromosome identifier
    
    Returns:
    - tuple: (gwas_files dict, pheno_files dict) or (empty, empty) if not all found
    """
    gwas_files = {}
    pheno_files = {}
    
    for pop in populations:
        pop_gwas_dir = os.path.join(gwas_base_dir, pop)
        pop_pheno_dir = os.path.join(pheno_base_dir, pop)
        
        if pop == central_pop:
            # Central population uses standard files
            gwas_filename = f"{pheno_causal}.{heritability}.{trait_type}.{chromosome}.step2.txt"
            pheno_filename = f"{pheno_causal}.{heritability}.{trait_type}.txt"
        else:
            # Other populations use drop_causal variant
            gwas_filename = f"{pheno_causal}.{heritability}.{trait_type}.{chromosome}.drop_causal.step2.txt"
            pheno_filename = f"{pheno_causal}.{heritability}.{trait_type}.txt"  # Phenotype is same
        
        # Find GWAS file
        gwas_path = os.path.join(pop_gwas_dir, gwas_filename)
        if os.path.exists(gwas_path):
            gwas_files[pop] = gwas_path
        else:
            if pop == central_pop:
                pattern = f"{pheno_causal}.*.{trait_type}.{chromosome}.step2.txt"
            else:
                pattern = f"{pheno_causal}.*.{trait_type}.{chromosome}.drop_causal.step2.txt"
            matching = glob.glob(os.path.join(pop_gwas_dir, pattern))
            if matching:
                gwas_files[pop] = matching[0]
        
        # Find phenotype file
        pheno_path = os.path.join(pop_pheno_dir, pheno_filename)
        if os.path.exists(pheno_path):
            pheno_files[pop] = pheno_path
        else:
            broad_pattern = f"{pheno_causal}.*.{trait_type}.txt"
            matching = glob.glob(os.path.join(pop_pheno_dir, broad_pattern))
            if matching:
                pheno_files[pop] = matching[0]
    
    return (gwas_files if len(gwas_files) == len(populations) else {}, 
            pheno_files if len(pheno_files) == len(populations) else {})

def find_pheno_and_gwas_files_matched_effects(gwas_base_dir, pheno_base_dir, populations, central_pop,
                                              pheno_causal, heritability, trait_type, chromosome):
    """
    Find GWAS and phenotype files for matched_effects analysis.
    Central population uses standard files. All other populations use matched_to_[central_pop] variant.
    
    Parameters:
    - gwas_base_dir: Base directory containing population subdirectories
    - pheno_base_dir: Base directory containing population subdirectories
    - populations: List of population names
    - central_pop: Central population (uses standard files, others matched to this)
    - pheno_causal: Phenotype pattern
    - heritability: Heritability string
    - trait_type: Trait type
    - chromosome: Chromosome identifier
    
    Returns:
    - tuple: (gwas_files dict, pheno_files dict) or (empty, empty) if not all found
    """
    gwas_files = {}
    pheno_files = {}
    
    for pop in populations:
        pop_gwas_dir = os.path.join(gwas_base_dir, pop)
        pop_pheno_dir = os.path.join(pheno_base_dir, pop)
        
        if pop == central_pop:
            # Central population uses standard files
            gwas_filename = f"{pheno_causal}.{heritability}.{trait_type}.{chromosome}.step2.txt"
            pheno_filename = f"{pheno_causal}.{heritability}.{trait_type}.txt"
        else:
            # Other populations use matched_to_[central_pop] variant
            gwas_filename = f"{pheno_causal}.{heritability}.{trait_type}.matched_to_{central_pop}.{chromosome}.step2.txt"
            pheno_filename = f"{pheno_causal}.{heritability}.{trait_type}.txt"  # Phenotype is same
        
        # Find GWAS file
        gwas_path = os.path.join(pop_gwas_dir, gwas_filename)
        if os.path.exists(gwas_path):
            gwas_files[pop] = gwas_path
        else:
            if pop == central_pop:
                pattern = f"{pheno_causal}.*.{trait_type}.{chromosome}.step2.txt"
            else:
                pattern = f"{pheno_causal}.*.{trait_type}.matched_to_{central_pop}.{chromosome}.step2.txt"
            matching = glob.glob(os.path.join(pop_gwas_dir, pattern))
            if matching:
                gwas_files[pop] = matching[0]
        
        # Find phenotype file
        pheno_path = os.path.join(pop_pheno_dir, pheno_filename)
        if os.path.exists(pheno_path):
            pheno_files[pop] = pheno_path
        else:
            broad_pattern = f"{pheno_causal}.*.{trait_type}.txt"
            matching = glob.glob(os.path.join(pop_pheno_dir, broad_pattern))
            if matching:
                pheno_files[pop] = matching[0]
    
    return (gwas_files if len(gwas_files) == len(populations) else {}, 
            pheno_files if len(pheno_files) == len(populations) else {})

# ---------------------------------------------------------------------------

# Core processing

# ---------------------------------------------------------------------------

def process_one_pattern(args_tuple) -> Tuple[List[str], List[str]]:

    """

    Build Multi-SuSiE commands for a single phenotype/locus pattern.

    Returns (cmd_list, warnings)

    """

    (gwas_example_file, gwas_base, pheno_base, ld_base, matrix_types, wrapper, out_base, 
     populations, skip_existing, target_heritability, analysis_type, central_pop) = args_tuple

    warnings: List[str] = []

    cmds: List[str] = []

    meta = parse_gwas_filename(os.path.basename(gwas_example_file))

    if meta is None:

        warnings.append(f"Unparsable GWAS file name: {gwas_example_file}")

        return cmds, warnings

    pheno_causal = meta["pheno_causal"]

    trait_type = meta["trait_type"]

    chrom = meta["chromosome"]

    heritability = meta["heritability"]

    # Find GWAS and phenotype files based on analysis type

    if analysis_type == "standard":

        gwas_files, pheno_files = find_pheno_and_gwas_files_standard(

            gwas_base, pheno_base, populations, pheno_causal, heritability, trait_type, chrom)

    elif analysis_type == "dropped_causal":

        gwas_files, pheno_files = find_pheno_and_gwas_files_dropped_causal(

            gwas_base, pheno_base, populations, central_pop, pheno_causal, heritability, trait_type, chrom)

    elif analysis_type == "matched_effects":

        gwas_files, pheno_files = find_pheno_and_gwas_files_matched_effects(

            gwas_base, pheno_base, populations, central_pop, pheno_causal, heritability, trait_type, chrom)

    else:

        warnings.append(f"Unknown analysis type: {analysis_type}")

        return cmds, warnings

    # Verify we have files for all populations

    if len(gwas_files) != len(populations):

        missing_pops = [p for p in populations if p not in gwas_files]

        warnings.append(f"Missing GWAS for {missing_pops}: {pheno_causal}, {chrom}")

        return cmds, warnings

    if len(pheno_files) != len(populations):

        missing_pops = [p for p in populations if p not in pheno_files]

        warnings.append(f"Missing PHENO for {missing_pops}: {pheno_causal}, {chrom}")

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

        # Construct file lists for Multi-SuSiE

        gwas_list = ",".join([gwas_files[pop] for pop in populations])

        pheno_list = ",".join([pheno_files[pop] for pop in populations])

        ld_list = ",".join([ld_files[pop] for pop in populations])

        ancestry_labels = ",".join(populations)

        # Determine wrapper launcher

        if wrapper.endswith(".py"):

            launcher = "python"

        elif wrapper.endswith(".R"):

            launcher = "Rscript"

        else:

            warnings.append(f"Unknown wrapper file type: {wrapper}")

            continue

        # Build Multi-SuSiE command

        # Rscript finemap_MultiSuSiE.R --gwas_files file1,file2 --pheno_files pheno1,pheno2 --ld_files ld1,ld2 --ancestry_labels AFR,EUR --output_file output.txt

        cmd = (f"{launcher} {wrapper} "

               f"--gwas_files {gwas_list} "

               f"--pheno_files {pheno_list} "

               f"--ld_files {ld_list} "

               f"--ancestry_labels {ancestry_labels} "

               f"--output_file {out_file}")

        cmds.append(cmd)

    return cmds, warnings

def parse_matrix_filter(matrix_string):

    """

    Parse comma-delimited matrix filter string.

    Parameters:

    - matrix_string: Comma-delimited string of matrix names

    Returns:

    - list: List of matrix names, or None if empty

    """

    if not matrix_string:

        return None

    # Split by comma and strip whitespace

    matrix_list = [name.strip() for name in matrix_string.split(',')]

    # Remove empty strings

    matrix_list = [name for name in matrix_list if name]

    return matrix_list if matrix_list else None

def parse_populations(populations_string):

    """

    Parse comma-delimited populations string.

    Parameters:

    - populations_string: Comma-delimited string of population names

    Returns:

    - list: List of population names

    """

    if not populations_string:

        return ["AFR", "EUR"] # Default populations

    # Split by comma and strip whitespace

    pop_list = [name.strip().upper() for name in populations_string.split(',')]

    # Remove empty strings

    pop_list = [name for name in pop_list if name]

    return pop_list if pop_list else ["AFR", "EUR"]

def process_causal_files(gwas_base_dir, populations, analysis_type, central_pop):

    """

    Filter GWAS causal files based on analysis type.

    For standard analysis: Include files without 'drop_causal' or 'matched_to_' in name.

    For dropped_causal analysis: Include standard files from central pop (others get drop_causal).

    For matched_effects analysis: Include standard files from central pop (others get matched_to_).

    Parameters:

    - gwas_base_dir: Base directory containing population subdirectories

    - populations: List of population names

    - analysis_type: Type of analysis ('standard', 'dropped_causal', or 'matched_effects')

    - central_pop: Central population for analysis

    Returns:

    - list: List of valid causal file patterns

    """

    # Use the central population as the source for discovering file patterns

    pop_gwas_dir = os.path.join(gwas_base_dir, central_pop)

    if not os.path.exists(pop_gwas_dir):

        return []

    all_causal_files = glob.glob(os.path.join(pop_gwas_dir, "*_causal*.step2.txt"))

    valid_files = []

    for causal_file in all_causal_files:

        filename = os.path.basename(causal_file)

        # All analysis types use standard files from central population

        if 'drop_causal' not in filename and 'matched_to_' not in filename:

            valid_files.append(causal_file)

    return valid_files

def create_multisusie_commands(gwas_base_dir, pheno_base_dir, ld_base_dir, output_base_dir, wrapper_script,

                               output_commands_file, populations, analysis_type, central_pop,

                               matrix_filter=None, n_processes=None, skip_existing=False, target_heritability=None):

    """

    Generate Multi-SuSiE multi-population fine-mapping commands for all loci using parallel processing.

    Parameters:

    - gwas_base_dir: Base directory containing GWAS summary statistics (with population subdirs)

    - pheno_base_dir: Base directory containing phenotype files (with population subdirs)

    - ld_base_dir: Base directory containing LD matrices

    - output_base_dir: Base directory for fine-mapping results

    - wrapper_script: Path to finemap_MultiSuSiE.py or finemap_MultiSuSiE.R script

    - output_commands_file: File to write commands to

    - populations: List of population names

    - analysis_type: Type of analysis ('standard', 'dropped_causal', or 'matched_effects')

    - central_pop: Central population for analysis

    - matrix_filter: Optional list of matrix subdirectory names to include

    - n_processes: Number of processes to use (None = use all available cores)

    - skip_existing: Whether to skip generating commands for existing complete output files

    - target_heritability: Target heritability value to filter for (None = no filter)

    """

    # Validate central_pop is in populations

    if central_pop not in populations:

        print(f"Error: central_pop '{central_pop}' not in populations {populations}")

        return []

    # Print analysis configuration

    print(f"\n{'='*60}")

    print(f"MULTISUSIE COMMAND GENERATION CONFIGURATION")

    print(f"{'='*60}")

    print(f"Analysis Type: {analysis_type}")

    print(f"Central Population: {central_pop}")

    print(f"Populations: {populations}")

    if target_heritability is not None:

        print(f"Heritability filter: {target_heritability}")

        print(f"Only processing traits with heritability = {target_heritability}")

    else:

        print(f"Heritability filter: DISABLED - processing all heritability values")

    # Automatically detect matrix types

    print(f"\nDetecting matrix subdirectories in: {ld_base_dir}")

    if matrix_filter:

        print(f"Filtering for matrix types: {matrix_filter}")

    matrix_types = detect_matrix_subdirs(ld_base_dir, matrix_filter)

    if not matrix_types:

        if matrix_filter:

            print(f"Error: No valid matrix configurations found matching filter: {matrix_filter}")

        else:

            print("Error: No valid matrix configurations found!")

        return []

    print(f"Matrix types to process: {matrix_types}")

    if skip_existing:

        print(f"File existence checking: ENABLED - will skip existing complete output files")

    else:

        print(f"File existence checking: DISABLED - will generate all commands")

    # Create output directories now

    print(f"\nCreating output directories...")

    for matrix_type in matrix_types:

        output_dir = os.path.join(output_base_dir, matrix_type)

        os.makedirs(output_dir, exist_ok=True)

        print(f" Created: {output_dir}")

    # Filter files based on analysis type

    causal_files = process_causal_files(gwas_base_dir, populations, analysis_type, central_pop)

    print(f"\nFound {len(causal_files)} causal files for {analysis_type} analysis")

    # Group files by phenotype pattern to avoid duplicates

    unique_patterns = set()

    gwas_file_patterns = []

    for causal_file in causal_files:

        filename = os.path.basename(causal_file)

        parsed = parse_gwas_filename(filename, target_heritability=target_heritability)

        if parsed:

            pattern_key = f"{parsed['pheno_causal']}.{parsed['trait_type']}.{parsed['chromosome']}"

            if pattern_key not in unique_patterns:

                unique_patterns.add(pattern_key)

                gwas_file_patterns.append(causal_file)

    print(f"Found {len(gwas_file_patterns)} unique phenotype patterns")

    # Determine number of processes

    if n_processes is None:

        n_processes = cpu_count()

    print(f"Using {n_processes} processes for parallel processing")

    # Prepare arguments for parallel processing

    process_args = []

    for gwas_file_pattern in gwas_file_patterns:

        args_tuple = (gwas_file_pattern, gwas_base_dir, pheno_base_dir, ld_base_dir, matrix_types,

                      wrapper_script, output_base_dir, populations, skip_existing, target_heritability,

                      analysis_type, central_pop)

        process_args.append(args_tuple)

    # Process file patterns in parallel

    print(f"\nProcessing {len(process_args)} file patterns in parallel...")

    all_commands = []

    all_warnings = []

    with Pool(processes=n_processes) as pool:

        results = pool.map(process_one_pattern, process_args)

    # Collect results

    for commands, warnings in results:

        all_commands.extend(commands)

        all_warnings.extend(warnings)

    # Print warnings

    if all_warnings:

        print(f"\nWarnings encountered ({len(all_warnings)}):") 

        for warning in all_warnings[:10]:

            print(f" {warning}")

        if len(all_warnings) > 10:

            print(f" ... and {len(all_warnings) - 10} more warnings")

    # Print processing summary

    print(f"\n{'='*60}")

    print(f"PROCESSING SUMMARY")

    print(f"{'='*60}")

    print(f"Total commands generated: {len(all_commands)}")

    # Show examples of commands

    if all_commands:

        print(f"\nCommand Examples (first 3):")

        for i, cmd in enumerate(all_commands[:3]):

            print(f" {i+1}: {cmd[:100]}...")

        if len(all_commands) > 3:

            print(f" ... and {len(all_commands) - 3} more commands")

    # Write commands to file

    if output_commands_file:

        print(f"\n{'='*60}")

        print(f"Writing {len(all_commands)} commands to: {output_commands_file}")

        print(f"{'='*60}")

        with open(output_commands_file, 'w') as f:

            for cmd in all_commands:

                f.write(cmd + '\n')

    return all_commands

def main():

    """Main function with command line argument parsing."""

    parser = argparse.ArgumentParser(

        description='Generate Multi-SuSiE multi-population fine-mapping commands',

        formatter_class=argparse.RawDescriptionHelpFormatter,

        epilog="""

ANALYSIS TYPES:

standard - Original analysis: paired identical filenames across populations

dropped_causal - One population has causal variant dropped, others use standard files

Pairing: central_pop uses standard files, others use drop_causal variant

matched_effects - Effect sizes matched to central population

Pairing: central_pop uses standard files, others use matched_to_[central_pop] variant

EXAMPLES:

Standard analysis (original):

python generate_multisusie_commands.py --populations EUR,AFR --analysis-type standard

Dropped causal with AFR as reference:

python generate_multisusie_commands.py --populations EUR,AFR --analysis-type dropped_causal --central-population AFR

Matched effects with EUR as reference:

python generate_multisusie_commands.py --populations EUR,AFR,EAS --analysis-type matched_effects --central-population EUR

""")

    # Define default paths

    default_wrapper = "/path/to/finemap_MultiSuSiE.R"

    default_gwas = "/path/to/gwas_summary_statistics"

    default_pheno = "/path/to/phenotypes"

    default_ld = "/path/to/ld_matrices"

    default_output = "/path/to/fine-mapping_results"

    default_commands = f"{default_output}/multisusie_commands.txt"

    parser.add_argument('--wrapper-script', '-w', default=default_wrapper,

                        help=f'Path to finemap_MultiSuSiE.py or .R script (default: {default_wrapper})')

    parser.add_argument('--gwas-dir', '-g', default=default_gwas,

                        help=f'Base directory containing GWAS summary statistics with population subdirectories (default: {default_gwas})')

    parser.add_argument('--pheno-dir', '-ph', default=default_pheno,

                        help=f'Base directory containing phenotype files with population subdirectories (default: {default_pheno})')

    parser.add_argument('--ld-dir', '-l', default=default_ld,

                        help=f'Base directory containing LD matrices (default: {default_ld})')

    parser.add_argument('--output-dir', '-o', default=default_output,

                        help=f'Base directory for fine-mapping results (default: {default_output})')

    parser.add_argument('--commands-file', '-c', default=default_commands,

                        help=f'Output file for commands (default: {default_commands})')

    parser.add_argument('--populations', '-pop', type=str, default="AFR,EUR",

                        help='Comma-delimited list of population names (default: "AFR,EUR")')

    parser.add_argument('--analysis-type', '-a', type=str, choices=['standard', 'dropped_causal', 'matched_effects'],

                        default='standard', help='Type of analysis to run (default: standard)')

    parser.add_argument('--central-population', '-cp', type=str, default=None,

                        help='Central/reference population for dropped_causal or matched_effects analyses. '

                             'For dropped_causal: this pop uses standard files, others use drop_causal. '

                             'For matched_effects: this pop uses standard files, others matched to this pop. '

                             'If not specified, uses first population in --populations list.')

    parser.add_argument('--matrix-types', '-m', type=str, default=None,

                        help='Comma-delimited list of matrix subdirectory names to process (e.g., "raw_in-sample,PC-corrected"). '

                             'If not specified, all detected matrix types will be processed.')

    parser.add_argument('--processes', '-p', type=int, default=None,

                        help='Number of processes to use for parallel processing (default: use all available cores)')

    parser.add_argument('--skip-existing', '-s', action='store_true',

                        help='Skip generating commands for output files that already exist')

    parser.add_argument('--heritability', type=float, default=None,

                        help='Target heritability value to filter traits (default: process all)')

    parser.add_argument('--dry-run', '-d', action='store_true',

                        help='Print first 10 commands without writing to file')

    parser.add_argument('--list-matrices', action='store_true',

                        help='List available matrix subdirectories and exit')

    args = parser.parse_args()

    # Parse populations

    populations = parse_populations(args.populations)

    print(f"Populations: {populations}")

    # Set central population

    central_pop = args.central_population if args.central_population else populations[0]

    print(f"Central population: {central_pop}")

    # Validate central population

    if central_pop not in populations:

        print(f"Error: central_population '{central_pop}' not in populations {populations}")

        return 1

    # Handle list-matrices option

    if args.list_matrices:

        print(f"Scanning for matrix subdirectories in: {args.ld_dir}")

        if not os.path.exists(args.ld_dir):

            print(f"Error: LD directory does not exist: {args.ld_dir}")

            return 1

        matrix_types = detect_matrix_subdirs(args.ld_dir)

        if matrix_types:

            print("\nAvailable matrix types:")

            for matrix_type in matrix_types:

                print(f" {matrix_type}")

        else:

            print("No valid matrix configurations found.")

        return 0

    # Check if directories exist

    for dir_path, name in [(args.gwas_dir, 'GWAS'), (args.pheno_dir, 'Phenotype'), (args.ld_dir, 'LD matrices')]:

        if not os.path.exists(dir_path):

            print(f"Error: {name} directory does not exist: {dir_path}")

            return 1

    # Check population directories exist

    for pop in populations:

        pop_dir = os.path.join(args.gwas_dir, pop)

        if not os.path.exists(pop_dir):

            print(f"Error: Population GWAS directory does not exist: {pop_dir}")

            return 1

        pop_dir = os.path.join(args.pheno_dir, pop)

        if not os.path.exists(pop_dir):

            print(f"Error: Population phenotype directory does not exist: {pop_dir}")

            return 1

    if not os.path.exists(args.wrapper_script):

        print(f"Error: Wrapper script does not exist: {args.wrapper_script}")

        return 1

    # Parse matrix filter

    matrix_filter = parse_matrix_filter(args.matrix_types)

    if matrix_filter:

        print(f"Matrix filter applied: {matrix_filter}")

    # Validate number of processes

    if args.processes is not None and args.processes < 1:

        print(f"Error: Number of processes must be >= 1, got {args.processes}")

        return 1

    # Generate commands

    commands = create_multisusie_commands(

        gwas_base_dir=args.gwas_dir,

        pheno_base_dir=args.pheno_dir,

        ld_base_dir=args.ld_dir,

        output_base_dir=args.output_dir,

        wrapper_script=args.wrapper_script,

        output_commands_file=args.commands_file if not args.dry_run else None,

        populations=populations,

        analysis_type=args.analysis_type,

        central_pop=central_pop,

        matrix_filter=matrix_filter,

        n_processes=args.processes,

        skip_existing=args.skip_existing,

        target_heritability=args.heritability

    )

    if args.dry_run:

        print("\nDry run - first 10 commands:")

        for i, cmd in enumerate(commands[:10]):

            print(f"{i+1:3d}: {cmd}")

        if len(commands) > 10:

            print(f"... and {len(commands) - 10} more commands")

    else:

        print(f"\nCommands written to: {args.commands_file}")

        print(f"Total commands: {len(commands)}")

        print("\nTo run all commands, you can use:")

        print(f" bash {args.commands_file}")

        print("\nOr submit as array job to SLURM")

    return 0

if __name__ == "__main__":

    exit(main())
