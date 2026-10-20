#!/usr/bin/env python3

"""

SuSiEx Multi-Population Fine-Mapping Script

Updated to auto-detect chromosome and BP range from *.vars files

and generate PLINK-format frequency files with A1 frequency

EXTENDED: Now includes:

- Automatic conversion of text LD matrices to PLINK bin4 format

- PLINK-based frequency file generation for better SuSiEx compatibility

- Generation of required _ref.bim files that SuSiEx expects

- Optional cleanup of generated bin4 files after processing

- Filtering and deduplication of GWAS summary statistics

"""

import argparse

import os

import sys

import subprocess

import pandas as pd

import numpy as np

import re

import struct

import time

def parse_arguments():

    parser = argparse.ArgumentParser(

        description='Multi-population fine-mapping using SuSiEx (auto chr/bp from LD matrix .vars, PLINK frequency generation)',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--gwas_files', '-g', required=True,

                        help='Comma-delimited list of GWAS summary statistics files')

    parser.add_argument('--ld_files', '-l', required=True,

                        help='Comma-delimited list of LD matrix file prefixes (each must have a corresponding .vars file)')

    parser.add_argument('--output_file', '-o', required=True,

                        help='Full path for final output file')

    parser.add_argument('--n_gwas', '-n', required=True,

                        help='Comma-delimited list of GWAS sample sizes, same order as GWAS files')

    parser.add_argument('--ref_panels', required=True,

                        help='Comma-delimited list of reference panel prefixes (.bed/.bim/.fam) for PLINK frequency calculation')

    parser.add_argument('--maf', type=float, default=0.005)

    parser.add_argument('--n_sig', type=int, default=10)

    parser.add_argument('--level', type=int, default=0.95)

    parser.add_argument('--min_purity', type=float, default=0.5)

    parser.add_argument('--threads', type=int, default=1)

    parser.add_argument('--keep_ambig', action='store_true')

    parser.add_argument('--mult_step', action='store_true')

    parser.add_argument('--pval_thresh', type=float, default=1e-5)

    parser.add_argument('--max_iter', type=int, default=100)

    parser.add_argument('--tol', type=float, default=1e-4)

    parser.add_argument('--susiex_bin', required=False,

                        default='/sw/frontier/mde/bif154_mde/tools/SuSiEx/bin/SuSiEx',

                        help='Path to SuSiEx executable')

    parser.add_argument('--plink_bin', required=False,

                        default='/sw/frontier/mde/bif154_mde/tools/plink',

                        help='Path to PLINK executable')

    return parser.parse_args()

def convert_text_ld_to_bin4(text_file, output_prefix):

    """

    Convert a text LD matrix to PLINK bin4 format.

    """

    print(f"[LD CONVERT] Reading LD matrix from: {text_file}")

    if not os.path.exists(text_file):

        print(f"[LD CONVERT] Error: Input file {text_file} not found!")

        return False

    try:

        ld_matrix = np.loadtxt(text_file, dtype=np.float64)

        print(f"[LD CONVERT] Successfully loaded matrix of size: {ld_matrix.shape}")

    except Exception as e:

        print(f"[LD CONVERT] Error reading file: {e}")

        return False

    if ld_matrix.shape[0] != ld_matrix.shape[1]:

        print(f"[LD CONVERT] Error: Matrix is not square! Shape: {ld_matrix.shape}")

        return False

    n_variants = ld_matrix.shape[0]

    diag_values = np.diag(ld_matrix)

    if not np.allclose(diag_values, 1.0, atol=0.01):

        print("[LD CONVERT] Warning: Diagonal values are not close to 1.0")

    if ld_matrix.min() < -1.1 or ld_matrix.max() > 1.1:

        print(f"[LD CONVERT] Warning: LD values outside expected range [-1, 1]")

    ld_matrix_f32 = ld_matrix.astype(np.float32)

    output_bin = f"{output_prefix}.ld.bin"

    print(f"[LD CONVERT] Writing binary file: {output_bin}")

    try:

        with open(output_bin, 'wb') as f:

            f.write(ld_matrix_f32.tobytes('C'))

        print(f"[LD CONVERT] Successfully converted to bin4 format!")

        return True

    except Exception as e:

        print(f"[LD CONVERT] Error writing binary file: {e}")

        return False

def delete_bin4_files(ld_bin_files):

    """Delete bin4 LD matrix files."""

    print("\n=== Cleaning up bin4 LD files ===")

    for binf in ld_bin_files:

        if os.path.exists(binf):

            os.remove(binf) # Uncomment to actually delete files

            print(f"Deleted: {binf}")

        else:

            print(f"Not found (skip): {binf}")

def delete_ref_bim_files(ref_bim_files):

    """Delete _ref.bim files created for SuSiEx."""

    print("\n=== Cleaning up _ref.bim files ===")

    for bim_file in ref_bim_files:

        if os.path.exists(bim_file):

            os.remove(bim_file) # Uncomment to actually delete files

            print(f"Deleted: {bim_file}")

        else:

            print(f"Not found (skip): {bim_file}")

def delete_susiex_output_files(output_file):

    """

    Delete SuSiEx's default output files after creating the SuSiE-format output.

    SuSiEx creates three main output files:

    - {output_name}.snp: Contains PIPs and LogBF values

    - {output_name}.cs: Contains credible set information

    - {output_name}.summary: Contains summary information

    """

    print("\n=== Cleaning up SuSiEx output files ===")

    # Determine output file paths

    output_base = output_file.replace('.txt', '').replace('.susiex', '')

    susiex_files = [

        f"{output_base}.snp",

        f"{output_base}.cs",

        f"{output_base}.summary"

    ]

    deleted_count = 0

    for susiex_file in susiex_files:

        if os.path.exists(susiex_file):

            try:

                os.remove(susiex_file)

                print(f"Deleted: {susiex_file}")

                deleted_count += 1

            except Exception as e:

                print(f"Error deleting {susiex_file}: {e}")

        else:

            print(f"Not found (skip): {susiex_file}")

    if deleted_count > 0:

        print(f"Successfully deleted {deleted_count} SuSiEx output files")

    else:

        print("No SuSiEx output files found to delete")

def remove_duplicates_from_frq_file(frq_file):

    """

    Read a PLINK .frq file, remove duplicate variants, and rewrite the file.

    Preserves the original file delimiter and header structure.

    """

    if not os.path.exists(frq_file):

        print(f"[DEDUP] ERROR: Frequency file not found: {frq_file}")

        return False

    print(f"[DEDUP] Checking for duplicates in: {frq_file}")

    try:

        # Read the entire file

        df = pd.read_csv(frq_file, sep="\t", dtype=str)

        original_count = len(df)

        if original_count == 0:

            print(f"[DEDUP] WARNING: Frequency file is empty")

            return True

        # Identify duplicate variants based on only column since reading in as a single column

        # Find duplicates

        duplicates = df[df.duplicated(keep=False)].iloc[:,0]

        n_duplicates = len(duplicates)

        if n_duplicates == 0:

            print(f"[DEDUP] No duplicates found in {frq_file}")

            return True

        print(f"[DEDUP] Found {n_duplicates} duplicate variants")

        duplicate_snps = duplicates.unique()

        print(f"[DEDUP] Duplicate SNP IDs: {duplicate_snps[:10]}{'...' if len(duplicate_snps) > 10 else ''}")

        # Remove duplicates, keeping the first occurrence

        df_dedup = df.drop_duplicates(keep='first')

        final_count = len(df_dedup)

        removed_count = original_count - final_count

        print(f"[DEDUP] Removed {removed_count} duplicate variants")

        print(f"[DEDUP] Final variant count: {final_count}")

        # Write back to the same file, preserving the delimiter

        df_dedup.to_csv(frq_file, sep="\n", index=False, na_rep='NA')

        print(f"[DEDUP] Successfully rewrote {frq_file} without duplicates")

        return True

    except Exception as e:

        print(f"[DEDUP] ERROR: Failed to process {frq_file}: {e}")

        return False

def remove_duplicates_from_bim_file(bim_file):

    """

    Read a PLINK .bim file, remove duplicate variants, and rewrite the file.

    Preserves the original file delimiter (typically tab) and structure.

    .bim files have columns: CHR SNP MORGANS BP A1 A2

    """

    if not os.path.exists(bim_file):

        print(f"[DEDUP BIM] ERROR: .bim file not found: {bim_file}")

        return False

    print(f"[DEDUP BIM] Checking for duplicates in: {bim_file}")

    try:

        # Read the file and detect delimiter

        with open(bim_file, 'r') as f:

            first_line = f.readline().strip()

            second_line = f.readline().strip() if f.tell() < os.path.getsize(bim_file) else ""

        # Detect delimiter (typically tab for .bim files)

        delimiter = '\t' # Default to tab

        if first_line.count(' ') > first_line.count('\t') and second_line.count(' ') > second_line.count('\t'):

            delimiter = ' '

        elif first_line.count(',') > first_line.count('\t'):

            delimiter = ','

        print(f"[DEDUP BIM] Detected delimiter: {'TAB' if delimiter == chr(9) else repr(delimiter)}")

        # Read the entire file - .bim files don't have headers

        column_names = ['CHR', 'SNP', 'MORGANS', 'BP', 'A1', 'A2']

        df = pd.read_csv(bim_file, sep=delimiter, header=None, names=column_names, dtype=str)

        original_count = len(df)

        if original_count == 0:

            print(f"[DEDUP BIM] WARNING: .bim file is empty")

            return True

        # Find duplicates based on SNP column (column 1, 0-indexed)

        duplicates = df[df.duplicated(subset=['SNP'], keep=False)]

        n_duplicates = len(duplicates)

        if n_duplicates == 0:

            print(f"[DEDUP BIM] No duplicates found in {bim_file}")

            return True

        print(f"[DEDUP BIM] Found {n_duplicates} duplicate variants")

        duplicate_snps = duplicates['SNP'].unique()

        print(f"[DEDUP BIM] Duplicate SNP IDs: {duplicate_snps[:10]}{'...' if len(duplicate_snps) > 10 else ''}")

        # Remove duplicates, keeping the first occurrence

        df_dedup = df.drop_duplicates(subset=['SNP'], keep='first')

        final_count = len(df_dedup)

        removed_count = original_count - final_count

        print(f"[DEDUP BIM] Removed {removed_count} duplicate variants")

        print(f"[DEDUP BIM] Final variant count: {final_count}")

        # Write back to the same file, preserving the delimiter and NO header

        df_dedup.to_csv(bim_file, sep=delimiter, index=False, header=False, na_rep='NA')

        print(f"[DEDUP BIM] Successfully rewrote {bim_file} without duplicates")

        return True

    except Exception as e:

        print(f"[DEDUP BIM] ERROR: Failed to process {bim_file}: {e}")

        return False

def filter_and_deduplicate_gwas(gwas_file, variant_ids, temp_dir="temp_gwas"):
    """
    Filter GWAS summary statistics to contain only variants in the locus
    and remove any duplicate variant IDs.
    
    Args:
        gwas_file: Path to original GWAS file
        variant_ids: List of variant IDs from .vars file
        temp_dir: Directory to store temporary files
    
    Returns:
        Path to filtered temporary GWAS file
    """
    print(f"[GWAS FILTER] Processing {os.path.basename(gwas_file)}")
    
    # Create temporary directory if it doesn't exist
    os.makedirs(temp_dir, exist_ok=True)
    
    # Create temporary filename
    base_name = os.path.basename(gwas_file)
    temp_file = os.path.join(temp_dir, f"filtered_{base_name}")
    
    try:
        # Try reading with tab separator first
        try:
            df = pd.read_csv(gwas_file, sep='\t')
        except:
            # Try space separator
            df = pd.read_csv(gwas_file, sep=' ')
        
        original_count = len(df)
        print(f"[GWAS FILTER] Original file has {original_count} variants")
        
        # Detect column structure
        detected_cols = detect_columns(gwas_file, verbose=False)
        snp_col_idx = detected_cols['snp'] - 1  # Convert to 0-based indexing
        
        # Get the SNP column name
        snp_col_name = df.columns[snp_col_idx]
        
        # Filter for variants in the locus
        print(f"[GWAS FILTER] Filtering for {len(variant_ids)} locus variants")
        df_filtered = df[df[snp_col_name].isin(variant_ids)]
        filtered_count = len(df_filtered)
        print(f"[GWAS FILTER] After locus filtering: {filtered_count} variants")
        
        # Remove duplicates based on SNP column
        duplicates_before = df_filtered.duplicated(subset=[snp_col_name], keep=False).sum()
        if duplicates_before > 0:
            print(f"[GWAS FILTER] Found {duplicates_before} duplicate variants")
            df_dedup = df_filtered.drop_duplicates(subset=[snp_col_name], keep='first')
            dedup_count = len(df_dedup)
            removed_count = filtered_count - dedup_count
            print(f"[GWAS FILTER] Removed {removed_count} duplicate variants")
        else:
            print(f"[GWAS FILTER] No duplicate variants found")
            df_dedup = df_filtered
            dedup_count = filtered_count
        
        # Write temporary file
        # Determine the original separator
        with open(gwas_file, 'r') as f:
            first_line = f.readline()
            if first_line.count('\t') > first_line.count(' '):
                sep = '\t'
            else:
                sep = ' '
        
        df_dedup.to_csv(temp_file, sep=sep, index=False, na_rep='NA')
        print(f"[GWAS FILTER] Created temporary file: {temp_file} ({dedup_count} variants)")
        
        return temp_file
        
    except Exception as e:
        print(f"[GWAS FILTER] ERROR: Failed to process {gwas_file}: {e}")
        return None

def create_filtered_gwas_files(gwas_files, all_variant_ids, temp_dir="temp_gwas"):
    """
    Create filtered and deduplicated temporary GWAS files for all populations.
    
    Args:
        gwas_files: List of original GWAS file paths
        all_variant_ids: List of variant ID lists (one per population)
        temp_dir: Directory to store temporary files
    
    Returns:
        List of temporary GWAS file paths
    """
    print("=== Creating filtered GWAS files ===")
    
    temp_gwas_files = []
    
    for i, gwas_file in enumerate(gwas_files):
        print(f"Population {i+1}:")
        
        # Use the variant IDs from the corresponding LD matrix
        variant_ids = all_variant_ids[i]
        
        temp_file = filter_and_deduplicate_gwas(gwas_file, variant_ids, temp_dir)
        
        if temp_file:
            temp_gwas_files.append(temp_file)
        else:
            print(f"[GWAS FILTER] ERROR: Failed to create filtered file for {gwas_file}")
            # Clean up any previously created temp files
            cleanup_temp_gwas_files(temp_gwas_files)
            sys.exit(1)
    
    return temp_gwas_files

def cleanup_temp_gwas_files(temp_gwas_files):
    """Delete temporary GWAS files."""
    print("\n=== Cleaning up temporary GWAS files ===")
    
    deleted_count = 0
    for temp_file in temp_gwas_files:
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
                print(f"Deleted: {temp_file}")
                deleted_count += 1
            except Exception as e:
                print(f"Error deleting {temp_file}: {e}")
        else:
            print(f"Not found (skip): {temp_file}")
    
    # Also try to remove the temporary directory if it's empty
    if temp_gwas_files:
        temp_dir = os.path.dirname(temp_gwas_files[0])
        try:
            if os.path.exists(temp_dir) and not os.listdir(temp_dir):
                os.rmdir(temp_dir)
                print(f"Removed empty temp directory: {temp_dir}")
        except:
            pass  # Directory not empty or other error, ignore
    
    if deleted_count > 0:
        print(f"Successfully deleted {deleted_count} temporary GWAS files")
    else:
        print("No temporary GWAS files found to delete")

def read_vars_file(vars_file):

    """Reads variant IDs from .vars file. Returns chromosome, min_bp, max_bp."""

    variant_ids = []

    with open(vars_file, 'r') as f:

        for line in f:

            v = line.strip()

            if v:

                variant_ids.append(v)

    chrom, positions = None, []

    for vid in variant_ids:

        m = re.match(r'^(\w+):(\d+):[^:]+:[^:]+$', vid)

        if m:

            c, pos = m.group(1), int(m.group(2))

            chrom = c if chrom is None else chrom

            positions.append(pos)

        else:

            print(f'WARNING: Variant ID does not match CHR:POS:REF:ALT: {vid}')

    if not positions or chrom is None:

        print(f'ERROR: Could not parse chromosome and positions from {vars_file}')

        sys.exit(1)

    return chrom, min(positions), max(positions), variant_ids

def detect_columns(gwas_file, verbose=True):

    """Detect column numbers in GWAS file."""

    try:

        df = pd.read_csv(gwas_file, nrows=0, sep='\t')

        headers = df.columns.tolist()

        if verbose:

            print(f"Headers in {os.path.basename(gwas_file)}: {headers}")

    except:

        try:

            df = pd.read_csv(gwas_file, nrows=0, sep=' ')

            headers = df.columns.tolist()

        except:

            print(f"ERROR: Could not read GWAS file {gwas_file}")

            sys.exit(1)

    column_patterns = {

        'chr': ['CHR', 'chr', 'chromosome', 'CHROM', 'chrom'],

        'snp': ['SNP', 'snp', 'rsid', 'RSID', 'ID', 'id', 'variant_id', 'VARIANT_ID', 'MarkerName', 'rs_number', 'MarkerID'],

        'bp': ['BP', 'bp', 'pos', 'POS', 'position', 'POSITION', 'base_pair_location'],

        'a1': ['A1', 'a1', 'effect_allele', 'EFFECT_ALLELE', 'ALT', 'alt', 'ea', 'Allele1'],

        'a2': ['A2', 'a2', 'other_allele', 'OTHER_ALLELE', 'REF', 'ref', 'nea', 'Allele2'],

        'eff': ['BETA', 'beta', 'Beta', 'OR', 'or', 'effect', 'EFFECT', 'b'],

        'se': ['SE', 'se', 'standard_error', 'STANDARD_ERROR', 'stderr', 'STDERR'],

        'pval': ['P', 'p', 'pval', 'PVAL', 'p_value', 'P_VALUE', 'p-value', 'P-VALUE', 'p.value'],

    }

    detected_cols = {}

    for col_type, patterns in column_patterns.items():

        col_num = None

        for i, header in enumerate(headers):

            if header in patterns:

                col_num = i + 1

                break

        if col_num is None:

            defaults = {'chr': 1, 'snp': 2, 'bp': 3, 'a1': 4, 'a2': 5, 'eff': 6, 'se': 7, 'pval': 8}

            col_num = defaults.get(col_type, None)

            if col_type != 'missing' and col_num is None:

                print(f"WARNING: Could not detect {col_type} column")

            elif col_type != 'missing':

                print(f"WARNING: Could not detect {col_type} column, using default {col_num}")

        detected_cols[col_type] = col_num

    if verbose:

        print(f"Detected columns: {detected_cols}")

    return detected_cols

def generate_plink_ref_bim_file(ld_prefix, vars_file, chrom, bp_start, bp_end, plink_bin, ref_panel=None):

    """

    Generate the _ref.bim file that SuSiEx expects using PLINK.

    This creates ld_prefix_ref.bim with the variants specified in the .vars file.

    After generation, removes duplicate variants from the .bim file.

    """

    print(f"[PLINK REF.BIM] Generating _ref.bim file for {ld_prefix}")

    # If no reference panel provided, look for common patterns

    if ref_panel is None:

        possible_refs = [

            ld_prefix, # Same prefix as LD matrix

            f"{ld_prefix}_ref", # Common SuSiEx pattern

            f"{os.path.dirname(ld_prefix)}/reference",

            f"{os.path.dirname(ld_prefix)}/ref"

        ]

        for ref_candidate in possible_refs:

            if os.path.exists(f"{ref_candidate}.bed"):

                ref_panel = ref_candidate

                print(f"[PLINK REF.BIM] Found reference panel: {ref_panel}")

                break

        if ref_panel is None:

            print(f"[PLINK REF.BIM] ERROR: No reference panel found for {ld_prefix}")

            print(f"[PLINK REF.BIM] Please provide --ref_panels argument or ensure reference files exist")

            return False

    # Run PLINK to generate reference subset that matches the LD matrix

    ref_output = f"{ld_prefix}_ref"

    plink_cmd = [

        plink_bin,

        '--bfile', ref_panel,

        '--extract', vars_file,

        '--chr', str(chrom),

        '--from-bp', str(bp_start),

        '--to-bp', str(bp_end),

        '--keep-allele-order',

        '--make-bed',

        '--out', ref_output

    ]

    print(f"[PLINK REF.BIM] Running: {' '.join(plink_cmd)}")

    try:

        result = subprocess.run(plink_cmd, check=True, capture_output=True, text=True)

        print(f"[PLINK REF.BIM] PLINK reference subset generation completed")

        # Check if _ref.bim file was created

        ref_bim_file = f"{ref_output}.bim"

        if os.path.exists(ref_bim_file):

            # Count variants before duplicate removal

            with open(ref_bim_file, 'r') as f:

                n_variants_before = sum(1 for line in f)

            print(f"[PLINK REF.BIM] Generated _ref.bim file: {ref_bim_file} ({n_variants_before} variants)")

            # Remove duplicates from the .bim file

            if remove_duplicates_from_bim_file(ref_bim_file):

                # Count variants after duplicate removal

                with open(ref_bim_file, 'r') as f:

                    n_variants_after = sum(1 for line in f)

                print(f"[PLINK REF.BIM] Final _ref.bim file: {ref_bim_file} ({n_variants_after} variants)")

            else:

                print(f"[PLINK REF.BIM] WARNING: Failed to remove duplicates from {ref_bim_file}")

            # Clean up temporary files, but keep the _ref.bim file

            for temp_file in [f"{ref_output}.log", f"{ref_output}.nosex"]:

                if os.path.exists(temp_file):

                    os.remove(temp_file)

            # Remove .bed and .fam files as SuSiEx only needs .bim

            for temp_file in [f"{ref_output}.bed", f"{ref_output}.fam"]:

                if os.path.exists(temp_file):

                    os.remove(temp_file)

            return ref_bim_file

        else:

            print(f"[PLINK REF.BIM] ERROR: Expected _ref.bim file {ref_bim_file} was not created")

            return False

    except subprocess.CalledProcessError as e:

        print(f"[PLINK REF.BIM] ERROR: PLINK failed with return code {e.returncode}")

        if e.stdout:

            print(f"[PLINK REF.BIM] STDOUT: {e.stdout}")

        if e.stderr:

            print(f"[PLINK REF.BIM] STDERR: {e.stderr}")

        return False

def generate_plink_frequency_file(ld_prefix, vars_file, chrom, bp_start, bp_end, plink_bin, ref_panel=None):

    """

    Generate PLINK-format frequency file using PLINK directly.

    If ref_panel is provided, use it. Otherwise, try to find a reference panel.

    After generation, removes duplicate variants from the .frq file.

    """

    print(f"[PLINK FREQ] Generating frequency file for {ld_prefix}")

    # If no reference panel provided, look for common patterns

    if ref_panel is None:

        # Try common reference panel naming patterns

        possible_refs = [

            ld_prefix, # Same prefix as LD matrix

            f"{ld_prefix}_ref", # Common SuSiEx pattern

            f"{os.path.dirname(ld_prefix)}/reference",

            f"{os.path.dirname(ld_prefix)}/ref"

        ]

        for ref_candidate in possible_refs:

            if os.path.exists(f"{ref_candidate}.bed"):

                ref_panel = ref_candidate

                print(f"[PLINK FREQ] Found reference panel: {ref_panel}")

                break

        if ref_panel is None:

            print(f"[PLINK FREQ] ERROR: No reference panel found for {ld_prefix}")

            print(f"[PLINK FREQ] Please provide --ref_panels argument or ensure reference files exist")

            return False

    # Run PLINK to generate frequency file

    freq_output = f"{ld_prefix}_frq"

    plink_cmd = [

        plink_bin,

        '--bfile', ref_panel,

        '--extract', vars_file,

        '--chr', str(chrom),

        '--from-bp', str(bp_start),

        '--to-bp', str(bp_end),

        '--freq', '--keep-allele-order',

        '--out', freq_output

    ]

    print(f"[PLINK FREQ] Running: {' '.join(plink_cmd)}")

    try:

        result = subprocess.run(plink_cmd, check=True, capture_output=True, text=True)

        print(f"[PLINK FREQ] PLINK frequency calculation completed")

        # Check if frequency file was created

        freq_file = f"{freq_output}.frq"

        if os.path.exists(freq_file):

            # Count variants before duplicate removal

            with open(freq_file, 'r') as f:

                n_variants_before = sum(1 for line in f) - 1 # Subtract header

            print(f"[PLINK FREQ] Generated frequency file: {freq_file} ({n_variants_before} variants)")

            # Remove duplicates from the .frq file

            if remove_duplicates_from_frq_file(freq_file):

                # Count variants after duplicate removal

                with open(freq_file, 'r') as f:

                    n_variants_after = sum(1 for line in f) - 1 # Subtract header

                print(f"[PLINK FREQ] Final frequency file: {freq_file} ({n_variants_after} variants)")

            else:

                print(f"[PLINK FREQ] WARNING: Failed to remove duplicates from {freq_file}")

            # Clean up temporary files

            for temp_file in [f"{freq_output}.log", f"{freq_output}.nosex"]:

                if os.path.exists(temp_file):

                    os.remove(temp_file)

            return freq_file

        else:

            print(f"[PLINK FREQ] ERROR: Expected frequency file {freq_file} was not created")

            return False

    except subprocess.CalledProcessError as e:

        print(f"[PLINK FREQ] ERROR: PLINK failed with return code {e.returncode}")

        if e.stdout:

            print(f"[PLINK FREQ] STDOUT: {e.stdout}")

        if e.stderr:

            print(f"[PLINK FREQ] STDERR: {e.stderr}")

        return False

def extract_population_names(gwas_files):

    """Extract population names from GWAS file paths."""

    pop_names = []

    for gwas_file in gwas_files:

        # Extract directory structure to find population name

        # Example: "/path/to/AFR/file.txt" -> "AFR"

        path_parts = gwas_file.split('/')

        for part in path_parts:

            # Look for common population codes

            if part.upper() in ['AFR', 'EUR', 'EAS', 'SAS', 'AMR', 'ASN', 'HISPANIC', 'LATINO'] or (len(part) == 3 and part.isupper()):

                pop_names.append(part.upper())

                break

        else:

            # If no clear population name found, use directory name containing the file

            dir_name = os.path.basename(os.path.dirname(gwas_file))

            pop_names.append(dir_name)

    return pop_names

def calculate_z_scores(gwas_files, all_detected_cols):

    """Calculate Z-scores for each population from GWAS files."""

    z_scores_by_pop = {}

    for i, gwas_file in enumerate(gwas_files):

        print(f"[Z-SCORE] Calculating Z-scores for population {i+1}: {gwas_file}")

        try:

            df = pd.read_csv(gwas_file, sep='\t')

        except:

            try:

                df = pd.read_csv(gwas_file, sep=' ')

            except Exception as e:

                print(f"[Z-SCORE] Error reading GWAS file {gwas_file}: {e}")

                continue

        detected = all_detected_cols[i]

        # Get column names (convert from 1-based to 0-based indexing)

        snp_col = detected['snp'] - 1

        eff_col = detected['eff'] - 1

        se_col = detected['se'] - 1

        # Calculate Z-score = effect / SE

        df['Z_SCORE'] = df.iloc[:, eff_col] / df.iloc[:, se_col]

        # Create variant ID key using SNP column

        variant_key = df.iloc[:, snp_col]

        # Store Z-scores with variant ID as key

        z_scores_by_pop[i] = dict(zip(variant_key, df['Z_SCORE']))

        print(f"[Z-SCORE] Calculated {len(z_scores_by_pop[i])} Z-scores for population {i+1}")

    return z_scores_by_pop

def read_native_credible_sets(cs_file):

    """Read SuSiEx's own 95% credible sets from its .cs file.

    Returns {variant_id: [cs_id, ...]} with the CS IDs sorted. SuSiEx writes 'NULL' when the

    model converged without a credible set and 'FAIL' when it did not converge (the .snp file

    then has no PIP(CS) columns either); both give no credible sets. SuSiEx credible sets can

    overlap, so a variant can belong to several of them."""

    with open(cs_file) as f:

        header = f.readline().strip()

    if not header.startswith('CS_ID'):

        print(f"[SUSIE FORMAT] SuSiEx .cs file reports '{header}'; no credible sets assigned")

        return {}

    cs_df = pd.read_csv(cs_file, sep='\t', usecols=['CS_ID', 'SNP'])

    return cs_df.groupby('SNP')['CS_ID'].apply(lambda ids: sorted(ids.astype(int))).to_dict()

def create_susie_format_output(output_file, gwas_files, all_detected_cols):

    """Create final output file in SuSiE format combining SuSiEx results with summary statistics."""

    print("=== Creating SuSiE-format output file ===")

    # Determine output file paths

    output_base = output_file.replace('.txt', '').replace('.susiex', '')

    susiex_snp_file = f"{output_base}.snp"

    susiex_cs_file = f"{output_base}.cs"

    susiex_summary_file = f"{output_base}.summary"

    final_output_file = f"{output_base}.txt"

    # Check if SuSiEx output files exist

    required_files = [susiex_snp_file, susiex_cs_file, susiex_summary_file]

    for req_file in required_files:

        if not os.path.exists(req_file):

            print(f"[SUSIE FORMAT] ERROR: Required SuSiEx output file not found: {req_file}")

            return False

    # Extract population names

    pop_names = extract_population_names(gwas_files)

    print(f"[SUSIE FORMAT] Population names: {pop_names}")

    # Calculate Z-scores from GWAS files

    z_scores_by_pop = calculate_z_scores(gwas_files, all_detected_cols)

    # Read SuSiEx .snp file

    print(f"[SUSIE FORMAT] Reading {susiex_snp_file}")

    try:

        snp_df = pd.read_csv(susiex_snp_file, sep='\t')

        print(f"[SUSIE FORMAT] Loaded {len(snp_df)} variants from .snp file")

    except Exception as e:

        print(f"[SUSIE FORMAT] Error reading .snp file: {e}")

        return False

    # Identify CS columns and LogBF columns

    pip_cols = [col for col in snp_df.columns if col.startswith('PIP(CS')]

    logbf_cols = [col for col in snp_df.columns if col.startswith('LogBF(CS')]

    print(f"[SUSIE FORMAT] Found {len(pip_cols)} PIP columns and {len(logbf_cols)} LogBF columns")

    # Calculate overall PIP (sum of CS-specific PIPs)

    if pip_cols:

        snp_df['PIP'] = snp_df[pip_cols].sum(axis=1)

    else:

        snp_df['PIP'] = 0.0

    # Create output dataframe

    output_df = pd.DataFrame()

    output_df['VARIANT_ID'] = snp_df['SNP']

    # Add Z-scores for each population (ALWAYS separate columns regardless of number of populations)

    n_pops = len(gwas_files)

    for i in range(n_pops):

        pop_name = pop_names[i] if i < len(pop_names) else f"Pop{i+1}"

        z_col_name = f"Z_{pop_name}"

        z_col_data = []

        for variant in output_df['VARIANT_ID']:

            z_score = z_scores_by_pop.get(i, {}).get(variant, 'NA')

            z_col_data.append(z_score)

        output_df[z_col_name] = z_col_data

    # Add PIP column

    output_df['PIP'] = snp_df['PIP']

    # Add converged column (all NA as specified)

    output_df['converged'] = 'NA'

    # Determine credible set membership from SuSiEx's own 95% credible sets (.cs file)

    # cs = lowest-numbered credible set containing the variant (0 = none)

    # cs_all = all credible sets containing the variant, comma-separated ('0' = none)

    native_cs = read_native_credible_sets(susiex_cs_file)

    memberships = output_df['VARIANT_ID'].map(native_cs)

    output_df['cs'] = [m[0] if isinstance(m, list) else 0 for m in memberships]

    output_df['cs_all'] = [','.join(map(str, m)) if isinstance(m, list) else '0' for m in memberships]

    # Add LogBF columns with population names

    for logbf_col in logbf_cols:

        # Parse column name like 'LogBF(CS1,Pop1)' to 'LogBF_CS1_PopName'

        match = re.match(r'LogBF\(CS(\d+),Pop(\d+)\)', logbf_col)

        if match:

            cs_num = match.group(1)

            pop_idx = int(match.group(2)) - 1 # Convert to 0-based

            pop_name = pop_names[pop_idx] if pop_idx < len(pop_names) else f"Pop{pop_idx+1}"

            new_col_name = f"LogBF_CS{cs_num}_{pop_name}"

            output_df[new_col_name] = snp_df[logbf_col]

    # Write output file

    print(f"[SUSIE FORMAT] Writing output to: {final_output_file}")

    try:

        output_df.to_csv(final_output_file, sep='\t', index=False, na_rep='NA')

        print(f"[SUSIE FORMAT] Successfully created output file with {len(output_df)} variants")

        print(f"[SUSIE FORMAT] Output columns: {list(output_df.columns)}")

        return final_output_file

    except Exception as e:

        print(f"[SUSIE FORMAT] Error writing output file: {e}")

        return False

def validate_inputs(args):

    """Validate input files and arguments."""

    print("=== Validating inputs ===")

    if not os.path.isfile(args.susiex_bin):

        print(f"ERROR: SuSiEx executable not found: {args.susiex_bin}")

        sys.exit(1)

    if not os.path.isfile(args.plink_bin):

        print(f"ERROR: PLINK executable not found: {args.plink_bin}")

        sys.exit(1)

    gwas_files = [f.strip() for f in args.gwas_files.split(',')]

    ld_files = [f.strip() for f in args.ld_files.split(',')]

    n_gwas_list = [int(n.strip()) for n in args.n_gwas.split(',')]

    ref_panels = None

    if args.ref_panels:

        ref_panels = [f.strip() for f in args.ref_panels.split(',')]

        if len(ref_panels) != len(ld_files):

            print("ERROR: Number of reference panels must match number of LD files")

            sys.exit(1)

    n_pops = len(gwas_files)

    if len(ld_files) != n_pops or len(n_gwas_list) != n_pops:

        print("ERROR: Number of populations must match across all inputs")

        sys.exit(1)

    for i, gwas_file in enumerate(gwas_files):

        if not os.path.isfile(gwas_file):

            print(f"ERROR: GWAS file not found: {gwas_file}")

            sys.exit(1)

        print(f"Population {i+1}: GWAS file {gwas_file}")

    for i, ld_prefix in enumerate(ld_files):

        vars_file = f"{ld_prefix}.vars"

        if not os.path.isfile(vars_file):

            print(f"ERROR: .vars file not found: {vars_file}")

            sys.exit(1)

        print(f"Population {i+1}: LD matrix {ld_prefix}, .vars ✓")

    output_dir = os.path.dirname(args.output_file)

    if output_dir:

        os.makedirs(output_dir, exist_ok=True)

    print(f"Validated input files for {n_pops} populations.")

    return gwas_files, ld_files, n_gwas_list, n_pops, ref_panels

def detect_all_columns(gwas_files):

    """Detect columns for all GWAS files."""

    print("=== Detecting column numbers ===")

    all_cols = {'chr': [], 'snp': [], 'bp': [], 'a1': [], 'a2': [], 'eff': [], 'se': [], 'pval': []}

    all_detected_cols = []

    for i, gwas_file in enumerate(gwas_files):

        print(f"Population {i+1}:")

        detected_cols = detect_columns(gwas_file)

        all_detected_cols.append(detected_cols)

        for col_type in all_cols.keys():

            all_cols[col_type].append(str(detected_cols[col_type]))

    col_strings = {k: ','.join(v) for k, v in all_cols.items()}

    return col_strings, all_detected_cols

def get_chr_and_bp_range(ld_files):

    """Get chromosome and BP range from .vars files."""

    chrs, min_pos, max_pos = [], [], []

    all_variant_ids = []

    for ld_prefix in ld_files:

        vars_file = f"{ld_prefix}.vars"

        c, mn, mx, variant_ids = read_vars_file(vars_file)

        chrs.append(c)

        min_pos.append(mn)

        max_pos.append(mx)

        all_variant_ids.append(variant_ids)

    if len(set(chrs)) > 1:

        print(f"ERROR: Chromosomes differ across populations: {chrs}")

        sys.exit(1)

    return chrs[0], str(min(min_pos)), str(max(max_pos)), all_variant_ids

def process_ld_matrices(ld_files):

    """Convert text LD matrices to bin4 format if needed."""

    print("=== Processing LD matrices ===")

    bin4_created = []

    for ld_prefix in ld_files:

        ld_bin = f"{ld_prefix}.ld.bin"

        ld_txt = f"{ld_prefix}"

        if not os.path.exists(ld_bin):

            if os.path.exists(ld_txt):

                print(f"Converting {ld_txt} to {ld_bin}")

                if convert_text_ld_to_bin4(ld_txt, ld_prefix):

                    bin4_created.append(ld_bin)

                else:

                    print(f"ERROR: Failed to convert {ld_txt} to bin4 format")

                    sys.exit(1)

            else:

                print(f"ERROR: Missing both {ld_bin} and {ld_txt} for LD matrix prefix {ld_prefix}")

                sys.exit(1)

        else:

            print(f"Found existing bin4 LD matrix: {ld_bin}")

    return bin4_created

def generate_all_frequency_and_ref_files(ld_files, chrom, bp_start, bp_end, plink_bin, ref_panels):

    """Generate PLINK frequency files and _ref.bim files for all populations."""

    print("=== Generating frequency files and _ref.bim files ===")

    frq_files = []

    ref_bim_files = []

    for i, ld_prefix in enumerate(ld_files):

        vars_file = f"{ld_prefix}.vars"

        ref_panel = ref_panels[i] if ref_panels else None

        print(f"Population {i+1}:")

        # Generate frequency file

        frq_file = generate_plink_frequency_file(

            ld_prefix, vars_file, chrom, bp_start, bp_end, plink_bin, ref_panel

        )

        if frq_file:

            frq_files.append(frq_file)

        else:

            print(f"ERROR: Failed to generate frequency file for {ld_prefix}")

            sys.exit(1)

        # Generate _ref.bim file

        ref_bim_file = generate_plink_ref_bim_file(

            ld_prefix, vars_file, chrom, bp_start, bp_end, plink_bin, ref_panel

        )

        if ref_bim_file:

            ref_bim_files.append(ref_bim_file)

        else:

            print(f"ERROR: Failed to generate _ref.bim file for {ld_prefix}")

            sys.exit(1)

    return frq_files, ref_bim_files

def run_susiex(args, temp_gwas_files, ld_files, n_gwas_list, col_strings, chrom, bp_start, bp_end):
    """Run SuSiEx fine-mapping with filtered GWAS files."""
    print("=== Running SuSiEx fine-mapping ===")
    
    output_dir = os.path.dirname(args.output_file)
    output_name = os.path.basename(args.output_file).replace('.susiex', '').replace('.txt', '')
    bp_range = f"{bp_start},{bp_end}"
    
    # Convert temp_gwas_files list to comma-separated string
    temp_gwas_files_str = ','.join(temp_gwas_files)
    
    cmd = [
        args.susiex_bin,
        '--sst_file', temp_gwas_files_str,  # Use filtered files instead of original
        '--n_gwas', args.n_gwas,
        '--ld_file', args.ld_files,
        '--out_dir', output_dir if output_dir else '.',
        '--out_name', output_name,
        '--chr', chrom,
        '--bp', bp_range,
        '--chr_col', col_strings['chr'],
        '--snp_col', col_strings['snp'],
        '--bp_col', col_strings['bp'],
        '--a1_col', col_strings['a1'],
        '--a2_col', col_strings['a2'],
        '--eff_col', col_strings['eff'],
        '--se_col', col_strings['se'],
        '--pval_col', col_strings['pval'],
        '--maf', str(args.maf),
        '--n_sig', str(args.n_sig),
        '--level', str(args.level),
        '--min_purity', str(args.min_purity),
        '--pval_thresh', str(args.pval_thresh),
        '--max_iter', str(args.max_iter),
        '--tol', str(args.tol),
        '--threads', str(args.threads),
        '--plink', str(args.plink_bin)
    ]
    
    if args.keep_ambig:
        cmd.append('--keep-ambig=True')
    if args.mult_step:
        cmd.append('--mult-step=True')
    
    print(' '.join(cmd))
    
    try:
        # Time the SuSiEx execution
        start_time = time.time()
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        end_time = time.time()
        susiex_execution_time = end_time - start_time
        
        print("SuSiEx completed successfully!")
        print(f"SuSiEx execution time: {susiex_execution_time:.2f} seconds")
        
        # Create timing directory and file
        output_dir = os.path.dirname(args.output_file) if args.output_file else '.'
        timing_dir = os.path.join(output_dir, 'timing')
        os.makedirs(timing_dir, exist_ok=True)
        
        # Create timing filename based on output file
        output_basename = os.path.basename(args.output_file)
        timing_file = os.path.join(timing_dir, f"{output_basename}")
        
        # Write timing file
        with open(timing_file, 'w') as f:
            f.write("susiex_time_seconds\n")
            f.write(f"{susiex_execution_time}\n")
        
        print(f"Timing information written to: {timing_file}")
        if result.stdout:
            print("STDOUT:", result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: SuSiEx failed with return code {e.returncode}")
        if e.stdout:
            print("STDOUT:", e.stdout)
        if e.stderr:
            print("STDERR:", e.stderr)
        sys.exit(1)

def cleanup_frequency_files(frq_files):

    """Delete frequency files after SuSiEx completion."""

    print("=== Cleaning up frequency files ===")

    for frq_file in frq_files:

        if os.path.exists(frq_file):

            os.remove(frq_file)

            print(f"Deleted: {frq_file}")

def main():

    print("=== SuSiEx Multi-Population Fine-Mapping (PLINK frequency + _ref.bim generation + GWAS filtering) ===")

    args = parse_arguments()

    # Validate inputs

    gwas_files, ld_files, n_gwas_list, n_pops, ref_panels = validate_inputs(args)

    # Detect columns

    col_strings, all_detected_cols = detect_all_columns(gwas_files)

    # Get chromosome and BP range

    chrom, bp_start, bp_end, all_variant_ids = get_chr_and_bp_range(ld_files)

    print(f"Chromosome: {chrom}, BP range: {bp_start}-{bp_end}")

    # Process LD matrices (convert text to bin4 if needed)

    bin4_created = process_ld_matrices(ld_files)

    # Initialize variables for cleanup in finally block
    temp_gwas_files = []
    frq_files = []
    ref_bim_files = []

    try:

        # Create filtered and deduplicated GWAS files
        temp_gwas_files = create_filtered_gwas_files(gwas_files, all_variant_ids)

        # Generate frequency files and _ref.bim files using PLINK

        frq_files, ref_bim_files = generate_all_frequency_and_ref_files(

            ld_files, chrom, bp_start, bp_end, args.plink_bin, ref_panels

        )

        # Run SuSiEx with filtered GWAS files

        run_susiex(args, temp_gwas_files, ld_files, n_gwas_list, col_strings, chrom, bp_start, bp_end)

        # CREATE SUSIE FORMAT OUTPUT (using original GWAS files for Z-score calculation)

        susie_output = create_susie_format_output(args.output_file, gwas_files, all_detected_cols)

        # DELETE SUSIEX OUTPUT FILES (only if SuSiE format was created successfully)

        if susie_output:

            delete_susiex_output_files(args.output_file)

    finally:

        # Clean up temporary GWAS files
        if temp_gwas_files:
            cleanup_temp_gwas_files(temp_gwas_files)

        # Clean up frequency files

        if frq_files:

            cleanup_frequency_files(frq_files)

        # Clean up _ref.bim files (commented out for safety)

        if ref_bim_files:

            delete_ref_bim_files(ref_bim_files)

        # Clean up bin4 files (commented out for safety)

        if bin4_created:

            delete_bin4_files(bin4_created)

if __name__ == '__main__':

    main()