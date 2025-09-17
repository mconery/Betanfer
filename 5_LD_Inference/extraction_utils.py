"""
extraction_utils.py

Utility functions for reading and processing GWAS summary statistics.
"""

import pandas as pd
import numpy as np
import os
import logging
from typing import List, Dict, Set, Tuple, Optional
import concurrent.futures
import pysam
import re
from pathlib import Path
import shutil
import torch
import random


# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_variant_list(file_path: str, variant_col: str = 'variant_id', sep: str = '\t') -> Set[str]:
    """
    Load the common set of variants to extract from a file.
    
    Args:
        file_path: Path to file containing variant IDs
        variant_col: Name of the column containing variant IDs
        sep: Column separator (default: tab)
        
    Returns:
        Set of variant IDs to extract
    """
    try:
        df = pd.read_csv(file_path, sep=sep)
        
        if variant_col not in df.columns:
            raise ValueError(f"ERROR: Column '{variant_col}' not found in {file_path}")
        
        variant_set = set(df[variant_col].astype(str))
        logger.info(f"NOTE: Loaded {len(variant_set)} variants from {file_path}")
        return variant_set
    except Exception as e:
        logger.error(f"ERROR: Error loading variant list from {file_path}: {e}")
        raise

def make_var_dict_for_positional_ids(
    variants_to_extract: Set[str],
    id_conv: str = "CHR:POS:REF:ALT"
) -> Dict[str, Tuple[str, str]]:
    """
    Creates a dictionary relating variant IDs to positions, inferring field positions
    for chromosome and position from the id_conv string.
    The function allows for multiple possible delimiters (":", "_", "-", etc.)
    and uses the same parsing logic for `variants_to_extract`.
    Args:
        variants_to_extract: Set of variant IDs we ultimately want to extract (required)
        id_conv: Naming convention for variant IDs, e.g. "CHR:POS:REF:ALT".
                 Must contain "CHR" and "POS". Delimiters can be any non-alphanumeric char.
    Returns:
        dict: Keys = variant IDs (str), Values = (CHR, POS) tuples (both as strings)
    """
    
    if "CHR" not in id_conv or "POS" not in id_conv:
        raise ValueError("id_conv must contain 'CHR' and 'POS' fields.")
    
    # Identify delimiter(s) — simply split using any non-alphanumeric
    conv_parts = re.split(r'\W+', id_conv)
    
    # Map the token name to its index
    field_to_index = {field: idx for idx, field in enumerate(conv_parts)}
    
    if "CHR" not in field_to_index or "POS" not in field_to_index:
        raise ValueError("CHR and POS placeholders not found in parsed id_conv.")
    
    chr_index = field_to_index["CHR"]
    pos_index = field_to_index["POS"]
    
    var_dict = {}
    for var_id in variants_to_extract:
        parts = re.split(r'\W+', var_id)
        try:
            chr_val = parts[chr_index]
            pos_val = parts[pos_index]
        except IndexError:
            raise ValueError(f"Variant ID '{var_id}' does not match the format: {id_conv}")
        var_dict[var_id] = (chr_val, pos_val)
    return var_dict

def extract_betas_batch(
    file_list: List[str],
    variant_col: str,
    beta_col: str,
    variants_to_extract: Set[str],
    variant_map: Dict[str, tuple],
    chunk_size: int = 100000
) -> pd.DataFrame:
    """
    Extracts beta values using tabix-indexed files for fast genomic position lookup
    
    Args:
        file_list: List of bgzipped, tabix-indexed file paths
        variant_col: Column name for variant identifiers
        beta_col: Column name for beta values
        variants_to_extract: Set of variant IDs to extract (required)
        variant_map: Dictionary mapping variant ID to (chromosome, position)
        chunk_size: Not used (retained for backward compatibility)
        
    Returns:
        pd.DataFrame: Variants as rows, files as columns
    """
    variants_list = sorted(variants_to_extract)
    results = []
    
    def process_single_file(file_path: str) -> pd.Series:
        """Process a single tabix-indexed file"""
        betas = pd.Series(np.nan, index=variants_list, name=file_path)
        try:
            with pysam.TabixFile(file_path) as tabix_file:
                header = tabix_file.header[0].strip().split('\t')
                vid_idx = header.index(variant_col)
                beta_idx = header.index(beta_col)
                
                for variant in variants_list:
                    chrom, pos = variant_map.get(variant, (None, None))
                    if not chrom or not pos:
                        continue
                    
                    try:
                        records = tabix_file.fetch(reference=chrom, start=int(pos)-1, end=int(pos)+1)
                        for record in records:
                            fields = record.split('\t')
                            if fields[vid_idx] == variant:
                                betas[variant] = float(fields[beta_idx])
                                break
                    except ValueError:
                        continue
        except Exception as e:
            logging.error(f"Error processing {file_path}: {str(e)}")
        return betas
    
    num_threads = min(len(file_list), os.cpu_count() or 1)
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_single_file, file) for file in file_list]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    
    return pd.concat(results, axis=1)

def bgzip_and_tbi_sumstats_with_label_cols(
    input_path: str,
    output_dir: str,
    chrom_col_label: str,
    pos_col_label: str,
    zero_based: bool = False,
    force: bool = True,
) -> tuple[str, str]:
    """
    BGZF-compress a GWAS summary stats file with a header line and tabix-index it,
    using column names (labels) to identify chromosome and position columns.
    
    Assumes the first line contains column labels, no other header lines, and
    no end column. The file must be sorted by chrom and position.
    
    Parameters
    ----------
    input_path : str
        Path to the input tab-delimited file.
    output_dir : str
        Directory to write the BGZF file and the .tbi index.
    chrom_col_label : str
        Column label for chromosome.
    pos_col_label : str
        Column label for position.
    zero_based : bool
        Whether the position column is zero-based (typically False for GWAS).
    force : bool
        Whether to overwrite outputs if they exist.
    Returns
    -------
    (bgzf_path, tbi_path) : tuple[str, str]
        Paths to BGZF-compressed file and .tbi index.
    """
    input_path = str(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Read the header line to find column indices (1-based)
    with open(input_path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
    try:
        chrom_col = header.index(chrom_col_label) + 1
    except ValueError:
        raise ValueError(f"Chromosome label '{chrom_col_label}' not found in header columns.")
    try:
        pos_col = header.index(pos_col_label) + 1
    except ValueError:
        raise ValueError(f"Position label '{pos_col_label}' not found in header columns.")
    # File naming for output
    in_name = Path(input_path).name
    out_name = in_name + ".gz" if not in_name.endswith(".gz") else in_name
    bgzf_path = output_dir / out_name
    tbi_path = str(bgzf_path) + ".tbi"
    if bgzf_path.exists():
        if force:
            bgzf_path.unlink()
        else:
            raise FileExistsError(f"Output file already exists: {bgzf_path}")
    # Write BGZF compressed output preserving all lines including header
    with open(input_path, "rb") as fin, pysam.BGZFile(bgzf_path, "w") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)
    if force and os.path.exists(tbi_path):
        os.remove(tbi_path)
    # Index, skipping the first line (header)
    pysam.tabix_index(
        str(bgzf_path),
        seq_col=chrom_col - 1,
        start_col=pos_col - 1,
        end_col=pos_col - 1,
        line_skip=1,
        zerobased=zero_based,
    )
    return str(bgzf_path), tbi_path

def merge_bim_files(bim_file_paths, output_path=None):
    """
    Merge multiple PLINK .bim files into a dataframe with the union of variants.
    
    Parameters:
    - bim_file_paths: list of str, file paths to .bim files
    - output_path: str or None, if provided the merged .bim file will be saved here
    
    Returns:
    - pandas.DataFrame with columns ["CHR", "SNP", "CM", "POS", "A1", "A2"]
    """
    
    column_names = ["CHR", "SNP", "CM", "POS", "A1", "A2"]
    dataframes = []
    for file_path in bim_file_paths:
        df = pd.read_csv(file_path, sep="\t", header=None, names=column_names, dtype=str)
        dataframes.append(df)
    # Concatenate and drop duplicates
    merged_df = pd.concat(dataframes, ignore_index=True).drop_duplicates()
    # Ensure order by chromosome and position if possible
    merged_df["POS"] = pd.to_numeric(merged_df["POS"], errors="coerce")
    merged_df["CHR"] = pd.to_numeric(merged_df["CHR"], errors="coerce")
    merged_df = merged_df.sort_values(by=["CHR", "POS"]).reset_index(drop=True)
    # Save to file if requested
    if output_path:
        merged_df.to_csv(output_path, sep="\t", header=False, index=False)
    return merged_df

def index_gwama_file(args):
    """
    Wrapper for bgzip_and_tbi_sumstats_with_bim_merge so it can be called in parallel.
    Expects a tuple:
      (gwama_file, out_dir, merged_bim, variant_col_gwama, chrom_col_label, pos_col_label)
    """
    gwama_file, out_dir, merged_bim, variant_col_gwama, chrom_col_label, pos_col_label = args
    return bgzip_and_tbi_sumstats_with_bim_merge(
        input_path=gwama_file,
        output_dir=out_dir,
        bim_df=merged_bim,
        variant_col_gwama=variant_col_gwama,
        chrom_col_label=chrom_col_label,
        pos_col_label=pos_col_label,
        zero_based=False,
        force=True
    )

def index_gwama_file_with_variant_ids(args):
    """
    Parallel wrapper for bgzip_and_tbi_sumstats_with_variant_ids.
    Expects: (input_path, output_dir, variant_col, chrom_col_label, pos_col_label, variant_id_format)
    """
    input_path, output_dir, variant_col, chrom_col_label, pos_col_label, variant_id_format = args
    return bgzip_and_tbi_sumstats_with_variant_ids(
        input_path=input_path,
        output_dir=output_dir,
        variant_col=variant_col,
        chrom_col_label=chrom_col_label,
        pos_col_label=pos_col_label,
        variant_id_format=variant_id_format,
        zero_based=False,
        force=True
    )

def bgzip_and_tbi_sumstats_with_bim_merge(
    input_path: str,
    output_dir: str,
    bim_df: pd.DataFrame,
    variant_col_gwama: str,
    chrom_col_label: str = "CHR",
    pos_col_label: str = "POS",
    zero_based: bool = False,
    force: bool = True,
) -> tuple[str, str]:
    """
    HEADER-PRESERVING VERSION: For GWAMA (or similar) summary stats missing chrom/pos,
    merge in chrom/pos from corresponding .bim file, then bgzip and tabix-index.
    
    KEY FIX: Adds '#' prefix to header so tabix will preserve it properly.
    """
    input_path = str(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Use SNP as merge key
    bim_df = bim_df.loc[:, ["CHR", "POS", "SNP"]].rename(columns={"SNP": variant_col_gwama})
    
    # Ensure proper data types for bim_df
    bim_df["CHR"] = bim_df["CHR"].astype(str)
    bim_df["POS"] = pd.to_numeric(bim_df["POS"], errors="coerce")
    
    # Load GWAMA file, infer delimiter
    with open(input_path) as f:
        first_line = f.readline()
    delim = "\t" if ("\t" in first_line) else r"\s+"
    gwama_df = pd.read_csv(
        input_path,
        sep=delim,
        dtype={variant_col_gwama: str}
    )
    
    # Merge: inner keeps only matching variants
    merged = pd.merge(bim_df, gwama_df, on=variant_col_gwama, how="inner")
    
    # Rename columns for output
    merged = merged.rename(columns={"CHR": chrom_col_label, "POS": pos_col_label})
    
    # Ensure proper sorting and data types
    merged[chrom_col_label] = merged[chrom_col_label].astype(str)
    merged[pos_col_label] = pd.to_numeric(merged[pos_col_label], errors="coerce")
    
    # Remove any rows with missing chromosome or position
    merged = merged.dropna(subset=[chrom_col_label, pos_col_label])
    
    # Sort by chromosome and position (ESSENTIAL for tabix)
    def sort_key(chrom):
        try:
            return (0, int(chrom))
        except ValueError:
            return (1, chrom)
    
    merged['_sort_key'] = merged[chrom_col_label].apply(sort_key)
    merged = merged.sort_values(['_sort_key', pos_col_label]).drop('_sort_key', axis=1)
    merged = merged.reset_index(drop=True)
    
    # Convert position back to int
    merged[pos_col_label] = merged[pos_col_label].astype(int)
    
    # Prepare output file
    in_name = Path(input_path).name
    out_name = in_name + ".merged" if ".merged" not in in_name else in_name
    merged_txt = output_dir / out_name
    
    # CRITICAL FIX: Write with '#' prefixed header for tabix compatibility
    with open(merged_txt, 'w') as f:
        # Write header with '#' prefix so tabix preserves it
        header_line = '#' + '\t'.join(merged.columns) + '\n'
        f.write(header_line)
        # Write data
        merged.to_csv(f, sep='\t', index=False, header=False)
    
    # Get header cols for indexing (without the '#' prefix)
    header_cols = list(merged.columns)
    
    # Clean-up memory
    del gwama_df, merged
    
    # Prepare BGZF and TBI output filenames
    bgzf_path = output_dir / (out_name + ".gz")
    tbi_path = str(bgzf_path) + ".tbi"
    
    # Clean up existing output if needed
    if bgzf_path.exists():
        if force:
            bgzf_path.unlink()
        else:
            raise FileExistsError(f"Output file already exists: {bgzf_path}")
    
    if force and os.path.exists(tbi_path):
        os.remove(tbi_path)
    
    # Write BGZF compressed output
    with open(merged_txt, "rb") as fin, pysam.BGZFile(bgzf_path, "w") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)
    
    # Delete the intermediate merged_txt file
    if os.path.exists(merged_txt):
        os.remove(merged_txt)
    
    # Prepare for indexing - find column indices
    try:
        chrom_col = header_cols.index(chrom_col_label)
    except ValueError:
        raise ValueError(f"Chromosome label '{chrom_col_label}' not found in merged file.")
    
    try:
        pos_col = header_cols.index(pos_col_label)
    except ValueError:
        raise ValueError(f"Position label '{pos_col_label}' not found in merged file.")
    
    # Create tabix index with meta character '#' specified
    try:
        pysam.tabix_index(
            str(bgzf_path),
            seq_col=chrom_col,      # 0-based column index
            start_col=pos_col,      # 0-based column index
            end_col=pos_col,        # Same as start for point data
            line_skip=0,            # Don't skip any lines now (# lines are preserved)
            zerobased=zero_based,   # Whether coordinates are 0-based
            meta_char='#',          # Specify that '#' lines are metadata/header
        )
        
        # Verify the index was created successfully
        if not os.path.exists(tbi_path):
            raise RuntimeError(f"Tabix index file was not created: {tbi_path}")
            
        # Test that we can read the header (should work now!)
        try:
            with pysam.TabixFile(str(bgzf_path)) as test_tabix:
                header_lines = list(test_tabix.header)
                if not header_lines:
                    logging.warning("Warning: Tabix file header is still empty")
                else:
                    logging.info(f"Success: Found {len(header_lines)} header line(s)")
        except Exception as e:
            raise RuntimeError(f"Created tabix file cannot be read properly: {e}")
            
    except Exception as e:
        # Clean up on failure
        if os.path.exists(str(bgzf_path)):
            os.remove(str(bgzf_path))
        if os.path.exists(tbi_path):
            os.remove(tbi_path)
        raise RuntimeError(f"Failed to create tabix index: {e}")
    
    return str(bgzf_path), tbi_path

def bgzip_and_tbi_sumstats_with_variant_ids(
    input_path: str,
    output_dir: str,
    variant_col: str,
    chrom_col_label: str = "CHR",
    pos_col_label: str = "POS",
    variant_id_format: str = "CHR:POS:REF:ALT",
    zero_based: bool = False,
    force: bool = True,
) -> tuple[str, str]:
    """
    For GWAMA (or similar) summary stats: extract chromosome and position
    from variant IDs, then bgzip and tabix-index. Does NOT require a .bim file.

    Args:
      input_path: GWAS/summary file path (tsv or space-delimited)
      output_dir: Where to store output files
      variant_col: Column name in input file with variant IDs
      chrom_col_label: Output chromosome column name
      pos_col_label: Output position column name
      variant_id_format: Format string showing fields in variant ID (e.g., "CHR:POS:REF:ALT")
      zero_based: Set True if positions are 0-based
      force: If True, overwrite outputs

    Returns:
      (bgzf_path, tbi_path): output file paths
    """
    input_path = str(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Identify the field indices
    conv_parts = re.split(r'\W+', variant_id_format)
    field_to_index = {field: idx for idx, field in enumerate(conv_parts)}
    if "CHR" not in field_to_index or "POS" not in field_to_index:
        raise ValueError("variant_id_format must contain 'CHR' and 'POS' fields.")
    chr_index = field_to_index["CHR"]
    pos_index = field_to_index["POS"]
    # read file, infer delimiter
    with open(input_path) as f:
        first_line = f.readline()
    delim = "\t" if ("\t" in first_line) else r"\s+"
    df = pd.read_csv(input_path, sep=delim, dtype={variant_col: str})
    # Parse each variant ID
    def parse_variant_id(var_id):
        if pd.isna(var_id) or var_id == '':
            return pd.NA, pd.NA
        try:
            parts = re.split(r'\W+', str(var_id))
            chr_val = parts[chr_index]
            pos_val = int(parts[pos_index])
            return chr_val, pos_val
        except (IndexError, ValueError):
            return pd.NA, pd.NA
    parsed = df[variant_col].apply(parse_variant_id)
    df[chrom_col_label] = [c[0] for c in parsed]
    df[pos_col_label] = [c[1] for c in parsed]
    df = df.dropna(subset=[chrom_col_label, pos_col_label])
    df[chrom_col_label] = df[chrom_col_label].astype(str)
    df[pos_col_label] = df[pos_col_label].astype(int)
    # Sort by chromosome and position
    def sort_key(chrom):
        try:
            return (0, int(chrom))
        except ValueError:
            return (1, chrom)
    df['_sort_key'] = df[chrom_col_label].apply(sort_key)
    df = df.sort_values(['_sort_key', pos_col_label]).drop('_sort_key', axis=1).reset_index(drop=True)
    # Output header: add prepended "#" for tabix compatibility
    in_name = Path(input_path).name
    out_name = in_name
    merged_txt = output_dir / out_name
    with open(merged_txt, 'w') as f:
        f.write('#' + '\t'.join(df.columns) + '\n')
        df.to_csv(f, sep='\t', index=False, header=False)
    header_cols = list(df.columns)
    del df
    bgzf_path = output_dir / (out_name + ".gz")
    tbi_path = str(bgzf_path) + ".tbi"
    if bgzf_path.exists():
        if force:
            bgzf_path.unlink()
        else:
            raise FileExistsError(f"Output file already exists: {bgzf_path}")
    if force and os.path.exists(tbi_path):
        os.remove(tbi_path)
    with open(merged_txt, "rb") as fin, pysam.BGZFile(bgzf_path, "w") as fout:
        shutil.copyfileobj(fin, fout, length=1024*1024)
    if os.path.exists(merged_txt):
        os.remove(merged_txt)
    chrom_col = header_cols.index(chrom_col_label)
    pos_col = header_cols.index(pos_col_label)
    pysam.tabix_index(
        str(bgzf_path),
        seq_col=chrom_col,
        start_col=pos_col,
        end_col=pos_col,
        line_skip=0,
        zerobased=zero_based,
        meta_char='#'
    )
    if not os.path.exists(tbi_path):
        raise RuntimeError(f"Tabix index file was not created: {tbi_path}")
    return str(bgzf_path), tbi_path

def save_beta_matrix(beta_matrix: pd.DataFrame, output_path: str, 
                    compress: bool = True) -> None:
    """
    Save the beta matrix to file.
    
    Args:
        beta_matrix: DataFrame with beta values
        output_path: Path to save the file
        compress: Whether to compress the output
    """
    try:
        if compress and not output_path.endswith('.gz'):
            output_path += '.gz'
        
        if output_path.endswith('.gz'):
            beta_matrix.to_csv(output_path, sep='\t', compression='gzip', header=False, index=False)
        else:
            beta_matrix.to_csv(output_path, sep='\t', header=False, index=False)
        
        logger.info(f"Saved beta matrix to {output_path}")
        
    except Exception as e:
        logger.error(f"Error saving beta matrix: {e}")
        raise

def compute_beta_correlation_matrix(
    gwas_files: List[str],
    variants_to_extract: Set[str],
    variant_map: Dict[str, tuple],
    variant_col: str = "rs_number",
    beta_col: str = "beta",
    chunk_size: int = 50000,
    device: str = 'cuda',
    save_prefix: Optional[str] = None
) -> np.ndarray:
    """
    Extract beta values from tabix-indexed GWAS files and compute correlation matrix.
    Optionally save as PLINK `.vcor1` binary matrix and `.vcor1.vars` variant list.
    
    Args:
        gwas_files: List of paths to tabix-indexed GWAS summary statistics files.
        variants_to_extract: Set of variant IDs to extract.
        variant_map: Dictionary mapping variant ID to (chromosome, position).
        variant_col: Column name for variant identifiers (default: "rs_number").
        beta_col: Column name for beta values (default: "beta").
        chunk_size: Chunk size for processing (default: 50000).
        device: Device to use for computation ('cuda' or 'cpu', default: 'cuda').
        save_prefix: If provided, base filename (without extension) to save:
                     - `{save_prefix}.vcor1` (float32 square matrix binary)
                     - `{save_prefix}.vcor1.vars` (newline-delimited variant IDs)
    
    Returns:
        np.ndarray: Correlation matrix as a numpy array.
    """
    # Extract beta matrix
    betas_df = extract_betas_batch(
        file_list=gwas_files,
        variant_col=variant_col,
        beta_col=beta_col,
        variants_to_extract=variants_to_extract,
        chunk_size=chunk_size,
        variant_map=variant_map
    )
    
    # Convert to tensor on device
    betas = torch.tensor(betas_df.values, device=device, dtype=torch.float32)
    
    # Center rows
    betas_centered = betas - betas.mean(dim=1, keepdim=True)
    
    # Normalize rows
    norms = torch.norm(betas_centered, dim=1, keepdim=True)
    betas_normed = betas_centered / norms
    
    # Compute correlation matrix
    corr = betas_normed @ betas_normed.T
    
    # Move to CPU and numpy
    corr_np = corr.cpu().numpy().astype(np.float32)
    
    # Optionally save PLINK-compatible files
    if save_prefix:
        # Save binary matrix
        vcor1_path = f"{save_prefix}.vcor1"
        corr_np.tofile(vcor1_path)
        
        # Save variant list in order
        vars_path = f"{save_prefix}.vcor1.vars"
        with open(vars_path, "w") as vf:
            for variant in betas_df.index:
                vf.write(f"{variant}\n")
    
    return corr_np

def sample_traits(indexed_files, n_traits, random_seed, active_chromosome=None):
    """
    Sample a specified number of traits from available indexed files.
    
    Args:
        indexed_files: List of (bgzf_file, tbi_file) tuples
        n_traits: Number of traits to sample
        random_seed: Random seed for reproducible sampling
        
    Returns:
        List of selected bgzf file paths
    """
    if n_traits > len(indexed_files):
        raise ValueError(f"Requested {n_traits} traits, but only {len(indexed_files)} indexed files available")
    
    # Set random seed for reproducible sampling
    random.seed(random_seed)
    
    #Filter for chromosome if need be
    if active_chromosome:
        indexed_files = [x for x in indexed_files if f'chr{active_chromosome}' in x[0]]
    
    # Sample the specified number of files
    sampled_files = random.sample(indexed_files, n_traits)
    
    # Return just the bgzf file paths
    return [bgzf_file for bgzf_file, tbi_file in sampled_files]

