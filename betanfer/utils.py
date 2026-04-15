"""
betanfer/utils.py

Utility functions for indexing GWAS summary statistics and computing LD matrices
from beta correlations. Core logic ported from paper_analyses/5_LD_Inference/.
"""

import io
import os
import re
import glob
import shutil
import logging
import tempfile
import concurrent.futures
from pathlib import Path
from collections import defaultdict
from functools import partial
from multiprocessing import cpu_count
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import pysam
import torch

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Column name lists for auto-detection
# ---------------------------------------------------------------------------

COMMON_CHROM_COLS = ["chr", "chrom", "chromosome", "hm_chrom", "contig"]
COMMON_POS_COLS = [
    "pos", "bp", "position", "base_pair_location", "hm_pos",
    "base_pair", "bp_hg38", "bp_hg19",
]
COMMON_VARIANT_COLS = [
    "rs_number", "rsid", "snp", "variant_id", "variant",
    "id", "markerid", "snpid", "rsnum",
]
COMMON_BETA_COLS = ["beta", "effect", "effect_size", "b", "beta_hat", "log_or"]
COMMON_PVAL_COLS = [
    "p_value", "pvalue", "p.value", "pval", "p", "p_value", "log10p",
]

# ---------------------------------------------------------------------------
# Header reading helpers
# ---------------------------------------------------------------------------

def read_header(file_path: str) -> List[str]:
    """
    Read the header line from a plain-text file, strip leading '#' from the
    first token, and return the list of column names.
    """
    with open(file_path, "r") as fh:
        line = fh.readline().rstrip("\n")
    delim = "\t" if "\t" in line else None
    cols = line.split("\t") if delim else line.split()
    if cols and cols[0].startswith("#"):
        cols[0] = cols[0][1:]
    return cols


def read_tabix_header(gz_file_path: str) -> List[str]:
    """Return the column list from a tabix-indexed BGZF file's header."""
    with pysam.TabixFile(gz_file_path) as tf:
        header_lines = list(tf.header)
    if not header_lines:
        raise ValueError(f"No header found in tabix file: {gz_file_path}")
    cols = header_lines[0].strip().split("\t")
    if cols and cols[0].startswith("#"):
        cols[0] = cols[0][1:]
    return cols


# ---------------------------------------------------------------------------
# Column auto-detection
# ---------------------------------------------------------------------------

def autodetect_chrom_pos_cols(file_path: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Try to find chromosome and position column names from a plain-text file
    by matching against common names (case-insensitive).

    Returns (chrom_col, pos_col) — either may be None if not found.
    """
    cols = read_header(file_path)
    cols_lower = [c.lower() for c in cols]

    chrom_col = None
    for common in COMMON_CHROM_COLS:
        if common.lower() in cols_lower:
            chrom_col = cols[cols_lower.index(common.lower())]
            break

    pos_col = None
    for common in COMMON_POS_COLS:
        if common.lower() in cols_lower:
            pos_col = cols[cols_lower.index(common.lower())]
            break

    return chrom_col, pos_col


def autodetect_variant_beta_cols(
    gz_file_path: str,
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Try to find variant ID, beta, and p-value column names from a tabix
    BGZF file's header, matching against common names (case-insensitive).

    Returns (variant_col, beta_col, pval_col) — any may be None if not found.
    """
    cols = read_tabix_header(gz_file_path)
    cols_lower = [c.lower() for c in cols]

    def _match(candidates):
        for name in candidates:
            if name.lower() in cols_lower:
                return cols[cols_lower.index(name.lower())]
        return None

    variant_col = _match(COMMON_VARIANT_COLS)
    beta_col = _match(COMMON_BETA_COLS)
    pval_col = _match(COMMON_PVAL_COLS)
    return variant_col, beta_col, pval_col


# ---------------------------------------------------------------------------
# Column validation
# ---------------------------------------------------------------------------

def validate_column_in_file(file_path: str, col_name: str) -> None:
    """
    Validate that *col_name* exists in the header of a plain-text file.
    Raises ValueError with a helpful message if not found.
    """
    cols = read_header(file_path)
    if col_name not in cols:
        raise ValueError(
            f"Column '{col_name}' not found in {file_path}.\n"
            f"Available columns: {', '.join(cols)}"
        )


def validate_column_in_tabix(gz_file_path: str, col_name: str) -> None:
    """
    Validate that *col_name* exists in the tabix header of a BGZF file.
    Raises ValueError with a helpful message if not found.
    """
    cols = read_tabix_header(gz_file_path)
    if col_name not in cols:
        raise ValueError(
            f"Column '{col_name}' not found in {gz_file_path}.\n"
            f"Available columns: {', '.join(cols)}"
        )


# ---------------------------------------------------------------------------
# Device detection
# ---------------------------------------------------------------------------

def autodetect_device(force_cpu: bool = False) -> str:
    """Return 'xpu', 'cuda', or 'cpu' depending on available hardware.

    Checks for Intel XPU (via IPEX) first, then NVIDIA CUDA, then falls
    back to CPU.  IPEX is provided by the ALCF frameworks module on Aurora.
    """
    if force_cpu:
        logger.info("CPU mode requested via --cpu flag.")
        return "cpu"
    try:
        import intel_extension_for_pytorch as ipex  # noqa: F401
    except ImportError:
        pass  # IPEX not required; native torch XPU support may still be available
    if torch.xpu.is_available():
        logger.info("Intel XPU detected — using XPU for correlation computation.")
        return "xpu"
    if torch.cuda.is_available():
        logger.info("CUDA GPU detected — using GPU for correlation computation.")
        return "cuda"
    logger.info("No GPU detected — falling back to CPU.")
    return "cpu"


# ---------------------------------------------------------------------------
# Indexed file discovery
# ---------------------------------------------------------------------------

def find_indexed_sumstats_files(
    indexed_dir: str, num_threads: int = None
) -> List[List[str]]:
    """
    Find all tabix-indexed summary statistics files in *indexed_dir*.

    Returns a list of [gz_file, tbi_file] pairs.
    """
    if num_threads is None:
        num_threads = os.cpu_count()

    gz_files = glob.glob(os.path.join(indexed_dir, "*.gz"))
    if not gz_files:
        logger.warning(f"No .gz files found in {indexed_dir}")
        return []

    logger.info(
        f"Found {len(gz_files)} .gz files. Checking for .tbi indices "
        f"(using {num_threads} threads)..."
    )

    def _check(gz_file):
        tbi = gz_file + ".tbi"
        return gz_file, os.path.exists(tbi)

    indexed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
        for gz_file, has_tbi in [f.result() for f in
                                   concurrent.futures.as_completed(
                                       [ex.submit(_check, f) for f in gz_files])]:
            if has_tbi:
                indexed.append([gz_file, gz_file + ".tbi"])

    logger.info(f"Found {len(indexed)} indexed files (with .tbi indices)")
    return indexed


def get_indexed_files_from_file_list(
    file_list_path: str,
) -> Tuple[List[List[str]], List[str]]:
    """
    Read a file list (one path per line) and return validated (gz, tbi) pairs
    plus a list of paths that are missing their index.

    Each line may be either the raw sumstats path (the .gz will be derived)
    or the .gz path directly.
    """
    with open(file_list_path) as fh:
        paths = [
            line.strip()
            for line in fh
            if line.strip() and not line.startswith("#")
        ]

    pairs: List[List[str]] = []
    missing: List[str] = []

    for path in paths:
        gz_path = path if path.endswith(".gz") else path + ".gz"
        tbi_path = gz_path + ".tbi"
        if not os.path.exists(gz_path):
            missing.append(f"{path} (.gz file not found)")
        elif not os.path.exists(tbi_path):
            missing.append(f"{gz_path} (.tbi index not found)")
        else:
            pairs.append([gz_path, tbi_path])

    return pairs, missing


# ---------------------------------------------------------------------------
# Variant loading & position mapping
# ---------------------------------------------------------------------------

def load_variant_list(file_path: str) -> Set[str]:
    """Load variant IDs from a single-column file (one ID per line, no header)."""
    try:
        df = pd.read_csv(file_path, sep="\t", header=None)
        variants = set(df.iloc[:, 0].astype(str))
        logger.info(f"Loaded {len(variants)} variants from {file_path}")
        return variants
    except Exception as e:
        logger.error(f"Error loading variant list from {file_path}: {e}")
        raise


def make_var_dict_for_positional_ids(
    variants_to_extract: Set[str],
    id_conv: str = "CHR:POS:REF:ALT",
) -> Dict[str, Tuple[str, str]]:
    """
    Build a dict mapping variant ID → (chromosome, position) by parsing each
    variant ID according to *id_conv* (e.g. "CHR:POS:REF:ALT").

    The delimiter is inferred as any non-alphanumeric character.
    """
    if "CHR" not in id_conv or "POS" not in id_conv:
        raise ValueError("id_conv must contain 'CHR' and 'POS' fields.")

    conv_parts = re.split(r"\W+", id_conv)
    field_to_index = {f: i for i, f in enumerate(conv_parts)}
    chr_idx = field_to_index["CHR"]
    pos_idx = field_to_index["POS"]

    var_dict: Dict[str, Tuple[str, str]] = {}
    for var_id in variants_to_extract:
        parts = re.split(r"\W+", var_id)
        try:
            var_dict[var_id] = (parts[chr_idx], parts[pos_idx])
        except IndexError:
            raise ValueError(
                f"Variant ID '{var_id}' does not match the format: {id_conv}"
            )
    return var_dict


def get_region_bounds_from_variant_map(
    variant_map: Dict[str, Tuple],
) -> Tuple[str, int, int]:
    """Return (chromosome, min_pos, max_pos) from a variant position map."""
    chroms: Set[str] = set()
    positions: List[int] = []
    for _var, (chrom, pos) in variant_map.items():
        chroms.add(chrom)
        try:
            positions.append(int(pos))
        except (ValueError, TypeError):
            pass

    if len(chroms) != 1:
        raise ValueError(f"Expected single chromosome, found: {chroms}")
    if not positions:
        raise ValueError("No valid positions found in variant map")

    return chroms.pop(), min(positions), max(positions)


# ---------------------------------------------------------------------------
# Beta extraction
# ---------------------------------------------------------------------------

def extract_betas_batch(
    file_list: List[str],
    variant_col: str,
    beta_col: str,
    variants_to_extract: Set[str],
    variant_map: Dict[str, tuple],
    num_processes: int = 1,
) -> pd.DataFrame:
    """
    Extract betas using per-variant tabix queries.

    Use this when variants are NOT from a contiguous genomic region
    (i.e. the --non-contiguous flag is set).
    """
    variants_list = sorted(variants_to_extract)

    def _process(file_path: str) -> pd.Series:
        betas = pd.Series(np.nan, index=variants_list, name=file_path)
        try:
            with pysam.TabixFile(file_path) as tf:
                header = tf.header[0].strip().split("\t") if tf.header else []
                if header and header[0].startswith("#"):
                    header[0] = header[0][1:]
                vid_idx = header.index(variant_col)
                beta_idx = header.index(beta_col)
                for variant in variants_list:
                    chrom, pos = variant_map.get(variant, (None, None))
                    if not chrom or not pos:
                        continue
                    try:
                        for record in tf.fetch(
                            reference=chrom, start=int(pos) - 1, end=int(pos) + 1
                        ):
                            fields = record.split("\t")
                            if fields[vid_idx] == variant:
                                betas[variant] = float(fields[beta_idx])
                                break
                    except ValueError:
                        continue
        except Exception as e:
            logger.error(f"Error processing {file_path}: {e}")
        return betas

    num_threads = min(len(file_list), num_processes)
    total_files = len(file_list)
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
        futures = [ex.submit(_process, f) for f in file_list]
        results = []
        for i, fut in enumerate(concurrent.futures.as_completed(futures)):
            results.append(fut.result())
            if (i + 1) % 100 == 0:
                logger.info(f"  Extracted betas from {i + 1}/{total_files} files...")
    return pd.concat(results, axis=1)


def extract_betas_batch_region(
    file_list: List[str],
    variant_col: str,
    beta_col: str,
    variants_to_extract: Set[str],
    variant_map: Dict[str, tuple],
    num_processes: int = 1,
) -> pd.DataFrame:
    """
    Extract betas using a single batch region query per file.

    Assumes all variants in *variants_to_extract* span a contiguous genomic
    region (default for --variant-file without --non-contiguous).
    """
    chrom, min_pos, max_pos = get_region_bounds_from_variant_map(variant_map)
    logger.info(f"Region bounds: {chrom}:{min_pos}-{max_pos}")
    variants_list = sorted(variants_to_extract)

    def _process(file_path: str) -> pd.Series:
        betas = pd.Series(np.nan, index=variants_list, name=file_path)
        try:
            with pysam.TabixFile(file_path) as tf:
                header = tf.header[0].strip().split("\t") if tf.header else []
                if header and header[0].startswith("#"):
                    header[0] = header[0][1:]
                if not header:
                    logger.warning(f"No header in {file_path}")
                    return betas
                try:
                    vid_idx = header.index(variant_col)
                    beta_idx = header.index(beta_col)
                except ValueError as e:
                    logger.warning(f"Column not found in {file_path}: {e}")
                    return betas

                records_by_variant: Dict[str, List] = defaultdict(list)
                for record in tf.fetch(
                    reference=chrom, start=int(min_pos) - 1, end=int(max_pos) + 1
                ):
                    fields = record.split("\t")
                    try:
                        records_by_variant[fields[vid_idx]].append(fields)
                    except IndexError:
                        continue

                for variant in variants_list:
                    if variant in records_by_variant:
                        try:
                            betas[variant] = float(
                                records_by_variant[variant][0][beta_idx]
                            )
                        except (ValueError, IndexError):
                            continue
        except Exception as e:
            logger.error(f"Error processing {file_path}: {e}")
        return betas

    num_threads = min(len(file_list), num_processes)
    total_files = len(file_list)
    logger.info(f"Using {num_threads} threads to process {total_files} files")
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
        futures = [ex.submit(_process, f) for f in file_list]
        results = []
        for i, fut in enumerate(concurrent.futures.as_completed(futures)):
            results.append(fut.result())
            if (i + 1) % 100 == 0:
                logger.info(f"  Extracted betas from {i + 1}/{total_files} files...")
    beta_df = pd.concat(results, axis=1)
    logger.info(
        f"Extracted betas: {beta_df.shape[0]} variants × {beta_df.shape[1]} files"
    )
    return beta_df


def extract_betas_from_range(
    file_list: List[str],
    variant_col: str,
    beta_col: str,
    chrom: str,
    start: int,
    end: int,
    num_processes: int = 1,
) -> pd.DataFrame:
    """
    Discover variants in a genomic range from the first indexed file, then
    extract betas for those variants from all files.

    Used when the user specifies --chrom / --start / --end instead of a
    variant file.
    """
    # Step 1: discover variants from the first file in the list that has any in the range
    variants_to_extract: Set[str] = set()
    for discovery_file in file_list:
        with pysam.TabixFile(discovery_file) as tf:
            header = tf.header[0].strip().split("\t") if tf.header else []
            if header and header[0].startswith("#"):
                header[0] = header[0][1:]
            if not header:
                continue
            try:
                vid_idx = header.index(variant_col)
            except ValueError:
                continue
            try:
                for record in tf.fetch(reference=chrom, start=start - 1, end=end):
                    fields = record.split("\t")
                    try:
                        variants_to_extract.add(fields[vid_idx])
                    except IndexError:
                        continue
            except Exception as e:
                logger.warning(
                    f"Error querying range {chrom}:{start}-{end} in {discovery_file}: {e}"
                )
                continue
        if variants_to_extract:
            logger.info(
                f"Discovered {len(variants_to_extract)} variants in range "
                f"{chrom}:{start}-{end} from: {discovery_file}"
            )
            break

    if not variants_to_extract:
        logger.warning(
            f"No variants found in range {chrom}:{start}-{end} in any input file. "
            "Returning empty beta matrix."
        )
        return pd.DataFrame()

    variants_list = sorted(variants_to_extract)

    # Step 2: extract betas from all files for those variants
    def _process(file_path: str) -> pd.Series:
        betas = pd.Series(np.nan, index=variants_list, name=file_path)
        try:
            with pysam.TabixFile(file_path) as tf:
                header = tf.header[0].strip().split("\t") if tf.header else []
                if header and header[0].startswith("#"):
                    header[0] = header[0][1:]
                if not header:
                    logger.warning(f"No header in {file_path}")
                    return betas
                try:
                    vid_idx = header.index(variant_col)
                    beta_idx = header.index(beta_col)
                except ValueError as e:
                    logger.warning(f"Column not found in {file_path}: {e}")
                    return betas

                records_by_variant: Dict[str, List] = defaultdict(list)
                try:
                    for record in tf.fetch(
                        reference=chrom, start=start - 1, end=end
                    ):
                        fields = record.split("\t")
                        try:
                            records_by_variant[fields[vid_idx]].append(fields)
                        except IndexError:
                            continue
                except Exception as e:
                    logger.warning(f"Error fetching range from {file_path}: {e}")
                    return betas

                for variant in variants_list:
                    if variant in records_by_variant:
                        try:
                            betas[variant] = float(
                                records_by_variant[variant][0][beta_idx]
                            )
                        except (ValueError, IndexError):
                            continue
        except Exception as e:
            logger.error(f"Error processing {file_path}: {e}")
        return betas

    num_threads = min(len(file_list), num_processes)
    total_files = len(file_list)
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) as ex:
        futures = [ex.submit(_process, f) for f in file_list]
        results = []
        for i, fut in enumerate(concurrent.futures.as_completed(futures)):
            results.append(fut.result())
            if (i + 1) % 100 == 0:
                logger.info(f"  Extracted betas from {i + 1}/{total_files} files...")
    beta_df = pd.concat(results, axis=1)
    logger.info(
        f"Extracted betas: {beta_df.shape[0]} variants × {beta_df.shape[1]} files"
    )
    return beta_df


# ---------------------------------------------------------------------------
# Correlation matrix computation
# ---------------------------------------------------------------------------

def compute_beta_correlation_matrix(
    betas_df: pd.DataFrame,
    device: str = "cuda",
) -> np.ndarray:
    """
    Compute a pairwise correlation matrix from a beta DataFrame.

    Uses pairwise deletion for missing (NaN) betas: each pair of variants
    is correlated only over the traits where both have non-missing values,
    avoiding the attenuation caused by zero-imputation.

    Args:
        betas_df: Variants × traits DataFrame.
        device: 'cuda', 'xpu', or 'cpu'.

    Returns:
        Correlation matrix as a float32 numpy array.
    """
    values = betas_df.values.astype(np.float32)  # (n_variants, n_traits)

    # Build mask: 1.0 where non-missing, 0.0 where NaN
    mask = (~np.isnan(values)).astype(np.float32)
    values_filled = np.where(mask.astype(bool), values, 0.0).astype(np.float32)

    # Move to device
    X = torch.tensor(values_filled, device=device, dtype=torch.float32)  # (n_variants, n_traits)
    M = torch.tensor(mask, device=device, dtype=torch.float32)            # (n_variants, n_traits)

    # Note: variants are in rows and traits are in columns, so the matmul
    # orientations below are transposed relative to the genotype-matrix case
    # where variants are in columns and samples are in rows.

    # Pairwise counts: N[i,j] = number of traits with both i and j non-missing
    N = M @ M.T                         # (n_variants, n_variants)

    # Pairwise sums over jointly non-missing traits.
    # X_filled[i,k]=0 when missing, so masking with M[j,:] selects only
    # jointly non-missing positions (0-filled positions contribute nothing).
    sum_x  = X @ M.T                    # sum of X[i] where X[j] non-missing
    sum_y  = M @ X.T                    # sum of X[j] where X[i] non-missing
    sum_x2 = (X * X) @ M.T             # sum of X[i]^2 where X[j] non-missing
    sum_y2 = M @ (X * X).T             # sum of X[j]^2 where X[i] non-missing
    sum_xy = X @ X.T                    # sum of X[i]*X[j] (0s cancel missing)
    del X, M                            # no longer needed; free cache slots for below

    # Pairwise variances and covariance
    # Each intermediate is deleted as soon as it is consumed so the caching
    # allocator can reuse that memory rather than growing peak XPU usage.
    var_x  = sum_x2 - (sum_x * sum_x) / N;   del sum_x2
    var_y  = sum_y2 - (sum_y * sum_y) / N;   del sum_y2
    cov_xy = sum_xy - (sum_x * sum_y) / N;   del sum_x, sum_y, sum_xy

    # Pearson r
    denom = torch.sqrt(var_x * var_y);        del var_x, var_y
    r     = cov_xy / denom;                   del cov_xy

    # Set to NaN where there are too few jointly non-missing traits or
    # where the denominator is invalid (zero variance in either variant).
    # Use in-place assignment to avoid allocating a full new n×n tensor.
    min_n = 2
    invalid = (N < min_n) | ~(denom > 0);    del N, denom
    r[invalid] = float("nan");               del invalid

    # Clamp to [-1, 1] in-place to absorb floating-point rounding error
    r.clamp_(-1.0, 1.0)

    corr_np = r.cpu().numpy().astype(np.float32)
    del r

    return corr_np


def compute_beta_correlation_matrix_chunked(
    betas_df: pd.DataFrame,
    device: str = "cuda",
    chunk_size: int = 15000,
) -> np.ndarray:
    """
    Block-wise chunked pairwise correlation matrix from a beta DataFrame.

    Computes the same pairwise-deletion Pearson correlation as
    compute_beta_correlation_matrix, but processes the variant dimension in
    blocks of `chunk_size` rows to keep peak GPU memory at O(chunk_size²)
    instead of O(n²). Recommended for loci with >~20 000 variants that would
    otherwise exhaust GPU memory.

    Args:
        betas_df:   Variants × traits DataFrame.
        device:     'cuda', 'xpu', or 'cpu'.
        chunk_size: Maximum number of variants per GPU block.

    Returns:
        Correlation matrix as a float32 numpy array.
    """
    values = betas_df.values.astype(np.float32)
    mask = (~np.isnan(values)).astype(np.float32)
    values_filled = np.where(mask.astype(bool), values, 0.0).astype(np.float32)

    n = len(betas_df)
    corr_np = np.empty((n, n), dtype=np.float32)

    n_chunks = (n + chunk_size - 1) // chunk_size
    n_blocks = n_chunks * (n_chunks + 1) // 2
    block_idx = 0

    for i_start in range(0, n, chunk_size):
        i_end = min(i_start + chunk_size, n)
        X_r = torch.tensor(values_filled[i_start:i_end, :], device=device, dtype=torch.float32)
        M_r = torch.tensor(mask[i_start:i_end, :], device=device, dtype=torch.float32)

        for j_start in range(i_start, n, chunk_size):
            j_end = min(j_start + chunk_size, n)
            block_idx += 1
            logger.info(f"  Block {block_idx}/{n_blocks}: rows {i_start}-{i_end} x cols {j_start}-{j_end}")

            if j_start == i_start:
                X_c, M_c = X_r, M_r
            else:
                X_c = torch.tensor(values_filled[j_start:j_end, :], device=device, dtype=torch.float32)
                M_c = torch.tensor(mask[j_start:j_end, :], device=device, dtype=torch.float32)

            N      = M_r @ M_c.T
            sum_x  = X_r @ M_c.T
            sum_y  = M_r @ X_c.T
            sum_x2 = (X_r * X_r) @ M_c.T
            sum_y2 = M_r @ (X_c * X_c).T
            sum_xy = X_r @ X_c.T

            var_x  = sum_x2 - (sum_x * sum_x) / N;  del sum_x2
            var_y  = sum_y2 - (sum_y * sum_y) / N;  del sum_y2
            cov_xy = sum_xy - (sum_x * sum_y) / N;  del sum_x, sum_y, sum_xy

            denom  = torch.sqrt(var_x * var_y);      del var_x, var_y
            r      = cov_xy / denom;                 del cov_xy

            min_n = 2
            invalid = (N < min_n) | ~(denom > 0);   del N, denom
            r[invalid] = float("nan");               del invalid
            r.clamp_(-1.0, 1.0)

            block_np = r.cpu().numpy().astype(np.float32)
            del r

            corr_np[i_start:i_end, j_start:j_end] = block_np
            if j_start != i_start:
                corr_np[j_start:j_end, i_start:i_end] = block_np.T
                del X_c, M_c

        del X_r, M_r

    return corr_np


def save_beta_matrix(
    beta_matrix: pd.DataFrame, output_path: str, compress: bool = False
) -> None:
    """Save a matrix (with index/column labels) to a TSV file."""
    if compress and not output_path.endswith(".gz"):
        output_path += ".gz"
    if output_path.endswith(".gz"):
        beta_matrix.to_csv(output_path, sep="\t", compression="gzip")
    else:
        beta_matrix.to_csv(output_path, sep="\t")
    logger.info(f"Saved matrix to {output_path}")


def save_plink2_matrix(
    corr_np: np.ndarray,
    variants_list: List[str],
    output_prefix: str,
) -> None:
    """Save a correlation matrix in plink2 format.

    Creates two files:
      {output_prefix}.vcor1       — tab-delimited matrix values, no row/column names
      {output_prefix}.vcor1.vars  — variant IDs, one per line (same order as matrix rows/cols)
    """
    vcor1_path = f"{output_prefix}.vcor1"
    vars_path  = f"{output_prefix}.vcor1.vars"
    # Write row-by-row in ≤64 KB chunks to avoid E2BIG on Lustre for wide matrices
    _CHUNK = 65536
    with open(vcor1_path, "wb") as fh:
        for i in range(corr_np.shape[0]):
            buf = io.BytesIO()
            np.savetxt(buf, corr_np[i:i+1], delimiter="\t", fmt="%.6g")
            line = buf.getvalue()
            for k in range(0, len(line), _CHUNK):
                fh.write(line[k:k+_CHUNK])
    with open(vars_path, "w") as vf:
        for variant in variants_list:
            vf.write(f"{variant}\n")
    logger.info(f"Saved plink2 matrix to {vcor1_path} and {vars_path}")


def save_binary_matrix(
    corr_np: np.ndarray,
    variants_list: List[str],
    output_prefix: str,
) -> None:
    """Save a correlation matrix in numpy binary format.

    Creates two files:
      {output_prefix}.bin   — raw float32 binary (row-major)
      {output_prefix}.vars  — variant IDs, one per line (same order as matrix rows/cols)
    """
    bin_path  = f"{output_prefix}.bin"
    vars_path = f"{output_prefix}.vars"
    corr_np.astype(np.float32).tofile(bin_path)
    with open(vars_path, "w") as vf:
        for variant in variants_list:
            vf.write(f"{variant}\n")
    logger.info(f"Saved binary matrix to {bin_path} and {vars_path}")


# ---------------------------------------------------------------------------
# P-value filtering — variant-map mode
# ---------------------------------------------------------------------------

def filter_files_by_chromosome(
    indexed_files: List[List[str]], active_chromosome: str
) -> List[str]:
    """
    Filter indexed file pairs to keep only those whose filename contains
    'chr{active_chromosome}.' (case-sensitive).

    Returns a list of gz file paths (not pairs).
    """
    pattern = f"chr{active_chromosome}."
    filtered = [
        pair[0]
        for pair in indexed_files
        if pattern in os.path.basename(pair[0])
    ]
    if not filtered:
        logger.warning(
            f"No files matched 'chr{active_chromosome}.' — using all files."
        )
        filtered = [pair[0] for pair in indexed_files]
    return filtered


def get_buffer_region_variants(
    variant_map: Dict[str, Tuple], buffer_bp: int
) -> Tuple[str, int, int]:
    """
    Return (chromosome, buffer_start, buffer_end) extending the variant
    region by *buffer_bp* base pairs on each side.
    """
    chrom, min_pos, max_pos = get_region_bounds_from_variant_map(variant_map)
    buffer_start = max(1, min_pos - buffer_bp)
    buffer_end = max_pos + buffer_bp
    logger.info(
        f"Buffer region: {chrom}:{buffer_start}-{buffer_end} "
        f"(original: {min_pos}-{max_pos}, buffer: {buffer_bp} bp)"
    )
    return chrom, buffer_start, buffer_end


def check_trait_pvalue_significance_optimized(
    trait_file: str,
    variants_to_check: Set[str],
    variant_map: Dict[str, Tuple],
    variant_col: str,
    pval_col: str,
    pval_threshold: float,
    buffer_region: Optional[Tuple] = None,
) -> Tuple[bool, Dict]:
    """
    Check whether a trait has any p-value below *pval_threshold* for the
    specified variants (and optional buffer region).

    Optimisations:
    - Single batch tabix fetch per region (not per variant)
    - Early exit on first significant p-value found
    """
    stats: Dict = {
        "trait_file": os.path.basename(trait_file),
        "variants_checked": 0,
        "variants_found": 0,
        "significant_variants": [],
        "significant_in_matrix": False,
        "significant_in_buffer": False,
        "min_pval": 1.0,
        "variants_with_data": [],
        "matrix_tabix_calls": 1,
        "buffer_tabix_calls": 1 if buffer_region else 0,
    }

    try:
        with pysam.TabixFile(trait_file) as tf:
            header = tf.header[0].strip().split("\t") if tf.header else []
            if header and header[0].startswith("#"):
                header[0] = header[0][1:]
            if not header:
                logger.warning(f"No header in {trait_file}")
                return False, stats
            try:
                variant_idx = header.index(variant_col)
                pval_idx = header.index(pval_col)
            except ValueError as e:
                logger.warning(f"Column not found in {trait_file}: {e}")
                return False, stats

            # Batch fetch for matrix region
            chrom, min_pos, max_pos = get_region_bounds_from_variant_map(variant_map)
            records_by_variant: Dict[str, List] = defaultdict(list)
            for record in tf.fetch(
                reference=chrom, start=min_pos - 1, end=max_pos + 1
            ):
                fields = record.split("\t")
                try:
                    records_by_variant[fields[variant_idx]].append(fields)
                except IndexError:
                    continue

            for variant in variants_to_check:
                stats["variants_checked"] += 1
                if variant in records_by_variant:
                    for fields in records_by_variant[variant]:
                        stats["variants_found"] += 1
                        stats["variants_with_data"].append(variant)
                        try:
                            pval = float(fields[pval_idx])
                            stats["min_pval"] = min(stats["min_pval"], pval)
                            if pval < pval_threshold:
                                stats["significant_variants"].append((variant, pval))
                                stats["significant_in_matrix"] = True
                                return True, stats  # Early exit
                        except (ValueError, IndexError):
                            continue
                        break

            # Buffer region check
            if buffer_region:
                buf_chrom, buf_start, buf_end = buffer_region
                try:
                    for record in tf.fetch(
                        reference=buf_chrom, start=buf_start, end=buf_end
                    ):
                        fields = record.split("\t")
                        try:
                            variant_id = fields[variant_idx]
                            pval = float(fields[pval_idx])
                            if variant_id in variants_to_check:
                                continue
                            stats["min_pval"] = min(stats["min_pval"], pval)
                            if pval < pval_threshold:
                                stats["significant_variants"].append(
                                    (variant_id, pval, "buffer")
                                )
                                stats["significant_in_buffer"] = True
                                return True, stats
                        except (ValueError, IndexError):
                            continue
                except Exception as e:
                    logger.debug(f"Buffer region query error in {trait_file}: {e}")

    except Exception as e:
        logger.warning(f"Error checking p-values in {trait_file}: {e}")
        return False, stats

    is_significant = stats["significant_in_matrix"] or stats["significant_in_buffer"]
    return is_significant, stats


def _check_trait_wrapper(
    trait_file: str,
    variants_to_check: Set[str],
    variant_map: Dict[str, Tuple],
    variant_col: str,
    pval_col: str,
    pval_threshold: float,
    buffer_region: Optional[Tuple] = None,
) -> Tuple[str, bool, Dict]:
    """Wrapper for parallel processing of p-value checking."""
    is_sig, stats = check_trait_pvalue_significance_optimized(
        trait_file=trait_file,
        variants_to_check=variants_to_check,
        variant_map=variant_map,
        variant_col=variant_col,
        pval_col=pval_col,
        pval_threshold=pval_threshold,
        buffer_region=buffer_region,
    )
    return trait_file, is_sig, stats


def filter_traits_by_pvalue(
    trait_files: List[str],
    variants_to_check: Set[str],
    variant_map: Dict[str, Tuple],
    variant_col: str,
    pval_col: str,
    pval_threshold: float,
    buffer_region: Optional[Tuple] = None,
    num_processes: Optional[int] = None,
) -> Tuple[List[str], Dict]:
    """
    Filter trait files keeping only those with no significant p-values in the
    variant region.  Parallelised across CPU cores.
    """
    if num_processes is None:
        num_processes = cpu_count()
    logger.info(f"Filtering {len(trait_files)} traits using {num_processes} threads")

    check_fn = partial(
        _check_trait_wrapper,
        variants_to_check=variants_to_check,
        variant_map=variant_map,
        variant_col=variant_col,
        pval_col=pval_col,
        pval_threshold=pval_threshold,
        buffer_region=buffer_region,
    )

    filtering_stats: Dict = {
        "total_traits": len(trait_files),
        "traits_checked": 0,
        "traits_retained": 0,
        "traits_filtered_out": 0,
        "min_pval_across_all": 1.0,
        "trait_details": [],
        "num_processes": num_processes,
        "mode": "variant_map",
    }
    filtered_files: List[str] = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_processes) as ex:
        futures = {ex.submit(check_fn, f): f for f in trait_files}
        for i, future in enumerate(concurrent.futures.as_completed(futures)):
            trait_file, is_sig, stats = future.result()
            filtering_stats["traits_checked"] += 1
            if (i + 1) % 100 == 0:
                logger.info(f"Filtered {i + 1}/{len(trait_files)} traits...")
            filtering_stats["min_pval_across_all"] = min(
                filtering_stats["min_pval_across_all"], stats["min_pval"]
            )
            filtering_stats["trait_details"].append(stats)
            if not is_sig:
                filtered_files.append(trait_file)
                filtering_stats["traits_retained"] += 1
            else:
                filtering_stats["traits_filtered_out"] += 1

    return filtered_files, filtering_stats


# ---------------------------------------------------------------------------
# P-value filtering — range mode
# ---------------------------------------------------------------------------

def _check_pvalue_in_range(args_tuple: tuple) -> Tuple[str, bool, Dict]:
    """
    Check whether any variant in a genomic range has p < threshold.
    Expects args_tuple = (trait_file, chrom, start, end, pval_col, pval_threshold).
    """
    trait_file, chrom, start, end, pval_col, pval_threshold = args_tuple
    details: Dict = {
        "trait_file": os.path.basename(trait_file),
        "min_pval": 1.0,
        "significant": False,
    }
    try:
        with pysam.TabixFile(trait_file) as tf:
            header = tf.header[0].strip().split("\t") if tf.header else []
            if not header:
                return trait_file, False, details
            try:
                pval_idx = header.index(pval_col)
            except ValueError:
                logger.warning(
                    f"Column '{pval_col}' not found in {trait_file}"
                )
                return trait_file, False, details
            try:
                for record in tf.fetch(
                    reference=chrom, start=start - 1, end=end
                ):
                    fields = record.split("\t")
                    try:
                        pval = float(fields[pval_idx])
                        details["min_pval"] = min(details["min_pval"], pval)
                        if pval < pval_threshold:
                            details["significant"] = True
                            return trait_file, True, details  # Early exit
                    except (ValueError, IndexError):
                        continue
            except Exception as e:
                logger.debug(f"Range query error in {trait_file}: {e}")
    except Exception as e:
        logger.warning(f"Error reading {trait_file}: {e}")
    return trait_file, False, details


def filter_traits_by_pvalue_in_range(
    file_list: List[str],
    chrom: str,
    start: int,
    end: int,
    pval_col: str,
    pval_threshold: float,
    buffer: int = 0,
    num_processes: Optional[int] = None,
) -> Tuple[List[str], Dict]:
    """
    Filter trait files by p-value significance in a genomic range.
    Used when variants are specified via --chrom/--start/--end.
    """
    if num_processes is None:
        num_processes = cpu_count()

    query_start = max(1, start - buffer)
    query_end = end + buffer
    logger.info(
        f"Filtering {len(file_list)} traits for p < {pval_threshold} "
        f"in {chrom}:{query_start}-{query_end} using {num_processes} threads"
    )

    args_list = [
        (f, chrom, query_start, query_end, pval_col, pval_threshold)
        for f in file_list
    ]
    filtering_stats: Dict = {
        "total_traits": len(file_list),
        "traits_retained": 0,
        "traits_filtered_out": 0,
        "trait_details": [],
        "mode": "range",
    }
    filtered_files: List[str] = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_processes) as ex:
        futures = {ex.submit(_check_pvalue_in_range, args): args[0] for args in args_list}
        for i, future in enumerate(concurrent.futures.as_completed(futures)):
            trait_file, is_sig, details = future.result()
            if (i + 1) % 100 == 0:
                logger.info(f"Filtered {i + 1}/{len(file_list)} traits...")
            filtering_stats["trait_details"].append(details)
            if not is_sig:
                filtered_files.append(trait_file)
                filtering_stats["traits_retained"] += 1
            else:
                filtering_stats["traits_filtered_out"] += 1

    return filtered_files, filtering_stats


def write_filtering_stats(
    stats: Dict, output_path: str, pval_threshold: float, buffer_bp: int = 0
) -> None:
    """Write per-trait filtering statistics to a TSV file."""
    mode = stats.get("mode", "variant_map")
    with open(output_path, "w") as fh:
        if mode == "range":
            fh.write("Trait\tMin_Pval\tStatus\n")
            for detail in stats["trait_details"]:
                status = "FILTERED" if detail["significant"] else "RETAINED"
                fh.write(
                    f"{detail['trait_file']}\t{detail['min_pval']:.2e}\t{status}\n"
                )
        else:
            fh.write("Trait\tChecked\tFound\tSignificant\tMin_Pval\tStatus\n")
            for detail in stats["trait_details"]:
                is_sig = detail.get("significant_in_matrix") or detail.get(
                    "significant_in_buffer"
                )
                status = "FILTERED" if is_sig else "RETAINED"
                fh.write(
                    f"{detail['trait_file']}\t"
                    f"{detail['variants_checked']}\t"
                    f"{detail['variants_found']}\t"
                    f"{len(detail['significant_variants'])}\t"
                    f"{detail['min_pval']:.2e}\t"
                    f"{status}\n"
                )


# ---------------------------------------------------------------------------
# Indexing functions
# ---------------------------------------------------------------------------

def _sort_df_by_chrom_pos(df: pd.DataFrame, chrom_col: str, pos_col: str) -> pd.DataFrame:
    """Sort a DataFrame by chromosome (numerically where possible) then position."""
    def _key(chrom):
        try:
            return (0, int(chrom))
        except (ValueError, TypeError):
            return (1, str(chrom))

    df = df.copy()
    df["_sort_key"] = df[chrom_col].apply(_key)
    df = (
        df.sort_values(["_sort_key", pos_col])
        .drop("_sort_key", axis=1)
        .reset_index(drop=True)
    )
    return df


def bgzip_and_tbi_with_chrom_pos_cols(
    input_path: str,
    output_dir: str,
    chrom_col: str,
    pos_col: str,
    zero_based: bool = False,
    force: bool = False,
) -> Tuple[str, str]:
    """
    BGZF-compress a GWAS summary stats file and tabix-index it, using
    pre-existing chromosome and position columns.

    Handles files whose header may begin with '#'. Sorts by chromosome and
    position (required by tabix) and writes a '#'-prefixed header so that
    the header line is preserved in the tabix index.
    """
    input_path = str(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read file
    with open(input_path) as fh:
        first_line = fh.readline()
    delim = "\t" if "\t" in first_line else r"\s+"
    df = pd.read_csv(input_path, sep=delim, dtype=str)

    # Strip leading '#' from first column name if present
    if df.columns[0].startswith("#"):
        df.columns = pd.Index([df.columns[0][1:]] + list(df.columns[1:]))

    # Validate requested columns
    for col in (chrom_col, pos_col):
        if col not in df.columns:
            raise ValueError(
                f"Column '{col}' not found in {input_path}. "
                f"Available: {', '.join(df.columns)}"
            )

    # Coerce types and drop bad rows
    df[pos_col] = pd.to_numeric(df[pos_col], errors="coerce")
    df = df.dropna(subset=[chrom_col, pos_col])
    df[chrom_col] = df[chrom_col].astype(str)
    df[pos_col] = df[pos_col].astype(int)

    # Sort
    df = _sort_df_by_chrom_pos(df, chrom_col, pos_col)
    header_cols = list(df.columns)

    # Write with '#' prefix using a temp file
    in_name = Path(input_path).name
    bgzf_path = output_dir / (in_name + ".gz")
    tbi_path = str(bgzf_path) + ".tbi"

    if bgzf_path.exists():
        if force:
            bgzf_path.unlink()
        else:
            raise FileExistsError(f"Output already exists: {bgzf_path}")
    if force and os.path.exists(tbi_path):
        os.remove(tbi_path)

    with tempfile.NamedTemporaryFile(
        mode="w", dir=str(output_dir), suffix=".betanfer_tmp", delete=False
    ) as tmp:
        tmp_path = tmp.name
        tmp.write("#" + "\t".join(header_cols) + "\n")
        df.to_csv(tmp, sep="\t", index=False, header=False)

    del df

    # BGZIP
    with open(tmp_path, "rb") as fin, pysam.BGZFile(bgzf_path, "w") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)
    os.remove(tmp_path)

    # Tabix index
    chrom_col_idx = header_cols.index(chrom_col)
    pos_col_idx = header_cols.index(pos_col)
    pysam.tabix_index(
        str(bgzf_path),
        seq_col=chrom_col_idx,
        start_col=pos_col_idx,
        end_col=pos_col_idx,
        line_skip=0,
        zerobased=zero_based,
        meta_char="#",
    )
    if not os.path.exists(tbi_path):
        raise RuntimeError(f"Tabix index was not created: {tbi_path}")

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
) -> Tuple[str, str]:
    """
    Extract chromosome and position from variant IDs, then BGZF-compress
    and tabix-index the file.  Does not require a .bim file.
    """
    input_path = str(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    conv_parts = re.split(r"\W+", variant_id_format)
    field_to_index = {f: i for i, f in enumerate(conv_parts)}
    if "CHR" not in field_to_index or "POS" not in field_to_index:
        raise ValueError("variant_id_format must contain 'CHR' and 'POS' fields.")
    chr_index = field_to_index["CHR"]
    pos_index = field_to_index["POS"]

    with open(input_path) as fh:
        first_line = fh.readline()
    delim = "\t" if "\t" in first_line else r"\s+"
    df = pd.read_csv(input_path, sep=delim, dtype={variant_col: str})

    # Strip leading '#' from first column name
    if df.columns[0].startswith("#"):
        df.columns = pd.Index([df.columns[0][1:]] + list(df.columns[1:]))

    def _parse(var_id):
        if pd.isna(var_id) or var_id == "":
            return pd.NA, pd.NA
        try:
            parts = re.split(r"\W+", str(var_id))
            return parts[chr_index], int(parts[pos_index])
        except (IndexError, ValueError):
            return pd.NA, pd.NA

    parsed = df[variant_col].apply(_parse)
    df[chrom_col_label] = [c[0] for c in parsed]
    df[pos_col_label] = [c[1] for c in parsed]
    df = df.dropna(subset=[chrom_col_label, pos_col_label])
    df[chrom_col_label] = df[chrom_col_label].astype(str)
    df[pos_col_label] = df[pos_col_label].astype(int)
    df = _sort_df_by_chrom_pos(df, chrom_col_label, pos_col_label)
    header_cols = list(df.columns)

    in_name = Path(input_path).name
    bgzf_path = output_dir / (in_name + ".gz")
    tbi_path = str(bgzf_path) + ".tbi"

    if bgzf_path.exists():
        if force:
            bgzf_path.unlink()
        else:
            raise FileExistsError(f"Output already exists: {bgzf_path}")
    if force and os.path.exists(tbi_path):
        os.remove(tbi_path)

    with tempfile.NamedTemporaryFile(
        mode="w", dir=str(output_dir), suffix=".betanfer_tmp", delete=False
    ) as tmp:
        tmp_path = tmp.name
        tmp.write("#" + "\t".join(header_cols) + "\n")
        df.to_csv(tmp, sep="\t", index=False, header=False)
    del df

    with open(tmp_path, "rb") as fin, pysam.BGZFile(bgzf_path, "w") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)
    os.remove(tmp_path)

    chrom_col_idx = header_cols.index(chrom_col_label)
    pos_col_idx = header_cols.index(pos_col_label)
    pysam.tabix_index(
        str(bgzf_path),
        seq_col=chrom_col_idx,
        start_col=pos_col_idx,
        end_col=pos_col_idx,
        line_skip=0,
        zerobased=zero_based,
        meta_char="#",
    )
    if not os.path.exists(tbi_path):
        raise RuntimeError(f"Tabix index was not created: {tbi_path}")
    return str(bgzf_path), tbi_path


def merge_bim_files(bim_file_paths: List[str], output_path: Optional[str] = None) -> pd.DataFrame:
    """Merge multiple PLINK .bim files into a single deduplicated DataFrame."""
    col_names = ["CHR", "SNP", "CM", "POS", "A1", "A2"]
    frames = [
        pd.read_csv(p, sep="\t", header=None, names=col_names, dtype=str)
        for p in bim_file_paths
    ]
    merged = pd.concat(frames, ignore_index=True).drop_duplicates()
    merged["POS"] = pd.to_numeric(merged["POS"], errors="coerce")
    merged["CHR"] = pd.to_numeric(merged["CHR"], errors="coerce")
    merged = merged.sort_values(["CHR", "POS"]).reset_index(drop=True)
    if output_path:
        merged.to_csv(output_path, sep="\t", header=False, index=False)
    return merged


def bgzip_and_tbi_sumstats_with_bim_merge(
    input_path: str,
    output_dir: str,
    bim_df: pd.DataFrame,
    variant_col_gwama: str,
    chrom_col_label: str = "CHR",
    pos_col_label: str = "POS",
    zero_based: bool = False,
    force: bool = True,
) -> Tuple[str, str]:
    """
    Merge chromosome/position from a .bim DataFrame into a summary stats file,
    then BGZF-compress and tabix-index.
    """
    input_path = str(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    bim_df = bim_df.loc[:, ["CHR", "POS", "SNP"]].rename(
        columns={"SNP": variant_col_gwama}
    )
    bim_df["CHR"] = bim_df["CHR"].astype(str)
    bim_df["POS"] = pd.to_numeric(bim_df["POS"], errors="coerce")

    with open(input_path) as fh:
        first_line = fh.readline()
    delim = "\t" if "\t" in first_line else r"\s+"
    gwama_df = pd.read_csv(input_path, sep=delim, dtype={variant_col_gwama: str})

    # Strip '#' from first column name
    if gwama_df.columns[0].startswith("#"):
        gwama_df.columns = pd.Index(
            [gwama_df.columns[0][1:]] + list(gwama_df.columns[1:])
        )

    merged = pd.merge(bim_df, gwama_df, on=variant_col_gwama, how="inner")
    merged = merged.rename(
        columns={"CHR": chrom_col_label, "POS": pos_col_label}
    )
    merged[chrom_col_label] = merged[chrom_col_label].astype(str)
    merged[pos_col_label] = pd.to_numeric(merged[pos_col_label], errors="coerce")
    merged = merged.dropna(subset=[chrom_col_label, pos_col_label])
    merged[pos_col_label] = merged[pos_col_label].astype(int)
    merged = _sort_df_by_chrom_pos(merged, chrom_col_label, pos_col_label)
    header_cols = list(merged.columns)

    in_name = Path(input_path).name
    out_name = in_name if ".merged" in in_name else in_name + ".merged"
    bgzf_path = output_dir / (out_name + ".gz")
    tbi_path = str(bgzf_path) + ".tbi"

    if bgzf_path.exists():
        if force:
            bgzf_path.unlink()
        else:
            raise FileExistsError(f"Output already exists: {bgzf_path}")
    if force and os.path.exists(tbi_path):
        os.remove(tbi_path)

    with tempfile.NamedTemporaryFile(
        mode="w", dir=str(output_dir), suffix=".betanfer_tmp", delete=False
    ) as tmp:
        tmp_path = tmp.name
        tmp.write("#" + "\t".join(header_cols) + "\n")
        merged.to_csv(tmp, sep="\t", index=False, header=False)
    del gwama_df, merged

    with open(tmp_path, "rb") as fin, pysam.BGZFile(bgzf_path, "w") as fout:
        shutil.copyfileobj(fin, fout, length=1024 * 1024)
    os.remove(tmp_path)

    chrom_col_idx = header_cols.index(chrom_col_label)
    pos_col_idx = header_cols.index(pos_col_label)
    pysam.tabix_index(
        str(bgzf_path),
        seq_col=chrom_col_idx,
        start_col=pos_col_idx,
        end_col=pos_col_idx,
        line_skip=0,
        zerobased=zero_based,
        meta_char="#",
    )
    if not os.path.exists(tbi_path):
        raise RuntimeError(f"Tabix index was not created: {tbi_path}")
    return str(bgzf_path), tbi_path
