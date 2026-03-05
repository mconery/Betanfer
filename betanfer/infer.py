"""
betanfer/infer.py

Infer mode: extract betas from tabix-indexed GWAS summary statistics and
compute pairwise LD correlations between variants.
"""

import os
import sys
import logging
from pathlib import Path
from typing import List, Optional, Set

import pandas as pd

from betanfer.utils import (
    # File discovery
    find_indexed_sumstats_files,
    get_indexed_files_from_file_list,
    # Column auto-detection / validation
    autodetect_variant_beta_cols,
    validate_column_in_tabix,
    # Variant loading
    load_variant_list,
    make_var_dict_for_positional_ids,
    get_region_bounds_from_variant_map,
    # Beta extraction
    extract_betas_batch,
    extract_betas_batch_region,
    extract_betas_from_range,
    # P-value filtering
    filter_files_by_chromosome,
    get_buffer_region_variants,
    filter_traits_by_pvalue,
    filter_traits_by_pvalue_in_range,
    write_filtering_stats,
    # Computation & output
    compute_beta_correlation_matrix,
    save_beta_matrix,
    save_plink2_matrix,
    save_binary_matrix,
    autodetect_device,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _resolve_indexed_files(args) -> List[List[str]]:
    """
    Return a list of [gz_path, tbi_path] pairs from --dir or --file-list.
    Exits with an error if any indexed file is missing its .tbi.
    """
    if args.dir and args.file_list:
        logger.error("Error: --dir and --file-list are mutually exclusive.")
        sys.exit(1)
    if not args.dir and not args.file_list:
        logger.error("Error: One of --dir or --file-list is required.")
        sys.exit(1)

    if args.dir:
        if not os.path.isdir(args.dir):
            logger.error(f"Error: Directory not found: {args.dir}")
            sys.exit(1)
        pairs = find_indexed_sumstats_files(args.dir)
        if not pairs:
            logger.error(
                f"Error: No tabix-indexed (.gz + .tbi) files found in {args.dir}. "
                "Run 'betanfer index' to create tabix indexes for your files."
            )
            sys.exit(1)
        return pairs
    else:
        if not os.path.exists(args.file_list):
            logger.error(f"Error: File list not found: {args.file_list}")
            sys.exit(1)
        pairs, missing = get_indexed_files_from_file_list(args.file_list)
        if missing:
            logger.error(
                "Error: The following files are missing a tabix index (.tbi):"
            )
            for m in missing:
                logger.error(f"  {m}")
            logger.error(
                "Run 'betanfer index' to create tabix indexes for your files."
            )
            sys.exit(1)
        if not pairs:
            logger.error("Error: No valid indexed files found in the provided list.")
            sys.exit(1)
        return pairs


def _resolve_columns(args, first_gz: str):
    """
    Resolve variant_col, beta_col, and pval_col from user args or by
    auto-detecting from the first indexed file.

    Returns (variant_col, beta_col, pval_col) where pval_col may be None
    if p-value filtering is not requested.
    """
    pval_filtering = bool(args.pval_threshold is not None or args.pval_col)

    # Auto-detect what we can from the file
    detected_variant, detected_beta, detected_pval = autodetect_variant_beta_cols(
        first_gz
    )

    # --- variant_col ---
    if args.variant_col:
        try:
            validate_column_in_tabix(first_gz, args.variant_col)
            variant_col = args.variant_col
        except ValueError as e:
            logger.error(f"Error: {e}")
            sys.exit(1)
    else:
        if detected_variant is None:
            logger.error(
                "Error: Could not auto-detect variant ID column. "
                "Please specify --variant-col."
            )
            sys.exit(1)
        variant_col = detected_variant
        logger.info(f"Auto-detected variant column: '{variant_col}'")

    # --- beta_col ---
    if args.beta_col:
        try:
            validate_column_in_tabix(first_gz, args.beta_col)
            beta_col = args.beta_col
        except ValueError as e:
            logger.error(f"Error: {e}")
            sys.exit(1)
    else:
        if detected_beta is None:
            logger.error(
                "Error: Could not auto-detect beta column. "
                "Please specify --beta-col."
            )
            sys.exit(1)
        beta_col = detected_beta
        logger.info(f"Auto-detected beta column: '{beta_col}'")

    # --- pval_col ---
    pval_col = None
    if pval_filtering:
        if args.pval_col:
            try:
                validate_column_in_tabix(first_gz, args.pval_col)
                pval_col = args.pval_col
            except ValueError as e:
                logger.error(f"Error: {e}")
                sys.exit(1)
        else:
            if detected_pval is None:
                logger.error(
                    "Error: Could not auto-detect p-value column. "
                    "Please specify --pval-col."
                )
                sys.exit(1)
            pval_col = detected_pval
            logger.info(f"Auto-detected p-value column: '{pval_col}'")

    return variant_col, beta_col, pval_col


# ---------------------------------------------------------------------------
# Helpers (continued)
# ---------------------------------------------------------------------------

def _create_empty_outputs(args) -> None:
    """Create empty sentinel output files when LD inference cannot proceed."""
    if args.save_binary:
        paths = [args.output + ".bin", args.output + ".vars"]
    elif args.plink2_format:
        paths = [args.output + ".vcor1", args.output + ".vcor1.vars"]
    else:
        paths = [args.output + ".txt"]
    for path in paths:
        open(path, "w").close()
        logger.info(f"Created empty output file: {path}")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_infer(args) -> None:
    """Entry point for 'betanfer infer'."""

    # Validate variant specification
    range_mode = bool(
        args.chrom is not None and args.start is not None and args.end is not None
    )
    variant_file_mode = bool(args.variant_file)

    if range_mode and variant_file_mode:
        logger.error(
            "Error: --variant-file and --chrom/--start/--end are mutually exclusive."
        )
        sys.exit(1)
    if not range_mode and not variant_file_mode:
        logger.error(
            "Error: Specify variants using --variant-file or "
            "--chrom / --start / --end."
        )
        sys.exit(1)
    if args.non_contiguous and range_mode:
        logger.error(
            "Error: --non-contiguous is only valid when using --variant-file."
        )
        sys.exit(1)

    # Check p-value filtering activation
    pval_filtering = bool(args.pval_threshold is not None or args.pval_col)
    pval_threshold = args.pval_threshold if args.pval_threshold is not None else 5e-8

    # Validate output path
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1: Get indexed file list
    # ------------------------------------------------------------------
    logger.info("Step 1: Discovering indexed summary statistics files...")
    indexed_pairs = _resolve_indexed_files(args)
    logger.info(f"Found {len(indexed_pairs)} indexed file(s)")

    # ------------------------------------------------------------------
    # Step 2: Determine compute device
    # ------------------------------------------------------------------
    device = autodetect_device(force_cpu=args.cpu)

    # ------------------------------------------------------------------
    # Step 3: Determine parallelism
    # ------------------------------------------------------------------
    num_processes = args.num_processes if args.num_processes else os.cpu_count()

    # ------------------------------------------------------------------
    # Step 4: Resolve column names
    # ------------------------------------------------------------------
    first_gz = indexed_pairs[0][0]
    logger.info("Step 4: Resolving column names...")
    variant_col, beta_col, pval_col = _resolve_columns(args, first_gz)
    logger.info(
        f"Using: variant_col='{variant_col}', beta_col='{beta_col}'"
        + (f", pval_col='{pval_col}'" if pval_col else "")
    )

    # ------------------------------------------------------------------
    # Step 5: Load / discover variants
    # ------------------------------------------------------------------
    if variant_file_mode:
        logger.info("Step 5: Loading variants from file...")
        if not os.path.exists(args.variant_file):
            logger.error(f"Error: Variant file not found: {args.variant_file}")
            sys.exit(1)
        variants_to_extract = load_variant_list(args.variant_file)
        if not variants_to_extract:
            logger.error(f"Error: No variants loaded from {args.variant_file}")
            sys.exit(1)

        # Build position map (needed for contiguous region extraction and pval filtering)
        logger.info("Building variant position map...")
        try:
            variant_map = make_var_dict_for_positional_ids(
                variants_to_extract, args.variant_id_format
            )
        except ValueError as e:
            logger.error(f"Error building variant position map: {e}")
            sys.exit(1)

        # Determine the active chromosome (must be single)
        chroms = set(v[0] for v in variant_map.values())
        if len(chroms) != 1:
            logger.error(
                f"Error: Variant list spans multiple chromosomes: {chroms}. "
                "All variants must be on the same chromosome."
            )
            sys.exit(1)
        active_chromosome = chroms.pop()
        logger.info(f"Active chromosome: {active_chromosome}")

    else:
        logger.info(
            f"Step 5: Using genomic range: {args.chrom}:{args.start}-{args.end}"
        )
        active_chromosome = args.chrom
        # variant_map / variants_to_extract not needed for range mode
        variant_map = None
        variants_to_extract = None

    # ------------------------------------------------------------------
    # Step 6: P-value filtering (if requested)
    # ------------------------------------------------------------------
    # Start with all gz paths
    gz_files = [pair[0] for pair in indexed_pairs]

    if pval_filtering:
        logger.info(
            f"Step 6: Filtering traits by p-value "
            f"(threshold={pval_threshold:.2e}, column='{pval_col}')..."
        )
        if args.buffer:
            logger.info(f"  Genomic buffer: {args.buffer} bp")

        if variant_file_mode:
            # Filter by chromosome first (uses filename convention chr{N}.)
            gz_files = filter_files_by_chromosome(indexed_pairs, active_chromosome)
            logger.info(
                f"  Chromosome-filtered files: {len(gz_files)} "
                f"(chr{active_chromosome})"
            )

            buffer_region = None
            if args.buffer:
                buffer_region = get_buffer_region_variants(variant_map, args.buffer)

            gz_files, filtering_stats = filter_traits_by_pvalue(
                trait_files=gz_files,
                variants_to_check=variants_to_extract,
                variant_map=variant_map,
                variant_col=variant_col,
                pval_col=pval_col,
                pval_threshold=pval_threshold,
                buffer_region=buffer_region,
                num_processes=num_processes,
            )
        else:
            # Range mode — no chromosome-based filename filter
            gz_files, filtering_stats = filter_traits_by_pvalue_in_range(
                file_list=gz_files,
                chrom=args.chrom,
                start=args.start,
                end=args.end,
                pval_col=pval_col,
                pval_threshold=pval_threshold,
                buffer=args.buffer or 0,
                num_processes=num_processes,
            )

        retained = filtering_stats.get("traits_retained", len(gz_files))
        filtered_out = filtering_stats.get("traits_filtered_out", 0)
        logger.info(
            f"  Retained {retained} traits; filtered out {filtered_out} traits"
        )

        if not gz_files:
            logger.error(
                "Error: No traits remaining after p-value filtering. "
                "Consider relaxing --pval-threshold."
            )
            sys.exit(1)
    else:
        filtering_stats = None

    # ------------------------------------------------------------------
    # Step 7: Extract betas
    # ------------------------------------------------------------------
    logger.info(f"Step 7: Extracting betas from {len(gz_files)} file(s)...")

    if range_mode:
        beta_df = extract_betas_from_range(
            file_list=gz_files,
            variant_col=variant_col,
            beta_col=beta_col,
            chrom=args.chrom,
            start=args.start,
            end=args.end,
        )
    elif args.non_contiguous:
        beta_df = extract_betas_batch(
            file_list=gz_files,
            variant_col=variant_col,
            beta_col=beta_col,
            variants_to_extract=variants_to_extract,
            variant_map=variant_map,
        )
    else:
        beta_df = extract_betas_batch_region(
            file_list=gz_files,
            variant_col=variant_col,
            beta_col=beta_col,
            variants_to_extract=variants_to_extract,
            variant_map=variant_map,
        )

    logger.info(
        f"Beta matrix: {beta_df.shape[0]} variants × {beta_df.shape[1]} traits"
    )
    if beta_df.size > 0:
        missing_frac = beta_df.isna().sum().sum() / beta_df.size
        if missing_frac > 0:
            logger.info(f"Missing data fraction: {missing_frac:.3f} (filled with 0)")

    # ------------------------------------------------------------------
    # Step 7.5: Filter traits with fewer than 2 variants in the locus
    # ------------------------------------------------------------------
    min_traits = args.min_traits

    # Count non-NaN betas per trait (column = one file/trait)
    trait_variant_counts = beta_df.notna().sum(axis=0)

    # Warn and drop traits with < 2 variants
    low_variant_traits = trait_variant_counts[trait_variant_counts < 2].index.tolist()
    if low_variant_traits:
        logger.warning(
            f"Dropping {len(low_variant_traits)} trait(s) with fewer than 2 variants "
            f"in the locus:"
        )
        for t in low_variant_traits:
            logger.warning(f"  {t}")
        beta_df = beta_df.drop(columns=low_variant_traits)

    n_valid_traits = beta_df.shape[1]
    logger.info(f"Traits with >= 2 variants: {n_valid_traits}")

    if n_valid_traits < min_traits:
        logger.warning(
            f"Insufficient traits for LD inference: {n_valid_traits} trait(s) with "
            f">= 2 variants in the locus (minimum required: {min_traits}). "
            "Creating empty output files and skipping inference."
        )
        _create_empty_outputs(args)
        return

    # ------------------------------------------------------------------
    # Step 8: Compute correlation matrix
    # ------------------------------------------------------------------
    logger.info(f"Step 8: Computing LD matrix on device='{device}'...")
    corr_np = compute_beta_correlation_matrix(
        betas_df=beta_df,
        device=device,
    )
    logger.info(f"LD matrix shape: {corr_np.shape[0]} × {corr_np.shape[1]}")

    # ------------------------------------------------------------------
    # Step 9: Save output
    # ------------------------------------------------------------------
    variants_list = list(beta_df.index)

    if args.save_binary:
        save_binary_matrix(corr_np, variants_list, args.output)
        logger.info(
            f"Step 9: Saved binary LD matrix:\n"
            f"  Matrix  : {args.output}.bin\n"
            f"  Variants: {args.output}.vars"
        )
    elif args.plink2_format:
        save_plink2_matrix(corr_np, variants_list, args.output)
        logger.info(
            f"Step 9: Saved plink2 LD matrix:\n"
            f"  Matrix  : {args.output}.vcor1\n"
            f"  Variants: {args.output}.vcor1.vars"
        )
    else:
        output_path = args.output + ".txt"
        logger.info(f"Step 9: Saving labeled LD matrix to {output_path}...")
        ld_df = pd.DataFrame(corr_np, index=variants_list, columns=variants_list)
        save_beta_matrix(ld_df, output_path, compress=False)

    # ------------------------------------------------------------------
    # Step 10: Save filtering statistics (optional)
    # ------------------------------------------------------------------
    if args.stats_output and filtering_stats is not None:
        logger.info(f"Step 10: Saving filtering statistics to {args.stats_output}...")
        write_filtering_stats(
            stats=filtering_stats,
            output_path=args.stats_output,
            pval_threshold=pval_threshold,
            buffer_bp=args.buffer or 0,
        )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    logger.info("=" * 60)
    logger.info("LD INFERENCE COMPLETE")
    logger.info("=" * 60)
    logger.info(f"  Traits used         : {beta_df.shape[1]}")
    logger.info(f"  Variants in matrix  : {corr_np.shape[0]}")
    logger.info(f"  Device              : {device}")
    if args.save_binary:
        logger.info(f"  Output (binary)     : {args.output}.bin / .vars")
    elif args.plink2_format:
        logger.info(f"  Output (plink2)     : {args.output}.vcor1 / .vcor1.vars")
    else:
        logger.info(f"  Output (labeled)    : {args.output}.txt")
    if args.stats_output and filtering_stats:
        logger.info(f"  Filtering stats     : {args.stats_output}")
    logger.info("=" * 60)
