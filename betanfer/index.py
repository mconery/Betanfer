"""
betanfer/index.py

Index mode: BGZF-compress and tabix-index a single GWAS summary statistics file.
"""

import os
import sys
import logging
from pathlib import Path
from typing import List, Optional

from betanfer.utils import (
    autodetect_chrom_pos_cols,
    validate_column_in_file,
    bgzip_and_tbi_with_chrom_pos_cols,
    bgzip_and_tbi_sumstats_with_variant_ids,
    bgzip_and_tbi_sumstats_with_bim_merge,
    merge_bim_files,
)

logger = logging.getLogger(__name__)


def _detect_column_spec_mode(args) -> str:
    """
    Determine which column-specification mode to use based on provided args.

    Returns one of: 'autodetect', 'chrom_pos', 'variant_id', 'bim_file'
    Raises SystemExit on conflicting arguments.
    """
    has_chrom_pos = bool(args.chrom_col or args.pos_col)
    has_variant_id_fmt = bool(args.variant_id_format)
    has_bim = bool(args.bim_file)
    has_variant_col = bool(args.variant_col)

    # Detect conflicts
    if has_chrom_pos and (has_variant_id_fmt or has_bim):
        logger.error(
            "Error: --chrom-col/--pos-col cannot be combined with "
            "--variant-id-format or --bim-file."
        )
        sys.exit(1)
    if has_variant_id_fmt and has_bim:
        logger.error(
            "Error: --variant-id-format and --bim-file are mutually exclusive."
        )
        sys.exit(1)
    if has_bim and not has_variant_col:
        logger.error("Error: --bim-file requires --variant-col.")
        sys.exit(1)
    if has_chrom_pos and (not args.chrom_col or not args.pos_col):
        logger.error(
            "Error: --chrom-col and --pos-col must be provided together."
        )
        sys.exit(1)

    if has_bim:
        return "bim_file"
    if has_variant_id_fmt or (has_variant_col and not has_chrom_pos):
        return "variant_id"
    if has_chrom_pos:
        return "chrom_pos"
    return "autodetect"


def run_index(args) -> None:
    """Entry point for 'betanfer index'."""
    input_file = args.file

    # Validate input file
    if not os.path.exists(input_file):
        logger.error(f"Error: Input file not found: {input_file}")
        sys.exit(1)

    logger.info(f"Indexing file: {input_file}")

    # Determine output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = str(Path(input_file).parent)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Check if already indexed (and --force not set)
    in_name = Path(input_file).name
    mode_for_skip = _detect_column_spec_mode(args)

    # Determine the expected .tbi path
    # For bim_merge mode, the output name has '.merged' inserted
    if mode_for_skip == "bim_file":
        out_name = in_name if ".merged" in in_name else in_name + ".merged"
        expected_tbi = os.path.join(output_dir, out_name + ".gz.tbi")
    else:
        expected_tbi = os.path.join(output_dir, in_name + ".gz.tbi")

    if os.path.exists(expected_tbi) and not args.force:
        logger.info(
            f"Skipping {input_file}: index already exists at {expected_tbi}. "
            f"Use --force to re-index."
        )
        sys.exit(0)

    # Determine mode (again, in case we fell through the skip)
    mode = _detect_column_spec_mode(args)

    # -----------------------------------------------------------------------
    # Auto-detect mode
    # -----------------------------------------------------------------------
    if mode == "autodetect":
        logger.info("No column specification provided — attempting auto-detection...")
        chrom_col, pos_col = autodetect_chrom_pos_cols(input_file)
        if not chrom_col:
            logger.error(
                "Error: Could not auto-detect chromosome column. "
                "Please specify --chrom-col or use --variant-col with "
                "--variant-id-format or --bim-file."
            )
            sys.exit(1)
        if not pos_col:
            logger.error(
                "Error: Could not auto-detect position column. "
                "Please specify --pos-col or use --variant-col with "
                "--variant-id-format or --bim-file."
            )
            sys.exit(1)
        logger.info(f"Auto-detected: chrom_col='{chrom_col}', pos_col='{pos_col}'")
        mode = "chrom_pos"
        args.chrom_col = chrom_col
        args.pos_col = pos_col

    # -----------------------------------------------------------------------
    # chrom_pos mode — validate and index
    # -----------------------------------------------------------------------
    if mode == "chrom_pos":
        try:
            validate_column_in_file(input_file, args.chrom_col)
            validate_column_in_file(input_file, args.pos_col)
        except ValueError as e:
            logger.error(f"Column validation failed: {e}")
            sys.exit(1)
        logger.info(
            f"Using chrom_col='{args.chrom_col}', pos_col='{args.pos_col}'"
        )
        try:
            gz_path, tbi_path = bgzip_and_tbi_with_chrom_pos_cols(
                input_path=input_file,
                output_dir=output_dir,
                chrom_col=args.chrom_col,
                pos_col=args.pos_col,
                force=args.force,
            )
        except Exception as e:
            logger.error(f"Error indexing {input_file}: {e}")
            sys.exit(1)

    # -----------------------------------------------------------------------
    # variant_id mode — validate and index
    # -----------------------------------------------------------------------
    elif mode == "variant_id":
        variant_col = args.variant_col
        variant_id_format = getattr(args, "variant_id_format", "CHR:POS:REF:ALT") or "CHR:POS:REF:ALT"
        try:
            validate_column_in_file(input_file, variant_col)
        except ValueError as e:
            logger.error(f"Column validation failed: {e}")
            sys.exit(1)
        logger.info(
            f"Using variant_col='{variant_col}', "
            f"variant_id_format='{variant_id_format}'"
        )
        try:
            gz_path, tbi_path = bgzip_and_tbi_sumstats_with_variant_ids(
                input_path=input_file,
                output_dir=output_dir,
                variant_col=variant_col,
                variant_id_format=variant_id_format,
                force=args.force,
            )
        except Exception as e:
            logger.error(f"Error indexing {input_file}: {e}")
            sys.exit(1)

    # -----------------------------------------------------------------------
    # bim_file mode — merge BIM and index
    # -----------------------------------------------------------------------
    elif mode == "bim_file":
        variant_col = args.variant_col
        try:
            validate_column_in_file(input_file, variant_col)
        except ValueError as e:
            logger.error(f"Column validation failed: {e}")
            sys.exit(1)
        logger.info(
            f"Using variant_col='{variant_col}', "
            f"bim_file(s)={args.bim_file}"
        )
        try:
            bim_df = merge_bim_files(args.bim_file)
            gz_path, tbi_path = bgzip_and_tbi_sumstats_with_bim_merge(
                input_path=input_file,
                output_dir=output_dir,
                bim_df=bim_df,
                variant_col_gwama=variant_col,
                force=args.force,
            )
        except Exception as e:
            logger.error(f"Error indexing {input_file}: {e}")
            sys.exit(1)

    # -----------------------------------------------------------------------
    # Success
    # -----------------------------------------------------------------------
    logger.info("Indexing complete.")
    logger.info(f"  BGZF file : {gz_path}")
    logger.info(f"  TBI index : {tbi_path}")
