"""
betanfer/cli.py

Command-line entry point for Betanfer.

Usage:
    betanfer index --file SUMSTATS [options]
    betanfer infer (--dir DIR | --file-list FILE) [options]
"""

import argparse
import sys


# ---------------------------------------------------------------------------
# Index subcommand
# ---------------------------------------------------------------------------

def _build_index_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "index",
        help="BGZF-compress and tabix-index a single GWAS summary statistics file.",
        description=(
            "Index a single GWAS summary statistics file so it can be queried "
            "efficiently by 'betanfer infer'. Run this command once per file; "
            "scale to many files using a shell loop or an HPC scheduler."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Column specification (one group is used; auto-detection if none given):
  --chrom-col / --pos-col     Use pre-existing chromosome and position columns.
  --variant-col / --variant-id-format
                              Parse CHR and POS from a variant ID column.
  --variant-col / --bim-file  Map variant IDs to positions via a PLINK .bim file.

Examples:
  # Auto-detect column names
  betanfer index --file trait1.txt --output-dir indexed/

  # Explicit chrom/pos columns
  betanfer index --file trait1.txt --chrom-col CHR --pos-col POS

  # Parse positions from variant IDs
  betanfer index --file trait1.txt --variant-col SNP --variant-id-format CHR:POS:REF:ALT

  # Use a .bim file for position mapping
  betanfer index --file trait1.txt --variant-col SNP --bim-file ref.bim

  # Re-index (overwrite existing index)
  betanfer index --file trait1.txt --chrom-col CHR --pos-col POS --force
""",
    )

    p.add_argument(
        "--file",
        metavar="SUMSTATS",
        required=True,
        help="Path to the input GWAS summary statistics file.",
    )
    p.add_argument(
        "--output-dir",
        metavar="DIR",
        default=None,
        help=(
            "Directory for output .gz and .tbi files. "
            "Defaults to the same directory as --file."
        ),
    )
    p.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Re-create the index even if it already exists.",
    )

    # Column specification — chrom/pos mode
    col_group = p.add_argument_group("chromosome/position columns (auto-detected if omitted)")
    col_group.add_argument(
        "--chrom-col",
        metavar="COL",
        default=None,
        help="Column name for chromosome (e.g. CHR, CHROM, chromosome).",
    )
    col_group.add_argument(
        "--pos-col",
        metavar="COL",
        default=None,
        help="Column name for base-pair position (e.g. POS, BP, position).",
    )

    # Column specification — variant ID mode
    vid_group = p.add_argument_group("variant ID-based position extraction")
    vid_group.add_argument(
        "--variant-col",
        metavar="COL",
        default=None,
        help="Column name containing variant IDs (required for --variant-id-format and --bim-file).",
    )
    vid_group.add_argument(
        "--variant-id-format",
        metavar="FMT",
        default=None,
        help=(
            "Format of variant IDs indicating where CHR and POS are encoded, "
            "e.g. 'CHR:POS:REF:ALT'. The delimiter is inferred automatically."
        ),
    )
    vid_group.add_argument(
        "--bim-file",
        metavar="BIM",
        nargs="+",
        default=None,
        help=(
            "One or more PLINK .bim files used to map variant IDs to genomic "
            "positions. Multiple files are merged."
        ),
    )


# ---------------------------------------------------------------------------
# Infer subcommand
# ---------------------------------------------------------------------------

def _build_infer_parser(subparsers: argparse._SubParsersAction) -> None:
    p = subparsers.add_parser(
        "infer",
        help="Infer LD correlations from beta effect sizes across indexed GWAS traits.",
        description=(
            "Compute pairwise LD correlations between variants by correlating "
            "their effect sizes (betas) across many GWAS summary statistics files."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Variant file (contiguous region, GPU auto-detected)
  betanfer infer --dir indexed/ --variant-file variants.txt --output ld.tsv

  # Genomic range
  betanfer infer --dir indexed/ --chrom 1 --start 1000000 --end 2000000 --output ld.tsv

  # P-value filtering (auto-detect pval column; default threshold 5e-8)
  betanfer infer --dir indexed/ --variant-file variants.txt --output ld.tsv \\
      --pval-threshold 5e-8 --buffer 50000 --stats-output filter_stats.tsv

  # Save in binary PLINK-like format (.vcor1 + .vcor1.vars)
  betanfer infer --dir indexed/ --variant-file variants.txt --output prefix --save-binary

  # Non-contiguous variants; force CPU; explicit column names
  betanfer infer --dir indexed/ --variant-file scattered_variants.txt \\
      --output ld.tsv --variant-col SNP --beta-col BETA --non-contiguous --cpu
""",
    )

    # Input source
    src_group = p.add_argument_group("input source (one required)")
    src_ex = src_group.add_mutually_exclusive_group(required=True)
    src_ex.add_argument(
        "--dir",
        metavar="DIR",
        default=None,
        help="Directory containing tabix-indexed (.gz + .tbi) summary statistics files.",
    )
    src_ex.add_argument(
        "--file-list",
        metavar="FILE",
        default=None,
        help=(
            "Text file listing paths to summary statistics files (one per line). "
            "Lines starting with '#' are ignored."
        ),
    )

    # Variant specification
    # Note: mutual exclusivity between --variant-file and --chrom/--start/--end
    # is validated in infer.py because argparse cannot natively express "group A
    # XOR group B" when group B has multiple arguments.
    var_group = p.add_argument_group(
        "variant specification (one required: --variant-file OR --chrom/--start/--end)"
    )
    var_group.add_argument(
        "--variant-file",
        metavar="FILE",
        default=None,
        help=(
            "File containing variant IDs (one per line). "
            "Assumed to span a contiguous genomic region unless --non-contiguous is set."
        ),
    )

    var_group.add_argument(
        "--chrom",
        metavar="CHR",
        default=None,
        help="Chromosome for genomic-range mode (e.g. '1', 'chrX').",
    )
    var_group.add_argument(
        "--start",
        metavar="POS",
        type=int,
        default=None,
        help="Start position (1-based, inclusive) for genomic-range mode.",
    )
    var_group.add_argument(
        "--end",
        metavar="POS",
        type=int,
        default=None,
        help="End position (1-based, inclusive) for genomic-range mode.",
    )

    # Output
    out_group = p.add_argument_group("output")
    out_group.add_argument(
        "--output",
        "-o",
        metavar="FILE",
        required=True,
        help=(
            "Output path. Without --save-binary: a TSV file with variants as "
            "row/column labels. With --save-binary: used as the prefix for "
            ".vcor1 (binary matrix) and .vcor1.vars (variant list) files."
        ),
    )
    out_group.add_argument(
        "--save-binary",
        action="store_true",
        default=False,
        help=(
            "Save output as binary .vcor1 + .vcor1.vars files (using --output as "
            "the prefix) instead of a labeled TSV."
        ),
    )
    out_group.add_argument(
        "--stats-output",
        metavar="FILE",
        default=None,
        help=(
            "Optional: write per-trait p-value filtering statistics to this file "
            "(only relevant when --pval-threshold or --pval-col is set)."
        ),
    )

    # Column names
    col_group = p.add_argument_group("column names (auto-detected from first file if omitted)")
    col_group.add_argument(
        "--variant-col",
        metavar="COL",
        default=None,
        help="Column name for variant IDs in the summary statistics.",
    )
    col_group.add_argument(
        "--variant-id-format",
        metavar="FMT",
        default="CHR:POS:REF:ALT",
        help=(
            "Format of variant IDs for extracting CHR and POS "
            "(e.g. 'CHR:POS:REF:ALT'). [default: CHR:POS:REF:ALT]"
        ),
    )
    col_group.add_argument(
        "--beta-col",
        metavar="COL",
        default=None,
        help="Column name for effect size (beta) in the summary statistics.",
    )

    # P-value filtering
    pval_group = p.add_argument_group("p-value filtering (optional)")
    pval_group.add_argument(
        "--pval-col",
        metavar="COL",
        default=None,
        help=(
            "Column name for p-values. Auto-detected if omitted and "
            "--pval-threshold is set."
        ),
    )
    pval_group.add_argument(
        "--pval-threshold",
        metavar="FLOAT",
        type=float,
        default=None,
        help=(
            "P-value threshold: traits with any variant at p < FLOAT in the "
            "target region are excluded. Setting this flag (or --pval-col) "
            "enables p-value filtering. [default: 5e-8 when filtering is active]"
        ),
    )
    pval_group.add_argument(
        "--buffer",
        metavar="BP",
        type=int,
        default=0,
        help=(
            "Extend the p-value check region by this many base pairs on each "
            "side of the target variants. [default: 0]"
        ),
    )

    # Computation options
    comp_group = p.add_argument_group("computation options")
    comp_group.add_argument(
        "--cpu",
        action="store_true",
        default=False,
        help="Force CPU computation even if a CUDA GPU is available.",
    )
    comp_group.add_argument(
        "--num-processes",
        metavar="N",
        type=int,
        default=None,
        help=(
            "Number of parallel processes for p-value filtering. "
            "[default: all available CPU cores]"
        ),
    )
    comp_group.add_argument(
        "--non-contiguous",
        action="store_true",
        default=False,
        help=(
            "Treat variants in --variant-file as non-contiguous (scattered "
            "across the genome). Uses per-variant tabix queries instead of a "
            "single batch region fetch. Slower but correct for non-contiguous "
            "variant sets."
        ),
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        prog="betanfer",
        description=(
            "Betanfer: infer linkage disequilibrium (LD) from GWAS summary "
            "statistics beta correlations."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Run 'betanfer index --help' or 'betanfer infer --help' for "
            "subcommand-specific options."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.0",
    )

    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")
    subparsers.required = True

    _build_index_parser(subparsers)
    _build_infer_parser(subparsers)

    args = parser.parse_args()

    if args.command == "index":
        from betanfer.index import run_index
        run_index(args)
    elif args.command == "infer":
        # Validate genomic range: all three or none
        if any([args.chrom, args.start, args.end]) and not all(
            [args.chrom, args.start, args.end]
        ):
            parser.error(
                "--chrom, --start, and --end must all be provided together "
                "for genomic-range mode."
            )
        from betanfer.infer import run_infer
        run_infer(args)


if __name__ == "__main__":
    main()
