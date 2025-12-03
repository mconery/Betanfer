#!/usr/bin/env python3
"""
Generate filtering commands for GWAS summary statistics.

This script scans GWAS summary statistics files, matches them to corresponding
causal variant files, checks if output files already exist, and generates
filtering commands for new files only.

Output: A command file that can be read by mapfile and submitted to flux.
"""

import argparse
import re
from pathlib import Path
from typing import List, Tuple, Optional
import sys


def extract_pheno_number(filename: str) -> Optional[str]:
    """
    Extract phenotype number from summary statistics filename.
    
    Parameters
    ----------
    filename : str
        Summary statistics filename
        
    Returns
    -------
    str or None
        Phenotype number (e.g., "1.2") or None if not found
        
    Examples
    --------
    >>> extract_pheno_number("Pheno1.2_100_causal.heritability_0.5.binary.chr1.step2.txt")
    '1.2'
    """
    pattern = r"Pheno([0-9]+\.[0-9]+)_[0-9]+_causal"
    match = re.search(pattern, filename)
    return match.group(1) if match else None


def find_causal_file(pheno_num: str, causal_dir: Path) -> Optional[Path]:
    """
    Find causal variants file matching phenotype number.
    
    Parameters
    ----------
    pheno_num : str
        Phenotype number (e.g., "1.2")
    causal_dir : Path
        Directory containing causal variant files
        
    Returns
    -------
    Path or None
        Path to matching causal file or None if not found
    """
    pattern = f"Pheno{pheno_num}_*_causal.variant_set.txt"
    matching_files = list(causal_dir.glob(pattern))
    
    if not matching_files:
        return None
    
    # Return first match
    return matching_files[0]


def generate_output_filename(input_filename: str) -> str:
    """
    Generate output filename from input filename.
    
    Parameters
    ----------
    input_filename : str
        Input summary statistics filename
        
    Returns
    -------
    str
        Output filename with "drop_causal" inserted before ".step2.txt"
        
    Examples
    --------
    >>> generate_output_filename("Pheno1.2_100_causal.heritability_0.5.binary.chr1.step2.txt")
    'Pheno1.2_100_causal.heritability_0.5.binary.chr1.drop_causal.step2.txt'
    """
    return re.sub(r"\.step2\.txt$", ".drop_causal.step2.txt", input_filename)


def build_filter_command(
    input_file: Path,
    output_file: Path,
    causal_file: Path,
    python_script: Path,
    marker_col: str = "MarkerID",
    chr_col: str = "CHR",
    pos_col: str = "POS"
) -> str:
    """
    Build filtering command.
    
    Parameters
    ----------
    input_file : Path
        Path to input summary statistics file
    output_file : Path
        Path to output filtered file
    causal_file : Path
        Path to causal variants file
    python_script : Path
        Path to drop_variants.py script
    marker_col : str
        Variant ID column name (default: MarkerID)
    chr_col : str
        Chromosome column name (default: CHR)
    pos_col : str
        Position column name (default: POS)
        
    Returns
    -------
    str
        Complete filtering command
    """
    cmd = (
        f"python {python_script} "
        f"-i {input_file} "
        f"-o {output_file} "
        f"-c "
        f"-C {causal_file} "
        f"-m {marker_col} "
        f"--chr-col {chr_col} "
        f"--pos-col {pos_col}"
    )
    return cmd


def generate_commands(
    results_dir: Path,
    causal_dir: Path,
    python_script: Path,
    marker_col: str = "MarkerID",
    chr_col: str = "CHR",
    pos_col: str = "POS",
    verbose: bool = False
) -> Tuple[List[str], int, int, int]:
    """
    Generate filtering commands for all summary statistics files.
    
    Parameters
    ----------
    results_dir : Path
        Root directory containing AFR and EUR subdirectories
    causal_dir : Path
        Directory containing causal variant files
    python_script : Path
        Path to drop_variants.py script
    marker_col : str
        Variant ID column name (default: MarkerID)
    chr_col : str
        Chromosome column name (default: CHR)
    pos_col : str
        Position column name (default: POS)
    verbose : bool
        Print detailed logging (default: False)
        
    Returns
    -------
    tuple
        (commands list, total files, skipped files, generated commands)
    """
    commands = []
    total_files = 0
    skipped_files = 0
    
    # Find all .step2.txt files (excluding already processed ones)
    results_dir_path = Path(results_dir)
    
    for sumstats_file in sorted(results_dir_path.rglob("*.step2.txt")):
        # Skip already processed files
        if "drop_causal" in sumstats_file.name:
            continue
        
        total_files += 1
        
        filename = sumstats_file.name
        dir_path = sumstats_file.parent
        
        # Extract phenotype number
        pheno_num = extract_pheno_number(filename)
        if not pheno_num:
            if verbose:
                print(f"[WARN] Could not extract phenotype from: {filename}", file=sys.stderr)
            continue
        
        # Find causal file
        causal_file = find_causal_file(pheno_num, Path(causal_dir))
        if not causal_file:
            if verbose:
                print(f"[WARN] Causal file not found for phenotype {pheno_num}", file=sys.stderr)
            continue
        
        # Generate output filename
        output_filename = generate_output_filename(filename)
        output_file = dir_path / output_filename
        
        # Check if output file already exists
        if output_file.exists():
            if verbose:
                print(f"[SKIP] Output already exists: {output_filename}")
            skipped_files += 1
            continue
        
        # Build command
        cmd = build_filter_command(
            input_file=sumstats_file,
            output_file=output_file,
            causal_file=causal_file,
            python_script=python_script,
            marker_col=marker_col,
            chr_col=chr_col,
            pos_col=pos_col
        )
        
        commands.append(cmd)
        
        if verbose:
            print(f"[ADD] {filename} -> {output_filename}")
    
    return commands, total_files, skipped_files, len(commands)


def main():
    parser = argparse.ArgumentParser(
        description="Generate filtering commands for GWAS summary statistics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate commands with default settings
  %(prog)s -r /path/to/results -c /path/to/causal -p /path/to/drop_variants.py -o commands.txt
  
  # Generate commands with custom column names
  %(prog)s -r /path/to/results -c /path/to/causal -p /path/to/drop_variants.py \\
    -o commands.txt -m SNP_ID --chr-col CHR --pos-col BP
  
  # Generate with verbose output
  %(prog)s -r /path/to/results -c /path/to/causal -p /path/to/drop_variants.py \\
    -o commands.txt -v
        """
    )
    
    parser.add_argument(
        "-r", "--results-dir",
        required=True,
        help="Root directory containing AFR and EUR subdirectories with summary statistics"
    )
    parser.add_argument(
        "-c", "--causal-dir",
        required=True,
        help="Directory containing causal variant files"
    )
    parser.add_argument(
        "-p", "--python-script",
        required=True,
        help="Path to drop_variants.py script"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file to write commands to"
    )
    parser.add_argument(
        "-m", "--marker-col",
        default="MarkerID",
        help="Variant ID column name (default: MarkerID)"
    )
    parser.add_argument(
        "--chr-col",
        default="CHR",
        help="Chromosome column name (default: CHR)"
    )
    parser.add_argument(
        "--pos-col",
        default="POS",
        help="Position column name (default: POS)"
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print verbose logging output"
    )
    
    args = parser.parse_args()
    
    # Validate paths
    results_dir = Path(args.results_dir)
    if not results_dir.is_dir():
        print(f"Error: Results directory not found: {results_dir}", file=sys.stderr)
        sys.exit(1)
    
    causal_dir = Path(args.causal_dir)
    if not causal_dir.is_dir():
        print(f"Error: Causal directory not found: {causal_dir}", file=sys.stderr)
        sys.exit(1)
    
    python_script = Path(args.python_script)
    if not python_script.is_file():
        print(f"Error: Python script not found: {python_script}", file=sys.stderr)
        sys.exit(1)
    
    # Generate commands
    if args.verbose:
        print(f"[INFO] Scanning for summary statistics files...")
        print(f"[INFO] Results directory: {results_dir}")
        print(f"[INFO] Causal directory: {causal_dir}")
        print(f"[INFO] Python script: {python_script}")
    
    commands, total, skipped, generated = generate_commands(
        results_dir=results_dir,
        causal_dir=causal_dir,
        python_script=python_script,
        marker_col=args.marker_col,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        verbose=args.verbose
    )
    
    if args.verbose:
        print(f"[INFO] Total files found: {total}")
        print(f"[INFO] Already processed: {skipped}")
        print(f"[INFO] New commands generated: {generated}")
    
    # Write output file
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        for cmd in commands:
            f.write(cmd + '\n')
    
    print(f"Commands written to: {output_path}")
    print(f"Total: {total} | Skipped: {skipped} | Commands: {generated}")
    
    if generated == 0:
        print("No new commands to process.", file=sys.stderr)
        sys.exit(0)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()
