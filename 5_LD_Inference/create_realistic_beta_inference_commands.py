#!/usr/bin/env python3

"""

create_realistic_beta_inference_commands.py

Generate commands for parallel LD matrix inference using realistic traits with p-value filtering.

This script creates a commands file for infer_matrix_with_pval_filter.py that filters traits
based on p-value significance thresholds, allowing the use of actual GWAS traits instead of
simulated null traits for LD matrix inference.

Usage:

python create_realistic_beta_inference_commands.py \
    --variant_dir /path/to/variant/files \
    --output_dir /path/to/output \
    --commands_file /path/to/commands.txt \
    --matrix_script /path/to/infer_matrix_script.py \
    --indexed_dir /path/to/indexed/sumstats \
    --pval_threshold 0.05 \
    [--buffer 100000] \
    [--pval_col pvalue] \
    [--variant_col '#rs_number'] \
    [--beta_col beta] \
    [--device cuda]

"""

import argparse
import glob
import os
import sys

def find_variant_files(variant_dir, pattern="Pheno*raw_in-sample*.vars"):
    """
    Find all variant files matching the specified pattern.

    Args:
        variant_dir: Directory containing variant files
        pattern: Glob pattern to match variant files

    Returns:
        List of variant file paths
    """
    search_pattern = os.path.join(variant_dir, pattern)
    variant_files = glob.glob(search_pattern)
    return sorted(variant_files)

def generate_output_filename(variant_file, output_dir):
    """
    Generate output filename based on input variant file.

    Args:
        variant_file: Path to input variant file
        output_dir: Output directory

    Returns:
        Output file path
    """
    # Extract base filename without path and extension
    base_name = os.path.basename(variant_file)

    # Remove .vars extension
    if base_name.endswith('.vars'):
        base_name = base_name[:-5]

    # Transform filename: remove raw_in-sample.unphased.vcor1 and add realistic.unphased.vcor1
    if '.raw_in-sample.unphased.vcor1' in base_name:
        base_name = base_name.replace('.raw_in-sample.unphased.vcor1', '.realistic.unphased.vcor1')
    elif '.raw_in-sample' in base_name:
        base_name = base_name.replace('.raw_in-sample', '.realistic')

    if not base_name.endswith('.unphased.vcor1'):
        base_name += '.unphased.vcor1'

    return os.path.join(output_dir, base_name)

def create_matrix_inference_command(matrix_script, indexed_dir, variant_file, output_file,
                                   pval_threshold, pval_col='pvalue', variant_col='#rs_number',
                                   beta_col='beta', buffer=None, device='cuda'):
    """
    Create a command string for running the realistic LD matrix inference script.

    Args:
        matrix_script: Path to the LD matrix inference script
        indexed_dir: Directory containing indexed summary statistics
        variant_file: Path to variant file
        output_file: Path to output file
        pval_threshold: P-value threshold for filtering traits
        pval_col: Column name for p-values
        variant_col: Column name for variant IDs
        beta_col: Column name for beta values
        buffer: Optional buffer region size in base pairs
        device: Computation device (cuda/cpu)

    Returns:
        Command string
    """
    cmd_parts = [
        'python', matrix_script,
        '--indexed_dir', indexed_dir,
        '--variant_file', variant_file,
        '--output_file', output_file,
        '--pval_col', pval_col,
        '--pval_threshold', str(pval_threshold),
        '--variant_col', f'{variant_col}',
        '--beta_col', beta_col,
        '--stats_output', output_file.replace("unphased.vcor1", "filtering.txt"),
        '--device', device
    ]

    # Add optional buffer argument if provided
    if buffer is not None:
        cmd_parts.extend(['--buffer', str(buffer)])

    return ' '.join(cmd_parts)

def create_copy_command(variant_file, output_file):
    """
    Create a command string for copying the variant file.

    Args:
        variant_file: Path to source variant file
        output_file: Path to output matrix file (variant file will be copied as output_file.vars)

    Returns:
        Command string
    """
    return f'cp {variant_file} {output_file}.vars'

def main():
    parser = argparse.ArgumentParser(
        description="Generate commands for parallel LD matrix inference using realistic traits with p-value filtering"
    )

    # Required arguments
    parser.add_argument("--variant_dir", type=str, required=True,
                        help="Directory containing variant files")
    parser.add_argument("--output_dir", type=str, required=True,
                        help="Output directory for LD matrices")
    parser.add_argument("--commands_file", type=str, required=True,
                        help="Output file for commands")
    parser.add_argument("--matrix_script", type=str, required=True,
                        help="Path to LD matrix inference script")
    parser.add_argument("--indexed_dir", type=str, required=True,
                        help="Directory containing indexed summary statistics")
    parser.add_argument("--pval_threshold", type=float, required=True,
                        help="P-value threshold for filtering traits")

    # Optional arguments
    parser.add_argument("--buffer", type=int, default=None,
                        help="Genomic distance (bp) for buffer region around matrix variants [default: None]")
    parser.add_argument("--pval_col", type=str, default='pvalue',
                        help="Column name for p-values [default: pvalue]")
    parser.add_argument("--variant_col", type=str, default='#rs_number',
                        help="Column name for variant IDs [default: #rs_number]")
    parser.add_argument("--beta_col", type=str, default='beta',
                        help="Column name for beta values [default: beta]")
    parser.add_argument("--device", type=str, default='cuda', choices=['cuda', 'cpu'],
                        help="Computation device [default: cuda]")
    parser.add_argument("--pattern", type=str, default="Pheno*raw_in-sample*.vars",
                        help="Glob pattern for variant files [default: Pheno*raw_in-sample*.vars]")

    args = parser.parse_args()

    # Validate inputs
    if not os.path.exists(args.variant_dir):
        print(f"Error: Variant directory does not exist: {args.variant_dir}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.indexed_dir):
        print(f"Error: Indexed directory does not exist: {args.indexed_dir}", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.matrix_script):
        print(f"Error: Matrix script does not exist: {args.matrix_script}", file=sys.stderr)
        sys.exit(1)

    if args.pval_threshold < 0 or args.pval_threshold > 1:
        print(f"Error: P-value threshold must be between 0 and 1, got: {args.pval_threshold}", file=sys.stderr)
        sys.exit(1)

    if args.buffer is not None and args.buffer < 0:
        print(f"Error: Buffer size must be non-negative, got: {args.buffer}", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Find variant files
    print(f"Searching for variant files in: {args.variant_dir}")
    print(f"Pattern: {args.pattern}")
    variant_files = find_variant_files(args.variant_dir, args.pattern)

    if not variant_files:
        print(f"Error: No variant files found matching pattern: {args.pattern}", file=sys.stderr)
        print(f"In directory: {args.variant_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(variant_files)} variant files")

    # Check variant col
    if args.variant_col[0] == '#':
        true_variant_col = "\"" + args.variant_col + "\""
    else:
        true_variant_col = args.variant_col

    # Generate commands
    print("Generating commands...")
    commands = []
    skipped_matrix = 0
    skipped_variants = 0

    for variant_file in variant_files:
        output_file = generate_output_filename(variant_file, args.output_dir)
        variant_copy_path = f"{output_file}.vars"

        # Check if matrix file already exists
        if not os.path.exists(output_file):
            matrix_command = create_matrix_inference_command(
                matrix_script=args.matrix_script,
                indexed_dir=args.indexed_dir,
                variant_file=variant_file,
                output_file=output_file,
                pval_threshold=args.pval_threshold,
                pval_col=args.pval_col,
                variant_col=true_variant_col,
                beta_col=args.beta_col,
                buffer=args.buffer,
                device=args.device
            )

            commands.append(matrix_command)
            print(f" Matrix: {os.path.basename(variant_file)} -> {os.path.basename(output_file)}")
        else:
            skipped_matrix += 1
            print(f" Matrix EXISTS: {os.path.basename(output_file)} (skipping)")

        # Check if variant copy already exists
        if not os.path.exists(variant_copy_path):
            copy_command = create_copy_command(variant_file, output_file)
            commands.append(copy_command)
            print(f" Variant: {os.path.basename(variant_file)} -> {os.path.basename(variant_copy_path)}")
        else:
            skipped_variants += 1
            print(f" Variant EXISTS: {os.path.basename(variant_copy_path)} (skipping)")

    # Write commands to file
    if commands:
        print(f"\nWriting {len(commands)} commands to: {args.commands_file}")
        try:
            with open(args.commands_file, 'w') as f:
                for command in commands:
                    f.write(command + '\n')
            print("Commands file created successfully")
        except Exception as e:
            print(f"Error writing commands file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("\nNo commands to write - all files already exist")
        # Still create an empty commands file
        try:
            with open(args.commands_file, 'w') as f:
                f.write("# All files already exist - no commands generated\n")
            print(f"Empty commands file created: {args.commands_file}")
        except Exception as e:
            print(f"Error writing commands file: {e}", file=sys.stderr)
            sys.exit(1)

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total variant files found: {len(variant_files)}")
    print(f"Commands to generate: {len(commands)}")
    print(f"P-value threshold: {args.pval_threshold}")
    if args.buffer is not None:
        print(f"Buffer region: {args.buffer} bp")
    else:
        print(f"Buffer region: None")
    print(f"P-value column: {args.pval_col}")
    print(f"Variant column: {args.variant_col}")
    print(f"Beta column: {args.beta_col}")
    print(f"Device: {args.device}")
    print("=" * 60)

if __name__ == "__main__":
    main()
