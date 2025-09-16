#!/usr/bin/env python3
"""
create_beta_inference_commands.py

Generate commands for parallel LD matrix inference across multiple variant files.
This script creates a commands file that can be executed in parallel using flux.

Usage:
    python create_ld_matrix_commands.py \
        --variant_dir /path/to/variant/files \
        --output_dir /path/to/output \
        --commands_file /path/to/commands.txt \
        --matrix_script /path/to/infer_matrix_script.py \
        --indexed_dir /path/to/indexed/sumstats \
        [--n_traits 2000] \
        [--variant_col '#rs_number'] \
        [--beta_col beta] \
        [--device cuda]
"""

import argparse
import glob
import os
import sys
from pathlib import Path


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
    
    # Transform filename: remove raw_in-sample.unphased.vcor1 and add inferred_2000-trait.unphased.vcor1
    if '.raw_in-sample.unphased.vcor1' in base_name:
        base_name = base_name.replace('.raw_in-sample.unphased.vcor1', '.inferred_2000-trait.unphased.vcor1')
    elif '.raw_in-sample' in base_name:
        base_name = base_name.replace('.raw_in-sample', '.inferred_2000-trait')
        if not base_name.endswith('.unphased.vcor1'):
            base_name += '.unphased.vcor1'
    else:
        # Fallback: just append the suffix
        base_name += '.inferred_2000-trait.unphased.vcor1'
    
    return os.path.join(output_dir, base_name)


def create_command(matrix_script, indexed_dir, variant_file, output_file, 
                  n_traits=2000, variant_col='#rs_number', beta_col='beta', device='cuda'):
    """
    Create a command string for running the LD matrix inference script.
    
    Args:
        matrix_script: Path to the LD matrix inference script
        indexed_dir: Directory containing indexed summary statistics
        variant_file: Path to variant file
        output_file: Path to output file
        n_traits: Number of traits to sample
        variant_col: Column name for variant IDs
        beta_col: Column name for beta values
        device: Computation device (cuda/cpu)
        
    Returns:
        Command string
    """
    cmd_parts = [
        'python', matrix_script,
        '--indexed_dir', indexed_dir,
        '--variant_file', variant_file,
        '--output_file', output_file,
        '--n_traits', str(n_traits),
        '--variant_col', f"'{variant_col}'",
        '--beta_col', beta_col,
        '--device', device,
        '&&',
        'cp', variant_file, f"{output_file}.vars"
    ]
    
    return ' '.join(cmd_parts)


def main():
    parser = argparse.ArgumentParser(
        description="Generate commands for parallel LD matrix inference"
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
    
    # Optional arguments
    parser.add_argument("--n_traits", type=int, default=2000,
                       help="Number of traits to sample [default: 2000]")
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
    
    # Generate commands
    print(f"Generating commands...")
    
    commands = []
    for variant_file in variant_files:
        output_file = generate_output_filename(variant_file, args.output_dir)
        
        command = create_command(
            matrix_script=args.matrix_script,
            indexed_dir=args.indexed_dir,
            variant_file=variant_file,
            output_file=output_file,
            n_traits=args.n_traits,
            variant_col=args.variant_col,
            beta_col=args.beta_col,
            device=args.device
        )
        
        commands.append(command)
        
        print(f"  {os.path.basename(variant_file)} -> {os.path.basename(output_file)}")
    
    # Write commands to file
    print(f"Writing {len(commands)} commands to: {args.commands_file}")
    
    try:
        with open(args.commands_file, 'w') as f:
            for command in commands:
                f.write(command + '\n')
        
        print(f"Commands file created successfully")
        
    except Exception as e:
        print(f"Error writing commands file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Summary
    print(f"\nSUMMARY:")
    print(f"  Variant files found: {len(variant_files)}")
    print(f"  Commands generated: {len(commands)}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Commands file: {args.commands_file}")
    print(f"  N traits: {args.n_traits}")
    print(f"  Device: {args.device}")


if __name__ == "__main__":
    main()