#!/usr/bin/env python3
"""
Script to randomly sample individuals from a plink dataset.

This script reads a plink .bim file, randomly samples a specified number of individuals,
and outputs a list that can be used to filter plink files.

Usage:
    python sample_individuals.py <bim_file> <num_samples> <output_file> [--seed SEED]

Arguments:
    bim_file: Path to the plink .bim file
    num_samples: Number of individuals to randomly sample
    output_file: Path to write the list of sampled individuals
    --seed: Random seed (default: 5)
"""

import argparse
import random
import sys
import os

def read_bim_file(bim_file):
    """
    Read a plink .bim file and extract unique individual IDs.

    Note: .bim files contain SNP information, not individual information.
    For individual information, we typically need .fam files.
    However, if the requirement is specifically for .bim files,
    this function can be adapted accordingly.

    Args:
        bim_file (str): Path to the .bim file

    Returns:
        list: List of individual IDs
    """
    individuals = set()

    try:
        with open(bim_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    # .bim format: chromosome, variant_id, genetic_distance, base_pair_position, allele1, allele2
                    # Note: .bim files don't typically contain individual information
                    # This might need to be a .fam file instead
                    parts = line.split()
                    if len(parts) >= 2:
                        # For demonstration, using variant_id as identifier
                        # In practice, you'd likely want to use a .fam file for individuals
                        individuals.add(parts[1])

    except FileNotFoundError:
        print(f"Error: File '{bim_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{bim_file}': {e}", file=sys.stderr)
        sys.exit(1)

    return list(individuals)

def read_fam_file(fam_file):
    """
    Read a plink .fam file and extract individual IDs.

    Args:
        fam_file (str): Path to the .fam file

    Returns:
        list: List of tuples (family_id, individual_id)
    """
    individuals = []

    try:
        with open(fam_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    # .fam format: family_id, individual_id, paternal_id, maternal_id, sex, phenotype
                    parts = line.split()
                    if len(parts) >= 2:
                        family_id = parts[0]
                        individual_id = parts[1]
                        individuals.append((family_id, individual_id))

    except FileNotFoundError:
        print(f"Error: File '{fam_file}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{fam_file}': {e}", file=sys.stderr)
        sys.exit(1)

    return individuals

def sample_individuals(individuals, num_samples, seed=5):
    """
    Randomly sample individuals from the list.

    Args:
        individuals (list): List of individuals
        num_samples (int): Number of individuals to sample
        seed (int): Random seed

    Returns:
        list: Randomly sampled individuals
    """
    if num_samples > len(individuals):
        print(f"Warning: Requested {num_samples} samples, but only {len(individuals)} individuals available.", file=sys.stderr)
        print(f"Returning all {len(individuals)} individuals.", file=sys.stderr)
        return individuals

    random.seed(seed)
    return random.sample(individuals, num_samples)

def write_sample_file(sampled_individuals, output_file, is_fam_format=False):
    """
    Write sampled individuals to output file in plink-compatible format.

    Args:
        sampled_individuals (list): List of sampled individuals
        output_file (str): Path to output file
        is_fam_format (bool): Whether the input was from a .fam file
    """
    try:
        with open(output_file, 'w') as f:
            if is_fam_format:
                # Write family_id and individual_id for plink --keep format
                for family_id, individual_id in sampled_individuals:
                    f.write(f"{family_id}\t{individual_id}\n")
            else:
                # Write individual IDs only
                for individual in sampled_individuals:
                    f.write(f"{individual}\n")

        print(f"Successfully wrote {len(sampled_individuals)} individuals to '{output_file}'")

    except Exception as e:
        print(f"Error writing to file '{output_file}': {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description="Randomly sample individuals from a plink dataset",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python sample_individuals.py data.bim 100 sampled_individuals.txt
    python sample_individuals.py data.fam 50 sampled_individuals.txt --seed 42

Note: This script works with both .bim and .fam files. For individual sampling,
.fam files are more appropriate as they contain individual information.
The output format is compatible with plink's --keep option.
        """
    )

    parser.add_argument('input_file', help='Path to the plink .bim or .fam file')
    parser.add_argument('num_samples', type=int, help='Number of individuals to randomly sample')
    parser.add_argument('output_file', help='Path to write the list of sampled individuals')
    parser.add_argument('--seed', type=int, default=5, help='Random seed (default: 5)')

    args = parser.parse_args()

    # Validate arguments
    if args.num_samples <= 0:
        print("Error: Number of samples must be positive.", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Determine file type and read accordingly
    is_fam_format = args.input_file.lower().endswith('.fam')
    is_bim_format = args.input_file.lower().endswith('.bim')

    if is_fam_format:
        print(f"Reading individuals from .fam file: {args.input_file}")
        individuals = read_fam_file(args.input_file)
        print(f"Found {len(individuals)} individuals")
    elif is_bim_format:
        print(f"Warning: .bim files contain SNP information, not individual information.")
        print(f"Consider using a .fam file instead for individual sampling.")
        print(f"Proceeding with .bim file: {args.input_file}")
        individuals = read_bim_file(args.input_file)
        print(f"Found {len(individuals)} entries")
        is_fam_format = False
    else:
        # Try to determine format by content
        print(f"File extension not recognized. Attempting to read as .fam file first...")
        try:
            individuals = read_fam_file(args.input_file)
            is_fam_format = True
            print(f"Successfully read as .fam format. Found {len(individuals)} individuals")
        except:
            print(f"Failed to read as .fam. Attempting to read as .bim format...")
            individuals = read_bim_file(args.input_file)
            is_fam_format = False
            print(f"Read as .bim format. Found {len(individuals)} entries")

    # Sample individuals
    print(f"Sampling {args.num_samples} individuals with seed {args.seed}")
    sampled = sample_individuals(individuals, args.num_samples, args.seed)

    # Write output file
    write_sample_file(sampled, args.output_file, is_fam_format)

    print(f"Output file '{args.output_file}' can be used with plink --keep option")

if __name__ == "__main__":
    main()
