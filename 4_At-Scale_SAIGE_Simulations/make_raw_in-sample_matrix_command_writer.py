#!/usr/bin/env python3

"""
make_raw_matrix_command_writer.py

For each text file whose first column contains variants in the form
chrom:pos:ref:alt

this script determines a locus window (min-pos – padding, max-pos + padding)
and writes individual shell scripts that extract those windows from one or more plink
datasets and compute LD matrices with `--r square --unphased`.

Usage
-----
python make_ld_commands.py \
--script_dir /path/to/script_directory \
--script_list_file /path/to/script_list.txt \
--matrix_dir /path/to/ld_matrices \
--plink_roots /path/CEU,/path/YRI \
--populations CEU,YRI \
--variant_lists_file list_of_variant_files.txt \
[--padding 250000]

The file *list_of_variant_files.txt* should contain one path per line.
"""

import argparse
from pathlib import Path
import re
from textwrap import dedent

VAR_RE = re.compile(r"^(\w+):(\d+):")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate individual shell scripts to compute LD matrices per locus."
    )

    parser.add_argument("--script_dir", required=True,
                        help="Directory to write individual shell scripts.")
    parser.add_argument("--script_list_file", required=True,
                        help="Path to write the list of shell script paths.")
    parser.add_argument("--matrix_dir", required=True,
                        help="Directory where LD matrices should be stored.")
    parser.add_argument("--plink_dirs", required=True,
                        help="Comma-delimited directories for each population's PLINK files.")
    parser.add_argument("--populations", required=True,
                        help="Comma-delimited population names in the same order as plink_roots.")
    parser.add_argument("--variant_lists_file", required=True,
                        help="File containing paths to variant text files, one per line.")
    parser.add_argument("--naming_conventions", type=str, default="chr*.plink.maf",
                        help="Comma-delimited patterns for each population's PLINK files.")
    parser.add_argument("--padding", type=int, default=250_000,
                        help="Padding (bp) added upstream & downstream of each locus (default 250 kb).")
    parser.add_argument("--variant-file-suffix", type=str, default="",
                        help="Optional suffix to remove from variant file name when making matrix file (default [empty string]).")

    return parser.parse_args()

def find_plink_root(plink_dir, naming_convention, chrom):
    # Replace * in naming convention with the chromosome number
    pattern = naming_convention.replace("*", str(chrom))
    bed_files = list(Path(plink_dir).glob(pattern + ".bed"))

    if not bed_files:
        raise FileNotFoundError(f"No .bed file found for chromosome {chrom} in {plink_dir} with pattern {pattern}")

    # Remove .bed extension to get the root
    return str(bed_files[0])[:-4]

def get_window(variant_file, padding):
    chrom = None
    min_pos, max_pos = None, None

    with open(variant_file) as fh:
        for line in fh:
            fields = line.strip().split()
            match = VAR_RE.match(fields[0])
            if not match:
                continue

            v_chrom, v_pos = match.group(1), int(match.group(2))

            if chrom is None:
                chrom = v_chrom
            if chrom != v_chrom:
                raise ValueError(f"Multiple chromosomes in {variant_file}")

            min_pos = v_pos if min_pos is None else min(min_pos, v_pos)
            max_pos = v_pos if max_pos is None else max(max_pos, v_pos)

    if chrom is None:
        raise ValueError(f"No variants found in {variant_file}")

    return chrom, max(min_pos - padding, 0), max_pos + padding

def main():
    # Set plink directories
    plink = "/lustre/orion/bif154/proj-shared/mconery/tools/bin/plink"

    # Parse arguments
    args = parse_args()

    plink_dirs = [p.strip() for p in args.plink_dirs.split(",") if p.strip()]
    populations = [p.strip() for p in args.populations.split(",") if p.strip()]

    if len(plink_dirs) != len(populations):
        raise ValueError("plink_roots and populations must have the same length")

    # Parse naming conventions
    temp = [p.strip() for p in args.naming_conventions.split(",") if p.strip()]
    if len(populations) != len(temp) and len(temp) == 1:
        naming_conventions = [temp[0] for x in populations]
    elif len(populations) != len(temp) and len(temp) != 1:
        raise ValueError("There must be either 1 naming convention or 1 per population.")
    else:
        naming_conventions = temp[:]

    # Ensure output dirs exist
    matrix_root = Path(args.matrix_dir)
    script_dir = Path(args.script_dir)
    script_dir.mkdir(parents=True, exist_ok=True)

    for pop in populations:
        (matrix_root / pop).mkdir(parents=True, exist_ok=True)

    script_paths = []
    script_counter = 0

    with open(args.variant_lists_file) as list_fh:
        for variant_file in map(str.strip, list_fh):
            if not variant_file:
                continue

            variant_path = Path(variant_file)
            chrom, start, end = get_window(variant_path, args.padding)
            locus_name = variant_path.stem.rstrip(args.variant_file_suffix)  # base filename without extension

            for plink_dir, pop, naming_conv in zip(plink_dirs, populations, naming_conventions):
                plink_root = find_plink_root(plink_dir, naming_conv, chrom)
                out_prefix = (matrix_root / pop / f"{locus_name}.raw_in-sample")

                cmd = f"{plink} --bfile {plink_root} --chr {chrom} --from-bp {start} --to-bp {end} --keep-allele-order --r square --make-just-bim --out {out_prefix}"

                # Create individual shell script
                script_name = f"make_matrix_{script_counter:06d}_{locus_name}_{pop}.sh"
                script_path = script_dir / script_name

                script_content = dedent(f"""#!/bin/bash
                set -e

                # Generated script for locus {locus_name}, population {pop}
                # Command: {cmd}

                {cmd};
                rm {script_name}
                """)

                with open(script_path, 'w') as script_fh:
                    script_fh.write(script_content)

                # Make script executable
                script_path.chmod(0o755)

                script_paths.append(str(script_path))
                script_counter += 1

    # Write list of all script paths
    with open(args.script_list_file, 'w') as list_fh:
        for script_path in script_paths:
            list_fh.write(script_path + "\n")

    print(f"Created {len(script_paths)} shell scripts in {args.script_dir}")
    print(f"Script list written to {args.script_list_file}")

if __name__ == "__main__":
    main()
