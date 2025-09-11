#!/usr/bin/env python3

"""
make_raw_matrix_command_writer_parallel.py

Parallelized version that creates individual shell scripts for LD matrix computation.
Uses multiprocessing to speed up file writing operations across all available CPUs.

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
[--padding 250000] \
[--n_processes 0]  # 0 = use all available CPUs

The file *list_of_variant_files.txt* should contain one path per line.
"""

import argparse
from pathlib import Path
import re
from textwrap import dedent
import multiprocessing as mp
from functools import partial
import os

VAR_RE = re.compile(r"^(\w+):(\d+):")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate individual shell scripts to compute LD matrices per locus (parallelized)."
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
    parser.add_argument("--n_processes", type=int, default=0,
                        help="Number of processes to use (0 = use all available CPUs).")

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

def write_script_batch(script_batch_info):
    """
    Write a batch of scripts. This function will be called by each worker process.

    script_batch_info is a tuple containing:
    (script_tasks, plink_path, matrix_root)

    where script_tasks is a list of tuples:
    (script_counter, locus_name, pop, cmd, script_dir)
    """
    script_tasks, plink_path, matrix_root = script_batch_info
    written_scripts = []

    for script_counter, locus_name, pop, cmd, script_dir in script_tasks:
        try:
            # Create individual shell script
            script_name = f"make_matrix_{script_counter:06d}_{locus_name}_{pop}.sh"
            script_path = Path(script_dir) / script_name

            script_content = dedent(f"""#!/bin/bash
            set -e

            # Generated script for locus {locus_name}, population {pop}
            # Command: {cmd}

            {cmd}
            """)

            with open(script_path, 'w') as script_fh:
                script_fh.write(script_content)

            # Make script executable
            script_path.chmod(0o755)

            written_scripts.append(str(script_path))

        except Exception as e:
            print(f"Error writing script {script_counter}: {e}")
            continue

    return written_scripts

def prepare_script_tasks(variant_files, plink_dirs, populations, naming_conventions, 
                        matrix_root, script_dir, plink_path, padding, variant_file_suffix):
    """
    Prepare all script writing tasks without actually writing them.
    Returns a list of tasks that can be distributed to worker processes.
    """
    all_tasks = []
    script_counter = 0

    for variant_file in variant_files:
        if not variant_file.strip():
            continue

        try:
            variant_path = Path(variant_file.strip())
            chrom, start, end = get_window(variant_path, padding)
            locus_name = variant_path.stem.rstrip(variant_file_suffix)

            for plink_dir, pop, naming_conv in zip(plink_dirs, populations, naming_conventions):
                plink_root = find_plink_root(plink_dir, naming_conv, chrom)
                out_prefix = (matrix_root / pop / f"{locus_name}.raw_in-sample")
                #Check if output files exist and skip if it does
                if os.path.exists(f"{out_prefix}.unphased.vcor1") and os.path.exists(f"{out_prefix}.unphased.vcor1.vars"):
                    continue
                cmd = f"{plink_path} --bfile {plink_root} --chr {chrom} --from-bp {start} --to-bp {end} --keep-allele-order --r square --make-just-bim --out {out_prefix}; "
                cmd = cmd + f'mv {out_prefix}.ld {out_prefix}.unphased.vcor1; awk ' + "\'{print $1}\' " + f'{out_prefix}.bim > {out_prefix}.unphased.vcor1.vars; rm {out_prefix}.bim' 

                # Add task to list
                task = (script_counter, locus_name, pop, cmd, str(script_dir))
                all_tasks.append(task)
                script_counter += 1

        except Exception as e:
            print(f"Error processing variant file {variant_file}: {e}")
            continue

    return all_tasks

def chunk_tasks(tasks, n_chunks):
    """Split tasks into roughly equal chunks for parallel processing."""
    chunk_size = max(1, len(tasks) // n_chunks)
    chunks = []

    for i in range(0, len(tasks), chunk_size):
        chunk = tasks[i:i + chunk_size]
        if chunk:  # Only add non-empty chunks
            chunks.append(chunk)

    return chunks

def main():
    # Set plink path
    plink = "/lustre/orion/bif154/proj-shared/mconery/tools/bin/plink"

    # Parse arguments
    args = parse_args()

    # Determine number of processes
    if args.n_processes <= 0:
        n_processes = mp.cpu_count()
    else:
        n_processes = min(args.n_processes, mp.cpu_count())

    print(f"Using {n_processes} processes for parallel script generation")

    plink_dirs = [p.strip() for p in args.plink_dirs.split(",") if p.strip()]
    populations = [p.strip() for p in args.populations.split(",") if p.strip()]

    if len(plink_dirs) != len(populations):
        raise ValueError("plink_dirs and populations must have the same length")

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

    # Read all variant files
    print("Reading variant files list...")
    with open(args.variant_lists_file) as list_fh:
        variant_files = [line.strip() for line in list_fh if line.strip()]

    print(f"Processing {len(variant_files)} variant files...")

    # Prepare all script writing tasks
    print("Preparing script tasks...")
    all_tasks = prepare_script_tasks(
        variant_files, plink_dirs, populations, naming_conventions,
        matrix_root, script_dir, plink, args.padding, args.variant_file_suffix
    )

    print(f"Prepared {len(all_tasks)} script writing tasks")

    if not all_tasks:
        print("No tasks to process!")
        return

    # Split tasks into chunks for parallel processing
    task_chunks = chunk_tasks(all_tasks, n_processes)
    print(f"Split tasks into {len(task_chunks)} chunks")

    # Create batch info for each chunk (include shared data)
    batch_infos = [(chunk, plink, matrix_root) for chunk in task_chunks]

    # Use multiprocessing to write scripts in parallel
    print("Writing scripts in parallel...")
    all_script_paths = []

    if n_processes == 1:
        # Single process - no multiprocessing overhead
        for batch_info in batch_infos:
            script_paths = write_script_batch(batch_info)
            all_script_paths.extend(script_paths)
    else:
        # Multiple processes
        with mp.Pool(processes=n_processes) as pool:
            results = pool.map(write_script_batch, batch_infos)

            # Collect all script paths
            for script_paths in results:
                all_script_paths.extend(script_paths)

    # Write list of all script paths
    print(f"Writing script list to {args.script_list_file}...")
    with open(args.script_list_file, 'w') as list_fh:
        for script_path in sorted(all_script_paths):  # Sort for consistent ordering
            list_fh.write(script_path + "\n")

    print(f"SUCCESS: Created {len(all_script_paths)} shell scripts in {args.script_dir}")
    print(f"Script list written to {args.script_list_file}")
    print(f"Used {n_processes} processes for parallel generation")

if __name__ == "__main__":
    main()
