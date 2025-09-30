#!/usr/bin/env python3

"""
make_raw_matrix_command_writer_parallel.py
MODIFIED FOR MULTIPLE KEEP FILES & SAMPLE COUNT IN SCRIPT NAMES

- Accepts comma-delimited list of --keep_individuals files (one per population)
- Filenames for the scripts include number of samples in each filter file
"""

import argparse
from pathlib import Path
import re
from textwrap import dedent
import multiprocessing as mp
import os

VAR_RE = re.compile(r"^(\w+):(\d+):")

def count_individuals(path):
    """Return the number of lines in a keep file or None if missing."""
    try:
        with open(path) as fh:
            return sum(1 for _ in fh)
    except Exception:
        return None

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate individual shell scripts to compute LD matrices per locus (parallelized)."
    )
    parser.add_argument("--script_dir", required=True, help="Directory to write individual shell scripts.")
    parser.add_argument("--script_list_file", required=True, help="Path to write list of shell script paths.")
    parser.add_argument("--matrix_dir", required=True, help="Directory where LD matrices should be stored.")
    parser.add_argument("--plink_dirs", required=True, help="Comma-delimited directories for each population's PLINK files.")
    parser.add_argument("--populations", required=True, help="Comma-delimited population names in the same order as plink_dirs.")
    parser.add_argument("--variant_lists_file", required=True, help="File containing paths to variant text files, one per line.")
    parser.add_argument("--naming_conventions", type=str, default="chr*.plink.maf", help="Comma-delimited patterns for PLINK files.")
    parser.add_argument("--padding", type=int, default=250_000, help="Padding (bp) upstream & downstream of each locus.")
    parser.add_argument("--variant-file-suffix", type=str, default="", help="Suffix to remove from variant file name.")
    parser.add_argument("--n_processes", type=int, default=0, help="Number of processes (0 = all).")
    parser.add_argument("--dedup_script", type=str, default="remove_matrix_duplicates.py", help="Deduplication helper script.")
    parser.add_argument("--keep_individuals", type=str, default=None,
        help="Comma-separated paths of individuals to keep files, one per population.")
    return parser.parse_args()

def find_plink_root(plink_dir, naming_convention, chrom):
    pattern = naming_convention.replace("*", str(chrom))
    bed_files = list(Path(plink_dir).glob(pattern + ".bed"))
    if not bed_files:
        raise FileNotFoundError(f"No .bed file for chromosome {chrom} in {plink_dir} with pattern {pattern}")
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
    script_tasks, plink_path, matrix_root = script_batch_info
    written_scripts = []
    for script_counter, locus_name, pop, cmd, script_dir, n_samples in script_tasks:
        try:
            script_name = f"make_matrix_{script_counter:06d}_{locus_name}_{pop}_N{safestr(n_samples)}.sh"
            script_path = Path(script_dir) / script_name
            script_content = dedent(f"""#!/bin/bash\n\nset -e\n\n# Generated script for locus {locus_name}, population {pop}, n_individuals={n_samples}\n{cmd}\n""")
            with open(script_path, 'w') as script_fh:
                script_fh.write(script_content)
            script_path.chmod(0o755)
            written_scripts.append(str(script_path))
        except Exception as e:
            print(f"Error writing script {script_counter}: {e}")
            continue
    return written_scripts

def safestr(val):
    return str(val) if val is not None else "NA"

def prepare_script_tasks(variant_files, plink_dirs, populations, naming_conventions, matrix_root, script_dir, plink_path, padding, variant_file_suffix, dedup_script_path, keep_individuals_files=None):
    all_tasks = []
    script_counter = 0
    for variant_file in variant_files:
        if not variant_file.strip():
            continue
        try:
            variant_path = Path(variant_file.strip())
            chrom, start, end = get_window(variant_path, padding)
            locus_name = variant_path.stem.rstrip(variant_file_suffix)
            for idx, (plink_dir, pop, naming_conv) in enumerate(zip(plink_dirs, populations, naming_conventions)):
                plink_root = find_plink_root(plink_dir, naming_conv, chrom)
                out_prefix = (matrix_root / pop / f"{locus_name}.raw_in-sample")
                if os.path.exists(f"{out_prefix}.unphased.vcor1") and os.path.exists(f"{out_prefix}.unphased.vcor1.vars"):
                    continue
                cmd_parts = []
                keep_file = keep_individuals_files[idx] if (keep_individuals_files and idx < len(keep_individuals_files)) else None
                n_samples = count_individuals(keep_file) if keep_file else None
                if keep_file:
                    plink_cmd = f"{plink_path} --bfile {plink_root} --keep {keep_file} --chr {chrom} --from-bp {start} --to-bp {end} --keep-allele-order --r square --make-just-bim --out {out_prefix}"
                else:
                    plink_cmd = f"{plink_path} --bfile {plink_root} --chr {chrom} --from-bp {start} --to-bp {end} --keep-allele-order --r square --make-just-bim --out {out_prefix}"
                cmd_parts.append(plink_cmd)
                cmd_parts.append(f"mv {out_prefix}.ld {out_prefix}.unphased.vcor1")
                cmd_parts.append(f"awk '{{print $2}}' {out_prefix}.bim > {out_prefix}.unphased.vcor1.vars")
                cmd_parts.append(f"rm {out_prefix}.bim")
                cmd_parts.append(f"python {dedup_script_path} {out_prefix}.unphased.vcor1 {out_prefix}.unphased.vcor1.vars")
                cmd = "; ".join(cmd_parts)
                task = (script_counter, locus_name, pop, cmd, str(script_dir), n_samples)
                all_tasks.append(task)
                script_counter += 1
        except Exception as e:
            print(f"Error processing variant file {variant_file}: {e}")
            continue
    return all_tasks

def chunk_tasks(tasks, n_chunks):
    chunk_size = max(1, len(tasks) // n_chunks)
    chunks = []
    for i in range(0, len(tasks), chunk_size):
        chunk = tasks[i:i + chunk_size]
        if chunk:
            chunks.append(chunk)
    return chunks

def main():
    plink = "/lustre/orion/bif154/proj-shared/mconery/tools/bin/plink"
    args = parse_args()
    dedup_script_path = Path(args.dedup_script)
    if not dedup_script_path.is_absolute():
        dedup_script_path = Path.cwd() / dedup_script_path
    if not dedup_script_path.exists():
        raise FileNotFoundError(f"Deduplication script not found: {dedup_script_path}")
    print(f"Using deduplication script: {dedup_script_path}")
    # Multiple keep files are now supported
    keep_files = [None]*len(args.populations.split(","))
    if args.keep_individuals:
        keep_files = [p.strip() for p in args.keep_individuals.split(",") if p.strip()]
        if len(keep_files) != len(args.populations.split(",")):
            raise ValueError("If --keep_individuals is provided, must match population count.")
        for path in keep_files:
            if not Path(path).exists():
                raise FileNotFoundError(f"Individuals file not found: {path}")
        print(f"Using individuals filter files per population: {keep_files}")
    else:
        print("No individual filtering will be applied")
    if args.n_processes <= 0:
        n_processes = mp.cpu_count()
    else:
        n_processes = min(args.n_processes, mp.cpu_count())
    print(f"Using {n_processes} processes for parallel script generation")
    plink_dirs = [p.strip() for p in args.plink_dirs.split(",") if p.strip()]
    populations = [p.strip() for p in args.populations.split(",") if p.strip()]
    if len(plink_dirs) != len(populations):
        raise ValueError("plink_dirs and populations must match in length")
    temp = [p.strip() for p in args.naming_conventions.split(",") if p.strip()]
    if len(populations) != len(temp) and len(temp) == 1:
        naming_conventions = [temp[0] for x in populations]
    elif len(populations) != len(temp) and len(temp) != 1:
        raise ValueError("One naming convention or one per population required")
    else:
        naming_conventions = temp[:]
    matrix_root = Path(args.matrix_dir)
    script_dir = Path(args.script_dir)
    script_dir.mkdir(parents=True, exist_ok=True)
    for pop in populations:
        (matrix_root / pop).mkdir(parents=True, exist_ok=True)
    print("Reading variant files list...")
    with open(args.variant_lists_file) as list_fh:
        variant_files = [line.strip() for line in list_fh if line.strip()]
    print(f"Processing {len(variant_files)} variant files...")
    print("Preparing script tasks...")
    all_tasks = prepare_script_tasks(
        variant_files, plink_dirs, populations, naming_conventions,
        matrix_root, script_dir, plink, args.padding,
        args.variant_file_suffix, str(dedup_script_path), keep_files
    )
    print(f"Prepared {len(all_tasks)} tasks")
    if not all_tasks:
        print("No tasks to process!")
        return
    task_chunks = chunk_tasks(all_tasks, n_processes)
    print(f"Split tasks into {len(task_chunks)} chunks")
    batch_infos = [(chunk, plink, matrix_root) for chunk in task_chunks]
    print("Writing scripts in parallel...")
    all_script_paths = []
    if n_processes == 1:
        for batch_info in batch_infos:
            script_paths = write_script_batch(batch_info)
            all_script_paths.extend(script_paths)
    else:
        with mp.Pool(processes=n_processes) as pool:
            results = pool.map(write_script_batch, batch_infos)
            for script_paths in results:
                all_script_paths.extend(script_paths)
    print(f"Writing script list to {args.script_list_file}...")
    with open(args.script_list_file, 'w') as list_fh:
        for script_path in sorted(all_script_paths):  # Sorted order
            list_fh.write(script_path + "\n")
    print(f"SUCCESS: Created {len(all_script_paths)} shell scripts in {args.script_dir}")
    print(f"Script list written to {args.script_list_file}")
    print(f"Used {n_processes} processes for parallel generation")
    print(f"Each script uses deduplication helper: {dedup_script_path}")
    print(f"Individual filtering files: {keep_files}")

if __name__ == "__main__":
    main()
