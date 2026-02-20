#!/usr/bin/env python3
"""
create_gwama_index_commands.py

Command generator for batch tabix indexing of GWAMA/BOLT summary statistic files.
This script generates individual commands for each file that needs to be indexed,
allowing for flexible job submission and parallelization strategies on HPC systems.

This script identifies all files matching a user-supplied pattern and regex filter,
then creates a list of commands that can be executed individually or submitted
as an array job to SLURM or other job schedulers.

Features:
- Parallelized command generation using all available cores
- Skip files that are already processed (both .gz and .tbi exist)
- Progress reporting during command generation

Usage
-----
python create_gwama_index_commands.py \
  --scratch_dir /path/to/scratch_dir \
  --pattern "/path/to/GWAS/*/Pheno*.meta.out" \
  --commands_file /path/to/commands.txt \
  [--variant_col rs_number] \
  [--variant_id_format CHR:POS:REF:ALT] \
  [--chrom_col CHR] \
  [--pos_col POS] \
  [--filename_regex "^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.meta\\.out$"] \
  [--indexer_script index_single_gwama_file.py] \
  [--force]

Arguments
---------
--scratch_dir         Directory to write .gz and .tbi outputs [REQUIRED]
--pattern             Glob pattern to search for input summary statistics files [REQUIRED]
--commands_file       Output file to write the list of commands [REQUIRED]
--variant_col         Column name holding variant IDs [default: rs_number]
--variant_id_format   Format string for variant IDs [default: CHR:POS:REF:ALT]
--chrom_col           Output column name for chromosome [default: CHR]
--pos_col             Output column name for position [default: POS]
--filename_regex      Regex pattern to filter filenames [optional]
--indexer_script      Path to the single-file indexer script [default: index_single_gwama_file.py]
--force               Force regeneration of commands even if output files exist [optional]

Output
------
Creates a text file with one command per line, where each command calls the
single-file indexer script for one input file.

"""

import os
import glob
import argparse
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

def find_valid_meta_files(pattern, filename_regex=None):
    """
    Find files matching the glob pattern and optionally filter by regex on basename.
    
    Args:
        pattern: Glob pattern to search for files
        filename_regex: Optional regex string to filter filenames (applied to basename only)
        
    Returns:
        List of file paths that match both pattern and regex filter
    """
    # Get all files matching the glob pattern
    all_files = glob.glob(pattern)
    
    # If no regex filter specified, return all files
    if not filename_regex:
        return all_files
    
    # Compile the regex pattern
    try:
        regex_pattern = re.compile(filename_regex)
    except re.error as e:
        raise ValueError(f"Invalid regex pattern '{filename_regex}': {e}")
    
    # Filter files by regex on basename
    filtered_files = [
        f for f in all_files 
        if regex_pattern.match(os.path.basename(f))
    ]
    
    return filtered_files

def check_output_exists(file_path, scratch_dir):
    """
    Check if both the .gz and .tbi files already exist for a given input file.
    
    Args:
        file_path: Path to the input file
        scratch_dir: Directory where output files would be written
        
    Returns:
        bool: True if both output files exist, False otherwise
    """
    input_basename = os.path.basename(file_path)
    output_basename = input_basename + ".gz"
    
    bgzf_path = os.path.join(scratch_dir, output_basename)
    tbi_path = bgzf_path + ".tbi"
    
    return os.path.exists(bgzf_path) and os.path.exists(tbi_path)

def create_single_command(file_path, scratch_dir, variant_col, variant_id_format, 
                         chrom_col, pos_col, indexer_script, force):
    """
    Create indexing command for a single file, checking if output already exists.
    
    Args:
        file_path: Path to input file
        scratch_dir: Output directory for indexed files
        variant_col: Column name for variant IDs
        variant_id_format: Format string for variant IDs
        chrom_col: Output column name for chromosome
        pos_col: Output column name for position
        indexer_script: Path to the single-file indexer script
        force: If True, create command even if output files exist
        
    Returns:
        tuple: (file_path, command_string or None, skip_reason or None)
    """
    # Check if output already exists (unless forced)
    if not force and check_output_exists(file_path, scratch_dir):
        return (file_path, None, "output_exists")
    
    # Build the command for this file
    cmd = f"python {indexer_script}"
    cmd += f" --input_file {file_path}"
    cmd += f" --scratch_dir {scratch_dir}"
    cmd += f" --variant_col {variant_col}"
    cmd += f" --variant_id_format {variant_id_format}"
    cmd += f" --chrom_col {chrom_col}"
    cmd += f" --pos_col {pos_col}"
    
    return (file_path, cmd, None)

def create_indexing_commands_parallel(files, scratch_dir, variant_col, variant_id_format, 
                                    chrom_col, pos_col, indexer_script, force=False):
    """
    Generate indexing commands for each file using parallel processing.
    
    Args:
        files: List of input files to process
        scratch_dir: Output directory for indexed files
        variant_col: Column name for variant IDs
        variant_id_format: Format string for variant IDs
        chrom_col: Output column name for chromosome
        pos_col: Output column name for position
        indexer_script: Path to the single-file indexer script
        force: If True, create commands even if output files exist
        
    Returns:
        tuple: (list_of_commands, skipped_count, skip_details)
    """
    commands = []
    skipped_files = []
    
    # Determine number of worker processes
    max_workers = min(os.cpu_count() or 1, len(files))
    
    print(f"Generating commands for {len(files)} files using {max_workers} processes...")
    
    # Create a partial function with fixed arguments
    process_file = partial(
        create_single_command,
        scratch_dir=scratch_dir,
        variant_col=variant_col,
        variant_id_format=variant_id_format,
        chrom_col=chrom_col,
        pos_col=pos_col,
        indexer_script=indexer_script,
        force=force
    )
    
    # Process files in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_file = {executor.submit(process_file, file_path): file_path for file_path in files}
        
        # Collect results with progress reporting
        completed = 0
        for future in as_completed(future_to_file):
            completed += 1
            if completed % 10000 == 0:  # Progress reporting every 100 files
                print(f"  Processed {completed}/{len(files)} files...")
            
            try:
                file_path, cmd, skip_reason = future.result()
                
                if cmd is not None:
                    commands.append(cmd)
                else:
                    skipped_files.append((file_path, skip_reason))
                    
            except Exception as exc:
                print(f"Error processing file {future_to_file[future]}: {exc}")
    
    print(f"Command generation complete: {len(commands)} commands, {len(skipped_files)} skipped")
    
    return commands, len(skipped_files), skipped_files

def main():
    parser = argparse.ArgumentParser(description="Generate commands for batch tabix indexing of GWAMA summary statistics files.")
    parser.add_argument("--scratch_dir", type=str, required=True,
                        help="Directory where indexed files (.gz and .tbi) will be written.")
    parser.add_argument("--pattern", type=str, required=True,
                        help="Glob pattern to search for GWAMA summary stats files (e.g. '/path/to/Pheno*.meta.out').")
    parser.add_argument("--commands_file", type=str, required=True,
                        help="Output file to write the list of commands.")
    parser.add_argument("--variant_col", type=str, default="rs_number",
                        help="Column name for variant IDs in the files. [default: rs_number]")
    parser.add_argument("--variant_id_format", type=str, default="CHR:POS:REF:ALT",
                        help="Format string describing variant ID structure. [default: 'CHR:POS:REF:ALT']")
    parser.add_argument("--chrom_col", type=str, default="CHR",
                        help="Output column name for chromosome [default: CHR]")
    parser.add_argument("--pos_col", type=str, default="POS",
                        help="Output column name for position [default: POS]")
    parser.add_argument("--filename_regex", type=str, default=None,
                        help="Optional regex pattern to filter filenames (applied to basename only).")
    parser.add_argument("--indexer_script", type=str, default="index_single_gwama_file.py",
                        help="Path to the single-file indexer script [default: index_single_gwama_file.py]")
    parser.add_argument("--force", action="store_true",
                        help="Force regeneration of commands even if output files already exist.")
    
    args = parser.parse_args()
    
    # Find files matching pattern and optional regex filter
    gwama_files = find_valid_meta_files(args.pattern, args.filename_regex)
    
    if not gwama_files:
        print(f"No files found matching pattern: {args.pattern}")
        if args.filename_regex:
            print(f"with filename regex: {args.filename_regex}")
        return
    
    print(f"Found {len(gwama_files)} files matching criteria")
    if args.filename_regex:
        print(f"Pattern: {args.pattern}")
        print(f"Filename regex: {args.filename_regex}")
        print("Sample matched files:")
        for f in gwama_files[:5]:  # Show first 5 files
            print(f"  {f}")
        if len(gwama_files) > 5:
            print(f"  ... and {len(gwama_files) - 5} more files")
    
    # Generate commands in parallel
    commands, skipped_count, skipped_files = create_indexing_commands_parallel(
        gwama_files,
        args.scratch_dir,
        args.variant_col,
        args.variant_id_format,
        args.chrom_col,
        args.pos_col,
        args.indexer_script,
        args.force
    )
    
    # Write commands to file
    if commands:
        with open(args.commands_file, 'w') as f:
            for cmd in commands:
                f.write(cmd + '\n')
    
    # Summary report
    print(f"\n{'='*60}")
    print(f"COMMAND GENERATION SUMMARY")
    print(f"{'='*60}")
    print(f"Total files found: {len(gwama_files)}")
    print(f"Commands generated: {len(commands)}")
    print(f"Files skipped: {skipped_count}")
    
    if skipped_count > 0:
        print(f"\nSkipped files breakdown:")
        skip_reasons = {}
        for _, reason in skipped_files:
            skip_reasons[reason] = skip_reasons.get(reason, 0) + 1
        for reason, count in skip_reasons.items():
            if reason == "output_exists":
                print(f"  Already processed (output exists): {count}")
            else:
                print(f"  {reason}: {count}")
    
    if commands:
        print(f"\nCommands written to: {args.commands_file}")
        print(f"Output directory: {args.scratch_dir}")
        
        # Show example usage
        print(f"\nExample usage:")
        print(f"# Execute commands sequentially:")
        print(f"while IFS= read -r cmd; do eval \"$cmd\"; done < {args.commands_file}")
        print(f"\n# Or submit as SLURM array job:")
        print(f"# sbatch --array=1-{len(commands)}%50 --time=00:30:00 array_job.slurm")
        print(f"# Where array_job.slurm contains:")
        print(f"# COMMAND=$(sed -n \"${{SLURM_ARRAY_TASK_ID}}p\" {args.commands_file})")
        print(f"# eval $COMMAND")
    else:
        print(f"\nNo commands generated. All files may already be processed.")
        print(f"Use --force to regenerate commands for existing output files.")
    
    print(f"{'='*60}")

if __name__ == "__main__":
    main()