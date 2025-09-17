#!/usr/bin/env python3

"""
create_saige_index_commands.py

Command generator for batch tabix indexing of SAIGE summary statistic files.

This script generates individual commands for each file that needs to be indexed,
allowing for flexible job submission and parallelization strategies on HPC systems.

This script identifies all files matching a user-supplied pattern and regex filter,
then creates a list of commands that can be executed individually or submitted
as an array job to SLURM or other job schedulers.

Features:
- Parallelized command generation using all available cores
- Skip files that are already processed (both .gz and .tbi exist)
- Progress reporting during command generation
- Support for EUR and AFR population-specific directories

Usage
-----
python create_saige_index_commands.py \
--scratch_dir /path/to/scratch_dir \
--pattern "/path/to/SAIGE/results/*/Pheno*.txt" \
--commands_file /path/to/commands.txt \
[--chrom_col CHR] \
[--pos_col POS] \
[--filename_regex "^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.step2\\.txt$"] \
[--indexer_script index_single_saige_file.py] \
[--force]

Arguments
---------
--scratch_dir       Directory to write .gz and .tbi outputs [REQUIRED]
--pattern           Glob pattern to search for input summary statistics files [REQUIRED]
--commands_file     Output file to write the list of commands [REQUIRED]
--chrom_col         Column name for chromosome [default: CHR]
--pos_col           Column name for position [default: POS]
--filename_regex    Regex pattern to filter filenames [optional]
--indexer_script    Path to the single-file indexer script [default: index_single_saige_file.py]
--force             Force regeneration of commands even if output files exist [optional]
--zero_based        Add --zero_based flag to indexer commands if positions are 0-based [optional]

Output
------
Creates a text file with one command per line, where each command calls the
single-file indexer script for one input file.

Examples
--------
# For EUR population
python create_saige_index_commands.py \
--scratch_dir /lustre/orion/bif154/scratch/mconery/indexed_EUR_sumstats \
--pattern "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/EUR/Pheno*.txt" \
--commands_file EUR_indexing_commands.txt

# For AFR population
python create_saige_index_commands.py \
--scratch_dir /lustre/orion/bif154/scratch/mconery/indexed_AFR_sumstats \
--pattern "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/step2_results/AFR/Pheno*.txt" \
--commands_file AFR_indexing_commands.txt
"""

import os
import glob
import argparse
import re
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

def find_valid_saige_files(pattern, filename_regex=None):
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

def create_single_command(file_path, scratch_dir, chrom_col, pos_col, 
                         indexer_script, force, zero_based):
    """
    Create indexing command for a single file, checking if output already exists.
    
    Args:
        file_path: Path to input file
        scratch_dir: Output directory for indexed files
        chrom_col: Column name for chromosome
        pos_col: Column name for position
        indexer_script: Path to the single-file indexer script
        force: If True, create command even if output files exist
        zero_based: If True, add --zero_based flag to command
    
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
    cmd += f" --chrom_col {chrom_col}"
    cmd += f" --pos_col {pos_col}"
    
    if zero_based:
        cmd += " --zero_based"
    
    return (file_path, cmd, None)

def create_indexing_commands_parallel(files, scratch_dir, chrom_col, pos_col, 
                                    indexer_script, force=False, zero_based=False):
    """
    Generate indexing commands for each file using parallel processing.
    
    Args:
        files: List of input files to process
        scratch_dir: Output directory for indexed files
        chrom_col: Column name for chromosome
        pos_col: Column name for position
        indexer_script: Path to the single-file indexer script
        force: If True, create commands even if output files exist
        zero_based: If True, add --zero_based flag to commands
    
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
        chrom_col=chrom_col,
        pos_col=pos_col,
        indexer_script=indexer_script,
        force=force,
        zero_based=zero_based
    )
    
    # Process files in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_file = {executor.submit(process_file, file_path): file_path for file_path in files}
        
        # Collect results with progress reporting
        completed = 0
        for future in as_completed(future_to_file):
            completed += 1
            if completed % 100 == 0:  # Progress reporting every 100 files
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
    parser = argparse.ArgumentParser(description="Generate commands for batch tabix indexing of SAIGE summary statistics files.")
    parser.add_argument("--scratch_dir", type=str, required=True,
                        help="Directory where indexed files (.gz and .tbi) will be written.")
    parser.add_argument("--pattern", type=str, required=True,
                        help="Glob pattern to search for SAIGE summary stats files (e.g. '/path/to/Pheno*.txt').")
    parser.add_argument("--commands_file", type=str, required=True,
                        help="Output file to write the list of commands.")
    parser.add_argument("--chrom_col", type=str, default="CHR",
                        help="Column name for chromosome in the files. [default: CHR]")
    parser.add_argument("--pos_col", type=str, default="POS",
                        help="Column name for position in the files. [default: POS]")
    parser.add_argument("--filename_regex", type=str, default="^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.step2\\.txt$",
                        help="Optional regex pattern to filter filenames (applied to basename only). [default: '^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.step2\\.txt$']")
    parser.add_argument("--indexer_script", type=str, default="index_single_saige_file.py",
                        help="Path to the single-file indexer script [default: index_single_saige_file.py]")
    parser.add_argument("--force", action="store_true",
                        help="Force regeneration of commands even if output files already exist.")
    parser.add_argument("--zero_based", action="store_true",
                        help="Add --zero_based flag to indexer commands (for 0-based positions).")
    
    args = parser.parse_args()
    
    # Find files matching pattern and optional regex filter
    saige_files = find_valid_saige_files(args.pattern, args.filename_regex)
    
    if not saige_files:
        print(f"No files found matching pattern: {args.pattern}")
        if args.filename_regex:
            print(f"with filename regex: {args.filename_regex}")
        return
    
    print(f"Found {len(saige_files)} files matching criteria")
    if args.filename_regex:
        print(f"Pattern: {args.pattern}")
        print(f"Filename regex: {args.filename_regex}")
    
    print("Sample matched files:")
    for f in saige_files[:5]:  # Show first 5 files
        print(f"  {f}")
    if len(saige_files) > 5:
        print(f"  ... and {len(saige_files) - 5} more files")
    
    # Generate commands in parallel
    commands, skipped_count, skipped_files = create_indexing_commands_parallel(
        saige_files,
        args.scratch_dir,
        args.chrom_col,
        args.pos_col,
        args.indexer_script,
        args.force,
        args.zero_based
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
    print(f"Total files found: {len(saige_files)}")
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