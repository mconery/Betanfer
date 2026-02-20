#!/usr/bin/env python3
"""
Fast parallel batch file splitter for HPC systems.
Splits a large command file into smaller batch files using multiprocessing.

Usage:
    python split_commands_parallel.py -i commands.txt -o batch_dir/ -b 50 -j 8
    
    -i, --input       Input commands file
    -o, --output      Output directory for batch files
    -b, --batch-size  Number of commands per batch file (default: 50)
    -j, --jobs        Number of parallel processes (default: CPU count)
    -p, --prefix      Prefix for batch filenames (default: batch_)
    -s, --suffix      Suffix for batch filenames (default: .sh)
"""

import argparse
import os
import sys
from pathlib import Path
from multiprocessing import Pool, Manager
import time
from collections import defaultdict


def calculate_batch_size(total_commands, num_batches): 
    """Calculate batch size from number of batches.""" 
    if num_batches <= 0: 
        raise ValueError("Number of batches must be > 0") 
    if total_commands == 0: 
        raise ValueError("No commands to process") 
    import math 
    return math.ceil(total_commands / num_batches)

def read_commands_chunked(input_file, batch_size, num_processes):
    """
    Read commands file and group into batches.
    Returns list of (batch_num, commands_list) tuples.
    Automatically filters out comments and empty lines.
    """
    batches = []
    current_batch = []
    batch_num = 0
    
    try:
        with open(input_file, 'r') as f:
            for line in f:
                # Skip comments and empty lines
                line = line.rstrip('\n')
                if not line or line.startswith('#'):
                    continue
                
                current_batch.append(line)
                
                if len(current_batch) >= batch_size:
                    batches.append((batch_num, current_batch))
                    current_batch = []
                    batch_num += 1
        
        # Add remaining commands
        if current_batch:
            batches.append((batch_num, current_batch))
    
    except IOError as e:
        print(f"ERROR: Failed to read input file '{input_file}': {e}", file=sys.stderr)
        sys.exit(1)
    
    return batches


def write_batch_file(args):
    """
    Write a single batch file. Designed to run in parallel.
    args: (batch_num, commands, output_dir, prefix, suffix, padding)
    """
    batch_num, commands, output_dir, prefix, suffix, padding = args
    
    try:
        # Format batch filename with zero-padding
        batch_filename = f"{prefix}{batch_num}{suffix}"
        batch_path = os.path.join(output_dir, batch_filename)
        
        # Write batch file
        with open(batch_path, 'w') as f:
            for cmd in commands:
                f.write(cmd + '\n')
        
        return (batch_num, batch_path, True, None)
    
    except Exception as e:
        return (batch_num, None, False, str(e))


def main():
    parser = argparse.ArgumentParser(
        description='Fast parallel batch file splitter for HPC systems',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help='Input commands file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for batch files')
    parser.add_argument('-b', '--batch-size', type=int, default=None, 
                        help='Number of commands per batch (mutually exclusive with --num-batches)')
    parser.add_argument('-n', '--num-batches', type=int, default=None, 
                        help='Number of batches to create (mutually exclusive with --batch-size)')
    parser.add_argument('-j', '--jobs', type=int, default=None,
                        help='Number of parallel processes (default: CPU count)')
    parser.add_argument('-p', '--prefix', default='batch_',
                        help='Prefix for batch filenames (default: batch_)')
    parser.add_argument('-s', '--suffix', default='.sh',
                        help='Suffix for batch filenames (default: .sh)')
    parser.add_argument('--verbose', action='store_true',
                        help='Print verbose output')
    
    args = parser.parse_args()
    
    # Validate that exactly one of batch_size or num_batches is provided
    if args.batch_size is None and args.num_batches is None:
        parser.error("ERROR: Must specify either --batch-size or --num-batches")
    if args.batch_size is not None and args.num_batches is not None:
        parser.error("ERROR: Cannot specify both --batch-size and --num-batches")
    
    # Validate inputs
    if not os.path.isfile(args.input):
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output)
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(f"ERROR: Failed to create output directory '{args.output}': {e}", file=sys.stderr)
        sys.exit(1)

    # Determine number of processes
    num_processes = args.jobs if args.jobs else os.cpu_count()
    if args.verbose:
        print(f"Input file: {args.input}")
        print(f"Output directory: {output_dir}")
        print(f"Parallel processes: {num_processes}")
        print()

    # Read and group commands into batches
    print("Reading and grouping commands...", end='', flush=True)
    start_time = time.time()
    # First pass: read all commands to count them if using --num-batches
    if args.num_batches is not None:
        # Count total commands first
        total_commands = 0
        try:
            with open(args.input, 'r') as f:
                for line in f:
                    line = line.rstrip('\n')
                    if not line or line.startswith('#'):
                        continue
                    total_commands += 1
        except IOError as e:
            print(f"\nERROR: Failed to read input file '{args.input}': {e}", file=sys.stderr)
            sys.exit(1)
        
        # Calculate batch size from number of batches
        try:
            args.batch_size = calculate_batch_size(total_commands, args.num_batches)
            if args.verbose:
                print(f"\nCalculated batch size: {args.batch_size} commands per batch")
                print(f"Total commands to process: {total_commands}")
                print(f"Number of batches: {args.num_batches}\n")
        except ValueError as e:
            print(f"\nERROR: {e}", file=sys.stderr)
            sys.exit(1)
    
    batches = read_commands_chunked(args.input, args.batch_size, num_processes)
    read_time = time.time() - start_time
    print(f" Done ({read_time:.2f}s)")
    total_commands = sum(len(cmds) for _, cmds in batches)
    print(f"Total commands: {total_commands}")
    print(f"Total batches: {len(batches)}")
    print()
    
    # Calculate zero-padding for batch numbers
    padding = len(str(len(batches)))
    # Prepare arguments for parallel processing
    write_args = [
        (batch_num, commands, str(output_dir), args.prefix, args.suffix, padding)
        for batch_num, commands in batches
    ]
    
    # Write batch files in parallel
    print("Writing batch files in parallel...", end='', flush=True)
    write_start = time.time()
    try:
        with Pool(processes=num_processes) as pool:
            results = pool.map(write_batch_file, write_args)
        write_time = time.time() - write_start
    except Exception as e:
        print(f"\nERROR: Parallel processing failed: {e}", file=sys.stderr)
        sys.exit(1)
    print(f" Done ({write_time:.2f}s)")
    
    # Check for errors
    failed = [r for r in results if not r[2]]
    if failed:
        print(f"\nERROR: {len(failed)} batch file(s) failed to write:", file=sys.stderr)
        for batch_num, _, _, error in failed:
            print(f" Batch {batch_num}: {error}", file=sys.stderr)
        sys.exit(1)
    
    # Verify batch creation
    batch_files = list(output_dir.glob(f"{args.prefix}*{args.suffix}"))
    if len(batch_files) == 0:
        print("ERROR: No batch files were created!", file=sys.stderr)
        sys.exit(1)
    if len(batch_files) != len(batches):
        print(f"WARNING: Expected {len(batches)} batch files, found {len(batch_files)}",
              file=sys.stderr)
    
    # Print summary
    total_time = time.time() - start_time
    print()
    print("=" * 60)
    print("BATCH FILE CREATION SUMMARY")
    print("=" * 60)
    print(f"Total commands processed: {total_commands}")
    print(f"Batch size: {args.batch_size}")
    print(f"Total batches created: {len(batch_files)}")
    print(f"Output directory: {output_dir}")
    print(f"Batch prefix: {args.prefix}")
    print(f"Batch suffix: {args.suffix}")
    print(f"Parallel processes used: {num_processes}")
    print(f"Total time: {total_time:.2f}s")
    print(f"Throughput: {total_commands/total_time:.0f} commands/sec")
    print("=" * 60)
    
    if args.verbose:
        print("\nBatch file listing:")
        for batch_file in sorted(batch_files):
            with open(batch_file) as f:
                num_lines = sum(1 for _ in f)
            print(f" {batch_file.name}: {num_lines} commands")

if __name__ == '__main__':
    main()
