#!/usr/bin/env python3
"""
index_gwama_files_mpi.py

Multi-node MPI parallelized utility to tabix-index GWAMA/BOLT summary statistic files,
extracting genomic coordinates directly from variant IDs using a flexible, user-specified format.

This script uses MPI to distribute work across multiple nodes on HPC systems like OLCF Frontier.
Each MPI process handles a subset of files and uses local multiprocessing for additional parallelization.

Usage
-----
# Submit with SLURM requesting multiple nodes
sbatch --nodes=4 --ntasks-per-node=1 --cpus-per-task=64 job_script.slurm

# Or run directly with srun:
srun --nodes=4 --ntasks-per-node=1 --cpus-per-task=64 python index_gwama_files_mpi.py \
  --scratch_dir /path/to/scratch_dir \
  --pattern "/path/to/GWAS/*/Pheno*.meta.out" \
  --filename_regex "^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.meta\\.out$"

Arguments
---------
--scratch_dir         Directory to write .gz and .tbi outputs [REQUIRED]
--pattern             Glob pattern to search for input summary statistics files [REQUIRED]
--variant_col         Column name holding variant IDs [default: rs_number]
--variant_id_format   Format string for variant IDs [default: CHR:POS:REF:ALT]
--chrom_col           Output column name for chromosome [default: CHR]
--pos_col             Output column name for position [default: POS]
--filename_regex      Regex pattern to filter filenames [optional]
--local_workers       Number of local processes per MPI rank [default: auto-detect]

MPI Strategy
------------
- Each MPI rank (node) gets a subset of files to process
- Each MPI rank uses multiprocessing locally for additional parallelization
- Results are gathered at the end for reporting

Example SLURM Script
-------------------
#!/bin/bash
#SBATCH --job-name=index_gwama_mpi
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00
#SBATCH --account=your_project

module load python/3.10
module load gcc/11.2.0

srun python index_gwama_files_mpi.py \
  --scratch_dir /lustre/orion/bif154/scratch/mconery/indexed_gwama_sumstats \
  --pattern "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/meta_analysis_results/*/Pheno*.meta.out" \
  --filename_regex "^Pheno\\d+\\.(binary|continuous)\\.chr\\d+\\.meta\\.out$"

"""

import os
import glob
import argparse
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    from mpi4py import MPI
    HAS_MPI = True
except ImportError:
    HAS_MPI = False
    print("Warning: mpi4py not available. Running in single-node mode.")

from extraction_utils import index_gwama_file_with_variant_ids

def find_valid_meta_files(pattern, filename_regex=None):
    """Find files matching the glob pattern and optionally filter by regex on basename."""
    all_files = glob.glob(pattern)
    
    if not filename_regex:
        return all_files
    
    try:
        regex_pattern = re.compile(filename_regex)
    except re.error as e:
        raise ValueError(f"Invalid regex pattern '{filename_regex}': {e}")
    
    filtered_files = [
        f for f in all_files 
        if regex_pattern.match(os.path.basename(f))
    ]
    
    return filtered_files

def distribute_files_across_ranks(files, comm):
    """Distribute files across MPI ranks as evenly as possible."""
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    # Calculate how many files each rank should get
    files_per_rank = len(files) // size
    remainder = len(files) % size
    
    # Determine start and end indices for this rank
    if rank < remainder:
        # First 'remainder' ranks get one extra file
        start_idx = rank * (files_per_rank + 1)
        end_idx = start_idx + files_per_rank + 1
    else:
        # Remaining ranks get the standard number of files
        start_idx = remainder * (files_per_rank + 1) + (rank - remainder) * files_per_rank
        end_idx = start_idx + files_per_rank
    
    # Return this rank's subset of files
    return files[start_idx:end_idx]

def process_files_on_rank(files, args, local_workers):
    """Process a list of files using local multiprocessing on this MPI rank."""
    if not files:
        return 0, 0
    
    rank = MPI.COMM_WORLD.Get_rank() if HAS_MPI else 0
    
    print(f"Rank {rank}: Processing {len(files)} files with {local_workers} local workers")
    
    task_args = [
        (fname, args.scratch_dir, args.variant_col, args.chrom_col, args.pos_col, args.variant_id_format)
        for fname in files
    ]
    
    success_count = 0
    error_count = 0
    
    with ThreadPoolExecutor(max_workers=local_workers) as executor:
        future_to_file = {
            executor.submit(index_gwama_file_with_variant_ids, task_arg): task_arg[0]
            for task_arg in task_args
        }
        
        for future in as_completed(future_to_file):
            fname = future_to_file[future]
            try:
                bgzf_path, tbi_path = future.result()
                print(f"Rank {rank}: Indexed {os.path.basename(fname)} -> {os.path.basename(bgzf_path)}")
                success_count += 1
            except Exception as exc:
                print(f"Rank {rank}: Error indexing {os.path.basename(fname)}: {exc}")
                error_count += 1
    
    return success_count, error_count

def main():
    parser = argparse.ArgumentParser(description="MPI parallelized tabix indexing of GWAMA summary statistics files.")
    parser.add_argument("--scratch_dir", type=str, required=True,
                        help="Directory where indexed files will be written.")
    parser.add_argument("--pattern", type=str, required=True,
                        help="Glob pattern to search for GWAMA summary stats files.")
    parser.add_argument("--variant_col", type=str, default="rs_number",
                        help="Column name for variant IDs [default: rs_number]")
    parser.add_argument("--variant_id_format", type=str, default="CHR:POS:REF:ALT",
                        help="Format string describing variant ID structure [default: CHR:POS:REF:ALT]")
    parser.add_argument("--chrom_col", type=str, default="CHR",
                        help="Output column name for chromosome [default: CHR]")
    parser.add_argument("--pos_col", type=str, default="POS",
                        help="Output column name for position [default: POS]")
    parser.add_argument("--filename_regex", type=str, default=None,
                        help="Optional regex pattern to filter filenames")
    parser.add_argument("--local_workers", type=int, default=None,
                        help="Number of local processes per MPI rank [default: auto-detect]")
    
    args = parser.parse_args()
    
    # Initialize MPI
    if HAS_MPI:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    else:
        # Fallback for non-MPI execution
        comm = None
        rank = 0
        size = 1
    
    # Determine local worker count
    if args.local_workers is None:
        # Use available CPU count, but be conservative on HPC systems
        local_workers = min(os.cpu_count() or 1, 64)  # Cap at 64 for HPC nodes
    else:
        local_workers = args.local_workers
    
    # Only rank 0 does file discovery to avoid race conditions
    if rank == 0:
        print(f"Running with {size} MPI ranks, {local_workers} local workers per rank")
        
        # Find all files matching criteria
        all_files = find_valid_meta_files(args.pattern, args.filename_regex)
        
        if not all_files:
            print(f"No files found matching pattern: {args.pattern}")
            if args.filename_regex:
                print(f"with filename regex: {args.filename_regex}")
            all_files = []
        else:
            print(f"Found {len(all_files)} files matching criteria")
            if args.filename_regex:
                print(f"Pattern: {args.pattern}")
                print(f"Filename regex: {args.filename_regex}")
        
        # Create scratch directory
        os.makedirs(args.scratch_dir, exist_ok=True)
    else:
        all_files = None
    
    # Broadcast file list to all ranks
    if HAS_MPI:
        all_files = comm.bcast(all_files, root=0)
    
    if not all_files:
        if rank == 0:
            print("No files to process. Exiting.")
        return
    
    # Distribute files across ranks
    if HAS_MPI:
        my_files = distribute_files_across_ranks(all_files, comm)
    else:
        my_files = all_files
    
    # Process files on this rank
    start_time = time.time()
    success_count, error_count = process_files_on_rank(my_files, args, local_workers)
    end_time = time.time()
    
    # Gather results from all ranks
    if HAS_MPI:
        all_success = comm.gather(success_count, root=0)
        all_errors = comm.gather(error_count, root=0)
        all_times = comm.gather(end_time - start_time, root=0)
    else:
        all_success = [success_count]
        all_errors = [error_count]
        all_times = [end_time - start_time]
    
    # Print final summary (only from rank 0)
    if rank == 0:
        total_success = sum(all_success)
        total_errors = sum(all_errors)
        max_time = max(all_times)
        
        print(f"\n{'='*60}")
        print(f"MPI INDEXING COMPLETE")
        print(f"{'='*60}")
        print(f"Total files processed: {total_success + total_errors}")
        print(f"Successful: {total_success}")
        print(f"Failed: {total_errors}")
        print(f"Wall time: {max_time:.2f} seconds")
        print(f"Used {size} MPI ranks with {local_workers} local workers each")
        print(f"{'='*60}")

if __name__ == "__main__":
    main()