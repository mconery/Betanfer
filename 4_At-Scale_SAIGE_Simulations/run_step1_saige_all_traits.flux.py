#!/usr/bin/env python3

import argparse
import sys
import os
from flux.job import JobspecV1
from flux import core
from flux.executor import FluxExecutor

# Get the number of nodes:
node_raw = os.getenv("SLURM_NODELIST")
node_list = [f"frontier{x}" for x  in node_raw.lstrip("frontier[").rstrip("]").split(",")]
num_nodes = len(node_list)

def parse_args():
    parser = argparse.ArgumentParser(description="Submit node-parallel tasks using Flux on OLCF Frontier")
    parser.add_argument("command_file", help="Path to text file containing commands (one per line)")
    parser.add_argument('--max-concurrent', type=int, default=num_nodes * 2, help="Max number of concurrent jobs")
    parser.add_argument('--timeout', type=int, default=0, help="Job timeout in seconds (0 for no timeout)")
    return parser.parse_args()

def read_commands(fname):
    with open(fname, "r") as f:
        return [line.strip() for line in f if line.strip()]

def main():
    args = parse_args()
    commands = read_commands(args.command_file)

    flux = core.Flux()  # Assumes you're running inside a Flux instance (e.g., started with srun)
    executor = FluxExecutor(flux)

    jobs = []
    for idx, cmd in enumerate(commands):
        spec = JobspecV1.from_command(
            command=["bash", "-c", cmd],
            num_tasks=1,
            cores_per_task=24,    # Changed from 64 to 24 cores
            gpus_per_task=3       # Changed from 8 to 3 GPUs
        )
        if args.timeout > 0:
            spec.timeout_secs = args.timeout
        jobs.append((idx, cmd, spec))

    # Submit and monitor jobs with a concurrency limit
    from concurrent.futures import ThreadPoolExecutor, as_completed

    results = []
    with ThreadPoolExecutor(max_workers=args.max_concurrent) as pool:
        future_map = {pool.submit(executor.submit, spec): (idx, cmd) for idx, cmd, spec in jobs}
        for future in as_completed(future_map):
            idx, cmd = future_map[future]
            try:
                jobid = future.result()
                print(f"Submitted job {idx} (cmd: {cmd}) as Flux jobid: {jobid}", flush=True)
                results.append((idx, cmd, jobid))
            except Exception as e:
                print(f"Failed to submit job {idx}: {e}", file=sys.stderr)

    print("All jobs submitted.")

if __name__ == "__main__":
    main()
