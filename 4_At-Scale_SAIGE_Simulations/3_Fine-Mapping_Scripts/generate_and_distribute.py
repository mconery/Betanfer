#!/usr/bin/env python3

"""
Master work queue generator and distributor for SuSiE fine-mapping.

This script:
1. Generates all SuSiE commands using your existing command generation logic
2. Creates a distributed work queue structure on Lustre
3. Submits one persistent manager job per node

Usage:
    python generate_and_distribute_commands.py \
        --all-commands all_commands.txt \
        --queue-dir /path/to/queue \
        --num-nodes 1000 \
        --commands-per-batch 41 \
        --node-manager node_manager.py

Prerequisites:
    - all_commands.txt file (output from generate_susie_commands.py)
    - node_manager.py script available
"""

import os
import json
import argparse
import sys
from pathlib import Path
from typing import List, Dict
import subprocess
import time

def create_queue_structure(queue_dir: str) -> None:
    """
    Create the distributed work queue directory structure on Lustre.
    
    Structure:
        queue_dir/
            ├── work/              (unclaimed batches)
            ├── in_progress/       (currently processing)
            ├── completed/         (finished batches)
            ├── logs/              (node execution logs)
            └── queue_state.json   (metadata)
    
    Args:
        queue_dir: Root path for queue on Lustre
    """
    queue_path = Path(queue_dir)
    
    required_dirs = ['work', 'in_progress', 'completed', 'logs']
    
    print(f"Creating queue structure at {queue_dir}...")
    
    for dirname in required_dirs:
        dirpath = queue_path / dirname
        dirpath.mkdir(parents=True, exist_ok=True)
        print(f"  ✓ Created {dirname}/")
    
    print("Queue structure created successfully")


def read_all_commands(commands_file: str) -> List[str]:
    """
    Read all commands from the master commands file.
    
    Args:
        commands_file: Path to file with one command per line
    
    Returns:
        List of command strings
    """
    commands = []
    
    try:
        with open(commands_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comments
                if line and not line.startswith('#'):
                    commands.append(line)
    except FileNotFoundError:
        print(f"Error: Commands file not found: {commands_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading commands file: {e}")
        sys.exit(1)
    
    print(f"Read {len(commands)} commands from {commands_file}")
    return commands


def create_batch_files(
    commands: List[str],
    queue_dir: str,
    commands_per_batch: int
) -> int:
    """
    Create batch JSON files in the work queue.
    
    Each batch is a JSON file containing commands_per_batch commands.
    
    Args:
        commands: List of all commands
        queue_dir: Root queue directory
        commands_per_batch: Number of commands per batch
    
    Returns:
        Number of batches created
    """
    queue_path = Path(queue_dir)
    work_dir = queue_path / "work"
    
    print(f"\nCreating {len(commands) // commands_per_batch + 1} batch files...")
    
    batch_num = 0
    current_batch = []
    total_batches = 0
    
    for cmd in commands:
        current_batch.append({
            'command': cmd,
            'metadata': {
                'index': len(current_batch),
            }
        })
        
        if len(current_batch) >= commands_per_batch:
            batch_file = work_dir / f"batch_{batch_num:06d}.json"
            
            try:
                with open(batch_file, 'w') as f:
                    json.dump(current_batch, f)
                total_batches += 1
                
                if batch_num % 100 == 0 and batch_num > 0:
                    print(f"  ✓ Created batch {batch_num}")
            
            except Exception as e:
                print(f"Error writing batch file {batch_file}: {e}")
                sys.exit(1)
            
            batch_num += 1
            current_batch = []
    
    # Write remaining commands
    if current_batch:
        batch_file = work_dir / f"batch_{batch_num:06d}.json"
        try:
            with open(batch_file, 'w') as f:
                json.dump(current_batch, f)
            total_batches += 1
            print(f"  ✓ Created batch {batch_num} (final, {len(current_batch)} commands)")
        except Exception as e:
            print(f"Error writing final batch file: {e}")
            sys.exit(1)
    
    return total_batches


def create_queue_state(queue_dir: str, total_commands: int, total_batches: int) -> None:
    """
    Create the queue state metadata file.
    
    Args:
        queue_dir: Root queue directory
        total_commands: Total number of commands
        total_batches: Total number of batches
    """
    queue_path = Path(queue_dir)
    state_file = queue_path / "queue_state.json"
    
    state = {
        'total_batches': total_batches,
        'total_commands': total_commands,
        'created_at': time.time(),
        'status': 'ready',
        'metadata': {
            'description': 'SuSiE fine-mapping distributed work queue'
        }
    }
    
    try:
        with open(state_file, 'w') as f:
            json.dump(state, f, indent=2)
        print(f"\n✓ Queue state file created")
    except Exception as e:
        print(f"Error creating queue state file: {e}")
        sys.exit(1)


def submit_node_managers(
    num_nodes: int,
    queue_dir: str,
    node_manager_script: str,
    workers_per_node: int = 28,
    dry_run: bool = False
) -> None:
    """
    Submit one persistent manager job per node to Flux.
    
    Args:
        num_nodes: Number of compute nodes to use
        queue_dir: Root queue directory
        node_manager_script: Path to node_manager.py
        workers_per_node: Number of worker threads per node
        dry_run: If True, print commands instead of executing
    """
    # Verify node_manager_script exists
    if not Path(node_manager_script).exists():
        print(f"Error: node_manager_script not found: {node_manager_script}")
        sys.exit(1)
    
    print(f"\nSubmitting {num_nodes} node manager jobs to Flux...")
    print(f"(Each job gets {workers_per_node} worker threads)")
    
    submitted = 0
    failed = 0
    
    for node_id in range(num_nodes):
        job_name = f"susie_manager_{node_id:05d}"
        
        log_file = f"{queue_dir}/logs/{job_name}.log"
        
        flux_cmd = (
            f"flux submit "
            f"-N1 -n1 -c{workers_per_node} "
            f"--job-name={job_name} "
            f"--output={log_file} "
            f"python {node_manager_script} "
            f"--queue-dir {queue_dir} "
            f"--node-id {node_id} "
            f"--workers {workers_per_node}"
        )
        
        if dry_run:
            if node_id < 3 or node_id >= num_nodes - 3:
                print(f"  [{node_id:05d}] {flux_cmd}")
            elif node_id == 3:
                print(f"  ... ({num_nodes - 6} more jobs) ...")
        else:
            try:
                result = subprocess.run(
                    flux_cmd,
                    shell=True,
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                
                if result.returncode == 0:
                    submitted += 1
                    if (submitted) % 50 == 0:
                        print(f"  ✓ Submitted {submitted} node managers...")
                else:
                    failed += 1
                    print(f"  ✗ Failed to submit job {node_id}: {result.stderr}")
            
            except subprocess.TimeoutExpired:
                failed += 1
                print(f"  ✗ Timeout submitting job {node_id}")
            except Exception as e:
                failed += 1
                print(f"  ✗ Error submitting job {node_id}: {e}")
        
        # Stagger submissions to avoid overwhelming Flux
        if not dry_run and submitted % 100 == 0 and submitted > 0:
            time.sleep(0.5)
    
    if dry_run:
        print(f"\n(Dry run) Would submit {num_nodes} jobs")
    else:
        print(f"\n✓ Submitted {submitted} node managers")
        if failed > 0:
            print(f"✗ Failed to submit {failed} jobs")


def main():
    """Main entry point with argument parsing."""
    parser = argparse.ArgumentParser(
        description='Create distributed work queue and submit node managers',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full workflow
  python generate_and_distribute_commands.py \\
      --all-commands all_commands.txt \\
      --queue-dir /lustre/orion/bif154/work_queue \\
      --num-nodes 1000 \\
      --commands-per-batch 41 \\
      --node-manager ./node_manager.py

  # Dry run to see what would happen
  python generate_and_distribute_commands.py \\
      --all-commands all_commands.txt \\
      --queue-dir /lustre/orion/bif154/work_queue \\
      --num-nodes 1000 \\
      --commands-per-batch 41 \\
      --node-manager ./node_manager.py \\
      --dry-run
        """
    )
    
    parser.add_argument(
        '--all-commands',
        required=True,
        help='Path to file with all SuSiE commands (one per line)'
    )
    parser.add_argument(
        '--queue-dir',
        required=True,
        help='Root directory on Lustre for work queue'
    )
    parser.add_argument(
        '--num-nodes',
        type=int,
        required=True,
        help='Number of compute nodes to use'
    )
    parser.add_argument(
        '--commands-per-batch',
        type=int,
        default=41,
        help='Number of commands per batch (default: 41)'
    )
    parser.add_argument(
        '--node-manager',
        required=True,
        help='Path to node_manager.py script'
    )
    parser.add_argument(
        '--workers-per-node',
        type=int,
        default=28,
        help='Worker threads per node (default: 28)'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print commands without executing'
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.num_nodes < 1:
        print(f"Error: num-nodes must be >= 1, got {args.num_nodes}")
        return 1
    
    if args.commands_per_batch < 1:
        print(f"Error: commands-per-batch must be >= 1, got {args.commands_per_batch}")
        return 1
    
    if args.workers_per_node < 1:
        print(f"Error: workers-per-node must be >= 1, got {args.workers_per_node}")
        return 1
    
    if not Path(args.all_commands).exists():
        print(f"Error: Commands file not found: {args.all_commands}")
        return 1
    
    if not Path(args.node_manager).exists():
        print(f"Error: Node manager script not found: {args.node_manager}")
        return 1
    
    print("=" * 70)
    print("DISTRIBUTED SUSIE FINE-MAPPING WORK QUEUE SETUP")
    print("=" * 70)
    
    # Step 1: Create queue structure
    create_queue_structure(args.queue_dir)
    
    # Step 2: Read commands
    commands = read_all_commands(args.all_commands)
    
    # Step 3: Create batch files
    total_batches = create_batch_files(
        commands,
        args.queue_dir,
        args.commands_per_batch
    )
    
    # Step 4: Create queue state
    create_queue_state(args.queue_dir, len(commands), total_batches)
    
    # Step 5: Submit node managers
    submit_node_managers(
        args.num_nodes,
        args.queue_dir,
        args.node_manager,
        args.workers_per_node,
        args.dry_run
    )
    
    print("\n" + "=" * 70)
    print("SETUP COMPLETE")
    print("=" * 70)
    print(f"\nQueue directory: {args.queue_dir}")
    print(f"Total commands: {len(commands)}")
    print(f"Total batches: {total_batches}")
    print(f"Commands per batch: {args.commands_per_batch}")
    print(f"Num nodes: {args.num_nodes}")
    print(f"Workers per node: {args.workers_per_node}")
    
    if not args.dry_run:
        print(f"\n✓ Node managers submitted!")
        print(f"Monitor progress with: python monitor_queue.py --queue-dir {args.queue_dir}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
