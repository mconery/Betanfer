'''
run_commands_by_gpus.py

The purpose of this script is to run a set of commands that are meant to be tied
to gpus. This script uses mpi4py and should be called using the following approach:
    srun -N <number of nodes> -n <number of parallel tasks> --gpus-per-task=1 /
    --exclusive python run_commands_by_gpus.py --commands <command-file> /
    --gpus <gpus needed per command>

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

import argparse
import parsl
from parsl.app.app import bash_app
from parsl.config import Config
from parsl.providers import LocalProvider
from parsl.executors import HighThroughputExecutor
from parsl.launchers import SrunLauncher
from parsl.addresses import address_by_hostname
import os

# Get the number of nodes:
node_raw = os.getenv("SLURM_NODELIST")
node_list = [f"frontier{x}" for x  in node_raw.lstrip("frontier[").rstrip("]").split(",")]
num_nodes = len(node_list)

@bash_app
def run_command(cmd):
    return f"""
    {cmd}
    """

def main():
    #Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Parse commands file and number of GPUs required per command."
    )
    parser.add_argument("-c", "--commands-file", type=str, required=True,
        help="Path to the file containing commands, one per line."
    )
    parser.add_argument("-g", "--gpus", type=int, required=True,
        help="Number of GPUs needed to run each command (integer)."
    )
    args = parser.parse_args()
    print(f"Commands file: {args.commands_file}")
    print(f"GPUs per command: {args.gpus}")
    
    #Calculate number of workers per node
    workers_per_node = 8 // args.gpus
    nodes_per_block = args.gpus // 8 + 1
    
    # Read commands from file
    with open(args.commands_file, 'r') as f:
        commands = [line.strip() for line in f if line.strip()]
    
    # Parsl configuration for OLCF Frontier (56 cores per node)
    config = Config(
        executors=[
            HighThroughputExecutor(
                label="frontier_htex",
                address=address_by_hostname(),
                max_workers_per_node=args.gpus,
                available_acclerators=8, 
                provider=LocalProvider(
                    # Number of nodes job
                    nodes_per_block=nodes_per_block,
                    launcher=SrunLauncher(overrides='--gpus-per-task=1 --gpu-bind closest --exclusive'),
                    init_blocks=1,
                    max_blocks=1
                ),
            )
        ],
        usage_tracking=True,
    )

    parsl.clear()
    parsl.load(config)
    
    # Submit each command as a Parsl task
    tasks = []
    for cmd in commands:
        task = run_command(cmd)
        tasks.append(task)

    # Wait for all tasks to complete
    results = []
    for i, task in enumerate(tasks):
        try:
            result = task.result()
            print(f"SUCCESS: Task {i} completed successfully")
            results.append(result)
        except Exception as e:
            print(f"ERROR: Task {i} failed with error: {e}")
            results.append(None)

    parsl.clear()


if __name__ == "__main__":
    main()


