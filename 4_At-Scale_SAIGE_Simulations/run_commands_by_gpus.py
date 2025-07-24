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
#Get number of gpus per node
gpus_per_node=8

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
        
    # Read commands from file
    with open(args.commands_file, 'r') as f:
        commands = [line.strip() for line in f if line.strip()]
    
    #Check number of gpus per node and process accordingly
    if  args.gpus < gpus_per_node:
        #Calculate number of workers per node
        workers_per_node = 8 // args.gpus
        # Parsl configuration for OLCF Frontier (56 cores & 8 GPUs per node)
        config = Config(
            executors=[
            	HighThroughputExecutor(
                        label="frontier_htex",
                        address=address_by_hostname(),
                        # Set max_workers_per_node to the max number of tasks that can run concurrently
                        # on a single node. You can only accommodate 2 tasks each attaching to 3 GPUs
                        max_workers_per_node=workers_per_node,
                        # You could also try passing in explicit GPU ids like:
                        available_accelerators=["0,1,2", "3,4,5"],
                        provider=LocalProvider(
                    		# Number of nodes job
                    		nodes_per_block=num_nodes, #Variable denoting total number of nodes requested in batch job
                            # The manager is the component that is launched, in this case, we need all the GPUs attach
                            # to each task. So I would recommend removing --gpus-per-task=1 or setting it to 8
                    		launcher=SrunLauncher(overrides='--gpu-bind closest --exclusive'),
                    		init_blocks=1,
                    		max_blocks=1
                        ),
            	)
            ],
            usage_tracking=True,
        )
    else:
        pass

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


