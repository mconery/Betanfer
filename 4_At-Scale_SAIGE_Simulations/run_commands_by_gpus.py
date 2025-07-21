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


import os
import subprocess
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
import argparse
import socket

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
    
    def myfunction(command):
        try:
            print("banana")
            local_rank = int(os.environ.get("SLURM_LOCALID", 0))
            gpus_per_command = int(args.gpus)
            first_gpu = local_rank * gpus_per_command
            gpu_ids = ",".join(str(i) for i in range(first_gpu, first_gpu + gpus_per_command))
            os.environ["ROCR_VISIBLE_DEVICES"] = gpu_ids  # for ROCm
            output = subprocess.check_output(command.split(" "), env=os.environ)
            print(output)
        except Exception as e:
            print(e)
    
    with MPICommExecutor(MPI.COMM_WORLD, root=0) as executor:
        if executor is not None:
            print("test")
            executor.map(myfunction, commands)

if __name__ == "__main__":
    main()


