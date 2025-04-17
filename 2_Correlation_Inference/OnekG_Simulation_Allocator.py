'''
OnekG_Simulation_Allocator.py

The purpose of this script is to act as a wrapper for 1000G_Simulation_Standalone.py
and manage the multi-node distribution of jobs using mpi4py

***Run in Conda finemap environment***
'''

#Remainder of package imports
import sys
import getopt
import os
from mpi4py import MPI
import numpy as np

#Define help function
def help(exit_num=0):
    print("""Usage: python OnekG_Simulation_Allocator.py [OPTIONS]
Options:
  --command_file        File of commands to run (required)
  -h, --help            Show this help message""")
    sys.exit(exit_num)


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND MULTI-PROCESS LOCI #############
###############################################################################

def main():
    
    # Parse command-line arguments
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "command_file="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
        
    #Set parameter default values
    command_file = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
            sys.exit()
        elif o == "--command_file":
            command_file = a
        else:
            print(o)
            assert False, "WARNING: Unhandled input option"
    
    # Validate required arguments
    if None in [command_file]:
        print("ERROR: Missing required command file!")
        sys.exit(2)

    # Path Validation checks
    if not all(os.path.exists(f) for f in [command_file]):
        print("ERROR: Input command file does not exist.")
        sys.exit(1)
        
    #Read in the command files
    with open(command_file, 'r') as file:
        commands = file.readlines()
    
    #Now loop over the test_parameters and submit the jobs
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    num_nodes = comm.size
        
    # Split up commands into a set of chunks
    chunk_dict = {x : [] for x in range(num_nodes)}
    loop_rank = 0
    i = 0
    while i < len(commands): 
        chunk_dict[loop_rank].append(commands[i])
        i += 1
        if loop_rank < num_nodes - 1:
            loop_rank += 1
        else:
            loop_rank = 0
        
    # Process chunks on nodes using all local cores
    for command in chunk_dict[rank]:
        os.system(command)
    

if __name__ == "__main__":
    main()
