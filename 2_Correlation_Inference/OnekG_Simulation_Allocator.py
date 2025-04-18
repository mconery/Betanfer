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

#Define help function
def help(exit_num=0):
    print("""Usage: python OnekG_Simulation_Allocator.py [OPTIONS]
Options:
  --command_file        File of commands to run (required)
  -h, --help            Show this help message""")
    sys.exit(exit_num)


###############################################################################
###############  DEFINE A FUNCTION TO COLLECT OUTPUT FILE NAMES  ##############
###############################################################################

def get_output_names(command):
    # Parse command-line arguments
    try:
        opts, args = getopt.getopt(
            command.split(" ")[2:],
            "h",
            ["help", "plink=", "pairs=", "traits=", "keep=", "outdir=", 
             "max_dist=", "trait_type=", "case_prob=", "keep_conv=", "plink_conv=",
             "n_causal=", "min_distance=", "effect_size="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
    
    plink_path = None
    n_pairs = None
    n_traits = None
    keep_file = None
    outdir = None
    max_distance = 1e6
    case_prob = 0.5
    trait_type = 'continuous'
    n_causal = 0
    min_distance = 10000000  # 10Mb
    effect_size = 0.1
    keep_conv = "1000G.<POP>.txt"
    plink_conv = "<CHROMOSOME>_10k_sample.bed"
    
    for o, a in opts:
        if o in ("-h", "--help"):
            help()
            sys.exit()
        elif o == "--plink":
            plink_path = a
        elif o == "--pairs":
            n_pairs = int(a)
        elif o == "--traits":
            n_traits = int(a)
        elif o == "--keep":
            keep_file = a
        elif o == "--outdir":
            outdir = a
        elif  o == "--max_dist":
            max_distance = a
        elif o == "--trait_type":
            trait_type = a.lower()
        elif o == "--case_prob":
            case_prob = float(a)
            if not (0 < case_prob < 1):
                print("ERROR: --case_prob must be between 0 and 1")
                sys.exit(1)
        elif o == "--n_causal":
            n_causal = int(a)
        elif o == "--min_distance":
            min_distance = int(a)
        elif o == "--effect_size":
            effect_size = float(a)
            if not (effect_size > 0):
                print("ERROR: Specified effect size must be >0.")
                sys.exit(1)
        elif o == "--keep_conv":
            keep_conv = a
        elif o == "--plink_conv":
            plink_conv = a
        else:
            assert False, "WARNING: Unhandled input option"
    
    # Validate required arguments
    if None in [plink_path, n_pairs, n_traits, outdir]:
        print("ERROR: Missing required arguments!")
        sys.exit(2)

    # Path Validation checks
    if not all(os.path.exists(f) for f in [plink_path, outdir]):
        print("ERROR: Given output directory or plink file not found")
        sys.exit(1)
    
    #Numeric Value Validation Checks
    def check_int(check_num, description):
        try:
           check_num = int(check_num)
           return check_num
        except ValueError:
            print(f"ERROR: Invalid {description} given.")
            sys.exit(1)
    max_distance = check_int(max_distance, "maximum distance of variant pairs")
    n_traits = check_int(n_traits, "number of traits")
    n_pairs = check_int(n_pairs, "number of variant pairs")
    n_causal = check_int(n_causal, "number of causal variants")
    min_distance = check_int(min_distance, "minimum distance of causal variants")
    
    # Validate trait type
    if trait_type not in ['continuous', 'binary']:
        print("ERROR: Invalid trait type give. Select continuous/binary.")
        sys.exit(2)
    
    #Check naming conventions for needed parameters
    if '<POP>' not in keep_conv and keep_file:
        print("ERROR: Invalid keep file naming convention given. Convention must include <POP>.")
        sys.exit(2)
    if '<CHROMOSOME>' not in plink_conv and keep_file:
        print("ERROR: Invalid plink file naming convention given. Convention must include <CHROMOSOME>.")
        sys.exit(2)
    
    #Set validation on number of causal variants
    if n_causal < 0 or n_causal > 5:
        print("ERROR: --n_causal must be 0-5")
        sys.exit(1)
    
    #Extract chromosome
    plink_name_remnant = plink_conv.split("<CHROMOSOME>")
    chromo = os.path.basename(plink_path)
    for remnant in plink_name_remnant:
        chromo = chromo.replace(remnant, "")
    
    #Extract population name
    if keep_file: #Checks for whether keep file given
        #Collect population name if filtering
        keep_name_remnant = keep_conv.split("<POP>")
        pop_name = os.path.basename(keep_file)
        for remnant in keep_name_remnant:
            pop_name = pop_name.replace(remnant, "")
    else:
        pop_name='ALL'
       
    #Create output file name
    out_file_name = f"{chromo}.{pop_name}.{n_traits}_traits.{n_pairs}_pairs.{f'{trait_type}_{case_prob}' if trait_type == 'binary' else f'{trait_type}'}{f'.{n_causal}_causal_{effect_size}' if n_causal > 0 else f'.{n_causal}_causal'}.csv"
    output_path = os.path.join(outdir, out_file_name)
    return(output_path)

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
