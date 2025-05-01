'''
OnekG_Simulation_looping_Allocator.py

The purpose of this script is to act as a wrapper for 1000G_Simulation_Standalone.py
and manage the multi-node distribution of jobs simply by using the PBS scheduler

***Run in Conda finemap_old environment***
'''

#Remainder of package imports
import sys
import getopt
import os

#Define help function
def help(exit_num=0):
    print("""Usage: python OnekG_Simulation_Allocator.py [OPTIONS]
Options:
  --command_file        File of commands to run (required)
  --total_nodes         Total number of nodes allocated to job (required)
  --node_num            Number of current node (required)
  --log_dir             Directory for output logs (default: /flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs; optional)
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
            ["help", "command_file=", "log_dir=", "node_num=", "total_nodes="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)
        
    #Set parameter default values
    command_file = None
    node_num = None
    total_nodes = None
    log_dir = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
            sys.exit()
        elif o == "--command_file":
            command_file = a
        elif o == "--node_num":
            node_num = a
        elif o == "--total_nodes":
            total_nodes = a
        elif o == "--log_dir":
            log_dir = a
        else:
            print(o)
            assert False, "WARNING: Unhandled input option"
    
    # Validate required arguments
    if None in [command_file, node_num, total_nodes]:
        print("ERROR: Missing required command file!")
        sys.exit(2)
    
    # Path Validation checks
    if not all(os.path.exists(f) for f in [command_file, log_dir]):
        print("ERROR: Input command file or log directory does not exist.")
        sys.exit(1)
    
    #Check node number is integer
    def check_int(check_num, description):
        try:
           check_num = int(check_num)
           return check_num
        except ValueError:
            print(f"ERROR: Invalid {description} given.")
            sys.exit(1)
    node_num = check_int(node_num, "number of current node")
    total_nodes = check_int(total_nodes, "total number of nodes")
    
    #Read in the command files
    with open(command_file, 'r') as file:
        commands = file.readlines()
    
    #Get the names of the output files
    out_files = [get_output_names(x) for x in commands]
    #Get the commands for the files that don't exist yet
    needed_commands = [commands[i] for i in range(len(out_files)) if not os.path.exists(out_files[i])]
    
    #Rename variables because I borrowed this code and am too lazy to rewrite it
    rank = node_num
    num_nodes = total_nodes
        
    # Split up commands into a set of chunks
    chunk_dict = {x : [] for x in range(num_nodes)}
    loop_rank = 0
    i = 0
    while i < len(needed_commands): 
        chunk_dict[loop_rank].append(needed_commands[i])
        i += 1
        if loop_rank < num_nodes - 1:
            loop_rank += 1
        else:
            loop_rank = 0
    
    #Change working directory
    os.chdir(log_dir)
    print("NOTE: Starting Simulations")
    
    #Loop over commands and send to console
    for command in chunk_dict[rank]:
        os.system(command)
    
    print(f"SUCCESS: Finished running all simulations for command file: {os.path.basename(command_file)}")
    
if __name__ == "__main__":
    main()

