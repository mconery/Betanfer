'''
make_matrices.py

The purpose of this script is to be a wrapper that determines for each locus
which plink files are needed to map the locus. It submits separate srun 
commands to combine the plink files and then ultimately makes map 
files, or both map files and matrices depending upon the input parameters

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import sys #Taking command line parameters as inputs
import getopt
from os.path import exists
import os
import pandas as pd

#Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
    ARGUMENTS
    -l => Input file of loci with header (chr, start, end, length) REQUIRED
    -o => Output directory REQUIRED
    -p => Explicit path PLINK 2 file prefix (without .pgen/.pvar/.psam) REQUIRED
    -c => Command output file REQUIRED
    --prefix => Output file prefix REQUIRED
    """)
    sys.exit(exit_num)

###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    try: 
        opts, args = getopt.getopt(sys.argv[1:], "l:o:p:c:nh", ["help", "prefix="])
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)
    
    ## Check for help flag and return empty string to break function
    help_ind = options_dict.get('-h', False)
    help_ind = options_dict.get("help", False) or help_ind
    if help_ind != False:
        help()
    
    ## Required arguments
    try:
        locus_file = options_dict['-l']
        output_dir = options_dict['-o']
        output_prefix = options_dict['--prefix']
        command_file = options_dict['-c']
        plink_prefix = options_dict['-p']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()
        
    #Confirm that the file locs exist
    if exists(locus_file) == False:
        print("ERROR: Input file of merged significant loci does not exist.")
        sys.exit(1)
    if any([exists(plink_prefix + x) == False for x in [".pvar", ".psam", ".pgen"]]):
        print("ERROR: One or more of the plink 2 files do not exist at the specified location.")
        sys.exit(1)
    #Check that directories exist
    if exists(output_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)
    
    print("Acceptable Inputs Given")
    #Call driver function
    driver(locus_file, output_dir, output_prefix, command_file, plink_prefix)

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(locus_file, output_dir, output_prefix, command_file, plink_prefix):
    
    #Set plink directories
    plink2 = "/lustre/orion/bif154/proj-shared/mconery/tools/bin/plink2"
        
    #Add slashes to directories
    output_dir = output_dir.rstrip("/") + "/"
    
    #Read in the locus file
    temp = pd.read_csv(locus_file, sep=",")
    locus_raw = []  
    for _, row in temp.iterrows():
        locus_raw.append([
            row['CHR'],
            row['START'],
            row['END'],
            row['LENGTH']
        ])
    
    #Create a list of all the matrix/map files in the output directory
    out_files = os.listdir(output_dir) 
    var_files = [x.replace('.unphased.vcor1.vars','') for x in out_files if x.split(".")[-1] == 'vars']
    vcor_files = [x.replace('.unphased.vcor1','') for x in out_files if x.split(".")[-1] == 'vcor1']
    #Find subset where both are present
    existing_subset = [x for x in var_files if x in vcor_files]
    
    #Create a blank list that will hold the commands
    locus_coms = []
    #Loop over the locus_raw list and start making matrices
    for locus in locus_raw:
        #Get parts of locus
        chr_num=locus[0]
        start=locus[1]
        end=locus[2]
        #Check if file exists yet
        out_file = f"{output_prefix}.chr{chr_num}.{start}.{end}"
        if out_file in existing_subset:
           print(f"NOTE: Skipping {output_prefix}.chr{chr_num}.{start}.{end} since files already exist.")
           continue
        #Create base command for each locus
        base_cmd = f"{plink2} --pfile {plink_prefix} --chr {chr_num} --from-bp {start} --to-bp {end} "
        #Create command for making LD Matrix generation  
        locus_cmd = base_cmd + f"--r-unphased square --out {output_dir}{output_prefix}.chr{chr_num}.{start}.{end}; "
        #Delete the temporary files
        locus_cmd = locus_cmd + f'rm  {output_dir}{output_prefix}.chr{chr_num}.{start}.{end}.log;'
        #Replace all the spaces and semi-colons for bash interpretation
        locus_cmd = locus_cmd.replace(" ", "___")
        locus_cmd = locus_cmd.replace(";", "^^^")
        #Append the commands to the locus_coms list
        locus_coms.append(locus_cmd)
    
    #Verify that commands aren't empty
    if len(locus_coms) > 0:
        #Add a do-nothing line at the end of the file
        locus_coms.append("echo___\"empty\"")
        #Print the list of commands to a file
        with open(command_file, "w") as out_path:
            out_path.write('\n'.join(locus_coms))
            print(f"SUCCESS: Wrote commands for making {len(locus_coms)-1} matrices.")
    else:
        print("SUCCESS: All necessary matrix files previously generated.")
    

###############################################################################


#cProfile.run('main(sys.argv[1:])')
         
if __name__ == "__main__":
    main(sys.argv[1:])

