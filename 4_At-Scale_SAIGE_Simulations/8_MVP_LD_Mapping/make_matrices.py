'''

make_matrices.py

The purpose of this script is to read in lists of variants for fine-mapping
and create LD matrices for each set of variants across specified ancestral 
populations. It handles variant allele strand flipping and generates srun 
commands to process the variants in parallel.

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python. ***

'''

# import modules that will be needed
import sys  # Taking command line parameters as inputs
import getopt
from os.path import exists
import os
import pandas as pd

# Define help function to be called if there is a problem
def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
-v => Input directory of variant files for fine-mapping (variant format: chr:pos:ref:alt) REQUIRED
-o => Base folder in which to output matrix/map files REQUIRED
-a => Ancestry (comma-separated list, e.g., AFR,EUR) REQUIRED
-i => Directory of files specifying individuals for each ancestry (base directory with pop subfolders) REQUIRED
-r => Registry file that lists out the ranges found in each plink segment REQUIRED
-c => File of commands to be submitted to the console REQUIRED
-p => Directory that houses plink files (default is /gpfs/arx2/bif154_mde/proj-shared/RITM039091/pgen) OPTIONAL
-b => Directory containing plink2 binary (default is /gpfs/arx2/bif154_mde/proj-shared/tools/bin) OPTIONAL
-t => Output type (options are Matrix/Map/Both; default is Both) Optional

""")
    sys.exit(exit_num)

###############################################################################
########################### COLLECT AND CHECK ARGUMENTS #####################
###############################################################################

## GLOBAL VARS
## MAIN ##
def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:o:a:i:r:c:p:b:t:n:nh")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    ## Check for help flag and return empty string to break function
    help_ind = options_dict.get('-h', False)
    if help_ind != False:
        help()

    ## Required arguments
    try:
        variant_dir = options_dict['-v']
        output_dir = options_dict['-o']
        ancestry_str = options_dict['-a']
        individual_dir = options_dict['-i']
        register_file = options_dict['-r']
        command_file = options_dict['-c']
    except KeyError:
        print("ERROR: One of your required arguments does not exist.")
        help()

    # Optional arguments
    output_type = options_dict.get('-t', "BOTH")
    plink_dir = options_dict.get('-p', "/gpfs/arx2/bif154_mde/proj-shared/RITM0309091/pgen")
    plink_bin_dir = options_dict.get('-b', "/gpfs/arx2/bif154_mde/proj-shared/tools/bin")

    # Parse ancestry list
    ancestry_list = [x.strip().upper() for x in ancestry_str.split(',')]
    
    # Validate ancestries
    valid_ancestries = ['AFR', 'AMR', 'EAS', 'EUR']
    for anc in ancestry_list:
        if anc not in valid_ancestries:
            print(f"ERROR: Ancestry {anc} must be one of {valid_ancestries}.")
            sys.exit(1)

    # Confirm that output_type is a valid option
    if type(output_type) == str:
        output_type = output_type.upper()
        if output_type not in ['MATRIX', 'MAP', 'BOTH']:
            print("ERROR: Output Type must be Matrix, Map, or Both.")
            sys.exit(1)
    else:
        print("ERROR: Output Type must be a string equivalent to Matrix/Map/Both.")
        sys.exit(1)

    # Confirm that the file locs exist
    if exists(variant_dir) == False:
        print("ERROR: Input directory of variant files does not exist.")
        sys.exit(1)

    if exists(individual_dir) == False:
        print("ERROR: Input directory of files listing individuals by ancestry does not exist.")
        sys.exit(1)

    if exists(register_file) == False:
        print("ERROR: Input file detailing breakdown of plink file ranges does not exist.")
        sys.exit(1)

    # Check that directories exist
    if exists(output_dir) == False:
        print("ERROR: Output directory does not exist.")
        sys.exit(1)

    if exists(plink_dir) == False:
        print("ERROR: Invalid directory given for plink files' location.")
        sys.exit(1)

    if exists(plink_bin_dir) == False:
        print("ERROR: Invalid directory given for plink binary location.")
        sys.exit(1)

    print("Acceptable Inputs Given")

    # Call driver function
    driver(variant_dir, output_dir, output_type, ancestry_list, individual_dir, 
           register_file, command_file, plink_dir, plink_bin_dir)

###############################################################################
############################# HELPFUL FUNCTIONS ##############################
###############################################################################

# Function that adds slashes to directory names
def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    else:
        return directory

# Function to parse variant ID in format chr:pos:ref:alt
def parse_variant_id(variant_id):
    """
    Parse variant ID in format chr:pos:ref:alt
    Returns: (chr, pos, ref, alt)
    """
    parts = variant_id.split(':')
    if len(parts) != 4:
        return None
    try:
        chr_num = int(parts[0])
        pos = int(parts[1])
        ref = parts[2].upper()
        alt = parts[3].upper()
        return (chr_num, pos, ref, alt)
    except:
        return None

# Function to get reverse complement of allele
def reverse_complement(allele):
    """Get reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(allele))

# Function to flip alleles (swap ref and alt)
def flip_alleles(ref, alt):
    """Swap reference and alternate alleles"""
    return alt, ref

###############################################################################
############################# DRIVER ########################################
###############################################################################

## drive the script ##
def driver(variant_dir, output_dir, output_type, ancestry_list, individual_dir, 
           register_file, command_file, plink_dir, plink_bin_dir):

    # Set plink executables
    plink2 = plink_bin_dir.rstrip('/') + "/plink2"

    # Add slashes to directories
    output_dir = add_slash(output_dir)
    plink_dir = add_slash(plink_dir)
    individual_dir = add_slash(individual_dir)
    variant_dir = add_slash(variant_dir)

    # Read in the registry file
    register_raw = pd.read_csv(register_file, sep=",")
    
    # Rename columns if they have spaces
    register_raw.columns = register_raw.columns.str.strip()
    
    # Add dose extension to chunk names if not already present
    if 'pgen' in register_raw.columns:
        register_raw["Chunk"] = register_raw["pgen"].astype(str) + '.dose'
    elif 'Chunk' in register_raw.columns:
        register_raw["Chunk"] = register_raw["Chunk"].astype(str) + '.dose'

    # Create a list to hold all commands across all ancestries
    all_commands = []

    # Loop over each ancestry
    for ancestry in ancestry_list:
        
        print(f"\nProcessing ancestry: {ancestry}")
        
        # Get the variant directory for this ancestry
        anc_variant_dir = variant_dir + ancestry + "/"
        if not exists(anc_variant_dir):
            print(f"WARNING: Variant directory not found for {ancestry}: {anc_variant_dir}")
            print(f"Skipping {ancestry}")
            continue
        
        # Get list of variant files for this ancestry (*.vars files)
        variant_files = [f for f in os.listdir(anc_variant_dir) if f.endswith('.vars')]
        if not variant_files:
            print(f"WARNING: No .vars files found in {anc_variant_dir}")
            print(f"Skipping {ancestry}")
            continue
        
        print(f"Found {len(variant_files)} variant files for {ancestry}")
        
        # Create ancestry-specific output subdirectory
        anc_output_dir = output_dir + ancestry + "/"
        if not exists(anc_output_dir):
            os.makedirs(anc_output_dir)
            print(f"Created directory: {anc_output_dir}")

        # Get the sample list for this ancestry
        sample_list_file = individual_dir + ancestry + "/sample_list.txt"
        if not exists(sample_list_file):
            print(f"WARNING: Sample list file not found for {ancestry}: {sample_list_file}")
            print(f"Skipping {ancestry}")
            continue

        # Create a blank list that will hold the commands for this ancestry
        anc_commands = []

        # Process each variant file for this ancestry
        for variant_file in sorted(variant_files):
            # Extract the file prefix (remove .vars extension)
            file_prefix = variant_file.replace('.vars', '')
            
            # Check if output files already exist
            vcor1_file = anc_output_dir + file_prefix + '.vcor1'
            vars_file = anc_output_dir + file_prefix + '.vcor1.vars'
            if exists(vcor1_file) and exists(vars_file):
                print(f"    NOTE: Output files already exist for {file_prefix}. Skipping.")
                continue
            
            # Full path to variant file
            variant_file_path = anc_variant_dir + variant_file
            print(f"  Processing variant file: {variant_file}")
            
            # Read in the variant file for this phenotype
            variants_raw = []
            with open(variant_file_path) as var_read:
                for line in var_read.readlines():
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parsed = parse_variant_id(line)
                    if parsed:
                        variants_raw.append(parsed)
            
            if not variants_raw:
                print(f"    WARNING: No valid variants found in {variant_file}")
                continue
            
            print(f"    Found {len(variants_raw)} variants")

            # Group variants by chromosome for efficient plink file selection
            variants_by_chr = {}
            for variant in variants_raw:
                chr_num = variant[0]
                if chr_num not in variants_by_chr:
                    variants_by_chr[chr_num] = []
                variants_by_chr[chr_num].append(variant)

            # Process each chromosome
            for chr_num in sorted(variants_by_chr.keys()):
                chr_variants = variants_by_chr[chr_num]
                
                # Find the range of variants for this chromosome
                positions = [v[1] for v in chr_variants]
                min_pos = min(positions)
                max_pos = max(positions)
                
                # Identify the plink files that are needed for this chromosome
                temp1 = register_raw[register_raw['Chr'] == chr_num]
                if temp1.empty:
                    print(f"WARNING: No plink files found for chromosome {chr_num}. Skipping.")
                    continue
                
                # Find the chunk that contains min_pos
                temp2 = temp1[temp1['Start'] <= min_pos]
                if temp2.empty:
                    temp2 = temp1.iloc[0, :]
                else:
                    temp2 = temp2.loc[temp2.index[-1], :]
                
                # Find the chunk that contains max_pos
                temp3 = temp1[temp1['End'] >= max_pos]
                if temp3.empty:
                    temp3 = temp1.iloc[len(temp1)-1, :]
                else:
                    temp3 = temp3.loc[temp3.index[0], :]
                
                # Create indices for plink files to use
                lower_index = temp1[temp1['Chunk'] == temp2['Chunk']].index[0] if isinstance(temp2, pd.Series) else temp2.name
                upper_index = temp1[temp1['Chunk'] == temp3['Chunk']].index[0] if isinstance(temp3, pd.Series) else temp3.name
                
                plink_count = upper_index - lower_index + 1
                                
                # Create command string for this chromosome's variants
                cmd = f'echo "START: {file_prefix} `date`"; '
                cmd += f'echo "HOST: {file_prefix} $HOSTNAME"; '
                
                # Build plink command
                if plink_count > 1:
                    print(f'NOTE: Merging {plink_count} plink files for chr{chr_num} ({file_prefix})')
                    
                    # Export files to a merge list
                    temp_chunks = register_raw.loc[lower_index:upper_index, 'Chunk']
                    merge_list_file = plink_dir + file_prefix + '.temp.merge_list.txt'
                    temp_chunks.to_csv(merge_list_file, sep="\t", header=False, index=False, na_rep='NA')
                    
                    # Add merge command
                    cmd += f'{plink2} --pmerge-list {merge_list_file} --keep {sample_list_file} --chr {chr_num} --from-bp {min_pos} --to-bp {max_pos} --make-bed --out {plink_dir}{file_prefix}.temp --memory 200000; '
                else:
                    print(f'NOTE: Single plink file needed for chr{chr_num} ({file_prefix})')
                    chunk_name = register_raw.loc[lower_index, 'Chunk']
                    cmd += f'{plink2} --pfile {plink_dir}{chunk_name} --keep {sample_list_file} --chr {chr_num} --from-bp {min_pos} --to-bp {max_pos} --make-bed --out {plink_dir}{file_prefix}.temp --memory 200000; '
                
                # Add variant extraction and matrix creation command
                # This will be handled by a python script that processes individual variants
                cmd += f'python3 {os.path.dirname(os.path.realpath(__file__))}/extract_matrices.py '
                cmd += f'-b {plink_dir}{file_prefix}.temp '
                cmd += f'-v {variant_file_path} '
                cmd += f'-c {chr_num} '
                cmd += f'-o {anc_output_dir} '
                cmd += f'-p {file_prefix}; '
                
                # Rename output file extensions
                cmd += f'mv {anc_output_dir}/{file_prefix}.ld {anc_output_dir}/{file_prefix}.vcor1;'
                cmd += f'mv {anc_output_dir}/{file_prefix}.variants.txt {anc_output_dir}/{file_prefix}.vcor1.vars;'
                # Clean up temporary files
                cmd += f'rm {plink_dir}{file_prefix}.temp.*; '
                cmd += f'rm {anc_output_dir}/{file_prefix}.log;'
                cmd += f'rm {anc_output_dir}/{file_prefix}.nosex'
    
                anc_commands.append(cmd)
        
        #Add ancestry commands to all commands
        all_commands.extend(anc_commands)

    # Write master command file
    master_command_file = output_dir + "make_matrices.master.txt"
    with open(master_command_file, 'w') as f:
        for cmd in all_commands:
            f.write(cmd + '\n')
    
    print(f"\nWrote {len(all_commands)} total commands to {master_command_file}")

###############################################################################
############################ CALL MAIN ######################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv)
