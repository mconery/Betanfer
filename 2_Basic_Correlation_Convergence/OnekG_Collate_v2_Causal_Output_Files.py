'''

OnekG_Collate_v2_Causal_Output_Files.py (Modified for LD-based naming convention)

The purpose of this script is to combine together the outputs of all the 1000G
Simulations that have been executed with the 2nd version of the causal variant-
selection based script.

***Run in Conda finemap environment***

Modified to handle new naming convention:
chr#.POP.#_traits.#_pairs.trait_type.subset_<trait_subset>.<fixed or variable>_<effect_size>.csv

'''

import sys
import os
import getopt
import pandas as pd
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

def help():
    print("""Usage: python combine_simulations.py -i <input_dir> -o <output_file>

Options:
    -i <path>       Path to input directory of simulation results (required)
    -o <path>       Name of output file (required)
    -m <name>       Name of missing data file (default is failed_simulations.tsv)
    -h, --help      Show this help message""")
    sys.exit(2)

###############################################################################
#################### PROCESS INDIVIDUAL SIMULATION FILE ######################
###############################################################################

def process_file(file_path):
    """Process individual file and extract metadata from filename"""
    try:
        # Extract base filename without extension
        fname = Path(file_path).stem
        
        # Split filename into components using the new naming convention
        # New format: chr#.POP.#_traits.trait_type.subset_<trait_subset>.<fixed or variable>_<effect_size>.csv
        parts = fname.split('.')
        
        # Basic components
        chrom, pop, n_traits, n_pairs, trait_type_part = parts[0], parts[1], parts[2], parts[3], parts[4]
        
        # Initialize metadata dictionary with new structure
        meta = {
            'chromosome': chrom,
            'population': pop,
            'n_traits': int(n_traits.replace('_traits', '')),
            'trait_type': None,
            'case_prob': None,
            'trait_subset': None,
            'effect_model': None,  # 'fixed' or 'variable'
            'effect_size': None
        }
        
        # Parse trait type and case probability (if binary)
        # e.g. "binary_0.5" or "quant"
        if '_' in trait_type_part:
            meta['trait_type'], case_prob = trait_type_part.split('_', 1)
            meta['case_prob'] = float(case_prob)
        else:
            meta['trait_type'] = trait_type_part
        
        # Handle subset_<trait_subset> where trait_subset can contain a decimal point,
        # which may have been split by '.' into multiple parts.
        # We know there is a literal "subset_" token to anchor on.
        subset_idx = None
        for i, p in enumerate(parts):
            if p.startswith('subset_'):
                subset_idx = i
                break
        if subset_idx is None:
            raise ValueError(f"No 'subset_' segment found in filename: {fname}")
        
        # Extract the trait_subset value, allowing for decimals split across dots.
        # Example segments: ["subset_0", "25", "fixed_0.1", ...]
        subset_rest = parts[subset_idx].replace('subset_', '')
        
        # Collect digits/decimal from subsequent segments until we reach one that
        # contains 'fixed_' or 'variable_' (the effect model marker).
        subset_pieces = [subset_rest]
        effect_segment_idx = None
        for j in range(subset_idx + 1, len(parts)):
            if parts[j].startswith('fixed_') or parts[j].startswith('variable_'):
                effect_segment_idx = j
                break
            subset_pieces.append(parts[j])
        
        trait_subset_str = '.'.join(subset_pieces)
        meta['trait_subset'] = float(trait_subset_str)
        
        # Now parse effect model and effect size from the effect segment
        # e.g. "fixed_0.1" or "variable_0.2"
        if effect_segment_idx is None:
            raise ValueError(f"No fixed_/variable_ segment found in filename: {fname}")
        
        effect_part = parts[effect_segment_idx]
        if effect_part.startswith('fixed_'):
            meta['effect_model'] = 'fixed'
            effect_val = effect_part.replace('fixed_', '')
        elif effect_part.startswith('variable_'):
            meta['effect_model'] = 'variable'
            effect_val = effect_part.replace('variable_', '')
        else:
            raise ValueError(f"Unrecognized effect model segment: {effect_part}")
        
        meta['effect_size'] = float(effect_val)
        
        # Read CSV and add metadata columns
        df = pd.read_csv(file_path, sep='\t')
        for k, v in meta.items():
            df[k] = v
        
        # Check if the beta correlations got processed for any instances
        if 'beta_corr' in df.columns and df['beta_corr'].isnull().all():
            missing_df = pd.DataFrame([meta])
        else:
            missing_df = pd.DataFrame([{k: None for k in meta.keys()}])
        
        return [df, missing_df]
        
    except Exception as e:
        print(f"ERROR: Error processing {file_path}: {str(e)}")
        return [pd.DataFrame(), pd.DataFrame()]

###############################################################################
############################ COLLATION FUNCTIONS #############################
###############################################################################

def combine_simulation_outputs(input_dir, output_path, missing_path, n_workers=cpu_count()):
    files = list(Path(input_dir).glob('*.csv'))
    print(f"NOTE: Found {len(files)} simulation files")
    
    with Pool(n_workers) as pool:
        results = list(tqdm(pool.imap(process_file, files), total=len(files)))
    
    combined_df = pd.concat([x[0] for x in results], ignore_index=True)
    combined_df.to_csv(output_path, index=False)
    
    missing_df = pd.concat([x[1] for x in results], ignore_index=True).dropna()
    missing_df.to_csv(missing_path, index=False)
    
    print(f"SUCCESS: Combined data saved to {output_path}")

def main(argv):
    input_dir = None
    output_path = None
    missing_name = "failed_simulations.tsv"
    
    try:
        opts, args = getopt.getopt(argv, "hi:o:m:", ["help", "input_dir=", "output_path=", "missing_file="])
    except getopt.GetoptError:
        help()
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            help()
        elif opt in ("-i", "--input_dir"):
            input_dir = arg
        elif opt in ("-o", "--output_path"):
            output_path = arg
        elif opt in ("-m", "--missing_file"):
            missing_name = arg
    
    if not input_dir or not output_path or not missing_name:
        help()
    
    if not os.path.exists(input_dir):
        print("ERROR: Input directory does not exist.")
        sys.exit(1)
    elif not os.path.exists(os.path.dirname(output_path)):
        print("ERROR: Directory of output file does not exist.")
        sys.exit(1)
    
    # Create full path for missing file
    missing_path = os.path.join(os.path.dirname(output_path), missing_name)
    
    # No errors detected so call function
    combine_simulation_outputs(input_dir, output_path, missing_path)

if __name__ == "__main__":
    main(sys.argv[1:])
