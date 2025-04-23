'''
OnekG_Collate_Output_Files.py

The purpose of this script is to combine together the outputs of all the 1000G
Simulations that have been executed.

***Run in Conda finemap environment***
'''

import sys
import os
import getopt
import pandas as pd
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

def help():
    print("""Usage: python combine_simulations.py -i <input_dir> -o <output_path>
      Options:
        -i              Path to input directory of simulation results (required)
        -o              Name of output file (required)
        -m              Name of missing data file (default is failed_simulations.tsv)
        -h, --help      Show this help message""")
    sys.exit(2)

###############################################################################
####################  PROCESS INDIVIDUAL SIMULATION FILE ######################
###############################################################################

def process_file(file_path):
    """Process individual file and extract metadata from filename"""
    try:
        # Extract base filename without extension
        fname = Path(file_path).stem
        
        # Split filename into components using the naming convention
        parts = fname.split('.')
        chrom, pop, n_traits, n_pairs = parts[0], parts[1], parts[2], parts[3]
        #Recombine parts as needed due to decimal splits
        if "causal" not in parts[5]:
            temp = parts.pop(5)
            parts[4] = parts[4] + '.' + temp
        if len(parts) > 6:
            temp = parts.pop(6)
            parts[5] = parts[5] + '.' + temp
        # Initialize metadata dictionary
        meta = {
            'chromosome': chrom,
            'population': pop,
            'n_traits': int(n_traits.replace('_traits', '')),
            'n_pairs': int(n_pairs.replace('_pairs', '')),
            'trait_type': None,
            'case_prob': None,
            'n_causal': 0,
            'effect_size': None
        }
        # Parse trait type and case probability (if binary)
        trait_part = parts[4]
        if '_' in trait_part:
            meta['trait_type'], case_prob = trait_part.split('_')
            meta['case_prob'] = float(case_prob)
        else:
            meta['trait_type'] = trait_part
        # Parse causal variants and effect size if present
        if len(parts) > 5:
            causal_part = parts[5]
            if causal_part[0] == '0':
                n_causal = 0
                meta['n_causal'] = int(n_causal)
            elif '_causal' in causal_part:
                n_causal, effect = causal_part.split('_causal_')
                if n_causal != 0 and effect == 0:
                    effect = '0.' + parts[-1]
                meta['n_causal'] = int(n_causal)
                meta['effect_size'] = float(effect)
            else:
                meta['n_causal'] = int(causal_part.replace('_causal', ''))
        # Read CSV and add metadata columns
        df = pd.read_csv(file_path, sep='\t')
        for k, v in meta.items():
            df[k] = v
        #Check if the beta correlations got processed for any instances
        if df['beta_corr'].isnull().all():
            missing_df = pd.DataFrame([meta])
        else:
            missing_df = pd.DataFrame([{k : None for k in meta.keys()}])
        return [df, missing_df]
    except Exception as e:
        print(f"ERROR: Error processing {file_path}: {str(e)}")
        return [pd.DataFrame(), pd.DataFrame()]

###############################################################################
############################  COLLATION FUNCTIONS #############################
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
    #Create full path for missing file
    missing_path = os.path.join(os.path.dirname(output_path), missing_name)
    #No errors detected so call function
    combine_simulation_outputs(input_dir, output_path, missing_path)

if __name__ == "__main__":
    main(sys.argv[1:])
