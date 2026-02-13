'''

OnekG_Collate_v3_Causal_Output_Files.py (Modified for LD-based naming convention)

The purpose of this script is to combine together the outputs of all the 1000G
Simulations that have been executed with the 3rd version of the causal variant-
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

        # Split filename into components (dot-separated)
        # New format:
        # chr#.POP.#_traits.trait_type.maf_<maf>.prop_<subset>.effect_<effect>.csv
        parts = fname.split('.')

        # Basic components
        chrom = parts[0]             # e.g. "chr1"
        pop = parts[1]               # e.g. "EUR"
        n_traits = parts[2]          # e.g. "10_traits"
        n_pairs = parts[3]

        # Initialize metadata dictionary
        meta = {
            'chromosome': chrom,
            'population': pop,
            'n_traits': int(n_traits.replace('_traits', '')),
            'trait_type': None,
            'case_prob': None,
            'trait_subset': None,
            'maf_threshold': None,
            'effect_size': None
        }

        # Parse trait type and case probability (if binary)
        # e.g. "binary_0.5" or "quant"
        if not parts[5].startswith("maf"):
            trait_type_part = '.'.join(parts[4:6])
        else:
            trait_type_part = parts[4]
        if '_' in trait_type_part:
            meta['trait_type'], case_prob = trait_type_part.split('_', 1)
            meta['case_prob'] = float(case_prob)
        else:
            meta['trait_type'] = trait_type_part

        # We now need to parse:
        #   maf_<maf_threshold>
        #   prop_<trait_subset>
        #   effect_<effect_size>
        #
        # All three numeric pieces can contain decimal points and therefore
        # may be split across multiple dot-separated parts.

        def collect_value(start_idx, prefix, stop_prefixes):
            """
            Collect a possibly decimal value that starts at parts[start_idx]
            with given prefix, concatenating subsequent parts with '.' until
            a part whose prefix is in stop_prefixes is seen.
            Returns (value_str, next_index).
            """
            if not parts[start_idx].startswith(prefix):
                raise ValueError(
                    f"Expected segment starting with '{prefix}' at index "
                    f"{start_idx} in filename: {fname}"
                )

            first = parts[start_idx].replace(prefix, '')
            pieces = [first]
            idx = start_idx + 1

            while idx < len(parts):
                # Check if this part starts the next named field
                if any(parts[idx].startswith(sp) for sp in stop_prefixes):
                    break
                pieces.append(parts[idx])
                idx += 1

            return '.'.join(pieces), idx

        # Find indices of maf_, prop_, effect_ segments
        maf_idx = prop_idx = effect_idx = None
        for i, p in enumerate(parts):
            if p.startswith('maf_'):
                maf_idx = i
            elif p.startswith('prop_'):
                prop_idx = i
            elif p.startswith('effect_'):
                effect_idx = i

        if maf_idx is None or prop_idx is None or effect_idx is None:
            raise ValueError(f"Missing maf_/prop_/effect_ segment in filename: {fname}")

        # Parse maf_threshold
        maf_str, next_after_maf = collect_value(
            maf_idx,
            prefix='maf_',
            stop_prefixes=['prop_', 'effect_']
        )
        meta['maf_threshold'] = float(maf_str)

        # Parse trait_subset
        subset_str, next_after_prop = collect_value(
            prop_idx,
            prefix='prop_',
            stop_prefixes=['effect_']
        )
        meta['trait_subset'] = float(subset_str)

        # Parse effect_size
        # effect_ should be the last named field, so we can collect to the end
        effect_first = parts[effect_idx].replace('effect_', '')
        effect_pieces = [effect_first]
        for j in range(effect_idx + 1, len(parts)):
            effect_pieces.append(parts[j])
        effect_str = '.'.join(effect_pieces)
        meta['effect_size'] = float(effect_str)

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
