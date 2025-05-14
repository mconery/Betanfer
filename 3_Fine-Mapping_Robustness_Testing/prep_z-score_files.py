'''
prep_z-score_files.py

This script filters a file of GWAS summary stats for the SNPs given in a 
ld-matrix variant file

***This needs to be run in the rsparsepro conda environment ***
'''

#!/usr/bin/env python3

import sys
import getopt
from os.path import exists
import pandas as pd

def help(exit_num=1):
    print("""-----------------------------------------------------------------
    ARGUMENTS
    -v => Input .vcor1.vars file (variant list) REQUIRED
    -g => GWAS summary statistics file (CSV/TSV) REQUIRED
    -o => Output file path REQUIRED
    --zscore => Name of Z-score column in GWAS (either this or beta+se required)
    --beta => Name of beta column in GWAS
    --se => Name of standard error column in GWAS
""")
    sys.exit(exit_num)

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:g:o:h", ["zscore=", "beta=", "se=", "help"])
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of flags")
        help()

    options = dict(opts)
    
    if '--help' in options or '-h' in options:
        help(0)
        
    # Validate required arguments
    try:
        vars_file = options['-v']
        gwas_file = options['-g']
        output_file = options['-o']
    except KeyError:
        print("ERROR: Missing required arguments")
        help()

    # Check input files exist
    if not exists(vars_file):
        print(f"ERROR: Variant file {vars_file} not found")
        sys.exit(1)
    if not exists(gwas_file):
        print(f"ERROR: GWAS file {gwas_file} not found")
        sys.exit(1)

    # Validate score calculation method
    zcol = options.get('--zscore')
    beta_col = options.get('--beta')
    se_col = options.get('--se')
    
    if not zcol and not (beta_col and se_col):
        print("ERROR: Must specify either --zscore or both --beta and --se")
        help()

    print("Acceptable inputs given")
    
    # Read variant IDs
    with open(vars_file) as f:
        target_variants = set(line.strip() for line in f)

    # Read GWAS data
    try:
        df = pd.read_csv(gwas_file, sep=None, engine='python')  # auto-detect delimiter
    except Exception as e:
        print(f"ERROR: Error reading GWAS file: {e}")
        sys.exit(1)

    # Validate columns exist
    required_cols = []
    if zcol:
        required_cols = [zcol, 'variant_id']  # Assumes variant ID column name
    else:
        required_cols = [beta_col, se_col, 'variant_id']
    
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        print(f"ERROR: Missing columns in GWAS file: {', '.join(missing)}")
        sys.exit(1)

    # Filter and calculate scores
    filtered = df[df['variant_id'].isin(target_variants)].copy()
    
    if zcol:
        filtered['zscore'] = filtered[zcol]
    else:
        filtered['zscore'] = filtered[beta_col] / filtered[se_col]

    # Write output
    try:
        filtered[['variant_id', 'zscore']].to_csv(
            output_file, 
            sep='\t', 
            index=False, 
            header=False
        )
    except Exception as e:
        print(f"ERROR: Error writing output: {e}")
        sys.exit(1)

    print(f"SUCCESS: Wrote {len(filtered)} variants to {output_file}")

if __name__ == "__main__":
    main(sys.argv[1:])
