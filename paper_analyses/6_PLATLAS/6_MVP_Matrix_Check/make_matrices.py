'''

make_matrices.py

The purpose of this script is to read a sampled locus list (CHR, START, END)
and generate plink commands to compute LD matrices for each (locus, population)
pair using the MVP PGEN files on Frontier.

Unlike the variant-list-based pipeline in 8_MVP_LD_Mapping/make_matrices.py,
this script extracts all variants within each locus region by coordinates
(--from-bp / --to-bp), so no pre-existing .vars files are required and
extract_matrices.py is not needed.

Per-job command flow (written as a single semicolon-chained line):
  1. plink2: extract region by coordinates -> temp BED/BIM/FAM
             (single chunk or pmerge-list for multi-chunk loci)
  2. plink1.9: --r square --keep-allele-order -> LD matrix (.ld file)
  3. Rename .ld -> .unphased.vcor1
  4. awk: parse temp .bim -> .unphased.vcor1.vars  (chr:bp:a1:a2 format)
  5. Cleanup temp files

Sentinel-based checkpointing: jobs whose .unphased.vcor1.vars already
exists are skipped, so resubmission is safe.

Commands are sorted ascending by locus size (smallest first).

Usage:
  python make_matrices.py -l <locus_csv> -o <output_dir> -i <individual_dir>
                          -r <registry_file> -c <command_file>
                         [-a <populations>] [-p <plink_dir>] [-b <plink_bin_dir>]

***Run on Frontier (Python/3.9, not Aurora)***

'''

import sys
import getopt
from os.path import exists
import os
import pandas as pd


# ── Defaults ────────────────────────────────────────────────────────────────
DEFAULT_PLINK_DIR     = "/gpfs/arx2/bif154_mde/proj-shared/RITM0309091/pgen"
DEFAULT_PLINK_BIN_DIR = "/gpfs/arx2/bif154_mde/proj-shared/tools/bin"
DEFAULT_POPULATIONS   = "AFR,EUR"
MAF_FILTER            = "0.01"
GENO_FILTER           = "0.1"
PLINK_MEMORY          = "200000"


def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
-l => Sampled locus CSV with columns CHR,START,END REQUIRED
-o => Base output directory (population subdirs created automatically) REQUIRED
-i => Individual directory: base dir containing {POP}/sample_list.txt REQUIRED
-r => Registry CSV file (hg19.chunk_register.csv) REQUIRED
-c => Output command file path REQUIRED
-a => Populations comma-separated (default: AFR,EUR) OPTIONAL
-p => Directory containing plink PGEN files
      (default: /gpfs/arx2/bif154_mde/proj-shared/RITM0309091/pgen) OPTIONAL
-b => Directory containing plink/plink2 binaries
      (default: /gpfs/arx2/bif154_mde/proj-shared/tools/bin) OPTIONAL
-h => Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "l:o:i:r:c:a:p:b:h")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    if '-h' in options_dict:
        help(0)

    # Required arguments
    try:
        locus_csv      = options_dict['-l']
        output_dir     = options_dict['-o']
        individual_dir = options_dict['-i']
        registry_file  = options_dict['-r']
        command_file   = options_dict['-c']
    except KeyError:
        print("ERROR: One or more required arguments are missing.")
        help()

    # Optional arguments
    populations   = [p.strip().upper()
                     for p in options_dict.get('-a', DEFAULT_POPULATIONS).split(',')]
    plink_dir     = options_dict.get('-p', DEFAULT_PLINK_DIR)
    plink_bin_dir = options_dict.get('-b', DEFAULT_PLINK_BIN_DIR)

    # Validate paths
    for path, label in [
        (locus_csv,      "locus CSV"),
        (individual_dir, "individual directory"),
        (registry_file,  "registry file"),
        (plink_dir,      "plink directory"),
        (plink_bin_dir,  "plink binary directory"),
    ]:
        if not exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    valid_populations = {'AFR', 'AMR', 'EAS', 'EUR'}
    for pop in populations:
        if pop not in valid_populations:
            print(f"ERROR: Population '{pop}' not in valid set {valid_populations}.")
            sys.exit(1)

    print("Acceptable Inputs Given")

    driver(locus_csv, output_dir, populations, individual_dir,
           registry_file, command_file, plink_dir, plink_bin_dir)


###############################################################################
############################# HELPFUL FUNCTIONS ##############################
###############################################################################

def add_slash(directory):
    if directory[-1] != '/':
        return directory + '/'
    return directory


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(locus_csv, output_dir, populations, individual_dir,
           registry_file, command_file, plink_dir, plink_bin_dir):

    # Normalise directory paths
    output_dir     = add_slash(output_dir)
    plink_dir      = add_slash(plink_dir)
    individual_dir = add_slash(individual_dir)

    plink2 = plink_bin_dir.rstrip('/') + "/plink2"
    plink  = plink_bin_dir.rstrip('/') + "/plink"

    # ── Read registry ──────────────────────────────────────────────────────
    print(f"NOTE: Reading registry file: {registry_file}")
    register_raw = pd.read_csv(registry_file, sep=",")
    register_raw.columns = register_raw.columns.str.strip()

    # Normalise chunk column name and append .dose extension
    if 'pgen' in register_raw.columns:
        register_raw["Chunk"] = register_raw["pgen"].astype(str) + '.dose'
    elif 'Chunk' in register_raw.columns:
        register_raw["Chunk"] = register_raw["Chunk"].astype(str) + '.dose'
    else:
        print("ERROR: Registry file has neither 'pgen' nor 'Chunk' column.")
        sys.exit(1)

    # ── Read sampled locus list ────────────────────────────────────────────
    print(f"NOTE: Reading sampled locus file: {locus_csv}")
    loci = pd.read_csv(locus_csv)

    required_cols = {'CHR', 'START', 'END'}
    if not required_cols.issubset(loci.columns):
        print(f"ERROR: Locus CSV missing required columns. Expected: {required_cols}")
        sys.exit(1)

    print(f"NOTE: {len(loci)} loci loaded from {locus_csv}")

    # ── Generate commands ──────────────────────────────────────────────────
    all_commands = []   # list of (locus_size, command_str)

    for _, locus in loci.iterrows():
        chr_num = int(locus['CHR'])
        start   = int(locus['START'])
        end     = int(locus['END'])
        size    = end - start

        for pop in populations:
            # Output paths
            pop_output_dir = output_dir + pop + "/"
            locus_prefix   = f"chr{chr_num}_{start}_{end}"
            output_prefix  = pop_output_dir + locus_prefix
            sentinel       = output_prefix + ".unphased.vcor1.vars"

            # Sentinel check — skip completed jobs
            if exists(sentinel):
                print(f"    NOTE: Sentinel exists, skipping: {locus_prefix}.{pop}")
                continue

            # Sample list for this population
            sample_list = individual_dir + pop + "/sample_list.txt"
            if not exists(sample_list):
                print(f"WARNING: Sample list not found for {pop}: {sample_list} — skipping")
                continue

            # ── Chunk discovery ────────────────────────────────────────────
            chr_chunks = register_raw[register_raw['Chr'] == chr_num]
            if chr_chunks.empty:
                print(f"WARNING: No plink files found for chr{chr_num} — skipping {locus_prefix}")
                continue

            # Chunk containing start position (or first chunk on chromosome)
            lower_candidates = chr_chunks[chr_chunks['Start'] <= start]
            if lower_candidates.empty:
                lower_row = chr_chunks.iloc[0]
            else:
                lower_row = lower_candidates.iloc[-1]

            # Chunk containing end position (or last chunk on chromosome)
            upper_candidates = chr_chunks[chr_chunks['End'] >= end]
            if upper_candidates.empty:
                upper_row = chr_chunks.iloc[-1]
            else:
                upper_row = upper_candidates.iloc[0]

            lower_index = lower_row.name
            upper_index = upper_row.name
            plink_count = upper_index - lower_index + 1

            # Temp prefix stored in plink_dir (cleaned up after each job)
            temp_prefix = plink_dir + f"{locus_prefix}.{pop}.temp"

            # ── Build command string ───────────────────────────────────────
            cmd  = f'echo "START: {locus_prefix}.{pop} `date`"; '
            cmd += f'echo "HOST: {locus_prefix}.{pop} $HOSTNAME"; '
            cmd += f'mkdir -p {pop_output_dir}; '

            if plink_count > 1:
                print(f"NOTE: {locus_prefix}.{pop}: merging {plink_count} plink chunks")
                # Write merge list inline via printf
                chunk_paths = [plink_dir + str(register_raw.loc[i, 'Chunk'])
                               for i in range(lower_index, upper_index + 1)]
                merge_list_file = plink_dir + f"{locus_prefix}.{pop}.temp.merge_list.txt"
                # Create merge list file as part of the command
                merge_lines = "\\n".join(chunk_paths)
                cmd += f'printf "{merge_lines}\\n" > {merge_list_file}; '
                cmd += (f'{plink2} --pmerge-list {merge_list_file} '
                        f'--chr {chr_num} --from-bp {start} --to-bp {end} '
                        f'--keep {sample_list} '
                        f'--maf {MAF_FILTER} --geno {GENO_FILTER} '
                        f'--make-bed --out {temp_prefix}; ')
                cmd += f'rm {merge_list_file}; '
            else:
                chunk_name = register_raw.loc[lower_index, 'Chunk']
                cmd += (f'{plink2} --pfile {plink_dir}{chunk_name} '
                        f'--chr {chr_num} --from-bp {start} --to-bp {end} '
                        f'--keep {sample_list} '
                        f'--maf {MAF_FILTER} --geno {GENO_FILTER} '
                        f'--make-bed --out {temp_prefix}; ')

            # LD computation
            cmd += (f'{plink} --bfile {temp_prefix} '
                    f'--r square --keep-allele-order '
                    f'--out {output_prefix} --memory {PLINK_MEMORY}; ')

            # Rename matrix file
            cmd += f'mv {output_prefix}.ld {output_prefix}.unphased.vcor1; '

            # Create vars file from BIM (format: chr:bp:a1:a2)
            cmd += (f"awk '{{print $1\":\"$4\":\"$5\":\"$6}}' "
                    f"{temp_prefix}.bim > {output_prefix}.unphased.vcor1.vars; ")

            # Cleanup
            cmd += f'rm -f {temp_prefix}.*; '
            cmd += f'rm -f {output_prefix}.log; '
            cmd += f'echo "DONE: {locus_prefix}.{pop} `date`"'

            all_commands.append((size, cmd))

    print(f"NOTE: {len(all_commands)} commands generated "
          f"({len(loci)} loci × {len(populations)} populations)")

    # Sort ascending by locus size (small loci first)
    all_commands.sort(key=lambda x: x[0])

    # ── Write command file ─────────────────────────────────────────────────
    cmd_dir = os.path.dirname(command_file)
    if cmd_dir:
        os.makedirs(cmd_dir, exist_ok=True)

    with open(command_file, 'w') as fh:
        for _, cmd_str in all_commands:
            fh.write(cmd_str + '\n')

    print(f"NOTE: Commands written to {command_file}")


###############################################################################
################################ CALL MAIN ###################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv)
