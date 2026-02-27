'''
locus_parsl_allocator.py

The purpose of this script is to manage the multi-node distribution of
locus-defining jobs using parsl, then combine per-trait output files into
per-population master loci lists.

Each command in the supplied command file corresponds to one trait/population
locus-defining run (define_sig_loci_for_trait.py).  The allocator:
  1. Reads the command file produced by create_locus_commands.py.
  2. Determines which jobs have already produced a .sig_loci.csv output file
     and skips them, allowing safe resubmission of the PBS job without
     duplication.
  3. Submits remaining commands to parsl workers via @bash_app.
  4. Blocks until all futures are resolved.
  5. Combines per-trait .sig_loci.csv files into per-population master files,
     each with an OVERLAP_COUNT column computed via IntervalTree.

The parsl configuration (locus_defining_config) is imported from the local
parsl_config.py, which uses a LocalProvider over pre-allocated nodes.

Usage:
  python locus_parsl_allocator.py --command_file <path> [--log_dir <path>]

***Run in Conda finemap_old environment***
'''

import sys
import getopt
import os
from collections import defaultdict
import pandas as pd
from intervaltree import Interval, IntervalTree
import parsl
from parsl import bash_app
from parsl_config import locus_defining_config

# Default path for parsl worker log output
DEFAULT_LOG_DIR = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"
# Default output directory (used when -o is not parseable from a command)
DEFAULT_OUTPUT_DIR = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci"


def help(exit_num=0):
    print("""Usage: python locus_parsl_allocator.py [OPTIONS]
Options:
  --command_file  File of define_sig_loci_for_trait.py commands (required)
  --log_dir       Directory for parsl worker log output
                  (default: standard script_logs path; optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
############  HELPER: DETERMINE EXPECTED SENTINEL PATH FROM A COMMAND  ########
###############################################################################

def parse_command(command):
    """Return (output_dir, base_name) parsed from a define_sig_loci_for_trait.py command.

    Parses -o and -b flags from the command string.
    Returns (None, None) if either cannot be found.
    """
    parts = command.split()
    output_dir = DEFAULT_OUTPUT_DIR
    base_name  = None

    i = 0
    while i < len(parts):
        if parts[i] == '-o' and i + 1 < len(parts):
            output_dir = parts[i + 1]
        elif parts[i] == '-b' and i + 1 < len(parts):
            base_name = parts[i + 1]
        i += 1

    return output_dir, base_name


def get_sentinel_path(command):
    """Return the expected .sig_loci.csv path for a command, or None on parse failure."""
    output_dir, base_name = parse_command(command)
    if base_name is None:
        return None
    return os.path.join(output_dir, f"{base_name}.sig_loci.csv")


###############################################################################
##########  HELPER: COMBINE PER-TRAIT FILES INTO POPULATION MASTER FILE  ######
###############################################################################

def combine_population(pop, base_names, output_dir):
    """Read all per-trait sig_loci.csv files for a population and write master file.

    Adds a TRAIT column and computes OVERLAP_COUNT using IntervalTree.
    Skips traits whose output file is missing or empty.
    """
    frames = []
    for base_name in base_names:
        csv_path = os.path.join(output_dir, f"{base_name}.sig_loci.csv")
        if not os.path.exists(csv_path):
            print(f"WARNING: Missing output file for {base_name}, skipping in master")
            continue
        try:
            df = pd.read_csv(csv_path, sep=",")
        except pd.errors.EmptyDataError:
            continue  # trait had no significant loci
        if df.empty:
            continue
        # Extract trait name: everything before the last '.' in base_name
        trait = base_name.rsplit('.', 1)[0]
        df['TRAIT'] = trait
        frames.append(df)

    if not frames:
        print(f"NOTE: No significant loci found for population {pop} — master file not written")
        return

    combined = pd.concat(frames, ignore_index=True)

    # Compute OVERLAP_COUNT: number of OTHER loci from any trait that overlap each locus
    combined['OVERLAP_COUNT'] = 0
    for chr_name, chr_group in combined.groupby('CHR'):
        tree = IntervalTree()
        for idx, row in chr_group.iterrows():
            tree.add(Interval(row['START'], row['END'], idx))
        for idx, row in chr_group.iterrows():
            overlaps = tree.overlap(row['START'], row['END'])
            combined.loc[idx, 'OVERLAP_COUNT'] = len(overlaps) - 1  # subtract self

    master_path = os.path.join(output_dir, f"sig_loci.{pop}.master.csv")
    combined.to_csv(master_path, index=False)
    print(f"SUCCESS: Wrote master file {master_path} ({len(combined)} loci from {len(frames)} traits)")


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND DISTRIBUTE JOBS  ###############
###############################################################################

def main():

    # ── Parse arguments ────────────────────────────────────────────────────────
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "command_file=", "log_dir="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    command_file = None
    log_dir      = DEFAULT_LOG_DIR

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--command_file":
            command_file = a
        elif o == "--log_dir":
            log_dir = a
        else:
            assert False, f"WARNING: Unhandled input option: {o}"

    # ── Validate required arguments ────────────────────────────────────────────
    if command_file is None:
        print("ERROR: Missing required --command_file")
        sys.exit(2)

    for path, label in [(command_file, "command file"), (log_dir, "log directory")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Read command file ──────────────────────────────────────────────────────
    with open(command_file, 'r') as fh:
        all_commands = [line.rstrip('\n') for line in fh if line.strip()]

    # ── Skip commands whose sentinel output already exists ─────────────────────
    pending  = []
    skipped  = 0
    unparsed = 0
    for cmd in all_commands:
        sentinel = get_sentinel_path(cmd)
        if sentinel is None:
            print(f"WARNING: Could not parse output path for command: {cmd[:80]}")
            unparsed += 1
            pending.append(cmd)  # run it anyway; script will handle gracefully
        elif os.path.exists(sentinel):
            skipped += 1
        else:
            pending.append(cmd)

    print(f"NOTE: {len(all_commands)} total commands  |  "
          f"{skipped} already complete  |  {len(pending)} to run")
    if unparsed:
        print(f"WARNING: {unparsed} command(s) could not be parsed for output check")

    # ── Dispatch pending jobs via parsl ────────────────────────────────────────
    if pending:
        @bash_app
        def run_locus_job(command, stdout=None, stderr=None):
            return command

        os.chdir(log_dir)
        print("NOTE: Starting locus-defining jobs via parsl")

        with parsl.load(locus_defining_config):
            futures = []
            for cmd in pending:
                _, base_name = parse_command(cmd)
                label = base_name or f"unknown_{len(futures)}"
                futures.append(run_locus_job(
                    cmd,
                    stdout=os.path.join(log_dir, f"{label}.locus_out.txt"),
                    stderr=os.path.join(log_dir, f"{label}.locus_err.txt"),
                ))
            failed_cmds = []
            for cmd, future in zip(pending, futures):
                try:
                    future.result()
                except Exception as e:
                    _, base_name = parse_command(cmd)
                    label = base_name or cmd[:60]
                    print(f"ERROR: Task {label} failed: {e}", flush=True)
                    failed_cmds.append(cmd)

        if failed_cmds:
            print(f"WARNING: {len(failed_cmds)} task(s) failed — "
                  f"see .locus_err.txt / .locus_out.txt files in {log_dir}")
        else:
            print(f"SUCCESS: All {len(pending)} locus-defining jobs complete")
    else:
        print("NOTE: All jobs already complete.  Skipping parsl dispatch.")

    # ── Combine per-population master files ────────────────────────────────────
    # Group all commands (complete + newly run) by population.
    # Population = last '.'-delimited token of the -b base_name.
    # e.g. "Phe_008_52.EUR"    → pop "EUR"
    #      "A1C_Max_INT.MVP_EUR" → pop "MVP_EUR"
    pop_bases   = defaultdict(list)   # pop → [base_name, ...]
    pop_outdirs = {}                  # pop → output_dir (should be the same for all)

    for cmd in all_commands:
        output_dir, base_name = parse_command(cmd)
        if base_name is None:
            continue
        pop = base_name.rsplit('.', 1)[-1]
        pop_bases[pop].append(base_name)
        pop_outdirs[pop] = output_dir

    print(f"NOTE: Combining results for {len(pop_bases)} population(s): {sorted(pop_bases.keys())}")
    for pop in sorted(pop_bases.keys()):
        combine_population(pop, pop_bases[pop], pop_outdirs[pop])

    print("NOTE: locus_parsl_allocator.py completed successfully.")


if __name__ == "__main__":
    main()
