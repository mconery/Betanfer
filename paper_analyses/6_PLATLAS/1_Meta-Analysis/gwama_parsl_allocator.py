'''
gwama_parsl_allocator.py

The purpose of this script is to act as a wrapper for run_single_gwama.py and
manage the multi-node distribution of GWAMA meta-analysis jobs using parsl.

Each command in the supplied command file corresponds to one trait/population
meta-analysis.  The allocator:
  1. Reads the command file produced by create_gwama_commands.py.
  2. Determines which jobs have already produced a .done sentinel file and skips
     them, allowing safe resubmission of the outer PBS job without duplication.
  3. Submits the remaining commands to parsl workers via @bash_app.
  4. Blocks until all futures are resolved.

The parsl configuration (gwama_meta_analysis_config) is imported from the
local parsl_config.py, which uses a LocalProvider over pre-allocated nodes.

Usage:
  python gwama_parsl_allocator.py --command_file <path> [--log_dir <path>]

***Run in Conda finemap_old environment***
'''

import sys
import getopt
import os
import parsl
from parsl import bash_app
from parsl_config import gwama_meta_analysis_config

# Default path for parsl worker log output
DEFAULT_LOG_DIR = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"
# Default base output directory (used when --output_dir is not in the command)
OUTPUT_BASE = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results"


def help(exit_num=0):
    print("""Usage: python gwama_parsl_allocator.py [OPTIONS]
Options:
  --command_file  File of run_single_gwama.py commands to execute (required)
  --log_dir       Directory for parsl worker log output
                  (default: standard script_logs path; optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
############  HELPER: DETERMINE EXPECTED OUTPUT PATH FROM A COMMAND  ##########
###############################################################################

def get_output_path(command):
    """Return the expected .done sentinel path for a run_single_gwama.py command.

    Parses --trait_pop and (optionally) --output_dir from the command string.
    Returns None if the command cannot be parsed.
    """
    parts = command.split()
    trait_pop  = None
    output_dir = OUTPUT_BASE

    for i, part in enumerate(parts):
        if part == "--trait_pop" and i + 1 < len(parts):
            trait_pop = parts[i + 1]
        elif part == "--output_dir" and i + 1 < len(parts):
            output_dir = parts[i + 1]

    if trait_pop is None:
        return None

    # trait_pop is e.g. "Phe_008.AFR" — split on the last '.' to get pop
    dot_parts = trait_pop.rsplit('.', 1)
    if len(dot_parts) != 2:
        return None
    trait_name, pop = dot_parts
    return os.path.join(output_dir, pop, trait_name + ".done")


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

    # ── Skip commands whose output already exists ──────────────────────────────
    pending  = []
    skipped  = 0
    unparsed = 0
    for cmd in all_commands:
        out_path = get_output_path(cmd)
        if out_path is None:
            print(f"WARNING: Could not parse output path for command: {cmd[:80]}")
            unparsed += 1
            pending.append(cmd)   # run it anyway; run_single_gwama will skip if done
        elif os.path.exists(out_path):
            skipped += 1
        else:
            pending.append(cmd)

    print(f"NOTE: {len(all_commands)} total commands  |  "
          f"{skipped} already complete  |  {len(pending)} to run")
    if unparsed:
        print(f"WARNING: {unparsed} command(s) could not be parsed for output check")

    if not pending:
        print("NOTE: All jobs already complete.  Exiting.")
        sys.exit(0)

    # ── Define parsl bash_app ──────────────────────────────────────────────────
    @bash_app
    def run_meta_analysis(command):
        return command

    # ── Change to log dir so parsl worker logs are co-located with other logs ──
    os.chdir(log_dir)
    print("NOTE: Starting GWAMA meta-analyses via parsl")

    # ── Load parsl config and dispatch all pending jobs ────────────────────────
    with parsl.load(gwama_meta_analysis_config):
        futures = [run_meta_analysis(cmd) for cmd in pending]
        # Block until every future resolves (success or retried failure)
        results = [f.result() for f in futures]

    print(f"SUCCESS: All {len(pending)} GWAMA meta-analysis jobs complete")


if __name__ == "__main__":
    main()
