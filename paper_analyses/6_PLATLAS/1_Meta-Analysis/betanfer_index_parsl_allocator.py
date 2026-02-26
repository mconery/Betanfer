'''
betanfer_index_parsl_allocator.py

The purpose of this script is to manage the multi-node distribution of
betanfer index jobs using parsl.

Each command in the supplied command file corresponds to one .out file.
The allocator:
  1. Reads the command file produced by create_index_commands.py.
  2. Determines which jobs have already produced a .gz.tbi index file and
     skips them, allowing safe resubmission of the outer PBS job without
     duplication.
  3. Submits the remaining commands to parsl workers via @bash_app.
  4. Blocks until all futures are resolved.

The parsl configuration (betanfer_index_config) is imported from the
local parsl_config_betanfer.py, which uses a LocalProvider over pre-
allocated nodes.

Usage:
  python betanfer_index_parsl_allocator.py --command_file <path> [--log_dir <path>]

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os
import parsl
from parsl import bash_app
from parsl_config_betanfer import betanfer_index_config

# Default path for parsl worker log output
DEFAULT_LOG_DIR = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"


def help(exit_num=0):
    print("""Usage: python betanfer_index_parsl_allocator.py [OPTIONS]
Options:
  --command_file  File of betanfer index commands to execute (required)
  --log_dir       Directory for parsl worker log output
                  (default: standard script_logs path; optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
############  HELPER: DETERMINE EXPECTED OUTPUT PATH FROM A COMMAND  ##########
###############################################################################

def get_tbi_path(command):
    """Return the expected .gz.tbi path for a betanfer index command.

    Parses --file from the command string.
    Returns None if the command cannot be parsed.
    """
    parts = command.split()
    for i, part in enumerate(parts):
        if part == "--file" and i + 1 < len(parts):
            return parts[i + 1] + ".gz.tbi"
    return None


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

    # ── Skip commands whose .gz.tbi already exists ─────────────────────────────
    pending  = []
    skipped  = 0
    unparsed = 0
    for cmd in all_commands:
        tbi_path = get_tbi_path(cmd)
        if tbi_path is None:
            print(f"WARNING: Could not parse --file from command: {cmd[:80]}")
            unparsed += 1
            pending.append(cmd)   # run it anyway; betanfer --force handles partials
        elif os.path.exists(tbi_path):
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
    def run_betanfer_index(command):
        return command

    # ── Change to log dir so parsl worker logs are co-located with other logs ──
    os.chdir(log_dir)
    print("NOTE: Starting betanfer index runs via parsl")

    # ── Load parsl config and dispatch all pending jobs ────────────────────────
    with parsl.load(betanfer_index_config):
        futures = [run_betanfer_index(cmd) for cmd in pending]
        # Block until every future resolves (success or retried failure)
        results = [f.result() for f in futures]

    print(f"SUCCESS: All {len(pending)} betanfer index jobs complete")


if __name__ == "__main__":
    main()
