'''
rsparsepro_parsl_allocator.py

Distribute GWAS RSparsePro fine-mapping commands across pre-allocated HPC nodes
using parsl.

Each command in the supplied file is a Python invocation of run_rsparsepro.py
(produced by create_rsparsepro_commands.py).  The allocator:
  1. Reads the command file.
  2. Identifies commands whose output file already exists and skips them,
     enabling safe resubmission of the outer PBS job without duplication.
  3. Submits remaining commands to parsl workers via @bash_app.
  4. Blocks until all futures resolve, collecting failures fault-tolerantly
     so that one failed locus does not abort the rest.

The parsl configuration (gwas_finemap_config) is imported from the local
parsl_config.py, which uses a LocalProvider over pre-allocated nodes.

The expected output file is extracted from the --output-file flag in each
command string.

Usage:
  python rsparsepro_parsl_allocator.py --command-file <path> [--log-dir <path>]

***Run with: conda activate finemap (parsl is installed there)***
'''

import sys
import getopt
import os
import parsl
from parsl import bash_app
from parsl_config import gwas_finemap_config

DEFAULT_LOG_DIR = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"


def help(exit_num=0):
    print("""Usage: python rsparsepro_parsl_allocator.py [OPTIONS]
Options:
  --command-file  File of RSparsePro commands to execute (required)
  --log-dir       Directory for parsl worker log output
                  (default: standard script_logs path; optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  HELPER: DETERMINE SENTINEL OUTPUT PATH FROM A COMMAND  ############
###############################################################################

def get_output_path(command):
    """Return the expected output file path from --output-file in a command.

    Scans whitespace-split tokens for the --output-file flag and returns the
    following token as the output path.  Returns None on parse failure.
    """
    parts = command.split()
    for i, part in enumerate(parts):
        if part == "--output-file" and i + 1 < len(parts):
            return parts[i + 1]
    return None


###############################################################################
##########  HELPER: DERIVE A SHORT LABEL FOR LOG FILE NAMING  #################
###############################################################################

def get_label(command):
    """Return a short, filesystem-safe label derived from the output path."""
    output_path = get_output_path(command)
    if output_path is None:
        return "unknown"
    # e.g. ".../rsparsepro/EUR/Phe_XXX.chr1_100_200.rsparsepro.txt"
    # → "EUR.Phe_XXX.chr1_100_200.rsparsepro"
    pop   = os.path.basename(os.path.dirname(output_path))
    label = os.path.splitext(os.path.basename(output_path))[0]
    return f"{pop}.{label}"


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND DISTRIBUTE JOBS  ###############
###############################################################################

def main():

    # ── Parse arguments ────────────────────────────────────────────────────────
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "command-file=", "log-dir="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    command_file = None
    log_dir      = DEFAULT_LOG_DIR

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--command-file":
            command_file = a
        elif o == "--log-dir":
            log_dir = a
        else:
            assert False, f"WARNING: Unhandled input option: {o}"

    # ── Validate required arguments ────────────────────────────────────────────
    if command_file is None:
        print("ERROR: Missing required --command-file")
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
            pending.append(cmd)
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
    def run_rsparsepro_job(command, stdout=None, stderr=None):
        return command

    # ── Change to log dir so parsl internal logs land there too ───────────────
    os.chdir(log_dir)
    print("NOTE: Starting RSparsePro jobs via parsl")

    # ── Load parsl config and dispatch all pending jobs ────────────────────────
    with parsl.load(gwas_finemap_config):
        futures = []
        for cmd in pending:
            label = get_label(cmd)
            futures.append(run_rsparsepro_job(
                cmd,
                stdout=os.path.join(log_dir, f"{label}.rsparsepro_out.txt"),
                stderr=os.path.join(log_dir, f"{label}.rsparsepro_err.txt"),
            ))

        # Block until every future resolves; collect failures without aborting
        failed_cmds = []
        for cmd, future in zip(pending, futures):
            try:
                future.result()
            except Exception as e:
                label = get_label(cmd)
                print(f"ERROR: Task {label} failed: {e}", flush=True)
                failed_cmds.append(cmd)

    if failed_cmds:
        print(f"\nWARNING: {len(failed_cmds)} task(s) failed — "
              f"see .rsparsepro_err.txt files in {log_dir}")
    else:
        print(f"\nSUCCESS: All {len(pending)} RSparsePro jobs complete")


if __name__ == "__main__":
    main()
