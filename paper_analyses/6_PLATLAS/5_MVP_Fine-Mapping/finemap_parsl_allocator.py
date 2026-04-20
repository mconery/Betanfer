'''
finemap_parsl_allocator.py

Distribute MVP fine-mapping commands across pre-allocated HPC nodes using
parsl.

Each command in the supplied file is a bash invocation of run_rsparsepro.sh
or run_susie.sh (produced by create_finemap_commands.py).  The allocator:
  1. Reads the command file.
  2. Identifies commands whose output file already exists and skips them,
     enabling safe resubmission of the outer PBS job without duplication.
  3. Submits remaining commands to parsl workers via @bash_app.
  4. Blocks until all futures resolve, collecting failures fault-tolerantly
     so that one failed locus does not abort the rest.

The parsl configuration (mvp_finemap_config) is imported from the local
parsl_config.py, which uses a LocalProvider over pre-allocated nodes.

The expected output file is the 5th whitespace-separated token (index 4) of
each command:
  bash {script} {sumstats_gz} {ld_file} {output_file} {CHR} {START} {END}
   [0]   [1]        [2]          [3]         [4]         [5]   [6]    [7]

Usage:
  python finemap_parsl_allocator.py --command_file <path> [--log_dir <path>]

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os
import parsl
from parsl import bash_app
from parsl_config import mvp_finemap_config

# Default path for parsl worker log output
DEFAULT_LOG_DIR = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"


def help(exit_num=0):
    print("""Usage: python finemap_parsl_allocator.py [OPTIONS]
Options:
  --command_file  File of fine-mapping commands to execute (required)
  --log_dir       Directory for parsl worker log output
                  (default: standard script_logs path; optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  HELPER: DETERMINE SENTINEL OUTPUT PATH FROM A COMMAND  ############
###############################################################################

def get_output_path(command):
    """Return the expected output file path for a fine-mapping command.

    Command format:
      bash {script} {sumstats_gz} {ld_file} {output_file} {CHR} {START} {END}

    The output_file is the 5th token (index 4).  Returns None on parse failure.
    """
    parts = command.split()
    if len(parts) < 5:
        return None
    return parts[4]


###############################################################################
##########  HELPER: DERIVE A SHORT LABEL FOR LOG FILE NAMING  #################
###############################################################################

def get_label(command):
    """Return a short, filesystem-safe label derived from the output path."""
    output_path = get_output_path(command)
    if output_path is None:
        return "unknown"
    # e.g. ".../rsparsepro/AFR/A1C_Max_INT.chr10_114250001_115250000.rsparsepro.txt"
    # Include the population subdirectory (AFR/EUR) to avoid collisions between
    # the two populations for the same trait+locus, which share the same filename.
    # → "AFR.A1C_Max_INT.chr10_114250001_115250000.rsparsepro"
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
            pending.append(cmd)   # run it anyway; wrapper script will handle gracefully
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
    def run_finemap_job(command, stdout=None, stderr=None):
        return command

    # ── Change to log dir so parsl internal logs land there too ───────────────
    os.chdir(log_dir)
    print("NOTE: Starting fine-mapping jobs via parsl")

    # ── Load parsl config and dispatch all pending jobs ────────────────────────
    with parsl.load(mvp_finemap_config):
        futures = []
        for cmd in pending:
            label = get_label(cmd)
            futures.append(run_finemap_job(
                cmd,
                stdout=os.path.join(log_dir, f"{label}.finemap_out.txt"),
                stderr=os.path.join(log_dir, f"{label}.finemap_err.txt"),
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
              f"see .finemap_err.txt files in {log_dir}")
    else:
        print(f"\nSUCCESS: All {len(pending)} fine-mapping jobs complete")


if __name__ == "__main__":
    main()
