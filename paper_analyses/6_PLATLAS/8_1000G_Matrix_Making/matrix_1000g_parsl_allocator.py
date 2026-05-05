'''
matrix_1000g_parsl_allocator.py

Distributes plink2 LD matrix jobs across pre-allocated Aurora nodes using parsl.

Each command in the supplied command file corresponds to one (locus, population)
plink2 --r-unphased run.  The allocator:
  1. Reads the command file produced by create_1000g_commands.py.
  2. Determines which jobs have already produced a .vcor1.vars sentinel file
     and skips them, allowing safe resubmission without duplication.
  3. Submits remaining commands to parsl workers via @bash_app, with per-task
     stdout/stderr captured to named log files.
  4. Collects results fault-tolerantly: failed tasks are logged but do not abort
     the remaining work.

The parsl configuration (matrix_1000g_config) is imported from the local
parsl_config.py, which uses a LocalProvider with 1 worker per node on Aurora.

Usage:
  python matrix_1000g_parsl_allocator.py --command_file <path> [--log_dir <path>]

***Run with: conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap***
'''

import sys
import getopt
import os
import parsl
from parsl import bash_app
from parsl_config import matrix_1000g_config

DEFAULT_LOG_DIR = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"


def help(exit_num=0):
    print("""Usage: python matrix_1000g_parsl_allocator.py [OPTIONS]
Options:
  --command_file  File of plink2 commands to execute (required)
  --log_dir       Directory for parsl worker log output
                  (default: standard script_logs path; optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
############  HELPER: DETERMINE EXPECTED SENTINEL PATH FROM A COMMAND  ########
###############################################################################

def get_sentinel_path(command):
    """Return the expected .vcor1.vars sentinel path for a plink2 command.

    Parses --out {prefix} from the command string and appends
    '.unphased.vcor1.vars' (plink2 --r-unphased appends 'unphased' to the
    vcor1 extension).
    Returns None if --out cannot be found.
    """
    parts = command.split()
    for i, part in enumerate(parts):
        if part == "--out" and i + 1 < len(parts):
            return parts[i + 1] + ".unphased.vcor1.vars"
    return None


def get_task_label(command):
    """Return a short label for per-task log filenames.

    Uses the basename of the --out prefix (e.g.
    'chr10_114250001_115250000.1000G.AFR'), prepended with the parent directory
    basename (e.g. 'A1C_Max_INT.MVP_AFR'), to give a unique, readable label.
    """
    sentinel = get_sentinel_path(command)
    if sentinel is None:
        return None
    # strip .vcor1.vars to get output prefix, then take the last two path components
    prefix   = sentinel[:-len(".vcor1.vars")]
    analysis = os.path.basename(os.path.dirname(prefix))   # e.g. A1C_Max_INT.MVP_AFR
    locus    = os.path.basename(prefix)                     # e.g. chr10_...1000G.AFR
    return f"{analysis}.{locus}"


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND DISTRIBUTE JOBS  ###############
###############################################################################

def main():

    # ── Parse arguments ──────────────────────────────────────────────────────────
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
            assert False, f"WARNING: Unhandled option: {o}"

    # ── Validate required arguments ───────────────────────────────────────────────
    if command_file is None:
        print("ERROR: Missing required --command_file")
        sys.exit(2)

    for path, label in [(command_file, "command file"), (log_dir, "log directory")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Read command file ─────────────────────────────────────────────────────────
    with open(command_file, 'r') as fh:
        all_commands = [line.rstrip('\n') for line in fh if line.strip()]

    # ── Skip commands whose .vcor1.vars sentinel already exists ──────────────────
    pending  = []
    skipped  = 0
    unparsed = 0
    for cmd in all_commands:
        sentinel = get_sentinel_path(cmd)
        if sentinel is None:
            print(f"WARNING: Could not parse --out from command: {cmd[:80]}")
            unparsed += 1
            pending.append(cmd)  # run it anyway; plink2 will create the output
        elif os.path.exists(sentinel):
            skipped += 1
        else:
            # Ensure the output directory exists before the worker tries to write
            os.makedirs(os.path.dirname(sentinel), exist_ok=True)
            pending.append(cmd)

    print(f"NOTE: {len(all_commands)} total commands  |  "
          f"{skipped} already complete  |  {len(pending)} to run")
    if unparsed:
        print(f"WARNING: {unparsed} command(s) could not be parsed for sentinel check")

    if not pending:
        print("NOTE: All jobs already complete.  Skipping parsl dispatch.")
        print("NOTE: matrix_1000g_parsl_allocator.py completed successfully.")
        sys.exit(0)

    # ── Define parsl bash_app ─────────────────────────────────────────────────────
    @bash_app
    def run_1000g_job(command, stdout=None, stderr=None):
        return command

    # ── Change to log dir so parsl worker logs are co-located with other logs ─────
    os.chdir(log_dir)
    print("NOTE: Starting plink2 LD matrix jobs via parsl")

    # ── Load parsl config and dispatch all pending jobs ───────────────────────────
    with parsl.load(matrix_1000g_config):
        futures = []
        for cmd in pending:
            label = get_task_label(cmd) or f"unknown_{len(futures)}"
            futures.append(run_1000g_job(
                cmd,
                stdout=os.path.join(log_dir, f"{label}.1000g_out.txt"),
                stderr=os.path.join(log_dir, f"{label}.1000g_err.txt"),
            ))

        # Collect results fault-tolerantly: a single failure must not abort the rest
        failed_cmds = []
        for cmd, future in zip(pending, futures):
            try:
                future.result()
            except Exception as e:
                label = get_task_label(cmd) or cmd[:60]
                print(f"ERROR: Task {label} failed: {e}", flush=True)
                failed_cmds.append(cmd)

    if failed_cmds:
        print(f"WARNING: {len(failed_cmds)} task(s) failed — "
              f"see .1000g_err.txt / .1000g_out.txt files in {log_dir}")
    else:
        print(f"SUCCESS: All {len(pending)} plink2 LD matrix jobs complete")

    print("NOTE: matrix_1000g_parsl_allocator.py completed successfully.")


if __name__ == "__main__":
    main()
