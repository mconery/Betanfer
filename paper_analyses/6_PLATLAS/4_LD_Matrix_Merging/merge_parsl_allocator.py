'''
merge_parsl_allocator.py

Manages the multi-node distribution of LD matrix merge jobs using parsl.

Each command in the supplied command file corresponds to one
(trait, analysis_pop, locus) merge_ld_matrices.py run.  The allocator:
  1. Reads the command file produced by create_merge_commands.py.
  2. Determines which jobs have already produced a .merged.vcor1 sentinel file
     and skips them, allowing safe resubmission of the PBS job without
     duplication.
  3. Submits remaining commands to parsl workers via @bash_app, with per-task
     stdout/stderr captured to named log files.
  4. Collects results fault-tolerantly: failed tasks are logged but do not abort
     the remaining work.

The parsl configuration (ld_merge_config) is imported from the local
parsl_config.py, which uses a LocalProvider over pre-allocated PBS nodes.

Usage:
  python merge_parsl_allocator.py --command_file <path> [--log_dir <path>]

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os
import parsl
from parsl import bash_app
from parsl_config import ld_merge_config

# Default path for parsl worker log output
DEFAULT_LOG_DIR = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/script_logs"


def help(exit_num=0):
    print("""Usage: python merge_parsl_allocator.py [OPTIONS]
Options:
  --command_file  File of merge_ld_matrices.py commands to execute (required)
  --log_dir       Directory for parsl worker log output
                  (default: standard script_logs path; optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
############  HELPER: DETERMINE EXPECTED SENTINEL PATH FROM A COMMAND  ########
###############################################################################

def get_sentinel_path(command):
    """Return the expected .merged.vcor1 sentinel path for a merge command.

    Parses --output-prefix {prefix} from the command string and appends
    '.merged.vcor1'.  Returns None if --output-prefix cannot be found.
    """
    parts = command.split()
    for i, part in enumerate(parts):
        if part == "--output-prefix" and i + 1 < len(parts):
            return parts[i + 1] + ".merged.vcor1"
    return None


def get_task_label(command):
    """Return a short label for per-task log filenames.

    Derives a label from the --output-prefix path:
      {trait}.{analysis_pop}.{locus_stem}
    e.g. Phe_008_52.ALL.chr1_100000_200000.merged
    """
    sentinel = get_sentinel_path(command)
    if sentinel is None:
        return None
    # sentinel = {LD_BASE}/{trait}.{analysis_pop}/chr...merged.vcor1
    prefix   = sentinel[:-len(".merged.vcor1")]          # strip .merged.vcor1
    trait_pop = os.path.basename(os.path.dirname(prefix)) # e.g. Phe_008_52.ALL
    locus_stem = os.path.basename(prefix)                  # e.g. chr1_100000_200000.merged
    return f"{trait_pop}.{locus_stem}"


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
            assert False, f"WARNING: Unhandled option: {o}"

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

    # ── Skip commands whose .merged.vcor1 sentinel already exists ─────────────
    pending  = []
    skipped  = 0
    unparsed = 0
    for cmd in all_commands:
        sentinel = get_sentinel_path(cmd)
        if sentinel is None:
            print(f"WARNING: Could not parse --output-prefix from command: {cmd[:80]}")
            unparsed += 1
            pending.append(cmd)  # run it anyway; merge script handles gracefully
        elif os.path.exists(sentinel):
            skipped += 1
        else:
            # Ensure output directory exists before the worker tries to write
            out_dir = os.path.dirname(sentinel)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            pending.append(cmd)

    print(f"NOTE: {len(all_commands)} total commands  |  "
          f"{skipped} already complete  |  {len(pending)} to run")
    if unparsed:
        print(f"WARNING: {unparsed} command(s) could not be parsed for sentinel check")

    if not pending:
        print("NOTE: All jobs already complete.  Skipping parsl dispatch.")
        print("NOTE: merge_parsl_allocator.py completed successfully.")
        sys.exit(0)

    # ── Define parsl bash_app ─────────────────────────────────────────────────
    @bash_app
    def run_merge_job(command, stdout=None, stderr=None):
        return command

    # ── Change to log dir so parsl worker logs are co-located with other logs ──
    os.chdir(log_dir)
    print("NOTE: Starting LD matrix merge jobs via parsl")

    # ── Load parsl config and dispatch all pending jobs ───────────────────────
    with parsl.load(ld_merge_config):
        futures = []
        for cmd in pending:
            label = get_task_label(cmd) or f"unknown_{len(futures)}"
            futures.append(run_merge_job(
                cmd,
                stdout=os.path.join(log_dir, f"{label}.merge_out.txt"),
                stderr=os.path.join(log_dir, f"{label}.merge_err.txt"),
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
              f"see .merge_err.txt / .merge_out.txt files in {log_dir}")
    else:
        print(f"SUCCESS: All {len(pending)} merge jobs complete")

    print("NOTE: merge_parsl_allocator.py completed successfully.")


if __name__ == "__main__":
    main()
