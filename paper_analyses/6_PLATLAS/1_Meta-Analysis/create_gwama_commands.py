'''
create_gwama_commands.py

The purpose of this script is to read the meta-analysis linker file and produce
a text file of run_single_gwama.py commands — one command per unique
trait/population combination.

Commands are sorted ascending by the total compressed size of their input
summary statistics files.  This ensures that the smallest jobs complete first
and free worker slots quickly, while the largest jobs (e.g. ALL-population
meta-analyses combining MVP/AFR + UKB/EUR files) run when the queue is not
fully saturated.

Usage:
  python create_gwama_commands.py [OPTIONS]

  --linker_file     Path to linker TSV (default: standard path; optional)
  --sumstats_base   Base directory for sumstats files (default: standard; optional)
  --output_file     Path to write commands file
                    (default: <script_dir>/gwama_meta_analysis.commands.txt; optional)
  --output_dir      Base output directory (default: standard path; optional)
  --temp_dir        Base temp directory (default: standard path; optional)
  -h, --help        Show this help message

***Run in Conda finemap_old environment***
'''

import sys
import getopt
import os
from collections import defaultdict

# ── Default paths ──────────────────────────────────────────────────────────────
SUMSTATS_BASE = "/flare/GeomicVar/multi-biobank-meta/sumstats"
OUTPUT_BASE   = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results"
LINKER_FILE   = "/flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv"
TEMP_BASE     = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/gwama_temp"

SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CMDS  = os.path.join(SCRIPT_DIR, "gwama_meta_analysis.commands.txt")
STANDALONE    = os.path.join(SCRIPT_DIR, "run_single_gwama.py")


def help(exit_num=0):
    print("""Usage: python create_gwama_commands.py [OPTIONS]
Options:
  --linker_file     Path to linker TSV (default: standard path; optional)
  --sumstats_base   Base directory for sumstats files (default: standard path; optional)
  --output_file     Path to write the commands file (default: <script_dir>/gwama_meta_analysis.commands.txt; optional)
  --output_dir      Base output directory (default: standard path; optional)
  --temp_dir        Base temp directory for decompressed files (default: standard path; optional)
  -h, --help        Show this help message""")
    sys.exit(exit_num)


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND BUILD COMMAND LIST  ############
###############################################################################

def main():

    # ── Parse arguments ────────────────────────────────────────────────────────
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "linker_file=", "sumstats_base=", "output_file=",
             "output_dir=", "temp_dir="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    linker_file   = LINKER_FILE
    sumstats_base = SUMSTATS_BASE
    output_file   = DEFAULT_CMDS
    output_dir    = OUTPUT_BASE
    temp_dir      = TEMP_BASE

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--linker_file":
            linker_file = a
        elif o == "--sumstats_base":
            sumstats_base = a
        elif o == "--output_file":
            output_file = a
        elif o == "--output_dir":
            output_dir = a
        elif o == "--temp_dir":
            temp_dir = a
        else:
            assert False, f"WARNING: Unhandled input option: {o}"

    # ── Path validation ────────────────────────────────────────────────────────
    for path, label in [
        (linker_file,   "linker file"),
        (sumstats_base, "sumstats base directory"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    if not os.path.exists(STANDALONE):
        print(f"ERROR: run_single_gwama.py not found at: {STANDALONE}")
        sys.exit(1)

    # ── Parse linker file: group relative paths by trait_pop ─────────────────
    trait_pop_files = defaultdict(list)
    with open(linker_file, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 2:
                continue
            trait_pop = cols[0]
            full_path = os.path.join(sumstats_base, cols[1])
            trait_pop_files[trait_pop].append(full_path)

    # ── Build sorted command list ─────────────────────────────────────────────
    # Each entry: (total_compressed_bytes, trait_pop, command_string)
    job_list      = []
    skipped_count = 0

    for trait_pop, gz_files in trait_pop_files.items():
        # Check all input files exist; skip the combo if any are missing
        missing = [f for f in gz_files if not os.path.exists(f)]
        if missing:
            print(f"WARNING: Skipping {trait_pop} — missing input file(s):")
            for m in missing:
                print(f"  {m}")
            skipped_count += 1
            continue

        # Sum compressed sizes (bytes) for size-aware scheduling
        total_size = sum(os.path.getsize(f) for f in gz_files)

        cmd = (
            f"python {STANDALONE} "
            f"--trait_pop {trait_pop} "
            f"--linker_file {linker_file} "
            f"--sumstats_base {sumstats_base} "
            f"--output_dir {output_dir} "
            f"--temp_dir {temp_dir}"
        )
        job_list.append((total_size, trait_pop, cmd))

    # Sort ascending by total input size (smallest jobs first)
    job_list.sort(key=lambda x: x[0])

    # ── Write commands file ───────────────────────────────────────────────────
    with open(output_file, 'w') as fh:
        for _, _, cmd in job_list:
            fh.write(cmd + '\n')

    # ── Summary statistics ────────────────────────────────────────────────────
    print(f"NOTE: {len(job_list)} commands written to {output_file}")
    if skipped_count:
        print(f"WARNING: {skipped_count} trait/population combination(s) skipped "
              f"due to missing input files")

    if job_list:
        sizes = [s for s, _, _ in job_list]
        total_gb  = sum(sizes) / 1e9
        min_mb    = min(sizes) / 1e6
        max_gb    = max(sizes) / 1e9
        n_files_total = sum(
            len(trait_pop_files[tp]) for _, tp, _ in job_list
        )
        print(f"NOTE: Total compressed input data across all jobs: {total_gb:.1f} GB "
              f"({n_files_total} files)")
        print(f"NOTE: Smallest job: {min_mb:.1f} MB  |  "
              f"Largest job: {max_gb:.2f} GB (compressed)")


if __name__ == "__main__":
    main()
