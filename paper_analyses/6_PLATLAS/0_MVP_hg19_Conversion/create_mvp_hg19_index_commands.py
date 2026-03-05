'''
create_mvp_hg19_index_commands.py

The purpose of this script is to scan the MVP_hg19 output directory and
produce a text file of betanfer index commands — one command per
.sumstats.hg19.txt file across the AFR and EUR population subdirectories.

Files whose .gz.tbi index already exists are silently skipped.

Commands are sorted ascending by file size so the smallest jobs complete
first and free worker slots quickly.

Usage:
  python create_mvp_hg19_index_commands.py [OPTIONS]

  --hg19_dir     Base directory containing AFR/, EUR/ subdirs of hg19 sumstats
                 (default: standard MVP_hg19 path; optional)
  --output_file  Path to write commands file
                 (default: <script_dir>/mvp_hg19_index.commands.txt; optional)
  -h, --help     Show this help message

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os

# ── Default paths ──────────────────────────────────────────────────────────────
HG19_BASE  = "/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19"
POPS       = ["AFR", "EUR"]

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CMDS = os.path.join(SCRIPT_DIR, "mvp_hg19_index.commands.txt")


def help(exit_num=0):
    print("""Usage: python create_mvp_hg19_index_commands.py [OPTIONS]
Options:
  --hg19_dir     Base directory containing AFR/, EUR/ subdirs of hg19 sumstats
                 (default: standard MVP_hg19 path; optional)
  --output_file  Path to write the commands file
                 (default: <script_dir>/mvp_hg19_index.commands.txt; optional)
  -h, --help     Show this help message""")
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
            ["help", "hg19_dir=", "output_file="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    hg19_dir    = HG19_BASE
    output_file = DEFAULT_CMDS

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--hg19_dir":
            hg19_dir = a
        elif o == "--output_file":
            output_file = a
        else:
            assert False, f"WARNING: Unhandled input option: {o}"

    # ── Path validation ────────────────────────────────────────────────────────
    if not os.path.exists(hg19_dir):
        print(f"ERROR: hg19 directory not found: {hg19_dir}")
        sys.exit(1)

    # ── Scan AFR and EUR subdirectories ────────────────────────────────────────
    # Each entry: (file_size_bytes, abs_path)
    job_list      = []
    skipped_count = 0
    missing_pop   = []

    for pop in POPS:
        pop_dir = os.path.join(hg19_dir, pop)
        if not os.path.isdir(pop_dir):
            print(f"WARNING: Population subdirectory not found: {pop_dir}")
            missing_pop.append(pop)
            continue

        for fname in os.listdir(pop_dir):
            if not fname.endswith(".sumstats.hg19.txt"):
                continue

            abs_path = os.path.join(pop_dir, fname)

            # Skip if .gz.tbi already exists (job previously completed)
            if os.path.exists(abs_path + ".gz.tbi"):
                skipped_count += 1
                continue

            file_size = os.path.getsize(abs_path)
            job_list.append((file_size, abs_path))

    # Sort ascending by file size (smallest jobs first)
    job_list.sort(key=lambda x: x[0])

    # ── Build command strings ──────────────────────────────────────────────────
    # betanfer index with explicit chrom/pos columns since hg19 files have
    # standard CHR and POS columns in hg19 space
    commands = []
    for _, abs_path in job_list:
        cmd = (
            f"python -m betanfer index"
            f" --file {abs_path}"
            f" --chrom-col CHR"
            f" --pos-col POS"
            f" --force"
        )
        commands.append(cmd)

    # ── Write commands file ────────────────────────────────────────────────────
    with open(output_file, 'w') as fh:
        for cmd in commands:
            fh.write(cmd + '\n')

    # ── Summary statistics ─────────────────────────────────────────────────────
    print(f"NOTE: {len(commands)} commands written to {output_file}")
    print(f"NOTE: {skipped_count} file(s) already indexed (skipped)")
    if missing_pop:
        print(f"WARNING: Missing population directories: {missing_pop}")

    if job_list:
        sizes    = [s for s, _ in job_list]
        total_gb = sum(sizes) / 1e9
        min_gb   = min(sizes) / 1e9
        max_gb   = max(sizes) / 1e9
        print(f"NOTE: Total uncompressed input data to index: {total_gb:.1f} GB")
        print(f"NOTE: Smallest file: {min_gb:.2f} GB  |  "
              f"Largest file: {max_gb:.2f} GB")


if __name__ == "__main__":
    main()
