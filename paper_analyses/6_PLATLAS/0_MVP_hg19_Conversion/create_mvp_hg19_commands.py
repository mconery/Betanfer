'''
create_mvp_hg19_commands.py

The purpose of this script is to scan the MVP summary statistics directory and
produce a text file of convert_mvp_to_hg19.py commands — one command per
primary .sumstats.txt.gz file across the AFR and EUR population subdirectories.

Non-primary files (*.hm3_sumstats.txt.gz, *.adj.lead.tsv) are excluded.
Files whose output .sumstats.hg19.txt already exists are silently skipped.

Commands are sorted ascending by compressed input file size so the smallest
jobs complete first and free worker slots quickly.

Usage:
  python create_mvp_hg19_commands.py [OPTIONS]

  --sumstats_dir    Base directory containing AFP/, EUR/ subdirs of MVP sumstats
                    (default: standard MVP path; optional)
  --output_dir      Base directory for MVP_hg19 output (default: standard path; optional)
  --annotation_file Path to annotation .txt.gz file (default: standard path; optional)
  --output_file     Path to write commands file
                    (default: <script_dir>/mvp_hg19.commands.txt; optional)
  -h, --help        Show this help message

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os

# ── Default paths ──────────────────────────────────────────────────────────────
SUMSTATS_BASE   = "/flare/GeomicVar/multi-biobank-meta/sumstats/MVP"
OUTPUT_BASE     = "/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19"
ANNOTATION_FILE = "/flare/GeomicVar/multi-biobank-meta/Annotation_gwPheWAS_dbGaP_b38.noAF_R2.txt.gz"
POPS            = ["AFR", "EUR"]

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
CONVERTER    = os.path.join(SCRIPT_DIR, "convert_mvp_to_hg19.py")
DEFAULT_CMDS = os.path.join(SCRIPT_DIR, "mvp_hg19.commands.txt")

# Suffixes that identify non-primary MVP files — excluded from processing
EXCLUDE_SUFFIXES = (".hm3_sumstats.txt.gz", ".adj.lead.tsv")


def help(exit_num=0):
    print("""Usage: python create_mvp_hg19_commands.py [OPTIONS]
Options:
  --sumstats_dir    Base directory containing AFR/, EUR/ subdirs of MVP sumstats
                    (default: standard MVP path; optional)
  --output_dir      Base directory for MVP_hg19 output
                    (default: standard MVP_hg19 path; optional)
  --annotation_file Path to annotation .txt.gz file
                    (default: standard annotation path; optional)
  --output_file     Path to write the commands file
                    (default: <script_dir>/mvp_hg19.commands.txt; optional)
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
            ["help", "sumstats_dir=", "output_dir=", "annotation_file=", "output_file="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    sumstats_dir    = SUMSTATS_BASE
    output_dir      = OUTPUT_BASE
    annotation_file = ANNOTATION_FILE
    output_file     = DEFAULT_CMDS

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--sumstats_dir":
            sumstats_dir = a
        elif o == "--output_dir":
            output_dir = a
        elif o == "--annotation_file":
            annotation_file = a
        elif o == "--output_file":
            output_file = a
        else:
            assert False, f"WARNING: Unhandled input option: {o}"

    # ── Path validation ────────────────────────────────────────────────────────
    for path, label in [(sumstats_dir, "sumstats directory"),
                        (annotation_file, "annotation file")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Scan AFR and EUR subdirectories ────────────────────────────────────────
    # Each entry: (file_size_bytes, input_abs_path, output_abs_path)
    job_list      = []
    skipped_count = 0
    missing_pop   = []

    for pop in POPS:
        pop_in_dir  = os.path.join(sumstats_dir, pop)
        pop_out_dir = os.path.join(output_dir, pop)

        if not os.path.isdir(pop_in_dir):
            print(f"WARNING: Population subdirectory not found: {pop_in_dir}")
            missing_pop.append(pop)
            continue

        for fname in os.listdir(pop_in_dir):
            # Must end in .sumstats.txt.gz
            if not fname.endswith(".sumstats.txt.gz"):
                continue
            # Exclude non-primary files
            if any(fname.endswith(suf) for suf in EXCLUDE_SUFFIXES):
                continue

            input_path = os.path.join(pop_in_dir, fname)

            # Derive output filename: replace .sumstats.txt.gz → .sumstats.hg19.txt
            base      = fname[:-len(".sumstats.txt.gz")]
            out_fname = base + ".sumstats.hg19.txt"
            out_path  = os.path.join(pop_out_dir, out_fname)

            # Skip if output already exists
            if os.path.exists(out_path):
                skipped_count += 1
                continue

            file_size = os.path.getsize(input_path)
            job_list.append((file_size, input_path, out_path))

    # Sort ascending by input file size (smallest jobs first)
    job_list.sort(key=lambda x: x[0])

    # ── Build command strings ──────────────────────────────────────────────────
    commands = []
    for _, input_path, out_path in job_list:
        cmd = (
            f"python {CONVERTER}"
            f" --input_file {input_path}"
            f" --output_file {out_path}"
            f" --annotation_file {annotation_file}"
        )
        commands.append(cmd)

    # ── Write commands file ────────────────────────────────────────────────────
    with open(output_file, 'w') as fh:
        for cmd in commands:
            fh.write(cmd + '\n')

    # ── Summary statistics ─────────────────────────────────────────────────────
    print(f"NOTE: {len(commands)} commands written to {output_file}")
    print(f"NOTE: {skipped_count} file(s) already converted (skipped)")
    if missing_pop:
        print(f"WARNING: Missing population directories: {missing_pop}")

    if job_list:
        sizes    = [s for s, _, _ in job_list]
        total_gb = sum(sizes) / 1e9
        min_gb   = min(sizes) / 1e9
        max_gb   = max(sizes) / 1e9
        print(f"NOTE: Total compressed input data to process: {total_gb:.1f} GB")
        print(f"NOTE: Smallest file: {min_gb:.2f} GB  |  "
              f"Largest file: {max_gb:.2f} GB")


if __name__ == "__main__":
    main()
