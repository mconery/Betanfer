'''
create_locus_commands.py

The purpose of this script is to scan input directories for GWAS summary
statistics files and produce a text file of define_sig_loci_for_trait.py
commands — one command per trait/population combination.

Covers three meta-analysis population inputs (GWAMA output):
  EUR, AFR, ALL

MVP populations are handled separately by create_mvp_loci_from_s11.py,
which writes master locus files directly from Supplementary Table S11.

Commands are sorted ascending by compressed input file size (smallest first)
so that the smallest jobs complete early and free parsl worker slots, while
the largest jobs run when the queue is not fully saturated.

Usage:
  python create_locus_commands.py [OPTIONS]

  --meta_results_dir  Base directory of GWAMA meta-analysis outputs
                      (default: standard meta-analysis_results path; optional)
  --output_dir        Directory where per-trait loci CSVs will be written
                      (default: standard mapping_loci path; optional)
  --output_file       Path to write the commands file
                      (default: <script_dir>/locus_defining.commands.txt; optional)
  -h, --help          Show this help message

***Run in Conda finemap_old environment***
'''

import sys
import getopt
import os
import glob

# ── Default paths ──────────────────────────────────────────────────────────────
META_RESULTS_BASE = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results"
OUTPUT_BASE       = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci"

SCRIPT_DIR        = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CMDS      = os.path.join(SCRIPT_DIR, "locus_defining.commands.txt")
TRAIT_SCRIPT      = os.path.join(SCRIPT_DIR, "define_sig_loci_for_trait.py")

META_POPS = ["EUR", "AFR", "ALL"]


def help(exit_num=0):
    print("""Usage: python create_locus_commands.py [OPTIONS]
Options:
  --meta_results_dir  Base directory of GWAMA meta-analysis outputs
                      (default: standard path; optional)
  --output_dir        Directory for per-trait loci CSVs
                      (default: standard mapping_loci path; optional)
  --output_file       Path to write the commands file
                      (default: <script_dir>/locus_defining.commands.txt; optional)
  -h, --help          Show this help message""")
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
            ["help", "meta_results_dir=", "output_dir=", "output_file="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    meta_results_dir = META_RESULTS_BASE
    output_dir       = OUTPUT_BASE
    output_file      = DEFAULT_CMDS

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--meta_results_dir":
            meta_results_dir = a.rstrip('/')
        elif o == "--output_dir":
            output_dir = a.rstrip('/')
        elif o == "--output_file":
            output_file = a
        else:
            assert False, f"WARNING: Unhandled input option: {o}"

    # ── Path validation ────────────────────────────────────────────────────────
    for path, label in [
        (meta_results_dir, "meta-analysis results directory"),
        (TRAIT_SCRIPT,     "define_sig_loci_for_trait.py"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Build command list ─────────────────────────────────────────────────────
    # Each entry: (compressed_bytes, base_name, command_string)
    job_list = []

    # ── Meta-analysis populations (GWAMA output format) ───────────────────────
    # Files: Phe_*.out.gz
    # p-value column: "p-value"  (GWAMA header)
    # Index column: #rs_number (rsID) — NOT CHR:POS:REF:ALT
    for pop in META_POPS:
        pop_dir = os.path.join(meta_results_dir, pop)
        if not os.path.isdir(pop_dir):
            print(f"WARNING: Meta-analysis directory not found, skipping: {pop_dir}")
            continue

        gz_files = sorted(glob.glob(os.path.join(pop_dir, "Phe_*.out.gz")))
        if not gz_files:
            print(f"WARNING: No Phe_*.out.gz files found in {pop_dir}")
            continue

        for f in gz_files:
            trait = os.path.basename(f).replace('.out.gz', '')
            base_name = f"{trait}.{pop}"
            cmd = (
                f"python {TRAIT_SCRIPT} "
                f"-g {f} "
                f"-o {output_dir} "
                f"-p 5e-8 -s 1e-6 -t 100000 "
                f"-b {base_name} "
                f"--chr-col CHR --pos-col POS --p-col p-value"
            )
            job_list.append((os.path.getsize(f), base_name, cmd))

        print(f"NOTE: Meta-analysis {pop}: found {len(gz_files)} trait files")

    # ── Sort ascending by compressed file size ────────────────────────────────
    job_list.sort(key=lambda x: x[0])

    # ── Write commands file ───────────────────────────────────────────────────
    os.makedirs(os.path.dirname(output_file) if os.path.dirname(output_file) else '.', exist_ok=True)
    with open(output_file, 'w') as fh:
        for _, _, cmd in job_list:
            fh.write(cmd + '\n')

    # ── Summary statistics ────────────────────────────────────────────────────
    print(f"NOTE: {len(job_list)} total commands written to {output_file}")

    if job_list:
        sizes     = [s for s, _, _ in job_list]
        total_gb  = sum(sizes) / 1e9
        min_mb    = min(sizes) / 1e6
        max_gb    = max(sizes) / 1e9
        print(f"NOTE: Total compressed input: {total_gb:.1f} GB")
        print(f"NOTE: Smallest job: {min_mb:.1f} MB  |  Largest job: {max_gb:.2f} GB (compressed)")


if __name__ == "__main__":
    main()
