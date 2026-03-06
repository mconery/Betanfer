'''
create_matrix_commands.py

The purpose of this script is to build the per-cohort-population file lists and
generate betanfer infer commands for every significant locus in each of the five
population master CSVs produced by the locus-defining pipeline.

Each command runs betanfer infer for one (locus, cohort+population) pair:
  - File list:     all sumstats files from that cohort+population that participated
                   in the analysis (from the linker for EUR/AFR/ALL; from the MVP
                   sumstats directory for MVP_EUR/MVP_AFR)
  - Region:        the locus CHR/START/END
  - P-val filter:  --pval-threshold 5e-8 (exclude traits with no locus signal)
  - Output:        {LD_BASE}/{TRAIT}.{ANALYSIS_POP}/chr{CHR}_{START}_{END}.{COHORT}.{COHORT_POP}

File lists are written once to {LD_BASE}/file_lists/ and reused across all loci.
Commands are sorted ascending by locus size so small jobs finish first.

Usage:
  python create_matrix_commands.py --command_file <path> [--log_dir <path>]

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os
import glob
import csv
import multiprocessing
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

# ── Default paths ───────────────────────────────────────────────────────────────
LINKER_FILE   = "/flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv"
LOCI_DIR      = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci"
LD_BASE       = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices"
SUMSTATS_BASE = "/flare/GeomicVar/multi-biobank-meta/sumstats"
MVP_HG19_SUMSTATS = os.path.join(SUMSTATS_BASE, "MVP_hg19")

SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CMDS  = os.path.join(SCRIPT_DIR, "matrix_making.commands.txt")

# Populations served by the linker file (meta-analysis populations)
LINKER_POPS = {"EUR", "AFR", "ALL"}
# Populations served by the MVP single-cohort directory scan
MVP_POPS = {"MVP_EUR": "EUR", "MVP_AFR": "AFR"}  # analysis_pop → mvp_dir_pop

# betanfer infer options applied to every command
PVAL_THRESHOLD = "5e-8"
NUM_PROCESSES  = 8   # matches physical cores per GPU-tile worker (16 logical / 2)

# Parallel workers for filesystem existence checks (Lustre handles many concurrent stats well)
NUM_WORKERS = multiprocessing.cpu_count()


def help(exit_num=0):
    print("""Usage: python create_matrix_commands.py [OPTIONS]
Options:
  --command_file  Path to write the betanfer infer commands file (required)
  --log_dir       Directory for parsl worker log output (optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
##############  BUILD PER-(ANALYSIS_POP, COHORT, COHORT_POP) FILE LISTS  ######
###############################################################################

def build_file_lists(linker_file, mvp_sumstats, ld_base):
    """Read linker and scan MVP directories to produce per-group file lists.

    Returns a dict keyed by (analysis_pop, cohort, cohort_pop) → sorted list of
    absolute sumstats file paths.  Also creates the file_lists/ subdirectory and
    writes one .txt file per group (overwriting only if the set of files changed).
    """
    groups = defaultdict(set)   # (analysis_pop, cohort, cohort_pop) → set of paths

    # ── Linker-based groups (EUR, AFR, ALL meta-analyses) ──────────────────────
    # Phase 1: parse all entries without any filesystem calls
    raw_entries = []  # (analysis_pop, cohort, cohort_pop, full_path)
    with open(linker_file, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                print(f"WARNING: Unexpected linker line format, skipping: {line[:80]}")
                continue
            trait_pop, rel_path = parts
            # Extract analysis population: last '.' token of trait_pop
            dot_parts = trait_pop.rsplit('.', 1)
            if len(dot_parts) != 2:
                print(f"WARNING: Cannot parse analysis pop from '{trait_pop}', skipping")
                continue
            analysis_pop = dot_parts[1]
            if analysis_pop not in LINKER_POPS:
                print(f"WARNING: Unexpected population '{analysis_pop}' in linker, skipping")
                continue
            # Extract cohort and cohort_pop from relative path: COHORT/COHORT_POP/filename
            path_parts = rel_path.split('/')
            if len(path_parts) < 3:
                print(f"WARNING: Cannot parse cohort/pop from path '{rel_path}', skipping")
                continue
            cohort     = path_parts[0]
            cohort_pop = path_parts[1]
            full_path  = os.path.join(SUMSTATS_BASE, rel_path)
            raw_entries.append((analysis_pop, cohort, cohort_pop, full_path))

    # Phase 2: check all unique paths in parallel, then filter
    unique_paths = list({e[3] for e in raw_entries})
    print(f"NOTE: Checking {len(unique_paths)} linker paths in parallel ({NUM_WORKERS} workers)...")
    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as ex:
        path_exists = dict(zip(unique_paths, ex.map(os.path.exists, unique_paths)))

    for analysis_pop, cohort, cohort_pop, full_path in raw_entries:
        if path_exists[full_path]:
            groups[(analysis_pop, cohort, cohort_pop)].add(full_path)
        else:
            print(f"WARNING: Linker file not found, excluding: {full_path}")

    # ── MVP directory scans for MVP_EUR and MVP_AFR ────────────────────────────
    for analysis_pop, dir_pop in MVP_POPS.items():
        mvp_dir = os.path.join(mvp_sumstats, dir_pop)
        if not os.path.isdir(mvp_dir):
            print(f"WARNING: MVP directory not found, skipping: {mvp_dir}")
            continue
        all_gz = glob.glob(os.path.join(mvp_dir, "MVP_R4.1000G_AGR.*.sumstats.hg19.txt.gz"))
        gz_files = [f for f in all_gz if not f.endswith("hm3_sumstats.txt.gz")]
        if not gz_files:
            print(f"WARNING: No MVP sumstats files found in {mvp_dir}")
            continue
        for f in gz_files:
            groups[(analysis_pop, "MVP", dir_pop)].add(f)
        print(f"NOTE: {analysis_pop}: found {len(gz_files)} MVP/{dir_pop} files")

    # ── Write file list .txt files ─────────────────────────────────────────────
    file_lists_dir = os.path.join(ld_base, "file_lists")
    os.makedirs(file_lists_dir, exist_ok=True)

    written = 0
    for (analysis_pop, cohort, cohort_pop), paths in sorted(groups.items()):
        fname = f"{analysis_pop}.{cohort}.{cohort_pop}.txt"
        fpath = os.path.join(file_lists_dir, fname)
        sorted_paths = sorted(paths)
        with open(fpath, 'w') as fh:
            fh.write('\n'.join(sorted_paths) + '\n')
        written += 1

    print(f"NOTE: Wrote {written} file list(s) to {file_lists_dir}")

    # Return sorted lists for command generation
    return {k: sorted(v) for k, v in groups.items()}


###############################################################################
####################  GENERATE BETANFER INFER COMMANDS  ######################
###############################################################################

def generate_commands(groups, loci_dir, ld_base):
    """Scan the five population master CSVs and produce betanfer infer commands.

    Each command targets one (locus, cohort, cohort_pop) combination.  Jobs whose
    sentinel .vcor1 file already exists are skipped (checkpoint).

    Returns a list of (locus_size_bp, command_str) tuples.
    """
    file_lists_dir = os.path.join(ld_base, "file_lists")
    warnings       = 0

    # Master CSV files produced by locus_parsl_allocator.py
    master_files = sorted(glob.glob(os.path.join(loci_dir, "sig_loci.*.master.csv")))
    if not master_files:
        print(f"ERROR: No sig_loci.*.master.csv files found in {loci_dir}")
        sys.exit(1)

    # Phase 1: build all (locus_size, sentinel, cmd) work items without any
    # filesystem calls — just string construction from CSV data.
    all_work = []   # [(locus_size, sentinel_path, cmd_str), ...]

    for master_csv in master_files:
        # Derive analysis_pop from filename: sig_loci.{POP}.master.csv
        basename     = os.path.basename(master_csv)
        analysis_pop = basename.split('.')[1]

        # Find constituent (cohort, cohort_pop) pairs for this analysis population
        pairs = [(cohort, cohort_pop)
                 for (ap, cohort, cohort_pop) in groups
                 if ap == analysis_pop]
        if not pairs:
            print(f"WARNING: No file lists found for analysis population '{analysis_pop}', "
                  f"skipping {basename}")
            warnings += 1
            continue

        # Read the master CSV
        with open(master_csv, 'r') as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)

        loci_count = 0
        for row in rows:
            chrom = row['CHR']
            start = row['START']
            end   = row['END']
            trait = row['TRAIT']
            locus_size = int(end) - int(start)

            for cohort, cohort_pop in pairs:
                file_list_name = f"{analysis_pop}.{cohort}.{cohort_pop}.txt"
                file_list_path = os.path.join(file_lists_dir, file_list_name)
                output_prefix  = os.path.join(
                    ld_base,
                    f"{trait}.{analysis_pop}",
                    f"chr{chrom}_{start}_{end}.{cohort}.{cohort_pop}"
                )
                sentinel = output_prefix + ".vcor1"
                cmd = (
                    f"python -m betanfer infer"
                    f" --file-list {file_list_path}"
                    f" --chrom {chrom} --start {start} --end {end}"
                    f" --pval-threshold {PVAL_THRESHOLD}"
                    f" --num-processes {NUM_PROCESSES}"
                    f" --output {output_prefix}"
                    f" --plink2-format"
                )
                all_work.append((locus_size, sentinel, cmd))
            loci_count += 1

        print(f"NOTE: {analysis_pop}: processed {loci_count} loci × {len(pairs)} cohort(s)")

    # Phase 2: check all sentinel files in parallel, then filter
    total = len(all_work)
    print(f"NOTE: Checking {total} sentinel files in parallel ({NUM_WORKERS} workers)...")

    def _check(item):
        locus_size, sentinel, cmd = item
        return None if os.path.exists(sentinel) else (locus_size, cmd)

    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as ex:
        results = list(ex.map(_check, all_work))

    job_list = [r for r in results if r is not None]
    skipped  = total - len(job_list)

    print(f"NOTE: {total} total (locus, cohort+pop) pairs  |  "
          f"{skipped} already complete  |  {len(job_list)} to run  |  "
          f"{warnings} population(s) with no file list")

    return job_list


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND RUN  ###########################
###############################################################################

def main():

    # ── Parse arguments ─────────────────────────────────────────────────────────
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
    log_dir      = None   # unused in this script; accepted for CLI consistency

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--command_file":
            command_file = a
        elif o == "--log_dir":
            log_dir = a
        else:
            assert False, f"WARNING: Unhandled option: {o}"

    if command_file is None:
        print("ERROR: Missing required --command_file")
        sys.exit(2)

    # ── Validate required input paths ───────────────────────────────────────────
    for path, label in [
        (LINKER_FILE, "linker file"),
        (LOCI_DIR,    "loci directory"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # Ensure LD base directory exists
    os.makedirs(LD_BASE, exist_ok=True)

    # ── Build file lists ────────────────────────────────────────────────────────
    print("NOTE: Building constituent file lists...")
    groups = build_file_lists(LINKER_FILE, MVP_HG19_SUMSTATS, LD_BASE)

    # ── Generate commands ────────────────────────────────────────────────────────
    print("NOTE: Generating betanfer infer commands...")
    job_list = generate_commands(groups, LOCI_DIR, LD_BASE)

    # Sort ascending by locus size (small loci first)
    job_list.sort(key=lambda x: x[0])

    # ── Write commands file ──────────────────────────────────────────────────────
    cmd_dir = os.path.dirname(command_file)
    if cmd_dir:
        os.makedirs(cmd_dir, exist_ok=True)
    with open(command_file, 'w') as fh:
        for _, cmd in job_list:
            fh.write(cmd + '\n')

    total_size_gb = sum(s for s, _ in job_list) / 1e9
    print(f"NOTE: {len(job_list)} commands written to {command_file}")
    if job_list:
        min_kb = min(s for s, _ in job_list) / 1e3
        max_mb = max(s for s, _ in job_list) / 1e6
        print(f"NOTE: Locus sizes: smallest {min_kb:.1f} kb  |  largest {max_mb:.1f} Mb  |  "
              f"total {total_size_gb:.1f} Gb across all loci")


if __name__ == "__main__":
    main()
