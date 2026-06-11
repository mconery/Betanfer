'''
create_merge_commands.py

Generates per-(trait, analysis_pop, locus) merge commands for the LD matrix
merging pipeline.

For each locus in each population's master CSV, the script:
  1. Builds two cohort-pair sets from the linker:
       - all_pairs[analysis_pop]: union across all traits — matches what
         create_matrix_commands.py created, used for the completeness check.
       - trait_pairs[(trait, analysis_pop)]: only cohorts that contributed to
         that specific trait's meta-analysis — used for the merge command so
         that cohorts without sumstats for a given trait are excluded.
  2. Checks that ALL all_pairs constituent .vcor1 files exist for this
     (trait, locus) combination.  If any are missing, the locus is flagged as
     "incomplete" and skipped — the merge should not proceed until all
     matrix-making jobs have finished.
  3. Skips loci whose merged output (.merged.vcor1) already exists (checkpoint).
  4. From the trait-specific pairs, filters out empty constituent files
     (0-byte .vcor1 = betanfer found no variants for that cohort) and looks up
     sample sizes for the remainder.
  5. Generates a merge_ld_matrices.py command for every locus with at least
     one valid non-empty input.

Commands are written to --command_file, sorted ascending by total input bytes
(smallest merge jobs first).

Usage:
  python create_merge_commands.py --command_file <path>

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os
import csv
import glob
import multiprocessing
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

# ── Default paths ─────────────────────────────────────────────────────────────
LINKER_FILE      = "/flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv"
SUMSTATS_BASE    = "/flare/GeomicVar/multi-biobank-meta/sumstats"
LD_BASE          = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices"
LOCI_DIR         = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci"
SAMPLE_SIZES_FILE = os.path.join(LD_BASE, "sample_sizes.tsv")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Populations whose meta-analysis matrices need merging
MERGE_POPS = {"EUR", "AFR", "ALL"}

# Parallel workers for filesystem checks
NUM_WORKERS = multiprocessing.cpu_count()


def help(exit_num=0):
    print("""Usage: python create_merge_commands.py [OPTIONS]
Options:
  --command_file  Path to write the merge commands file (required)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
################  LINKER PARSING  #############################################
###############################################################################

def parse_linker(linker_file):
    """Parse linker and return two structures:

    all_pairs[analysis_pop] = union of (cohort, cohort_pop) across all traits
        for that analysis_pop.  Matches the scope of create_matrix_commands.py
        and is used for the per-locus completeness check.

    trait_pairs[(trait, analysis_pop)] = (cohort, cohort_pop) pairs that
        contributed to that specific trait's meta-analysis.  Used to build
        merge commands so that cohorts absent from a trait's linker entry are
        not included even if their .vcor1 file happens to be non-empty.

    Returns (all_pairs, trait_pairs).
    """
    all_pairs   = defaultdict(set)   # analysis_pop -> {(cohort, cohort_pop), ...}
    trait_pairs = defaultdict(set)   # (trait, analysis_pop) -> {(cohort, cohort_pop), ...}

    with open(linker_file, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                continue
            trait_pop, rel_path = parts
            dot_parts = trait_pop.rsplit('.', 1)
            if len(dot_parts) != 2:
                continue
            trait, analysis_pop = dot_parts
            if analysis_pop not in MERGE_POPS:
                continue
            path_parts = rel_path.split('/')
            if len(path_parts) < 3:
                continue
            cohort     = path_parts[0]
            cohort_pop = path_parts[1]
            all_pairs[analysis_pop].add((cohort, cohort_pop))
            trait_pairs[(trait, analysis_pop)].add((cohort, cohort_pop))

    for pop, pairs in sorted(all_pairs.items()):
        print(f"NOTE: {pop}: {len(pairs)} expected (cohort, cohort_pop) pairs: "
              f"{sorted(pairs)}")

    return all_pairs, trait_pairs


###############################################################################
################  SAMPLE SIZES  ###############################################
###############################################################################

def load_sample_sizes(sample_sizes_file):
    """Load sample_sizes.tsv into a dict: (trait, cohort, cohort_pop) -> N."""
    sample_sizes = {}
    if not os.path.exists(sample_sizes_file):
        print(f"ERROR: Sample sizes file not found: {sample_sizes_file}")
        print("ERROR: Run extract_sample_sizes.py first.")
        sys.exit(1)
    with open(sample_sizes_file, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            key = (row['TRAIT'], row['COHORT'], row['COHORT_POP'])
            sample_sizes[key] = int(row['N'])
    print(f"NOTE: Loaded {len(sample_sizes)} sample size entries")
    return sample_sizes


###############################################################################
################  COMMAND GENERATION  #########################################
###############################################################################

def generate_commands(all_pairs, trait_pairs, sample_sizes, loci_dir, ld_base):
    """Scan master loci CSVs and produce merge commands.

    all_pairs[(analysis_pop)]        – union cohort set; used for completeness check.
    trait_pairs[(trait,analysis_pop)] – per-trait cohort set; used for merge command.

    Returns a list of (sort_key, command_str) tuples.
    """
    merge_script = os.path.join(SCRIPT_DIR, "merge_ld_matrices.py")

    all_work = []

    master_files = sorted(glob.glob(os.path.join(loci_dir, "sig_loci.*.master.csv")))
    if not master_files:
        print(f"ERROR: No sig_loci.*.master.csv files found in {loci_dir}")
        sys.exit(1)

    for master_csv in master_files:
        basename     = os.path.basename(master_csv)
        analysis_pop = basename.split('.')[1]

        if analysis_pop not in MERGE_POPS:
            continue

        completeness_pairs = sorted(all_pairs.get(analysis_pop, set()))
        if not completeness_pairs:
            print(f"WARNING: No expected cohort pairs for '{analysis_pop}', skipping {basename}")
            continue

        with open(master_csv, 'r') as fh:
            reader = csv.DictReader(fh)
            rows   = list(reader)

        loci_count = 0
        for row in rows:
            chrom = row['CHR']
            start = row['START']
            end   = row['END']
            trait = row['TRAIT']

            # All files the matrix-making step should have created (completeness check)
            completeness_files = [
                os.path.join(
                    ld_base,
                    f"{trait}.{analysis_pop}",
                    f"chr{chrom}_{start}_{end}.{cohort}.{cohort_pop}.vcor1"
                )
                for cohort, cohort_pop in completeness_pairs
            ]

            # Only cohorts that contributed to this trait's meta-analysis
            merge_pairs = sorted(trait_pairs.get((trait, analysis_pop), set()))
            merge_files = [
                os.path.join(
                    ld_base,
                    f"{trait}.{analysis_pop}",
                    f"chr{chrom}_{start}_{end}.{cohort}.{cohort_pop}.vcor1"
                )
                for cohort, cohort_pop in merge_pairs
            ]

            merged_output = os.path.join(
                ld_base,
                f"{trait}.{analysis_pop}",
                f"chr{chrom}_{start}_{end}.merged"
            )
            locus_size = int(end) - int(start)
            all_work.append((analysis_pop, trait, chrom, start, end,
                             completeness_pairs, completeness_files,
                             merge_pairs, merge_files,
                             merged_output, locus_size))
            loci_count += 1

        print(f"NOTE: {analysis_pop}: {loci_count} loci × {len(completeness_pairs)} cohort(s) = "
              f"{loci_count * len(completeness_pairs)} constituent files to check")

    if not all_work:
        print("NOTE: No loci found to process.")
        return []

    # ── Parallel filesystem checks ────────────────────────────────────────────
    # Check completeness files and merged sentinels in one pass.
    # merge_files are a subset of completeness_files so no need to add separately.
    total_files_to_check = sum(len(item[6]) for item in all_work)
    print(f"NOTE: Checking {total_files_to_check} constituent files + "
          f"{len(all_work)} merged sentinels in parallel ({NUM_WORKERS} workers)...")

    unique_paths = set()
    for item in all_work:
        unique_paths.update(item[6])                # completeness_files
        unique_paths.add(item[9] + ".merged.vcor1")       # merged sentinel
        unique_paths.add(item[9] + ".merged.vcor1.vars")  # vars sentinel

    def _check_exists(path):
        return path, os.path.exists(path)

    def _get_size(path):
        try:
            return path, os.path.getsize(path)
        except OSError:
            return path, 0

    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as ex:
        exists_map = dict(ex.map(_check_exists, unique_paths))

    existing_paths = [p for p in unique_paths if exists_map.get(p, False)]
    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as ex:
        size_map = dict(ex.map(_get_size, existing_paths))

    # ── Build command list ────────────────────────────────────────────────────
    job_list     = []
    n_complete   = 0
    n_incomplete = 0
    n_no_valid   = 0

    for (analysis_pop, trait, chrom, start, end,
         completeness_pairs, completeness_files,
         merge_pairs, merge_files,
         merged_output, locus_size) in all_work:

        merged_sentinel = merged_output + ".merged.vcor1"
        vars_sentinel   = merged_output + ".merged.vcor1.vars"

        # Checkpoint: already merged (both matrix and variant list must exist)
        if exists_map.get(merged_sentinel, False) and exists_map.get(vars_sentinel, False):
            n_complete += 1
            continue

        # Completeness check: all matrix-making files must exist
        missing = [f for f in completeness_files if not exists_map.get(f, False)]
        if missing:
            n_incomplete += 1
            continue

        # Filter trait-specific files to non-empty only
        nonempty_files = [(f, cohort, cohort_pop)
                         for f, (cohort, cohort_pop)
                         in zip(merge_files, merge_pairs)
                         if size_map.get(f, 0) > 0]

        if not nonempty_files:
            n_no_valid += 1
            continue

        # Look up sample sizes (merge_pairs are linker-derived so should always resolve)
        ld_files_arg     = []
        sample_sizes_arg = []
        for f, cohort, cohort_pop in nonempty_files:
            n_val = sample_sizes.get((trait, cohort, cohort_pop))
            if n_val is None:
                print(f"WARNING: No sample size for ({trait}, {cohort}, {cohort_pop}), "
                      f"skipping this cohort for {trait}.{analysis_pop} "
                      f"chr{chrom}_{start}_{end}")
                continue
            ld_files_arg.append(f)
            sample_sizes_arg.append(str(n_val))

        if not ld_files_arg:
            n_no_valid += 1
            continue

        # Build command
        ld_files_str     = ' '.join(ld_files_arg)
        sample_sizes_str = ' '.join(sample_sizes_arg)
        cmd = (
            f"python {merge_script}"
            f" --ld-files '{ld_files_str}'"
            f" --sample-sizes '{sample_sizes_str}'"
            f" --output-prefix {merged_output}"
        )

        sort_key = locus_size * len(ld_files_arg)
        job_list.append((sort_key, cmd))

    print(f"NOTE: {len(all_work)} total loci  |  "
          f"{n_complete} already merged  |  "
          f"{n_incomplete} waiting for constituent matrices  |  "
          f"{n_no_valid} with no valid non-empty inputs  |  "
          f"{len(job_list)} to run")

    return job_list


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND RUN  ###########################
###############################################################################

def main():

    # ── Parse arguments ────────────────────────────────────────────────────────
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "command_file="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    command_file = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--command_file":
            command_file = a
        else:
            assert False, f"WARNING: Unhandled option: {o}"

    if command_file is None:
        print("ERROR: Missing required --command_file")
        sys.exit(2)

    # ── Validate input paths ──────────────────────────────────────────────────
    for path, label in [
        (LINKER_FILE,       "linker file"),
        (LOCI_DIR,          "loci directory"),
        (SAMPLE_SIZES_FILE, "sample sizes file"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            if path == SAMPLE_SIZES_FILE:
                print("ERROR: Run extract_sample_sizes.py first.")
            sys.exit(1)

    # ── Parse linker ──────────────────────────────────────────────────────────
    print("NOTE: Parsing linker file for expected cohort pairs...")
    all_pairs, trait_pairs = parse_linker(LINKER_FILE)

    # ── Load sample sizes ─────────────────────────────────────────────────────
    print("NOTE: Loading sample sizes...")
    sample_sizes = load_sample_sizes(SAMPLE_SIZES_FILE)

    # ── Generate commands ─────────────────────────────────────────────────────
    print("NOTE: Generating merge commands with completeness + checkpoint checks...")
    job_list = generate_commands(all_pairs, trait_pairs, sample_sizes, LOCI_DIR, LD_BASE)

    # Sort ascending by sort_key (smallest jobs first)
    job_list.sort(key=lambda x: x[0])

    # ── Write command file ─────────────────────────────────────────────────────
    cmd_dir = os.path.dirname(command_file)
    if cmd_dir:
        os.makedirs(cmd_dir, exist_ok=True)
    with open(command_file, 'w') as fh:
        for _, cmd in job_list:
            fh.write(cmd + '\n')

    print(f"NOTE: {len(job_list)} commands written to {command_file}")


if __name__ == "__main__":
    main()
