'''
extract_sample_sizes.py

One-time extraction of sample sizes (max N) from every sumstats file referenced
in the linker file.  Results are saved to a TSV so that subsequent runs of the
merging pipeline can look up sample sizes without re-reading any sumstats files.

For each unique (trait, cohort, cohort_pop) combination in the linker, the script
decompresses the corresponding bgzip sumstats file and reads the first non-header
data line to extract column 17 (the N column).  N is constant within a
trait-specific sumstats file, so reading one line is sufficient.

Parallelises over unique file paths using ProcessPoolExecutor.

Output TSV columns: TRAIT, COHORT, COHORT_POP, N

Usage:
  python extract_sample_sizes.py [--force-recompute]

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os
import subprocess
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

# ── Default paths ────────────────────────────────────────────────────────────
LINKER_FILE      = "/flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv"
SUMSTATS_BASE    = "/flare/GeomicVar/multi-biobank-meta/sumstats"
LD_BASE          = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices"
SAMPLE_SIZES_FILE = os.path.join(LD_BASE, "sample_sizes.tsv")

# Populations served by the linker (meta-analysis populations only)
LINKER_POPS = {"EUR", "AFR", "ALL"}

# N is column 17 (1-indexed) in the sumstats header:
# #ID CHR POS REF ALT BETA SE P LOG10P SE_LDSC P_LDSC LOG10P_LDSC AAF AAF_CASE AAC AAC_CASE N ...
N_COL_INDEX = 16  # 0-indexed


def help(exit_num=0):
    print("""Usage: python extract_sample_sizes.py [OPTIONS]
Options:
  --force-recompute   Overwrite existing sample_sizes.tsv even if it exists
  -h, --help          Show this help message""")
    sys.exit(exit_num)


###############################################################################
################  PARSE LINKER → UNIQUE (trait, cohort, cohort_pop, path)  ####
###############################################################################

def parse_linker(linker_file):
    """Parse linker file and return unique (trait, cohort, cohort_pop, abs_path) tuples.

    Deduplicates so that the same physical file is processed at most once per
    (trait, cohort, cohort_pop) key (a file can appear for both ALL and AFR).
    Returns a list of (trait, cohort, cohort_pop, abs_path) with no duplicates.
    """
    seen_keys = set()
    entries   = []

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

            # Extract analysis population and trait name
            dot_parts = trait_pop.rsplit('.', 1)
            if len(dot_parts) != 2:
                print(f"WARNING: Cannot parse trait.pop from '{trait_pop}', skipping")
                continue
            trait, analysis_pop = dot_parts
            if analysis_pop not in LINKER_POPS:
                continue  # skip MVP_EUR / MVP_AFR single-cohort entries

            # Extract cohort and cohort_pop from path: COHORT/COHORT_POP/filename
            path_parts = rel_path.split('/')
            if len(path_parts) < 3:
                print(f"WARNING: Cannot parse cohort/pop from path '{rel_path}', skipping")
                continue
            cohort     = path_parts[0]
            cohort_pop = path_parts[1]
            abs_path   = os.path.join(SUMSTATS_BASE, rel_path)

            key = (trait, cohort, cohort_pop)
            if key not in seen_keys:
                seen_keys.add(key)
                entries.append((trait, cohort, cohort_pop, abs_path))

    return entries


###############################################################################
############  WORKER: EXTRACT N FROM FIRST DATA LINE OF A SUMSTATS FILE  #####
###############################################################################

def extract_n_from_file(args):
    """Extract N from the first non-header data line of a bgzip sumstats file.

    args: (trait, cohort, cohort_pop, abs_path)
    Returns (trait, cohort, cohort_pop, N) or None on error.
    """
    trait, cohort, cohort_pop, abs_path = args

    if not os.path.exists(abs_path):
        print(f"WARNING: File not found, skipping: {abs_path}")
        return None

    try:
        result = subprocess.run(
            ["zcat", abs_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        lines = result.stdout.decode('utf-8', errors='replace').splitlines()

        # Skip header line (starts with '#ID' or '#')
        data_line = None
        for line in lines:
            if not line.startswith('#'):
                data_line = line
                break

        if data_line is None:
            print(f"WARNING: No data lines found in {abs_path}, skipping")
            return None

        fields = data_line.split('\t')
        if len(fields) <= N_COL_INDEX:
            print(f"WARNING: Line has fewer than {N_COL_INDEX + 1} columns in {abs_path}, skipping")
            return None

        n_val = fields[N_COL_INDEX].strip()
        # N may be an integer or float string
        n_int = int(float(n_val))
        return (trait, cohort, cohort_pop, n_int)

    except Exception as e:
        print(f"WARNING: Error reading {abs_path}: {e}")
        return None


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND RUN  ###########################
###############################################################################

def main():

    # ── Parse arguments ────────────────────────────────────────────────────────
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "force-recompute"]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    force_recompute = False

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--force-recompute":
            force_recompute = True
        else:
            assert False, f"WARNING: Unhandled option: {o}"

    # ── Checkpoint: skip if output already exists ─────────────────────────────
    if os.path.exists(SAMPLE_SIZES_FILE) and not force_recompute:
        print(f"NOTE: Sample sizes file already exists: {SAMPLE_SIZES_FILE}")
        print("NOTE: Use --force-recompute to overwrite.")
        sys.exit(0)

    # ── Validate paths ────────────────────────────────────────────────────────
    if not os.path.exists(LINKER_FILE):
        print(f"ERROR: Linker file not found: {LINKER_FILE}")
        sys.exit(1)

    os.makedirs(LD_BASE, exist_ok=True)

    # ── Parse linker ──────────────────────────────────────────────────────────
    print("NOTE: Parsing linker file...")
    entries = parse_linker(LINKER_FILE)
    print(f"NOTE: Found {len(entries)} unique (trait, cohort, cohort_pop) entries to process")

    # ── Extract N values in parallel ──────────────────────────────────────────
    num_workers = multiprocessing.cpu_count()
    print(f"NOTE: Extracting N values using {num_workers} workers...")

    results = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for result in executor.map(extract_n_from_file, entries):
            if result is not None:
                results.append(result)

    print(f"NOTE: Successfully extracted N for {len(results)} of {len(entries)} entries")

    # ── Write output TSV ──────────────────────────────────────────────────────
    with open(SAMPLE_SIZES_FILE, 'w') as fh:
        fh.write("TRAIT\tCOHORT\tCOHORT_POP\tN\n")
        for trait, cohort, cohort_pop, n_val in sorted(results):
            fh.write(f"{trait}\t{cohort}\t{cohort_pop}\t{n_val}\n")

    print(f"SUCCESS: Wrote {len(results)} sample size entries to {SAMPLE_SIZES_FILE}")


if __name__ == "__main__":
    main()
