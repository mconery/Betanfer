'''
create_1000g_commands.py

Generates plink2 --r-unphased commands for computing LD matrices from the
1000 Genomes Phase 3 hg19 reference panel for every significant locus in the
MVP_AFR and MVP_EUR population master CSVs.

Two steps:
  1. build_keep_files() — parse all_phase3.psam and write per-SuperPop keep
     files (1000G_phase3.AFR.txt, 1000G_phase3.EUR.txt) to the 1000G directory.
     Skipped if files already exist.
  2. generate_commands() — for each (locus, population) pair, construct a
     plink2 command that filters to the matching 1000G samples, restricts to
     the locus window, applies a MAF filter, and writes a square LD matrix in
     plink2 vcor1 format.

Commands are sorted ascending by locus size so small jobs finish first.
The sentinel file is {output_prefix}.vcor1.vars; existing sentinels are skipped.

Output naming convention (matches existing betanfer LD matrices):
  {LD_BASE}/{TRAIT}.MVP_{POP}/chr{CHR}_{START}_{END}.1000G.{POP}.vcor1

Usage:
  python create_1000g_commands.py --command_file <path> [--log_dir <path>]

***Run with: conda activate /flare/GeomicVar/mconery/tools/miniconda3/envs/finemap***
'''

import sys
import getopt
import os
import csv
import glob
import multiprocessing
from concurrent.futures import ThreadPoolExecutor

# ── Default paths ───────────────────────────────────────────────────────────────
LOCI_DIR    = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci"
LD_BASE     = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices"
PLINK2_BIN  = "/flare/GeomicVar/mconery/tools/bin/plink2"
KG_PFILE    = "/flare/GeomicVar/mconery/resources/1000G/all_phase3"
KG_KEEP_DIR = "/flare/GeomicVar/mconery/resources/1000G"
KG_PSAM     = "/flare/GeomicVar/mconery/resources/1000G/all_phase3.psam"
NUM_THREADS = 208
MAF_FILTER  = "0.01"

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))

# MVP analysis populations and their corresponding 1000G SuperPop codes
MVP_POPS = {
    "MVP_AFR": "AFR",
    "MVP_EUR": "EUR",
}

# Parallel workers for sentinel existence checks (Lustre handles concurrent stats well)
NUM_WORKERS = multiprocessing.cpu_count()


def help(exit_num=0):
    print("""Usage: python create_1000g_commands.py [OPTIONS]
Options:
  --command_file  Path to write the plink2 commands file (required)
  --log_dir       Directory for parsl worker log output (optional; unused here)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
####################  STEP 1: BUILD PER-POPULATION KEEP FILES  ################
###############################################################################

def build_keep_files():
    """Parse all_phase3.psam and write per-SuperPop keep files.

    Reads the SuperPop column from the psam header, extracts IIDs for AFR and
    EUR, and writes 1000G_phase3.{POP}.txt to KG_KEEP_DIR.  Each keep file has
    a '#IID' header followed by one sample ID per line, as required by plink2
    --keep.

    Files are only written if they do not already exist.  Returns a dict
    mapping SuperPop code → keep file path for the two target populations.
    """
    target_pops = set(MVP_POPS.values())   # {"AFR", "EUR"}
    keep_files  = {}                        # pop → path

    # Check which keep files are missing
    missing_pops = []
    for pop in target_pops:
        fpath = os.path.join(KG_KEEP_DIR, f"1000G_phase3.{pop}.txt")
        keep_files[pop] = fpath
        if os.path.exists(fpath):
            print(f"NOTE: Keep file already exists, skipping: {fpath}")
        else:
            missing_pops.append(pop)

    if not missing_pops:
        return keep_files

    # Parse psam to build {pop → [iid, ...]} mapping
    samples_by_pop = {pop: [] for pop in missing_pops}

    with open(KG_PSAM, 'r') as fh:
        header_line = fh.readline().rstrip('\n')
        # Remove leading '#' if present, then split on whitespace
        cols = header_line.lstrip('#').split()
        # First column is IID (psam may name it 'IID' or appear after a '#')
        iid_idx    = 0   # always first column in PLINK2 psam
        # Find SuperPop column
        if 'SuperPop' not in cols:
            print(f"ERROR: 'SuperPop' column not found in psam header: {cols}")
            sys.exit(1)
        superpop_idx = cols.index('SuperPop')

        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) <= superpop_idx:
                continue
            iid      = parts[iid_idx]
            superpop = parts[superpop_idx]
            if superpop in samples_by_pop:
                samples_by_pop[superpop].append(iid)

    # Write keep files for missing populations
    for pop in missing_pops:
        iids  = samples_by_pop[pop]
        fpath = keep_files[pop]
        with open(fpath, 'w') as fh:
            fh.write('#IID\n')
            for iid in iids:
                fh.write(iid + '\n')
        print(f"NOTE: Wrote {len(iids)} {pop} samples to {fpath}")

    return keep_files


###############################################################################
####################  STEP 2: GENERATE PLINK2 COMMANDS  ######################
###############################################################################

def generate_commands(keep_files, loci_dir, ld_base):
    """Scan MVP population master CSVs and produce plink2 LD commands.

    Each command targets one (locus, population) combination.  Jobs whose
    sentinel .vcor1.vars file already exists are skipped (checkpoint).

    Returns a list of (locus_size_bp, command_str) tuples.
    """
    warnings = 0

    # Phase 1: build all work items from CSV data — no filesystem calls.
    all_work = []   # [(locus_size, sentinel_path, cmd_str), ...]

    for analysis_pop, kg_pop in sorted(MVP_POPS.items()):
        master_csv = os.path.join(loci_dir, f"sig_loci.{analysis_pop}.master.csv")
        if not os.path.exists(master_csv):
            print(f"WARNING: Loci file not found, skipping {analysis_pop}: {master_csv}")
            warnings += 1
            continue

        keep_file = keep_files[kg_pop]

        with open(master_csv, 'r') as fh:
            reader = csv.DictReader(fh)
            rows   = list(reader)

        loci_count = 0
        for row in rows:
            chrom      = row['CHR']
            start      = row['START']
            end        = row['END']
            trait      = row['TRAIT']
            locus_size = int(end) - int(start)

            output_prefix = os.path.join(
                ld_base,
                f"{trait}.{analysis_pop}",
                f"chr{chrom}_{start}_{end}.1000G.{kg_pop}"
            )
            sentinel = output_prefix + ".unphased.vcor1.vars"

            # plink2 --set-all-var-ids uses $r/$a as REF/ALT tokens; wrap in
            # double quotes so bash does not expand them as shell variables.
            # --new-id-max-allele-len 1000 suppresses the error plink2 raises
            # when an allele code exceeds the default 23-character limit (common
            # in Phase 3 data for long indels/structural variants).  Such
            # variants will not match sumstats IDs anyway, so they are naturally
            # excluded in the downstream fine-mapping intersection step.
            cmd = (
                f"{PLINK2_BIN}"
                f" --pfile {KG_PFILE}"
                f" --keep {keep_file}"
                f" --chr {chrom}"
                f" --from-bp {start}"
                f" --to-bp {end}"
                f" --maf {MAF_FILTER}"
                f' --set-all-var-ids "@:#:\\$r:\\$a"'
                f" --new-id-max-allele-len 1000"
                f" --r-unphased square"
                f" --out {output_prefix}"
                f" --threads {NUM_THREADS}"
            )
            all_work.append((locus_size, sentinel, cmd))
            loci_count += 1

        print(f"NOTE: {analysis_pop}: processed {loci_count} loci")

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

    print(f"NOTE: {total} total (locus, population) pairs  |  "
          f"{skipped} already complete  |  {len(job_list)} to run  |  "
          f"{warnings} population(s) skipped")

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
    log_dir      = None   # accepted for CLI consistency; unused in this script

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
        (KG_PSAM,   "1000G psam file"),
        (LOCI_DIR,  "loci directory"),
        (PLINK2_BIN, "plink2 binary"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    for suffix in (".pgen", ".pvar", ".psam"):
        if not os.path.exists(KG_PFILE + suffix):
            print(f"ERROR: 1000G pfile missing: {KG_PFILE + suffix}")
            sys.exit(1)

    os.makedirs(LD_BASE, exist_ok=True)

    # ── Step 1: Build per-population keep files ─────────────────────────────────
    print("NOTE: Building 1000G population keep files from all_phase3.psam...")
    keep_files = build_keep_files()

    # Verify keep files exist before proceeding
    for pop, fpath in keep_files.items():
        if not os.path.exists(fpath):
            print(f"ERROR: Keep file missing after generation: {fpath}")
            sys.exit(1)

    # ── Step 2: Generate plink2 commands ────────────────────────────────────────
    print("NOTE: Generating plink2 LD matrix commands...")
    job_list = generate_commands(keep_files, LOCI_DIR, LD_BASE)

    # Sort ascending by locus size (small loci first)
    job_list.sort(key=lambda x: x[0])

    # Create output trait directories upfront so parsl workers can write immediately
    for _, cmd in job_list:
        parts = cmd.split()
        for i, part in enumerate(parts):
            if part == "--out" and i + 1 < len(parts):
                out_dir = os.path.dirname(parts[i + 1])
                if out_dir:
                    os.makedirs(out_dir, exist_ok=True)
                break

    # ── Write commands file ──────────────────────────────────────────────────────
    cmd_dir = os.path.dirname(command_file)
    if cmd_dir:
        os.makedirs(cmd_dir, exist_ok=True)
    with open(command_file, 'w') as fh:
        for _, cmd in job_list:
            fh.write(cmd + '\n')

    print(f"NOTE: {len(job_list)} commands written to {command_file}")
    if job_list:
        min_kb = min(s for s, _ in job_list) / 1e3
        max_mb = max(s for s, _ in job_list) / 1e6
        print(f"NOTE: Locus sizes: smallest {min_kb:.1f} kb  |  largest {max_mb:.1f} Mb")


if __name__ == "__main__":
    main()
