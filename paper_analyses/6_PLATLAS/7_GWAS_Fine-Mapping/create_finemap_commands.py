'''
create_finemap_commands.py

Generate fine-mapping commands for GWAS meta-analysis loci (EUR, AFR, ALL).

Reads the three meta-analysis master locus CSV files and, for each
(TRAIT, CHR, START, END, ANALYSIS_POP) locus, emits one SuSiEx command
and one MultiSuSiE command.  Each command calls the per-locus Python
wrapper (run_susiex.py or run_multisusie.py) with comma-separated lists
of per-cohort sumstats paths, LD matrix prefixes, sample sizes, and
ancestry labels.

Inputs:
  Locus files:   {locus_dir}/sig_loci.{EUR,AFR,ALL}.master.csv
                 Columns: CHR, START, END, LENGTH, TRAIT, OVERLAP_COUNT
  Linker file:   {linker_file}
                 Tab-delimited two-column: {TRAIT}.{ANALYSIS_POP}  {rel_path}
                 where rel_path = {COHORT}/{COHORT_POP}/filename.txt.gz
  Sample sizes:  {sample_sizes}
                 Columns: TRAIT, COHORT, COHORT_POP, N
  LD matrices:   {ld_dir}/{TRAIT}.{ANALYSIS_POP}/
                   chr{CHR}_{START}_{END}.{COHORT}.{COHORT_POP}.vcor1

Output file naming:
  SuSiEx:    {output_dir}/susiex/{ANALYSIS_POP}/{TRAIT}.chr{CHR}_{START}_{END}.susiex.txt
  MultiSuSiE:{output_dir}/multisusie/{ANALYSIS_POP}/{TRAIT}.chr{CHR}_{START}_{END}.multisusie.txt

Two-phase design:
  Phase 1 — pure string construction from CSV / linker data; zero filesystem calls.
  Phase 2 — parallel sentinel/input checks using multiprocessing.cpu_count()
             threads; all os.path.exists() calls happen here.

Commands are sorted by locus size (END - START) ascending so that the
smallest (fastest) jobs run first, keeping workers busy throughout.

Usage:
  python create_finemap_commands.py [OPTIONS]

***Run with: conda activate finemap (parsl is installed there)***
'''

import sys
import os
import getopt
import multiprocessing
import pandas as pd
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor


# Parallel workers for filesystem existence checks
NUM_WORKERS = multiprocessing.cpu_count()

# Full paths to each method's Python interpreter
FINEMAP_PYTHON    = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/python"
MULTISUSIE_PYTHON = "/flare/GeomicVar/mconery/tools/miniconda3/envs/MultiSuSiE/bin/python"

# Analysis populations to process
ANALYSIS_POPS = ["EUR", "AFR", "ALL"]


def help(exit_num=0):
    print("""Usage: python create_finemap_commands.py [OPTIONS]
Options:
  --linker-file   Path to meta_sumstats.linker.betanfer.tsv (required)
  --locus-dir     Directory containing sig_loci.{EUR,AFR,ALL}.master.csv (required)
  --ld-dir        Base directory for LD matrices (required)
                  Layout: {ld_dir}/{TRAIT}.{POP}/chr*_{COHORT}.{COHORT_POP}.vcor1
  --sumstats-base Base directory for cohort summary statistics (required)
                  Layout: {sumstats_base}/{COHORT}/{COHORT_POP}/file.txt.gz
  --sample-sizes  Path to sample_sizes.tsv (TRAIT, COHORT, COHORT_POP, N) (required)
  --output-dir    Base directory for fine-mapping results (required)
                  Subdirs susiex/{EUR,AFR,ALL}/timing/ and
                  multisusie/{EUR,AFR,ALL}/timing/ are created.
  --script-dir    Directory containing run_susiex.py and run_multisusie.py (required)
  --output-file   Path to write the commands file (required)
  -v, --verbose   Print all warnings (default: show first 20 only)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  LINKER PARSING  ###################################################
###############################################################################

def parse_linker(linker_file, sumstats_base):
    """Parse the linker file into a dict keyed by (trait, analysis_pop).

    Returns:
        dict: (trait, analysis_pop) -> list of (cohort, cohort_pop, abs_sumstats_path)
    """
    linker = defaultdict(list)
    with open(linker_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            key_col, rel_path = parts[0], parts[1]
            # key_col format: {TRAIT}.{ANALYSIS_POP}
            # Split on last '.' to separate trait (may contain dots) from pop
            dot_idx = key_col.rfind('.')
            if dot_idx < 0:
                continue
            trait        = key_col[:dot_idx]
            analysis_pop = key_col[dot_idx + 1:]
            # rel_path format: {COHORT}/{COHORT_POP}/filename.txt.gz
            path_parts = rel_path.strip().split('/')
            if len(path_parts) < 3:
                continue
            cohort     = path_parts[0]
            cohort_pop = path_parts[1]
            abs_path   = os.path.join(sumstats_base, rel_path.strip())
            linker[(trait, analysis_pop)].append((cohort, cohort_pop, abs_path))
    return linker


###############################################################################
##########  PATH HELPERS (pure string construction)  ##########################
###############################################################################

def ld_vcor1_path(ld_dir, trait, analysis_pop, cohort, cohort_pop, chrom, start, end):
    stem = f"chr{chrom}_{start}_{end}"
    return os.path.join(
        ld_dir, f"{trait}.{analysis_pop}",
        f"{stem}.{cohort}.{cohort_pop}.vcor1"
    )


def susiex_output_path(output_dir, analysis_pop, trait, chrom, start, end):
    stem = f"chr{chrom}_{start}_{end}"
    return os.path.join(output_dir, "susiex", analysis_pop, f"{trait}.{stem}.susiex.txt")


def multisusie_output_path(output_dir, analysis_pop, trait, chrom, start, end):
    stem = f"chr{chrom}_{start}_{end}"
    return os.path.join(output_dir, "multisusie", analysis_pop,
                        f"{trait}.{stem}.multisusie.txt")


def build_command(python_bin, script, out_path, gz_list, ld_list, n_list, label_list,
                  chrom, start, end):
    return (
        f"{python_bin} {script}"
        f" --gwas-files {','.join(gz_list)}"
        f" --ld-files {','.join(ld_list)}"
        f" --n-gwas {','.join(str(n) for n in n_list)}"
        f" --ancestry-labels {','.join(label_list)}"
        f" --output-file {out_path}"
        f" --chr {chrom} --start {start} --end {end}"
    )


###############################################################################
##########  PHASE 1: BUILD WORK ITEMS (no filesystem calls)  #################
###############################################################################

def build_work_item(row, analysis_pop, linker, sample_sizes, ld_dir, output_dir, script_dir):
    """Construct all paths and metadata for one locus.

    Returns a dict with per-cohort lists and locus metadata, or None if the
    (trait, analysis_pop) combination is not present in the linker.
    """
    trait      = row['TRAIT']
    chrom      = int(row['CHR'])
    start      = int(row['START'])
    end        = int(row['END'])
    locus_size = end - start

    cohort_entries = linker.get((trait, analysis_pop), [])
    if not cohort_entries:
        return None

    gz_list    = []
    ld_list    = []
    n_list     = []
    label_list = []

    for cohort, cohort_pop, abs_sumstats in cohort_entries:
        ld_path = ld_vcor1_path(ld_dir, trait, analysis_pop, cohort, cohort_pop,
                                chrom, start, end)
        n_val   = sample_sizes.get((trait, cohort, cohort_pop), 0)
        gz_list.append(abs_sumstats)
        ld_list.append(ld_path)
        n_list.append(n_val)
        label_list.append(f"{cohort}_{cohort_pop}")

    return {
        'locus_size':   locus_size,
        'trait':        trait,
        'analysis_pop': analysis_pop,
        'chrom':        chrom,
        'start':        start,
        'end':          end,
        'gz_list':      gz_list,
        'ld_list':      ld_list,
        'n_list':       n_list,
        'label_list':   label_list,
        'susiex_out':   susiex_output_path(output_dir, analysis_pop, trait, chrom, start, end),
        'multisusie_out': multisusie_output_path(output_dir, analysis_pop, trait,
                                                 chrom, start, end),
        'script_dir':   script_dir,
    }


###############################################################################
##########  PHASE 2: CHECK AND FILTER (all filesystem calls here)  ###########
###############################################################################

def check_work_item(item):
    """Check which commands are needed for this locus.

    Returns (commands_list, warnings_list, skipped_count).

    Fast path: both outputs exist → return immediately without checking inputs.
    Slow path: at least one output missing → validate inputs, emit commands.
    """
    susiex_done     = os.path.exists(item['susiex_out'])
    multisusie_done = os.path.exists(item['multisusie_out'])

    # Fast path: both complete
    if susiex_done and multisusie_done:
        return [], [], 2

    # Validate per-cohort inputs; build filtered lists excluding empty/missing LD files.
    valid_gz    = []
    valid_ld    = []
    valid_n     = []
    valid_label = []
    warnings    = []

    for gz, ld_path, n_val, label in zip(
        item['gz_list'], item['ld_list'], item['n_list'], item['label_list']
    ):
        vars_file = ld_path + ".vars"
        if not os.path.exists(gz):
            warnings.append(
                f"WARN: sumstats not found for "
                f"{item['trait']}.{item['analysis_pop']}: {gz}"
            )
            continue
        # Missing LD file → betanfer produced no output for this cohort; skip silently
        if not os.path.exists(ld_path):
            continue
        # 0-byte LD file → no variants passed betanfer p-value threshold; skip
        if os.path.getsize(ld_path) == 0:
            continue
        if not os.path.exists(vars_file):
            warnings.append(f"WARN: .vars file not found: {vars_file}")
            continue
        valid_gz.append(gz)
        valid_ld.append(ld_path)
        valid_n.append(n_val)
        valid_label.append(label)

    if not valid_ld:
        return [], warnings, 0

    script_dir  = item['script_dir']
    chrom, start, end = item['chrom'], item['start'], item['end']
    commands = []
    skipped  = 0

    if susiex_done:
        skipped += 1
    else:
        cmd = build_command(
            FINEMAP_PYTHON,
            f"{script_dir}/run_susiex.py",
            item['susiex_out'],
            valid_gz, valid_ld, valid_n, valid_label,
            chrom, start, end,
        )
        commands.append((item['locus_size'], cmd))

    if multisusie_done:
        skipped += 1
    else:
        cmd = build_command(
            MULTISUSIE_PYTHON,
            f"{script_dir}/run_multisusie.py",
            item['multisusie_out'],
            valid_gz, valid_ld, valid_n, valid_label,
            chrom, start, end,
        )
        commands.append((item['locus_size'], cmd))

    return commands, warnings, skipped


###############################################################################
##########  MAIN  #############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "hv",
            ["help", "verbose", "linker-file=", "locus-dir=", "ld-dir=",
             "sumstats-base=", "sample-sizes=", "output-dir=",
             "script-dir=", "output-file="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    linker_file       = None
    locus_dir         = None
    ld_dir            = None
    sumstats_base     = None
    sample_sizes_file = None
    output_dir        = None
    script_dir        = None
    output_file       = None
    verbose           = False

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o == "--linker-file":
            linker_file = a
        elif o == "--locus-dir":
            locus_dir = a.rstrip('/')
        elif o == "--ld-dir":
            ld_dir = a.rstrip('/')
        elif o == "--sumstats-base":
            sumstats_base = a.rstrip('/')
        elif o == "--sample-sizes":
            sample_sizes_file = a
        elif o == "--output-dir":
            output_dir = a.rstrip('/')
        elif o == "--script-dir":
            script_dir = a.rstrip('/')
        elif o == "--output-file":
            output_file = a
        else:
            assert False, f"Unhandled option: {o}"

    required = {
        "--linker-file":  linker_file,
        "--locus-dir":    locus_dir,
        "--ld-dir":       ld_dir,
        "--sumstats-base": sumstats_base,
        "--sample-sizes": sample_sizes_file,
        "--output-dir":   output_dir,
        "--script-dir":   script_dir,
        "--output-file":  output_file,
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        print(f"ERROR: Missing required arguments: {', '.join(missing)}")
        sys.exit(2)

    for path, label in [
        (linker_file,       "linker file"),
        (locus_dir,         "locus directory"),
        (ld_dir,            "LD matrices directory"),
        (sumstats_base,     "sumstats base directory"),
        (sample_sizes_file, "sample sizes file"),
        (script_dir,        "script directory"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Create output directories ──────────────────────────────────────────────
    print("NOTE: Creating output directories...")
    for method in ("susiex", "multisusie"):
        for pop in ANALYSIS_POPS:
            os.makedirs(os.path.join(output_dir, method, pop, "timing"), exist_ok=True)
    print(f"      {output_dir}/{{susiex,multisusie}}/{{EUR,AFR,ALL}}/timing/ ready")

    # ── Parse linker file ──────────────────────────────────────────────────────
    print("NOTE: Parsing linker file...")
    linker = parse_linker(linker_file, sumstats_base)
    print(f"      {len(linker)} (trait, analysis_pop) combinations loaded")

    # ── Load sample sizes ──────────────────────────────────────────────────────
    print("NOTE: Loading sample sizes...")
    ss_df = pd.read_csv(sample_sizes_file, sep='\t')
    sample_sizes = {
        (row['TRAIT'], row['COHORT'], row['COHORT_POP']): int(row['N'])
        for _, row in ss_df.iterrows()
    }
    print(f"      {len(sample_sizes)} (trait, cohort, cohort_pop) entries loaded")

    # ── Read master locus files and build all work items ───────────────────────
    all_items = []
    locus_counts = {}

    for pop in ANALYSIS_POPS:
        master_path = os.path.join(locus_dir, f"sig_loci.{pop}.master.csv")
        if not os.path.exists(master_path):
            print(f"WARNING: Master locus file not found, skipping {pop}: {master_path}")
            continue
        df = pd.read_csv(master_path)
        locus_counts[pop] = len(df)
        print(f"NOTE: {pop}: {len(df)} loci")

        for _, row in df.iterrows():
            item = build_work_item(
                row, pop, linker, sample_sizes, ld_dir, output_dir, script_dir
            )
            if item is not None:
                all_items.append(item)

    total_loci = sum(locus_counts.values())
    skipped_linker = total_loci - len(all_items)
    print(f"NOTE: Built {len(all_items)} work items "
          f"({skipped_linker} loci had no linker entries and were skipped)")

    # ── Phase 2: check sentinels and inputs in parallel ───────────────────────
    print(f"NOTE: Checking {len(all_items)} loci ({NUM_WORKERS} workers)...")
    all_commands  = []
    all_warnings  = []
    total_skipped = 0

    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as pool:
        for cmds, warns, skipped in pool.map(check_work_item, all_items):
            all_commands.extend(cmds)
            all_warnings.extend(warns)
            total_skipped += skipped

    # ── Sort by locus size ascending (smallest/fastest jobs first) ─────────────
    all_commands.sort(key=lambda x: x[0])
    commands_only = [cmd for _, cmd in all_commands]

    # ── Print warnings ─────────────────────────────────────────────────────────
    if all_warnings:
        print(f"\nWARNINGS ({len(all_warnings)}):")
        warnings_to_show = all_warnings if verbose else all_warnings[:20]
        for w in warnings_to_show:
            print(f"  {w}")
        if not verbose and len(all_warnings) > 20:
            print(f"  ... and {len(all_warnings) - 20} more (use --verbose to show all)")

    # ── Summary ────────────────────────────────────────────────────────────────
    susiex_count     = sum(1 for c in commands_only if "run_susiex.py"     in c)
    multisusie_count = sum(1 for c in commands_only if "run_multisusie.py" in c)
    print(f"\nNOTE: Processing summary:")
    print(f"      Total loci processed : {len(all_items)}")
    print(f"      Already complete     : {total_skipped}")
    print(f"      Commands to run      : {len(commands_only)}")
    print(f"        SuSiEx             : {susiex_count}")
    print(f"        MultiSuSiE         : {multisusie_count}")

    # ── Write commands file ────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(output_file, 'w') as fh:
        for cmd in commands_only:
            fh.write(cmd + '\n')
    print(f"\nSUCCESS: Wrote {len(commands_only)} commands to {output_file}")


if __name__ == "__main__":
    main()
