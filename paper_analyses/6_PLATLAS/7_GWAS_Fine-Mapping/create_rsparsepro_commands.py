'''
create_rsparsepro_commands.py

Generate RSparsePro fine-mapping commands for GWAS meta-analysis loci (EUR, AFR, ALL).

For each (TRAIT, CHR, START, END, ANALYSIS_POP) locus in the master locus files,
emits one RSparsePro command.  Unlike the SuSiEx/MultiSuSiE pipeline, RSparsePro
takes the population-level meta-analysis sumstats and the merged weighted-average LD
matrix (produced by step 4) rather than per-cohort files.

Inputs:
  Locus files:   {locus_dir}/sig_loci.{EUR,AFR,ALL}.master.csv
                 Columns: CHR, START, END, LENGTH, TRAIT, OVERLAP_COUNT
  Meta-analysis: {meta_dir}/{POP}/{TRAIT}.out.gz  (bgzip+tabix indexed)
  LD matrices:   {ld_dir}/{TRAIT}.{POP}/
                   chr{CHR}_{START}_{END}.merged.merged.vcor1
                   chr{CHR}_{START}_{END}.merged.merged.vcor1.vars

Output file naming:
  {output_dir}/rsparsepro/{ANALYSIS_POP}/{TRAIT}.chr{CHR}_{START}_{END}.rsparsepro.txt

Two-phase design:
  Phase 1 — pure string construction from CSV data; zero filesystem calls.
  Phase 2 — parallel sentinel/input checks using multiprocessing.cpu_count() threads.

Commands are sorted by locus size (END - START) ascending so that the smallest
(fastest) jobs run first, keeping workers busy throughout.

Usage:
  python create_rsparsepro_commands.py [OPTIONS]

***Run with: conda activate finemap (parsl is installed there)***
'''

import sys
import os
import getopt
import multiprocessing
import pandas as pd
from concurrent.futures import ThreadPoolExecutor


NUM_WORKERS    = multiprocessing.cpu_count()
FINEMAP_PYTHON = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/python"
ANALYSIS_POPS  = ["EUR", "AFR", "ALL"]


def help(exit_num=0):
    print("""Usage: python create_rsparsepro_commands.py [OPTIONS]
Options:
  --locus-dir   Directory containing sig_loci.{EUR,AFR,ALL}.master.csv (required)
  --meta-dir    Base directory for meta-analysis sumstats (required)
                Layout: {meta_dir}/{POP}/{TRAIT}.out.gz
  --ld-dir      Base directory for merged LD matrices (required)
                Layout: {ld_dir}/{TRAIT}.{POP}/chr*_*.merged.merged.vcor1
  --output-dir  Base directory for fine-mapping results (required)
                Subdir rsparsepro/{EUR,AFR,ALL}/timing/ will be created.
  --script-dir  Directory containing run_rsparsepro.py (required)
  --output-file Path to write the commands file (required)
  -v, --verbose Print all warnings (default: show first 20 only)
  -h, --help    Show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  PATH HELPERS (pure string construction)  ##########################
###############################################################################

def meta_sumstats_path(meta_dir, pop, trait):
    return os.path.join(meta_dir, pop, f"{trait}.out.gz")


def merged_ld_path(ld_dir, trait, pop, chrom, start, end):
    stem = f"chr{chrom}_{start}_{end}"
    return os.path.join(ld_dir, f"{trait}.{pop}", f"{stem}.merged.merged.vcor1")


def rsparsepro_output_path(output_dir, pop, trait, chrom, start, end):
    stem = f"chr{chrom}_{start}_{end}"
    return os.path.join(output_dir, "rsparsepro", pop,
                        f"{trait}.{stem}.rsparsepro.txt")


def build_command(script_dir, sumstats, ld_file, ld_vars, chrom, start, end, out_path):
    return (
        f"{FINEMAP_PYTHON} {script_dir}/run_rsparsepro.py"
        f" --sumstats {sumstats}"
        f" --ld-file {ld_file}"
        f" --ld-vars {ld_vars}"
        f" --chr {chrom} --start {start} --end {end}"
        f" --output-file {out_path}"
    )


###############################################################################
##########  PHASE 1: BUILD WORK ITEMS (no filesystem calls)  #################
###############################################################################

def build_work_item(row, pop, meta_dir, ld_dir, output_dir, script_dir):
    trait      = row['TRAIT']
    chrom      = int(row['CHR'])
    start      = int(row['START'])
    end        = int(row['END'])
    locus_size = end - start
    ld_file    = merged_ld_path(ld_dir, trait, pop, chrom, start, end)
    return {
        'locus_size':  locus_size,
        'trait':       trait,
        'pop':         pop,
        'chrom':       chrom,
        'start':       start,
        'end':         end,
        'sumstats':    meta_sumstats_path(meta_dir, pop, trait),
        'ld_file':     ld_file,
        'ld_vars':     ld_file + ".vars",
        'output_path': rsparsepro_output_path(output_dir, pop, trait, chrom, start, end),
        'script_dir':  script_dir,
    }


###############################################################################
##########  PHASE 2: CHECK AND FILTER (all filesystem calls here)  ###########
###############################################################################

def check_work_item(item):
    """Check whether this locus needs a command and validate its inputs.

    Returns (commands_list, warnings_list, skipped_count).
    commands_list contains (locus_size, command_string) tuples (0 or 1 entry).
    """
    # Fast path: output already exists
    if os.path.exists(item['output_path']):
        return [], [], 1

    warnings = []

    if not os.path.exists(item['sumstats']):
        warnings.append(
            f"WARN: meta-analysis sumstats not found for "
            f"{item['trait']}.{item['pop']}: {item['sumstats']}"
        )
        return [], warnings, 0

    # Missing merged LD → step 4 not yet complete for this locus; skip silently
    if not os.path.exists(item['ld_file']):
        return [], warnings, 0

    if os.path.getsize(item['ld_file']) == 0:
        return [], warnings, 0

    if not os.path.exists(item['ld_vars']):
        warnings.append(f"WARN: .vars file not found: {item['ld_vars']}")
        return [], warnings, 0

    if os.path.getsize(item['ld_vars']) == 0:
        warnings.append(f"WARN: .vars file is empty: {item['ld_vars']}")
        return [], warnings, 0

    cmd = build_command(
        item['script_dir'],
        item['sumstats'],
        item['ld_file'],
        item['ld_vars'],
        item['chrom'],
        item['start'],
        item['end'],
        item['output_path'],
    )
    return [(item['locus_size'], cmd)], warnings, 0


###############################################################################
##########  MAIN  #############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "hv",
            ["help", "verbose", "locus-dir=", "meta-dir=", "ld-dir=",
             "output-dir=", "script-dir=", "output-file="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    locus_dir   = None
    meta_dir    = None
    ld_dir      = None
    output_dir  = None
    script_dir  = None
    output_file = None
    verbose     = False

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o == "--locus-dir":
            locus_dir = a.rstrip('/')
        elif o == "--meta-dir":
            meta_dir = a.rstrip('/')
        elif o == "--ld-dir":
            ld_dir = a.rstrip('/')
        elif o == "--output-dir":
            output_dir = a.rstrip('/')
        elif o == "--script-dir":
            script_dir = a.rstrip('/')
        elif o == "--output-file":
            output_file = a
        else:
            assert False, f"Unhandled option: {o}"

    required = {
        "--locus-dir":   locus_dir,
        "--meta-dir":    meta_dir,
        "--ld-dir":      ld_dir,
        "--output-dir":  output_dir,
        "--script-dir":  script_dir,
        "--output-file": output_file,
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        print(f"ERROR: Missing required arguments: {', '.join(missing)}")
        sys.exit(2)

    for path, label in [
        (locus_dir,  "locus directory"),
        (meta_dir,   "meta-analysis directory"),
        (ld_dir,     "LD matrices directory"),
        (script_dir, "script directory"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Create output directories ──────────────────────────────────────────────
    print("NOTE: Creating output directories...")
    for pop in ANALYSIS_POPS:
        os.makedirs(os.path.join(output_dir, "rsparsepro", pop, "timing"), exist_ok=True)
    print(f"      {output_dir}/rsparsepro/{{EUR,AFR,ALL}}/timing/ ready")

    # ── Read master locus files and build work items ───────────────────────────
    all_items    = []
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
            all_items.append(
                build_work_item(row, pop, meta_dir, ld_dir, output_dir, script_dir)
            )

    total_loci = sum(locus_counts.values())
    print(f"NOTE: Built {total_loci} work items across all populations")

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
    print(f"\nNOTE: Processing summary:")
    print(f"      Total loci processed : {len(all_items)}")
    print(f"      Already complete     : {total_skipped}")
    print(f"      Commands to run      : {len(commands_only)}")

    # ── Write commands file ────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(output_file, 'w') as fh:
        for cmd in commands_only:
            fh.write(cmd + '\n')
    print(f"\nSUCCESS: Wrote {len(commands_only)} commands to {output_file}")


if __name__ == "__main__":
    main()
