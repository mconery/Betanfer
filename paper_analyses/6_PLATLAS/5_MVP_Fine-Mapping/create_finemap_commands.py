'''
create_finemap_commands.py

Generate fine-mapping commands for MVP population-specific loci.

Reads the two MVP master locus CSV files (sig_loci.MVP_AFR.master.csv and
sig_loci.MVP_EUR.master.csv) and, for each (TRAIT, CHR, START, END, POP)
row, emits one command for RSparsePro and one for SuSiE.  Each command is
a self-contained bash invocation of the per-locus wrapper script
(run_rsparsepro.sh or run_susie.sh), which handles tabix extraction from
the bgzip-compressed MVP summary statistics before running the tool.

Master CSV format (produced by create_mvp_loci_from_s11.py):
  CHR, START, END, LENGTH, TRAIT, OVERLAP_COUNT

Sumstats path convention:
  {sumstats_dir}/{POP}/MVP_R4.1000G_AGR.{TRAIT}.{POP}.GIA.dbGaP.sumstats.hg19.txt.gz

LD matrix path convention:
  {ld_dir}/{TRAIT}.MVP_{POP}/chr{CHR}_{START}_{END}.MVP.{POP}.vcor1

Output file naming:
  RSparsePro: {output_dir}/rsparsepro/{POP}/{TRAIT}.chr{CHR}_{START}_{END}.rsparsepro.txt
  SuSiE:      {output_dir}/susie/{POP}/{TRAIT}.chr{CHR}_{START}_{END}.susie.txt

Two-phase design (mirrors create_matrix_commands.py):
  Phase 1 — pure string construction from CSV data; zero filesystem calls.
  Phase 2 — parallel sentinel/input checks using multiprocessing.cpu_count()
             threads; all os.path.exists() calls happen here.

Output files are only written on successful completion (bash scripts use a
temp-file + rename pattern), so os.path.exists() is a reliable completion
sentinel.  Input validation is skipped for loci where both outputs exist.

Commands are sorted by locus size (END - START) ascending so that the
smallest (fastest) jobs run first, keeping workers busy throughout.

Usage:
  python create_finemap_commands.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import os
import getopt
import multiprocessing
import pandas as pd
from concurrent.futures import ThreadPoolExecutor


# Parallel workers for filesystem existence checks
NUM_WORKERS = multiprocessing.cpu_count()


def help(exit_num=0):
    print("""Usage: python create_finemap_commands.py [OPTIONS]
Options:
  --master_afr    Path to sig_loci.MVP_AFR.master.csv (required)
  --master_eur    Path to sig_loci.MVP_EUR.master.csv (required)
  --sumstats_dir  Base directory for MVP_hg19 summary statistics (required)
                  Expected layout: {sumstats_dir}/{AFR,EUR}/MVP_R4.1000G_AGR.*.txt.gz
  --ld_dir        Base directory for LD matrices (required)
                  Expected layout: {ld_dir}/{TRAIT}.MVP_{POP}/chr*_*_*.MVP.{POP}.vcor1
  --output_dir    Base directory for fine-mapping results (required)
                  Subdirectories rsparsepro/{AFR,EUR}/ and susie/{AFR,EUR}/ are created.
  --script_dir    Directory containing run_rsparsepro.sh and run_susie.sh (required)
  --output_file   Path to write the commands file (required)
  -v, --verbose   Print all warnings (default: show first 20 only)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
##########  PATH HELPERS (pure string construction)  ##########################
###############################################################################

def sumstats_path(sumstats_dir, pop, trait):
    return os.path.join(
        sumstats_dir, pop,
        f"MVP_R4.1000G_AGR.{trait}.{pop}.GIA.dbGaP.sumstats.hg19.txt.gz"
    )


def ld_matrix_path(ld_dir, trait, pop, chrom, start, end):
    stem = f"chr{chrom}_{start}_{end}"
    return os.path.join(ld_dir, f"{trait}.MVP_{pop}", f"{stem}.MVP.{pop}.vcor1")


def rsp_output_path(output_dir, pop, trait, chrom, start, end):
    stem = f"chr{chrom}_{start}_{end}"
    return os.path.join(output_dir, "rsparsepro", pop, f"{trait}.{stem}.rsparsepro.txt")


def susie_output_path(output_dir, pop, trait, chrom, start, end):
    stem = f"chr{chrom}_{start}_{end}"
    return os.path.join(output_dir, "susie", pop, f"{trait}.{stem}.susie.txt")


###############################################################################
##########  PHASE 1: BUILD WORK ITEMS (no filesystem calls)  #################
###############################################################################

def build_work_item(row, pop, sumstats_dir, ld_dir, output_dir, script_dir):
    """
    Construct all paths and command strings for one locus row.
    Pure string construction — no os.path.exists() calls.
    Returns a dict with all paths and pre-built command strings.
    """
    trait      = row['TRAIT']
    chrom      = int(row['CHR'])
    start      = int(row['START'])
    end        = int(row['END'])
    locus_size = end - start

    gz        = sumstats_path(sumstats_dir, pop, trait)
    ld        = ld_matrix_path(ld_dir, trait, pop, chrom, start, end)
    vars_file = ld + ".vars"
    rsp_out   = rsp_output_path(output_dir, pop, trait, chrom, start, end)
    susie_out = susie_output_path(output_dir, pop, trait, chrom, start, end)

    return {
        'locus_size': locus_size,
        'trait':      trait,
        'pop':        pop,
        'chrom':      chrom,
        'gz':         gz,
        'ld':         ld,
        'vars_file':  vars_file,
        'rsp_out':    rsp_out,
        'susie_out':  susie_out,
        'rsp_cmd':    (f"bash {script_dir}/run_rsparsepro.sh "
                       f"{gz} {ld} {rsp_out} {chrom} {start} {end}"),
        'susie_cmd':  (f"bash {script_dir}/run_susie.sh "
                       f"{gz} {ld} {susie_out} {chrom} {start} {end}"),
    }


###############################################################################
##########  PHASE 2: CHECK AND FILTER (all filesystem calls here)  ###########
###############################################################################

def check_work_item(item):
    """
    Check which commands are needed for this locus.
    Returns (commands_list, warnings_list, skipped_count).

    Fast path: both outputs exist → return immediately without checking inputs.
    Slow path: at least one output missing → validate inputs, then emit commands.
    """
    rsp_done   = os.path.exists(item['rsp_out'])
    susie_done = os.path.exists(item['susie_out'])

    # Fast path: both complete — skip all input validation
    if rsp_done and susie_done:
        return [], [], 2

    # At least one output missing; validate inputs before emitting commands
    if not os.path.exists(item['gz']):
        return [], [
            f"WARN: sumstats not found for {item['trait']}.{item['pop']}: {item['gz']}"
        ], 0
    if not os.path.exists(item['ld']):
        return [], [
            f"WARN: LD matrix not found for {item['trait']}.{item['pop']}: {item['ld']}"
        ], 0
    if not os.path.exists(item['vars_file']):
        return [], [
            f"WARN: .vars file not found: {item['vars_file']}"
        ], 0

    commands = []
    skipped  = 0

    if rsp_done:
        skipped += 1
    else:
        commands.append((item['locus_size'], item['rsp_cmd']))

    if susie_done:
        skipped += 1
    else:
        commands.append((item['locus_size'], item['susie_cmd']))

    return commands, [], skipped


###############################################################################
##########  MAIN  #############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "hv",
            ["help", "verbose", "master_afr=", "master_eur=", "sumstats_dir=",
             "ld_dir=", "output_dir=", "script_dir=", "output_file="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    master_afr   = None
    master_eur   = None
    sumstats_dir = None
    ld_dir       = None
    output_dir   = None
    script_dir   = None
    output_file  = None
    verbose      = False

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o == "--master_afr":
            master_afr = a
        elif o == "--master_eur":
            master_eur = a
        elif o == "--sumstats_dir":
            sumstats_dir = a.rstrip('/')
        elif o == "--ld_dir":
            ld_dir = a.rstrip('/')
        elif o == "--output_dir":
            output_dir = a.rstrip('/')
        elif o == "--script_dir":
            script_dir = a.rstrip('/')
        elif o == "--output_file":
            output_file = a
        else:
            assert False, f"Unhandled option: {o}"

    required = {
        "--master_afr": master_afr, "--master_eur": master_eur,
        "--sumstats_dir": sumstats_dir, "--ld_dir": ld_dir,
        "--output_dir": output_dir, "--script_dir": script_dir,
        "--output_file": output_file,
    }
    missing = [k for k, v in required.items() if v is None]
    if missing:
        print(f"ERROR: Missing required arguments: {', '.join(missing)}")
        sys.exit(2)

    for path, label in [
        (master_afr,   "AFR master loci file"),
        (master_eur,   "EUR master loci file"),
        (sumstats_dir, "sumstats directory"),
        (ld_dir,       "LD matrices directory"),
        (script_dir,   "script directory"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Create output directories ──────────────────────────────────────────────
    print("NOTE: Creating output directories...")
    for tool in ("rsparsepro", "susie"):
        for pop in ("AFR", "EUR"):
            os.makedirs(os.path.join(output_dir, tool, pop, "timing"), exist_ok=True)
    print(f"      {output_dir}/{{rsparsepro,susie}}/{{AFR,EUR}}/timing/ ready")

    # ── Read master locus files ────────────────────────────────────────────────
    print("NOTE: Reading master locus files...")
    afr_df = pd.read_csv(master_afr)
    eur_df = pd.read_csv(master_eur)
    afr_df['_pop'] = 'AFR'
    eur_df['_pop'] = 'EUR'
    all_loci = pd.concat([afr_df, eur_df], ignore_index=True)
    print(f"      {len(afr_df)} AFR loci + {len(eur_df)} EUR loci = {len(all_loci)} total")

    # ── Phase 1: build all work items (no filesystem calls) ───────────────────
    print("NOTE: Building work items...")
    all_items = [
        build_work_item(
            all_loci.iloc[i], all_loci.iloc[i]['_pop'],
            sumstats_dir, ld_dir, output_dir, script_dir
        )
        for i in range(len(all_loci))
    ]

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
    rsp_count   = sum(1 for c in commands_only if "run_rsparsepro.sh" in c)
    susie_count = sum(1 for c in commands_only if "run_susie.sh"       in c)
    print(f"\nNOTE: Processing summary:")
    print(f"      Total loci processed : {len(all_loci)}")
    print(f"      Already complete     : {total_skipped}")
    print(f"      Commands to run      : {len(commands_only)}")
    print(f"        RSparsePro         : {rsp_count}")
    print(f"        SuSiE              : {susie_count}")

    # ── Write commands file ────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(output_file, 'w') as fh:
        for cmd in commands_only:
            fh.write(cmd + '\n')
    print(f"\nSUCCESS: Wrote {len(commands_only)} commands to {output_file}")


if __name__ == "__main__":
    main()
