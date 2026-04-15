'''
create_restricted_commands.py

For (locus × cohort+pop) betanfer infer commands that could not complete using
the full genomic-range mode (--chrom/--start/--end), this script produces
alternative commands that restrict the variant set to only those present in the
trait-specific meta-analysis (or MVP) summary statistics for that locus.

For each failing command (missing .vcor1.vars sentinel):
  - Queries the appropriate tabix-indexed source file to get variant IDs in the
    locus region, writing them to a per-locus variant list file.
  - Replaces --chrom/--start/--end with --variant-file <variant_list_path>,
    keeping all other betanfer infer arguments unchanged.

Variant list sources:
  EUR / AFR / ALL  →  GWAMA meta-analysis output (.out.gz), #rs_number column
  MVP_EUR          →  MVP hg19 EUR sumstats, #ID column
  MVP_AFR          →  MVP hg19 AFR sumstats, #ID column

Usage:
  python create_restricted_commands.py
      --input_commands <matrix_making.commands.txt>
      --command_file   <matrix_making.restricted.commands.txt>
      [--log_dir <path>]

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os
import re
import glob
import multiprocessing
from concurrent.futures import ThreadPoolExecutor

import pysam

# ── Default paths ────────────────────────────────────────────────────────────────
LD_BASE       = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices"
META_RESULTS  = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results"
MVP_HG19_SUMSTATS = "/flare/GeomicVar/multi-biobank-meta/sumstats/MVP_hg19"

VARIANT_LISTS_DIR = os.path.join(LD_BASE, "variant_lists")

SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))

# Analysis populations served by GWAMA meta-analysis output
LINKER_POPS = {"EUR", "AFR", "ALL"}
# Analysis populations served by MVP hg19 single-cohort sumstats
MVP_POP_MAP = {"MVP_EUR": "EUR", "MVP_AFR": "AFR"}

NUM_WORKERS = multiprocessing.cpu_count()


def help(exit_num=0):
    print("""Usage: python create_restricted_commands.py [OPTIONS]
Options:
  --input_commands  Path to the existing matrix_making.commands.txt (required)
  --command_file    Path to write the restricted commands file (required)
  --log_dir         Directory for parsl worker log output (optional)
  -h, --help        Show this help message""")
    sys.exit(exit_num)


###############################################################################
################  COMMAND PARSING HELPERS  ####################################
###############################################################################

def _extract_arg(cmd, flag, is_int=False):
    """Extract the value of a named flag from a command string."""
    m = re.search(rf'{re.escape(flag)}\s+(\S+)', cmd)
    if not m:
        return None
    val = m.group(1)
    return int(val) if is_int else val


def parse_cmd(cmd):
    """Parse a betanfer infer command string into its key components.

    Returns (file_list, chrom, start, end, output_prefix, pval_threshold, num_processes)
    or raises ValueError if required fields are missing.
    """
    file_list      = _extract_arg(cmd, '--file-list')
    chrom          = _extract_arg(cmd, '--chrom')
    start          = _extract_arg(cmd, '--start', is_int=True)
    end            = _extract_arg(cmd, '--end',   is_int=True)
    output_prefix  = _extract_arg(cmd, '--output')
    pval_threshold = _extract_arg(cmd, '--pval-threshold') or "5e-8"
    num_processes  = _extract_arg(cmd, '--num-processes')  or "8"

    missing = [name for name, val in [
        ('--file-list', file_list),
        ('--chrom',     chrom),
        ('--start',     start),
        ('--end',       end),
        ('--output',    output_prefix),
    ] if val is None]

    if missing:
        raise ValueError(f"Could not parse {missing} from command")

    return file_list, chrom, start, end, output_prefix, pval_threshold, num_processes


def parse_trait_pop(output_prefix):
    """Extract (trait, analysis_pop) from an output prefix path.

    Output prefix pattern:
      {LD_BASE}/{TRAIT}.{ANALYSIS_POP}/chr{CHR}_{START}_{END}.{COHORT}.{COHORT_POP}
    """
    rel = os.path.relpath(output_prefix, LD_BASE)
    # rel = "{TRAIT}.{ANALYSIS_POP}/chr..."
    trait_pop_dir = rel.split('/')[0]
    analysis_pop = trait_pop_dir.rsplit('.', 1)[1]
    trait        = trait_pop_dir.rsplit('.', 1)[0]
    return trait, analysis_pop


###############################################################################
################  IDENTIFY FAILING COMMANDS  ##################################
###############################################################################

def find_failing_commands(input_commands_file):
    """Read the main command file; return commands whose sentinel is missing."""
    with open(input_commands_file) as fh:
        all_cmds = [line.rstrip('\n') for line in fh if line.strip()]

    # Build (cmd, sentinel) pairs
    cmd_sentinel = []
    for cmd in all_cmds:
        m = re.search(r'--output\s+(\S+)', cmd)
        if not m:
            print(f"WARNING: Could not parse --output from command: {cmd[:80]}")
            continue
        sentinel = m.group(1) + ".vcor1.vars"
        cmd_sentinel.append((cmd, sentinel))

    print(f"NOTE: Checking {len(cmd_sentinel)} sentinel files in parallel "
          f"({NUM_WORKERS} workers)...")
    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as ex:
        sentinel_exists = list(ex.map(lambda cs: os.path.exists(cs[1]), cmd_sentinel))

    failing = [cs[0] for cs, exists in zip(cmd_sentinel, sentinel_exists) if not exists]
    done    = sum(sentinel_exists)
    print(f"NOTE: {done} already complete  |  {len(failing)} missing sentinel (to restrict)")
    return failing


###############################################################################
################  VARIANT LIST BUILDING  ######################################
###############################################################################

def get_source_gz(trait, analysis_pop):
    """Return the tabix-indexed source file path for the given trait+population."""
    if analysis_pop in LINKER_POPS:
        return os.path.join(META_RESULTS, analysis_pop, f"{trait}.out.gz")
    elif analysis_pop in MVP_POP_MAP:
        mvp_pop = MVP_POP_MAP[analysis_pop]
        mvp_dir = os.path.join(MVP_HG19_SUMSTATS, mvp_pop)
        matches = glob.glob(os.path.join(mvp_dir, f"MVP_R4.1000G_AGR.{trait}.*.sumstats.hg19.txt.gz"))
        matches = [f for f in matches if not f.endswith("hm3_sumstats.txt.gz")]
        if len(matches) == 1:
            return matches[0]
        elif len(matches) > 1:
            print(f"WARNING: Multiple MVP files matched for trait '{trait}' in {mvp_dir}; using first")
            return matches[0]
        return None
    return None


def build_variant_list(trait, analysis_pop, chrom, start, end):
    """Query the appropriate tabix source and write a per-locus variant list.

    Returns the variant list file path, or None if the source is unavailable or
    yields no variants.  Skips building if a non-empty file already exists
    (checkpoint).
    """
    vlist_dir  = os.path.join(VARIANT_LISTS_DIR, f"{trait}.{analysis_pop}")
    vlist_file = os.path.join(vlist_dir, f"chr{chrom}_{start}_{end}.variants.txt")

    # Checkpoint: already built
    if os.path.exists(vlist_file) and os.path.getsize(vlist_file) > 0:
        return vlist_file

    source_gz = get_source_gz(trait, analysis_pop)
    if source_gz is None:
        print(f"WARNING: Unknown analysis_pop '{analysis_pop}' for trait '{trait}', skipping")
        return None

    if not os.path.exists(source_gz):
        print(f"WARNING: Source file not found: {source_gz}")
        return None

    # pysam.TabixFile.fetch() uses 0-based half-open coordinates:
    #   fetch(chrom, start-1, end) captures 1-based positions [start, end]
    try:
        seen = set()
        variants = []
        with pysam.TabixFile(source_gz) as tbx:
            for line in tbx.fetch(str(chrom), start - 1, end):
                vid = line.split('\t')[0]
                if vid and vid not in seen:
                    seen.add(vid)
                    variants.append(vid)
    except Exception as exc:
        print(f"WARNING: Failed to query {source_gz} for "
              f"chr{chrom}:{start}-{end}: {exc}")
        return None

    if not variants:
        print(f"WARNING: No variants in {source_gz} for chr{chrom}:{start}-{end}")
        return None

    os.makedirs(vlist_dir, exist_ok=True)
    with open(vlist_file, 'w') as fh:
        fh.write('\n'.join(variants) + '\n')

    return vlist_file


###############################################################################
################  COMMAND GENERATION  ########################################
###############################################################################

def make_restricted_cmd(file_list, variant_file, pval_threshold, num_processes, output_prefix, analysis_pop):
    """Build a betanfer infer command using --variant-file instead of --chrom/--start/--end."""
    cmd = (
        f"python -m betanfer infer"
        f" --file-list {file_list}"
        f" --variant-file {variant_file}"
        f" --pval-threshold {pval_threshold}"
        f" --num-processes {num_processes}"
        f" --output {output_prefix}"
        f" --plink2-format"
        f" --chunk-size 15000"
    )
    if analysis_pop in MVP_POP_MAP:
        cmd += " --variant-col ID"
    return cmd


def count_variants(vlist_file):
    """Count the number of non-empty lines in a variant list file."""
    with open(vlist_file) as fh:
        return sum(1 for line in fh if line.strip())


###############################################################################
################  MAIN  #######################################################
###############################################################################

def main():

    # ── Parse arguments ──────────────────────────────────────────────────────────
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "input_commands=", "command_file=", "log_dir="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    input_commands = None
    command_file   = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--input_commands":
            input_commands = a
        elif o == "--command_file":
            command_file = a
        elif o == "--log_dir":
            pass  # accepted for CLI consistency; unused here
        else:
            assert False, f"WARNING: Unhandled option: {o}"

    if input_commands is None:
        print("ERROR: Missing required --input_commands")
        sys.exit(2)
    if command_file is None:
        print("ERROR: Missing required --command_file")
        sys.exit(2)
    if not os.path.exists(input_commands):
        print(f"ERROR: Input commands file not found: {input_commands}")
        sys.exit(1)

    os.makedirs(LD_BASE, exist_ok=True)

    # ── Step 1: Find commands whose sentinel is still missing ─────────────────────
    print("NOTE: Identifying failing commands...")
    failing = find_failing_commands(input_commands)
    if not failing:
        print("NOTE: No failing commands (all sentinels present). Nothing to do.")
        sys.exit(0)

    # ── Step 2: Parse failing commands ────────────────────────────────────────────
    print(f"NOTE: Parsing {len(failing)} failing command(s)...")
    parsed = []    # [(locus_key, file_list, output_prefix, pval_threshold, num_processes)]
    locus_vlist = {}  # locus_key → variant_list_path (built in step 3)
    parse_errors = 0

    for cmd in failing:
        try:
            file_list, chrom, start, end, output_prefix, pval_threshold, num_processes \
                = parse_cmd(cmd)
        except ValueError as e:
            print(f"WARNING: Skipping unparseable command ({e}): {cmd[:80]}")
            parse_errors += 1
            continue

        trait, analysis_pop = parse_trait_pop(output_prefix)
        locus_key = (trait, analysis_pop, chrom, start, end)
        parsed.append((locus_key, file_list, output_prefix, pval_threshold, num_processes))
        if locus_key not in locus_vlist:
            locus_vlist[locus_key] = None  # populated in step 3

    if parse_errors:
        print(f"WARNING: {parse_errors} command(s) could not be parsed and were skipped")

    unique_loci = len(locus_vlist)
    print(f"NOTE: {len(parsed)} parseable failing command(s)  |  "
          f"{unique_loci} unique (trait, pop, locus) combination(s)")

    # ── Step 3: Build variant list files for each unique locus ───────────────────
    print(f"NOTE: Building variant lists for {unique_loci} unique loci...")
    cached  = 0
    missing = 0

    for locus_key in locus_vlist:
        trait, analysis_pop, chrom, start, end = locus_key
        was_cached = (
            os.path.exists(
                os.path.join(VARIANT_LISTS_DIR, f"{trait}.{analysis_pop}",
                             f"chr{chrom}_{start}_{end}.variants.txt")
            )
        )
        vlist = build_variant_list(trait, analysis_pop, chrom, start, end)
        locus_vlist[locus_key] = vlist
        if vlist is None:
            missing += 1
        elif was_cached:
            cached += 1

    built = unique_loci - missing - cached
    print(f"NOTE: Variant lists — built: {built}  |  cached: {cached}  |  "
          f"unavailable: {missing}")

    # ── Step 4: Generate restricted commands ─────────────────────────────────────
    print("NOTE: Generating restricted commands...")
    restricted = []
    skipped    = 0

    for (locus_key, file_list, output_prefix, pval_threshold, num_processes) in parsed:
        vlist = locus_vlist[locus_key]
        if vlist is None:
            trait, analysis_pop, chrom, start, end = locus_key
            print(f"WARNING: No variant list for {trait}.{analysis_pop} "
                  f"chr{chrom}:{start}-{end}, skipping command")
            skipped += 1
            continue
        n_variants = count_variants(vlist)
        _, analysis_pop, _, _, _ = locus_key
        cmd = make_restricted_cmd(file_list, vlist, pval_threshold, num_processes, output_prefix, analysis_pop)
        restricted.append((n_variants, cmd))

    # Sort ascending by variant count (small matrices first)
    restricted.sort(key=lambda x: x[0])

    # ── Step 5: Write restricted commands file ────────────────────────────────────
    cmd_dir = os.path.dirname(command_file)
    if cmd_dir:
        os.makedirs(cmd_dir, exist_ok=True)

    with open(command_file, 'w') as fh:
        for _, cmd in restricted:
            fh.write(cmd + '\n')

    print(f"NOTE: {len(restricted)} restricted commands written to {command_file}")
    if skipped:
        print(f"WARNING: {skipped} command(s) skipped (no variant list available)")
    if restricted:
        min_v = min(n for n, _ in restricted)
        max_v = max(n for n, _ in restricted)
        print(f"NOTE: Variant counts: min={min_v}  max={max_v}")


if __name__ == "__main__":
    main()
