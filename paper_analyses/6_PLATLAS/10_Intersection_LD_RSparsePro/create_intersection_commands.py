'''
create_intersection_commands.py

Plan the trait-aware intersection LD merge and RSparsePro run for every
(TRAIT, POP, LOCUS) combination in the master locus files, using the PLATLAS
linker to find each locus's contributing cohort/population strata.

Unlike the union merge in 4_LD_Matrix_Merging (which fills missing cells with
NaN and averages only over cohorts with coverage), this step restricts every
locus to the SNPs present in ALL contributing per-cohort LD matrices AND in
every stratum's trait sumstats with valid BETA/SE AND in the meta-analysis
sumstats with a finite z, then takes the sample-size-weighted LD average on
that strict intersection. merge_intersection_ld.py and
run_rsparsepro_intersection.py (both in this directory) carry out the merge
and fine-mapping; this script only plans which loci are ready and emits one
command per locus for each.

Inputs:
  Linker:        {linker} (two tab-separated columns: "{trait}.{pop}" and
                 "{cohort}/{cohort_pop}/{sumstats filename}")
  Master loci:   {loci-dir}/sig_loci.{EUR,AFR,ALL}.master.csv
                 Columns: CHR, START, END, LENGTH, TRAIT, OVERLAP_COUNT
  Sample sizes:  {sample-sizes} (TSV: TRAIT, COHORT, COHORT_POP, N), produced
                 by 4_LD_Matrix_Merging/extract_sample_sizes.py
  Per-cohort LD: {ld-dir}/{trait}.{pop}/{locus}.{cohort}.{cohort_pop}.vcor1[.vars]
  Meta-analysis: {meta-dir}/{pop}/{trait}.out.gz (bgzip+tabix indexed)

Output file naming:
  Manifest:            {manifest}
  Merge commands:       {merge-commands}
  RSparsePro commands:  {rsparsepro-commands}
  Merge output:          {output-root}/ld/{trait}.{pop}/{locus}.trait-intersection-maf005.vcor1[.vars]
  RSparsePro output:     {output-root}/rsparsepro/{pop}/{trait}.{locus}.rsparsepro.txt

Two-phase design (as in 7_GWAS_Fine-Mapping/create_rsparsepro_commands.py):
  Phase 1 — pure string construction from the linker, loci, and sample sizes;
            zero filesystem calls.
  Phase 2 — parallel input checks (ThreadPoolExecutor) that mark each locus
            "ready" or "missing_input" and record its total LD input size.

Commands for "ready" loci are sorted by total LD input size descending
(largest first) so that the slowest jobs start first and do not become a
tail that keeps idle workers waiting at the end of the run.

Usage:
  python create_intersection_commands.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import os
import getopt
import csv
import shlex
import multiprocessing
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor


SCRIPT_DIR        = os.path.dirname(os.path.abspath(__file__))
NUM_WORKERS        = multiprocessing.cpu_count()
FINEMAP_PYTHON     = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/python"
TABIX_BIN          = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/tabix"
LINKER_FILE        = "/flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv"
SUMSTATS_DIR       = "/flare/GeomicVar/multi-biobank-meta/sumstats"
LOCUS_DIR          = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mapping_loci"
META_DIR           = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results"
LD_DIR             = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices"
SAMPLE_SIZES_FILE  = os.path.join(LD_DIR, "sample_sizes.tsv")
OUTPUT_ROOT        = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/intersection_finemap"
MERGER_SCRIPT      = os.path.join(SCRIPT_DIR, "merge_intersection_ld.py")
RUNNER_SCRIPT      = os.path.join(SCRIPT_DIR, "run_rsparsepro_intersection.py")
ANALYSIS_POPS      = ["EUR", "AFR", "ALL"]
MAF_THRESHOLD      = 0.005
MAX_MEMORY_GB      = 8.0

# Resource-class byte thresholds, informational only (manifest column; no
# per-class command files are written — see module docstring).
HUGE_BYTES   = 50_000_000_000
LARGE_BYTES  = 10_000_000_000
MEDIUM_BYTES = 2_000_000_000


def help(exit_num=0):
    print("""Usage: python create_intersection_commands.py [OPTIONS]
Options:
  --linker              Path to PLATLAS linker file (default: {linker})
  --loci-dir            Directory with sig_loci.{{EUR,AFR,ALL}}.master.csv (default: {loci})
  --ld-dir              Base directory for per-cohort LD matrices (default: {ld})
  --sample-sizes        TSV of TRAIT/COHORT/COHORT_POP/N (default: {ss})
  --sumstats-dir        Base directory for per-cohort sumstats (default: {sd})
  --meta-dir            Base directory for meta-analysis sumstats (default: {md})
  --output-root         Base directory for merge + RSparsePro outputs (default: {out})
  --manifest            Path to write the production manifest TSV (required)
  --merge-commands      Path to write the merge command list (required)
  --rsparsepro-commands Path to write the RSparsePro command list (required)
  --python-bin          Python interpreter for both steps (default: {py})
  --tabix-bin           Path to tabix (default: {tabix})
  --maf-threshold       Minimum AAF kept on both sides (default: {maf})
  --max-memory-gb       In-memory vs streaming merge threshold (default: {mem})
  --stat-workers        Threads for parallel input checks (default: cpu_count)
  -v, --verbose         Print all warnings (default: show first 20 only)
  -h, --help            Show this help message""".format(
        linker=LINKER_FILE, loci=LOCUS_DIR, ld=LD_DIR, ss=SAMPLE_SIZES_FILE,
        sd=SUMSTATS_DIR, md=META_DIR, out=OUTPUT_ROOT, py=FINEMAP_PYTHON,
        tabix=TABIX_BIN, maf=MAF_THRESHOLD, mem=MAX_MEMORY_GB))
    sys.exit(exit_num)


###############################################################################
##########  PATH HELPERS (pure string construction)  ##########################
###############################################################################

def locus_name(chrom, start, end):
    return f"chr{chrom}_{start}_{end}"


def ld_file_path(ld_dir, trait, pop, locus, cohort, cohort_pop):
    return os.path.join(ld_dir, f"{trait}.{pop}", f"{locus}.{cohort}.{cohort_pop}.vcor1")


def merge_output_path(output_root, trait, pop, locus):
    return os.path.join(output_root, "ld", f"{trait}.{pop}", f"{locus}.trait-intersection-maf005.vcor1")


def rsparsepro_output_path(output_root, pop, trait, locus):
    return os.path.join(output_root, "rsparsepro", pop, f"{trait}.{locus}.rsparsepro.txt")


def resource_class(ld_input_bytes):
    if ld_input_bytes > HUGE_BYTES:
        return "huge"
    if ld_input_bytes > LARGE_BYTES:
        return "large"
    if ld_input_bytes > MEDIUM_BYTES:
        return "medium"
    return "standard"


###############################################################################
##########  LINKER / SAMPLE SIZE PARSING  ######################################
###############################################################################

def read_linker(linker_file, sumstats_dir, populations):
    """Read linker strata in source order, deduplicated per (trait, pop)."""
    strata = defaultdict(list)
    seen = defaultdict(set)
    with open(linker_file) as fh:
        for line_number, line in enumerate(fh, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Malformed linker row {line_number}: {linker_file}")
            trait_pop, relative_sumstats = fields
            trait, pop = trait_pop.rsplit(".", 1)
            if pop not in populations:
                continue
            parts = relative_sumstats.split(os.sep)
            if len(parts) < 3:
                raise ValueError(f"Malformed linker path at row {line_number}: {relative_sumstats}")
            cohort, cohort_pop = parts[0], parts[1]
            key = (trait, pop)
            identity = (cohort, cohort_pop)
            if identity in seen[key]:
                continue
            seen[key].add(identity)
            strata[key].append({
                "cohort": cohort,
                "cohort_pop": cohort_pop,
                "sumstats": os.path.join(sumstats_dir, relative_sumstats),
            })
    return dict(strata)


def read_sample_sizes(path):
    """Read positive per-trait cohort sample sizes."""
    sample_sizes = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["TRAIT"], row["COHORT"], row["COHORT_POP"])
            value = int(row["N"])
            if value <= 0:
                raise ValueError(f"Non-positive sample size for {key}: {value}")
            sample_sizes[key] = value
    return sample_sizes


###############################################################################
##########  PHASE 1: BUILD WORK ITEMS (no filesystem calls)  #################
###############################################################################

def build_work_item(row, pop, linker, sample_sizes, args):
    trait = row["TRAIT"]
    chrom = row["CHR"]
    start = int(row["START"])
    end = int(row["END"])
    locus = locus_name(chrom, start, end)

    trait_strata = linker.get((trait, pop), [])
    missing = [] if trait_strata else ["linker"]
    ld_files = []
    n_values = []
    for stratum in trait_strata:
        ld_files.append(ld_file_path(args['ld_dir'], trait, pop, locus, stratum["cohort"], stratum["cohort_pop"]))
        key = (trait, stratum["cohort"], stratum["cohort_pop"])
        if key in sample_sizes:
            n_values.append(sample_sizes[key])
        else:
            missing.append(f"sample_size:{stratum['cohort']}.{stratum['cohort_pop']}")

    return {
        "pop": pop,
        "trait": trait,
        "chrom": chrom,
        "start": start,
        "end": end,
        "locus": locus,
        "strata": trait_strata,
        "ld_files": ld_files,
        "sample_sizes": n_values,
        "meta_sumstats": os.path.join(args['meta_dir'], pop, f"{trait}.out.gz"),
        "ld_output": merge_output_path(args['output_root'], trait, pop, locus),
        "rsparsepro_output": rsparsepro_output_path(args['output_root'], pop, trait, locus),
        "missing": missing,
        "ld_input_bytes": 0,
    }


###############################################################################
##########  PHASE 2: CHECK INPUTS IN PARALLEL (all filesystem calls here)  ###
###############################################################################

def inspect_inputs(items, workers):
    """Check every unique required input once and annotate missing/size."""
    required = set()
    for item in items:
        for ld_file in item["ld_files"]:
            required.add(ld_file)
            required.add(ld_file + ".vars")
        for stratum in item["strata"]:
            required.add(stratum["sumstats"])
            required.add(stratum["sumstats"] + ".tbi")
        required.add(item["meta_sumstats"])
        required.add(item["meta_sumstats"] + ".tbi")

    def stat(path):
        try:
            return path, os.path.getsize(path)
        except OSError:
            return path, -1

    with ThreadPoolExecutor(max_workers=workers) as pool:
        sizes = dict(pool.map(stat, required))

    for item in items:
        checks = []
        for ld_file in item["ld_files"]:
            checks.append(("ld", ld_file))
            checks.append(("ld_vars", ld_file + ".vars"))
        for stratum in item["strata"]:
            checks.append(("sumstats", stratum["sumstats"]))
            checks.append(("sumstats_tbi", stratum["sumstats"] + ".tbi"))
        checks.append(("meta_sumstats", item["meta_sumstats"]))
        checks.append(("meta_sumstats_tbi", item["meta_sumstats"] + ".tbi"))
        for label, path in checks:
            if sizes.get(path, -1) <= 0:
                item["missing"].append(f"{label}:{path}")
        item["ld_input_bytes"] = sum(max(0, sizes.get(path, -1)) for path in item["ld_files"])


###############################################################################
##########  COMMAND BUILDERS  ##################################################
###############################################################################

def merge_command(item, args):
    cmd = [
        str(args['python_bin']), MERGER_SCRIPT,
        "--ld-files", " ".join(str(p) for p in item["ld_files"]),
        "--sample-sizes", " ".join(str(v) for v in item["sample_sizes"]),
        "--tabix-bin", str(args['tabix_bin']),
        "--chrom", str(item["chrom"]), "--start", str(item["start"]), "--end", str(item["end"]),
        "--maf-threshold", str(args['maf_threshold']),
        "--max-memory-gb", str(args['max_memory_gb']),
        "--output", str(item["ld_output"]),
    ]
    if item["strata"]:
        cmd += ["--sumstats-files", " ".join(str(s["sumstats"]) for s in item["strata"])]
    cmd += ["--meta-sumstats", str(item["meta_sumstats"])]
    return shlex.join(cmd)


def rsparsepro_command(item, args):
    cmd = [
        str(args['python_bin']), RUNNER_SCRIPT,
        "--sumstats", str(item["meta_sumstats"]),
        "--ld-file", str(item["ld_output"]),
        "--ld-vars", f"{item['ld_output']}.vars",
        "--chrom", str(item["chrom"]), "--start", str(item["start"]), "--end", str(item["end"]),
        "--tabix-bin", str(args['tabix_bin']),
        "--python-bin", str(args['python_bin']),
        "--output-file", str(item["rsparsepro_output"]),
    ]
    return shlex.join(cmd)


###############################################################################
##########  MANIFEST / COMMAND FILE WRITERS  ##################################
###############################################################################

def write_manifest(items, path):
    fields = ["pop", "trait", "locus", "chrom", "start", "end", "n_strata", "strata",
              "ld_input_bytes", "resource_class", "status", "missing", "ld_output", "rsparsepro_output"]
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for item in items:
            writer.writerow({
                "pop": item["pop"], "trait": item["trait"], "locus": item["locus"],
                "chrom": item["chrom"], "start": item["start"], "end": item["end"],
                "n_strata": len(item["strata"]),
                "strata": ",".join(f"{s['cohort']}.{s['cohort_pop']}" for s in item["strata"]),
                "ld_input_bytes": item["ld_input_bytes"],
                "resource_class": resource_class(item["ld_input_bytes"]),
                "status": "ready" if not item["missing"] else "missing_input",
                "missing": ";".join(item["missing"]),
                "ld_output": item["ld_output"], "rsparsepro_output": item["rsparsepro_output"],
            })


def write_commands(commands, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        for cmd in commands:
            fh.write(cmd + "\n")


###############################################################################
##########  MAIN  #############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "hv",
            ["help", "verbose", "linker=", "loci-dir=", "ld-dir=", "sample-sizes=",
             "sumstats-dir=", "meta-dir=", "output-root=", "manifest=", "merge-commands=",
             "rsparsepro-commands=", "python-bin=", "tabix-bin=", "maf-threshold=",
             "max-memory-gb=", "stat-workers="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    linker = LINKER_FILE
    loci_dir = LOCUS_DIR
    ld_dir = LD_DIR
    sample_sizes_file = SAMPLE_SIZES_FILE
    sumstats_dir = SUMSTATS_DIR
    meta_dir = META_DIR
    output_root = OUTPUT_ROOT
    manifest = merge_commands = rsparsepro_commands = None
    python_bin = FINEMAP_PYTHON
    tabix_bin = TABIX_BIN
    maf_threshold = MAF_THRESHOLD
    max_memory_gb = MAX_MEMORY_GB
    stat_workers = NUM_WORKERS
    verbose = False

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o == "--linker":
            linker = a
        elif o == "--loci-dir":
            loci_dir = a.rstrip('/')
        elif o == "--ld-dir":
            ld_dir = a.rstrip('/')
        elif o == "--sample-sizes":
            sample_sizes_file = a
        elif o == "--sumstats-dir":
            sumstats_dir = a.rstrip('/')
        elif o == "--meta-dir":
            meta_dir = a.rstrip('/')
        elif o == "--output-root":
            output_root = a.rstrip('/')
        elif o == "--manifest":
            manifest = a
        elif o == "--merge-commands":
            merge_commands = a
        elif o == "--rsparsepro-commands":
            rsparsepro_commands = a
        elif o == "--python-bin":
            python_bin = a
        elif o == "--tabix-bin":
            tabix_bin = a
        elif o == "--maf-threshold":
            maf_threshold = float(a)
        elif o == "--max-memory-gb":
            max_memory_gb = float(a)
        elif o == "--stat-workers":
            stat_workers = int(a)
        else:
            assert False, f"Unhandled option: {o}"

    required = {"--manifest": manifest, "--merge-commands": merge_commands,
                "--rsparsepro-commands": rsparsepro_commands}
    missing_args = [k for k, v in required.items() if v is None]
    if missing_args:
        print(f"ERROR: Missing required arguments: {', '.join(missing_args)}")
        sys.exit(2)
    if not 0 <= maf_threshold < 0.5:
        print("ERROR: --maf-threshold must be in [0, 0.5)")
        sys.exit(2)
    if stat_workers < 1:
        print("ERROR: --stat-workers must be positive")
        sys.exit(2)

    for path, label in [(linker, "linker file"), (loci_dir, "loci directory"),
                         (ld_dir, "LD directory"), (sample_sizes_file, "sample sizes file"),
                         (sumstats_dir, "sumstats directory"), (meta_dir, "meta-analysis directory"),
                         (python_bin, "python binary"), (tabix_bin, "tabix binary"),
                         (MERGER_SCRIPT, "merger script"), (RUNNER_SCRIPT, "RSparsePro runner script")]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    args = {'ld_dir': ld_dir, 'meta_dir': meta_dir, 'output_root': output_root,
            'python_bin': python_bin, 'tabix_bin': tabix_bin,
            'maf_threshold': maf_threshold, 'max_memory_gb': max_memory_gb}

    linker_map = read_linker(linker, sumstats_dir, set(ANALYSIS_POPS))
    sample_sizes = read_sample_sizes(sample_sizes_file)

    all_items = []
    pop_counts = {}
    for pop in ANALYSIS_POPS:
        master_path = os.path.join(loci_dir, f"sig_loci.{pop}.master.csv")
        if not os.path.exists(master_path):
            print(f"WARNING: Master locus file not found, skipping {pop}: {master_path}")
            continue
        with open(master_path, newline="") as fh:
            rows = list(csv.DictReader(fh))
        pop_counts[pop] = len(rows)
        print(f"NOTE: {pop}: {len(rows)} loci")
        for row in rows:
            all_items.append(build_work_item(row, pop, linker_map, sample_sizes, args))

    print(f"NOTE: Built {len(all_items)} work items across all populations")
    print(f"NOTE: Checking inputs ({stat_workers} workers)...")
    inspect_inputs(all_items, stat_workers)

    ready = [item for item in all_items if not item["missing"]]
    warnings = [f"WARN: missing input(s) for {item['trait']}.{item['pop']}.{item['locus']}: "
                f"{';'.join(item['missing'])}" for item in all_items if item["missing"]]
    if warnings:
        print(f"\nWARNINGS ({len(warnings)}):")
        shown = warnings if verbose else warnings[:20]
        for w in shown:
            print(f"  {w}")
        if not verbose and len(warnings) > 20:
            print(f"  ... and {len(warnings) - 20} more (use --verbose to show all)")

    # ASCII ------------------------------------------------------------------
    # Largest LD input first (LPT-style): keeps the slowest jobs from
    # becoming a tail that leaves workers idle at the end of the run.
    ready.sort(key=lambda item: item["ld_input_bytes"], reverse=True)

    write_manifest(all_items, manifest)
    write_commands([merge_command(item, args) for item in ready], merge_commands)
    write_commands([rsparsepro_command(item, args) for item in ready], rsparsepro_commands)

    class_counts = defaultdict(int)
    for item in ready:
        class_counts[resource_class(item["ld_input_bytes"])] += 1

    print(f"\nNOTE: Processing summary:")
    print(f"      Total loci processed : {len(all_items)}")
    print(f"      Ready                : {len(ready)}")
    print(f"      Missing input        : {len(all_items) - len(ready)}")
    print("TOTAL: " + " ".join(f"{pop}={pop_counts.get(pop, 0)}" for pop in ANALYSIS_POPS)
          + f" all={len(all_items)}")
    print(f"READY: {len(ready)} missing_input={len(all_items) - len(ready)}")
    print("RESOURCE_CLASSES: " + " ".join(
        f"{name}={class_counts[name]}" for name in ("standard", "medium", "large", "huge")))
    print(f"MANIFEST: {manifest}")
    print(f"MERGE_COMMANDS: {merge_commands}")
    print(f"RSPARSEPRO_COMMANDS: {rsparsepro_commands}")


if __name__ == "__main__":
    main()
