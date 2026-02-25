'''
run_single_gwama.py

The purpose of this script is to run a single GWAMA fixed-effects meta-analysis
for one trait/population combination (e.g. "Phe_008.AFR").

Steps performed:
  1. Read the linker file to identify input summary statistics files.
  2. Check whether output already exists (skip if so).
  3. For each input file: decompress with zcat and reformat column headers via
     awk, writing to a temporary plain-text file on flare.  The original bgzip
     files and their .tbi indices are never modified.
  4. Write a GWAMA input file (.in) listing the reformatted temp files.
  5. Run GWAMA with the --quantitative flag.
  6. Always clean up the per-job temp directory on exit (success or failure).

Column renames applied so GWAMA recognises the headers:
  #ID    -> MARKER     (GWAMA standard variant-ID column)
  ALT    -> EA         (effect allele)
  REF    -> NEA        (non-effect allele)
  AAF    -> EAF        (effect allele frequency)
  EFFECT -> EFFECT_DIR (avoid clash with GWAMA's own EFFECT keyword; BETA used)
BETA, SE, and N are recognised natively by GWAMA and are left unchanged.

***Run in Conda finemap_old environment***
'''

import sys
import getopt
import os
import subprocess
import shutil

# ── Default paths ──────────────────────────────────────────────────────────────
SUMSTATS_BASE = "/flare/GeomicVar/multi-biobank-meta/sumstats"
OUTPUT_BASE   = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/meta-analysis_results"
LINKER_FILE   = "/flare/GeomicVar/multi-biobank-meta/meta_sumstats.linker.betanfer.tsv"
TEMP_BASE     = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/gwama_temp"


def help(exit_num=0):
    print("""Usage: python run_single_gwama.py [OPTIONS]
Options:
  --trait_pop     Trait.Population label (e.g. Phe_008.AFR) (required)
  --linker_file   Path to linker TSV (default: standard path; optional)
  --sumstats_base Base directory for sumstats files (default: standard path; optional)
  --output_dir    Base output directory (default: standard path; optional)
  --temp_dir      Base directory for temporary decompressed files
                  (default: standard path; original files are never modified; optional)
  -h, --help      Show this help message""")
    sys.exit(exit_num)


###############################################################################
#################  HELPER: DECOMPRESS AND REFORMAT ONE FILE  ##################
###############################################################################

def reformat_sumstats(gz_file, output_file):
    """Stream-decompress gz_file through awk column renaming into output_file.

    Uses a zcat→awk pipeline to avoid loading the entire uncompressed file into
    Python memory.  awk renames only the header columns that GWAMA does not
    recognise natively; data rows pass through unchanged.
    """
    # awk programme: rename specific header fields; pass data rows as-is
    awk_prog = (
        'BEGIN{OFS="\\t"} '
        'NR==1{'
        'for(i=1;i<=NF;i++){'
        'if($i=="#ID")       $i="MARKER";'
        'else if($i=="ALT")  $i="EA";'
        'else if($i=="REF")  $i="NEA";'
        'else if($i=="AAF")  $i="EAF";'
        'else if($i=="EFFECT") $i="EFFECT_DIR"'
        '} print'
        '} '
        'NR>1{print}'
    )

    with open(output_file, 'w') as out_fh:
        zcat = subprocess.Popen(
            ["zcat", gz_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        awk = subprocess.Popen(
            ["awk", awk_prog],
            stdin=zcat.stdout,
            stdout=out_fh,
            stderr=subprocess.PIPE
        )
        # Allow zcat to receive SIGPIPE if awk exits early
        zcat.stdout.close()
        _, awk_stderr = awk.communicate()
        zcat_rc = zcat.wait()
        awk_rc  = awk.returncode

    if zcat_rc != 0:
        raise RuntimeError(f"zcat failed for {gz_file} (exit code {zcat_rc})")
    if awk_rc != 0:
        raise RuntimeError(
            f"awk column-rename failed for {gz_file} "
            f"(exit code {awk_rc}): {awk_stderr.decode().strip()}"
        )


###############################################################################
#################  COLLECT / CHECK ARGUMENTS AND RUN META-ANALYSIS  ###########
###############################################################################

def main():

    # ── Parse arguments ────────────────────────────────────────────────────────
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "trait_pop=", "linker_file=", "sumstats_base=",
             "output_dir=", "temp_dir="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    trait_pop     = None
    linker_file   = LINKER_FILE
    sumstats_base = SUMSTATS_BASE
    output_dir    = OUTPUT_BASE
    temp_dir      = TEMP_BASE

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--trait_pop":
            trait_pop = a
        elif o == "--linker_file":
            linker_file = a
        elif o == "--sumstats_base":
            sumstats_base = a
        elif o == "--output_dir":
            output_dir = a
        elif o == "--temp_dir":
            temp_dir = a
        else:
            assert False, f"WARNING: Unhandled input option: {o}"

    # ── Validate required arguments ────────────────────────────────────────────
    if trait_pop is None:
        print("ERROR: --trait_pop is required")
        sys.exit(2)

    # trait_pop format: TraitName.POP  (split on the last '.')
    parts = trait_pop.rsplit('.', 1)
    if len(parts) != 2:
        print(f"ERROR: --trait_pop must be in the form TraitName.POP, got: {trait_pop}")
        sys.exit(2)
    trait_name, pop = parts

    # ── Path validation ────────────────────────────────────────────────────────
    for path, label in [
        (linker_file,   "linker file"),
        (sumstats_base, "sumstats base directory"),
        (output_dir,    "output directory"),
        (temp_dir,      "temp base directory"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Skip if output already exists ─────────────────────────────────────────
    pop_output_dir = os.path.join(output_dir, pop)
    output_prefix  = os.path.join(pop_output_dir, trait_name)
    if os.path.exists(output_prefix + ".out"):
        print(f"NOTE: Output already exists for {trait_pop}, skipping.")
        sys.exit(0)

    # ── Read linker file to get input files for this trait_pop ────────────────
    input_gz_files = []
    with open(linker_file, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) >= 2 and cols[0] == trait_pop:
                input_gz_files.append(os.path.join(sumstats_base, cols[1]))

    if not input_gz_files:
        print(f"ERROR: No input files found in linker for trait_pop: {trait_pop}")
        sys.exit(1)

    # Confirm every input file exists before we start writing temp files
    for gz in input_gz_files:
        if not os.path.exists(gz):
            print(f"ERROR: Input file not found: {gz}")
            sys.exit(1)

    # ── Create unique per-job temp directory ──────────────────────────────────
    # Replace '.' with '_' to keep the path a valid directory name
    job_temp_dir = os.path.join(temp_dir, trait_pop.replace('.', '_'))
    os.makedirs(job_temp_dir, exist_ok=True)

    try:
        # ── Decompress + reformat each input file ─────────────────────────────
        temp_txt_files = []
        for i, gz_file in enumerate(input_gz_files):
            temp_txt = os.path.join(job_temp_dir, f"input_{i}.txt")
            print(f"NOTE: Decompressing and reformatting [{i+1}/{len(input_gz_files)}]: {gz_file}")
            reformat_sumstats(gz_file, temp_txt)
            temp_txt_files.append(temp_txt)

        # ── Write GWAMA input (.in) file ──────────────────────────────────────
        gwama_in_file = os.path.join(job_temp_dir, "gwama_input.in")
        with open(gwama_in_file, 'w') as fh:
            for txt in temp_txt_files:
                fh.write(txt + '\n')

        # ── Ensure population output directory exists ─────────────────────────
        os.makedirs(pop_output_dir, exist_ok=True)

        # ── Run GWAMA ─────────────────────────────────────────────────────────
        # --name_marker MARKER is required because GWAMA v2.2.2's internal
        # default for the marker column is not "MARKER".  BETA, SE, EA, NEA,
        # EAF, and N are all recognised natively and need no --name_* override.
        gwama_cmd = [
            "GWAMA",
            "--filelist",    gwama_in_file,
            "--quantitative",
            "--name_marker", "MARKER",
            "--output",      output_prefix,
        ]
        print(f"NOTE: Running GWAMA for {trait_pop} "
              f"({len(input_gz_files)} input file(s))")
        result = subprocess.run(
            gwama_cmd,
            capture_output=True,
            text=True,
            cwd=job_temp_dir,
        )

        # Echo GWAMA output for log capture
        if result.stdout.strip():
            print(result.stdout)
        if result.stderr.strip():
            print(result.stderr, file=sys.stderr)

        if result.returncode != 0:
            print(f"ERROR: GWAMA failed for {trait_pop} "
                  f"(exit code {result.returncode})")
            # Remove any partial output files GWAMA may have created before
            # failing, so the next run does not falsely skip this job.
            for ext in (".out", ".log.out"):
                partial = output_prefix + ext
                if os.path.exists(partial):
                    os.remove(partial)
                    print(f"NOTE: Removed partial output file: {partial}")
            sys.exit(1)

        print(f"SUCCESS: GWAMA meta-analysis complete for {trait_pop}")

    finally:
        # Always remove the per-job temp directory; original files are untouched
        shutil.rmtree(job_temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
