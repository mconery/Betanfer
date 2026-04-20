'''
create_mvp_recompute_ld_commands.py

Generates betanfer infer commands to recompute LD matrices for MVP loci that
failed fine-mapping due to either:
  1. CHR:POS:REF:ALT .vars format (produced by restricted-mode betanfer before
     --chunk-size support was added).
  2. Truncated .vcor1 files (GPU memory exhaustion during large-locus computation).

Reads finemap_commands.txt to identify the failing (TRAIT, POP, CHR, START, END)
loci, then emits one regular-mode betanfer infer command per unique locus using
--chunk-size 15000 to handle large loci within GPU memory.

Command format (mirrors matrix_making.commands.txt):
  python -m betanfer infer --file-list MVP_{POP}.MVP.{POP}.txt
      --chrom CHR --start START --end END
      --pval-threshold 5e-8 --num-processes 12
      --output {LD_BASE}/{TRAIT}.MVP_{POP}/chr{CHR}_{START}_{END}.MVP.{POP}
      --plink2-format --chunk-size 15000

Usage:
  python create_mvp_recompute_ld_commands.py [OPTIONS]

***Run in parsl_gpu conda environment***
'''

import sys
import getopt
import os

LD_BASE        = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices"
FILE_LIST_BASE = os.path.join(LD_BASE, "file_lists")
CHUNK_SIZE     = 15000
PVAL_THRESHOLD = "5e-8"
NUM_PROCESSES  = 12


def help(exit_num=0):
    print("""Usage: python create_mvp_recompute_ld_commands.py [OPTIONS]
Options:
  --finemap_commands  Path to finemap_commands.txt produced by
                      create_finemap_commands.py (required)
  --output_file       Path to write betanfer infer commands (required)
  -h, --help          Show this help message""")
    sys.exit(exit_num)


def build_command(trait, pop, chrom, start, end):
    file_list = os.path.join(FILE_LIST_BASE, f"MVP_{pop}.MVP.{pop}.txt")
    output    = os.path.join(LD_BASE, f"{trait}.MVP_{pop}",
                             f"chr{chrom}_{start}_{end}.MVP.{pop}")
    return (
        f"python -m betanfer infer"
        f" --file-list {file_list}"
        f" --chrom {chrom} --start {start} --end {end}"
        f" --pval-threshold {PVAL_THRESHOLD}"
        f" --num-processes {NUM_PROCESSES}"
        f" --output {output}"
        f" --plink2-format --chunk-size {CHUNK_SIZE}"
    )


def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "finemap_commands=", "output_file="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    finemap_commands = None
    output_file      = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--finemap_commands":
            finemap_commands = a
        elif o == "--output_file":
            output_file = a
        else:
            assert False, f"Unhandled option: {o}"

    if not finemap_commands or not output_file:
        print("ERROR: --finemap_commands and --output_file are required")
        sys.exit(2)

    if not os.path.exists(finemap_commands):
        print(f"ERROR: finemap_commands not found: {finemap_commands}")
        sys.exit(1)

    # Extract unique (TRAIT, POP, CHR, START, END) loci from finemap_commands.txt.
    # Each locus appears twice (one RSparsePro command, one SuSiE command).
    # Command format: bash SCRIPT SUMSTATS LD OUTPUT CHR START END
    seen  = set()
    loci  = []
    with open(finemap_commands) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 8:
                continue
            out   = parts[4]            # full output path
            pop   = out.split('/')[-2]  # AFR or EUR
            stem  = out.split('/')[-1].replace('.txt', '')
            trait = stem.split('.chr')[0]
            chrom = parts[5]
            start = parts[6]
            end   = parts[7]
            key = (trait, pop, chrom, start, end)
            if key not in seen:
                seen.add(key)
                loci.append(key)

    print(f"NOTE: Found {len(loci)} unique loci to recompute")

    # Validate file lists exist
    for pop in sorted(set(pop for _, pop, _, _, _ in loci)):
        fl = os.path.join(FILE_LIST_BASE, f"MVP_{pop}.MVP.{pop}.txt")
        if not os.path.exists(fl):
            print(f"ERROR: File list not found: {fl}")
            sys.exit(1)
        print(f"NOTE: Using file list for {pop}: {fl}")

    # Write commands
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    with open(output_file, 'w') as fh:
        for trait, pop, chrom, start, end in loci:
            fh.write(build_command(trait, pop, chrom, start, end) + '\n')

    print(f"SUCCESS: Wrote {len(loci)} commands to {output_file}")


if __name__ == "__main__":
    main()
