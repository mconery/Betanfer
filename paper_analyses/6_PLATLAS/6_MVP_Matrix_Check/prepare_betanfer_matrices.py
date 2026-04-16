'''

prepare_betanfer_matrices.py

The purpose of this script is to collect the Betanfer LD matrices for the
100 sampled PLATLAS MVP loci and copy them into a single staging folder so
they can be transferred to the Frontier secure environment for comparison
against the raw plink-based MVP matrices.

Betanfer LD matrix outputs are organised per-trait:
  {LD_BASE}/{TRAIT}.{ANALYSIS_POP}/chr{CHR}_{START}_{END}.MVP.{POP_SUFFIX}.vcor1
  {LD_BASE}/{TRAIT}.{ANALYSIS_POP}/chr{CHR}_{START}_{END}.MVP.{POP_SUFFIX}.vcor1.vars

Because all (locus, cohort, cohort_pop) betanfer runs use the same MVP sumstats
file list, the matrix content is identical across trait directories for the same
coordinate.  This script takes the first existing matrix found for any trait
associated with each (CHR, START, END) coordinate.

Staging output structure (mirrors mvp_platlas_loci/ on Frontier):
  {staging_dir}/
    EUR/
      chr{CHR}_{START}_{END}.vcor1
      chr{CHR}_{START}_{END}.vcor1.vars
    AFR/
      chr{CHR}_{START}_{END}.vcor1
      chr{CHR}_{START}_{END}.vcor1.vars

Usage:
  python prepare_betanfer_matrices.py [-l <locus_csv>] [-m <master_csv>]
                                      [-b <ld_base>] [-o <staging_dir>]

***Run in Conda finemap_old environment on Aurora***

'''

import sys
import getopt
import os
import shutil
import pandas as pd

# ── Defaults ──────────────────────────────────────────────────────────────────
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

DEFAULT_LOCUS_CSV  = os.path.join(_SCRIPT_DIR, "sampled_loci_750kb.csv")
DEFAULT_MASTER_CSV = ("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/"
                      "mapping_loci/sig_loci.MVP_EUR.master.csv")
DEFAULT_LD_BASE    = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/LD_matrices"
DEFAULT_STAGING    = os.path.join(DEFAULT_LD_BASE, "mvp_platlas_betanfer")

# analysis_pop → pop_suffix used in the betanfer per-trait file name
MVP_POPS = {
    "MVP_EUR": "EUR",
    "MVP_AFR": "AFR",
}


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
-l => Sampled locus CSV with columns CHR,START,END
      (default: sampled_loci_750kb.csv next to this script)
-m => MVP master locus CSV used to look up TRAIT associations
      (default: sig_loci.MVP_EUR.master.csv)
-b => LD matrices base directory containing per-trait subdirs
      (default: /flare/.../LD_matrices)
-o => Staging output directory for the transfer folder
      (default: {ld_base}/mvp_platlas_betanfer)
-h => Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "l:m:b:o:h")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help(1)

    options_dict = dict(opts)

    if '-h' in options_dict:
        help()

    locus_csv  = options_dict.get('-l', DEFAULT_LOCUS_CSV)
    master_csv = options_dict.get('-m', DEFAULT_MASTER_CSV)
    ld_base    = options_dict.get('-b', DEFAULT_LD_BASE)
    staging    = options_dict.get('-o', DEFAULT_STAGING)

    for path, label in [
        (locus_csv,  "locus CSV"),
        (master_csv, "master CSV"),
        (ld_base,    "LD matrices base directory"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    print("Acceptable Inputs Given")

    driver(locus_csv, master_csv, ld_base, staging)


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(locus_csv, master_csv, ld_base, staging):

    # ── Read sampled loci ──────────────────────────────────────────────────
    print(f"NOTE: Reading sampled locus file: {locus_csv}")
    loci = pd.read_csv(locus_csv)
    print(f"NOTE: {len(loci)} sampled loci to stage")

    # ── Read master CSV to map (CHR, START, END) → list of TRAITs ─────────
    print(f"NOTE: Reading master CSV: {master_csv}")
    master = pd.read_csv(master_csv)

    coord_to_traits = {}
    for _, row in master.iterrows():
        key = (int(row['CHR']), int(row['START']), int(row['END']))
        coord_to_traits.setdefault(key, []).append(str(row['TRAIT']))

    # ── Create staging directory structure ────────────────────────────────
    os.makedirs(staging, exist_ok=True)
    for pop_suffix in MVP_POPS.values():
        os.makedirs(os.path.join(staging, pop_suffix), exist_ok=True)

    # ── Copy matrix files ─────────────────────────────────────────────────
    found   = 0
    missing = 0
    total   = len(loci) * len(MVP_POPS)

    for _, locus in loci.iterrows():
        chr_num = int(locus['CHR'])
        start   = int(locus['START'])
        end     = int(locus['END'])
        key     = (chr_num, start, end)
        locus_label = f"chr{chr_num}_{start}_{end}"

        traits = coord_to_traits.get(key, [])
        if not traits:
            print(f"WARNING: No traits in master CSV for {locus_label} — skipping")
            missing += len(MVP_POPS)
            continue

        for analysis_pop, pop_suffix in MVP_POPS.items():
            # Search trait directories for the first existing matrix pair
            src_vcor1 = None
            src_vars  = None
            for trait in traits:
                trait_dir   = os.path.join(ld_base, f"{trait}.{analysis_pop}")
                file_prefix = f"chr{chr_num}_{start}_{end}.MVP.{pop_suffix}"
                cand_vcor1  = os.path.join(trait_dir, file_prefix + ".vcor1")
                cand_vars   = os.path.join(trait_dir, file_prefix + ".vcor1.vars")
                if os.path.exists(cand_vcor1) and os.path.exists(cand_vars):
                    src_vcor1 = cand_vcor1
                    src_vars  = cand_vars
                    break

            if src_vcor1 is None:
                print(f"WARNING: No betanfer matrix found for "
                      f"{locus_label}.{analysis_pop} "
                      f"(searched {len(traits)} trait dir(s))")
                missing += 1
                continue

            dst_dir   = os.path.join(staging, pop_suffix)
            dst_vcor1 = os.path.join(dst_dir, f"{locus_label}.vcor1")
            dst_vars  = os.path.join(dst_dir, f"{locus_label}.vcor1.vars")

            shutil.copy2(src_vcor1, dst_vcor1)
            shutil.copy2(src_vars,  dst_vars)
            print(f"  Staged: {locus_label}.{pop_suffix}  "
                  f"(source: {os.path.basename(os.path.dirname(src_vcor1))})")
            found += 1

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\nSUMMARY: {found}/{total} matrix pairs staged to {staging}")
    if missing > 0:
        print(f"WARNING: {missing} matrix pair(s) not found — "
              f"check betanfer pipeline completion for those loci")
    else:
        print("SUCCESS: All matrix pairs staged successfully")


###############################################################################
################################ CALL MAIN ###################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv)
