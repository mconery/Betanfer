'''

subset_matrices.py

For each of the 100 sampled PLATLAS loci and 2 populations (200 pairs total),
this script:
  1. Identifies variants present in both the betanfer LD matrix and the raw MVP
     LD matrix, using rsid_pos_lookup.json for rsID-to-position mapping
  2. Extracts the shared-variant submatrix from each source
  3. Applies sign correction to the MVP submatrix (ALT vs A1 coding)
  4. Writes both submatrices and the shared variant list to disk

Both output matrices use betanfer variant order among the shared variants.

Output structure:
  {output_dir}/{POP}/
    chr{CHR}_{START}_{END}.betanfer.vcor1       — tab-delimited submatrix
    chr{CHR}_{START}_{END}.betanfer.vcor1.vars  — shared rsIDs, one per line
    chr{CHR}_{START}_{END}.mvp.vcor1            — tab-delimited submatrix (sign-corrected)
    chr{CHR}_{START}_{END}.mvp.vcor1.vars       — shared rsIDs, one per line (same order)

Usage:
  python subset_matrices.py -b <betanfer_dir> -m <mvp_dir> -l <locus_csv>
                            -k <rsid_lookup_json> -o <output_dir>
                            [-a <populations>] [-n <n_workers>]

***Run on Frontier (Python/3.9)***

'''

import sys
import getopt
import os
import json
import numpy as np
import pandas as pd
from multiprocessing import Pool

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_POPULATIONS = "AFR,EUR"
DEFAULT_N_WORKERS   = 20

_STRAND_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

# Global lookup dict — populated once per worker via pool initializer
_RSID_TO_INFO = None


def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
-b => Betanfer matrices base directory (with EUR/ and AFR/ subdirs) REQUIRED
-m => Raw MVP matrices base directory (with EUR/ and AFR/ subdirs) REQUIRED
-l => Sampled locus CSV with columns CHR,START,END REQUIRED
-k => rsID position lookup JSON (output of build_rsid_lookup.py) REQUIRED
      Format: { rsid: [chr, pos, ref, alt], ... }
-o => Output base directory for subset matrices REQUIRED
-a => Populations comma-separated (default: AFR,EUR) OPTIONAL
-n => Number of parallel worker processes (default: 20) OPTIONAL
-h => Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "b:m:l:k:o:a:n:h")
    except getopt.GetoptError:
        print("ERROR: Incorrect usage of getopts flags!")
        help()

    options_dict = dict(opts)

    if '-h' in options_dict:
        help(0)

    # Required arguments
    try:
        betanfer_dir = options_dict['-b']
        mvp_dir      = options_dict['-m']
        locus_csv    = options_dict['-l']
        lookup_file  = options_dict['-k']
        output_dir   = options_dict['-o']
    except KeyError:
        print("ERROR: One or more required arguments are missing.")
        help()

    # Optional arguments
    populations = [p.strip().upper()
                   for p in options_dict.get('-a', DEFAULT_POPULATIONS).split(',')]
    n_workers   = int(options_dict.get('-n', DEFAULT_N_WORKERS))

    for path, label in [
        (betanfer_dir, "betanfer matrices directory"),
        (mvp_dir,      "raw MVP matrices directory"),
        (locus_csv,    "locus CSV"),
        (lookup_file,  "rsID lookup JSON"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    print("Acceptable Inputs Given")

    driver(betanfer_dir, mvp_dir, locus_csv, lookup_file,
           output_dir, populations, n_workers)


###############################################################################
############################# HELPFUL FUNCTIONS ##############################
###############################################################################

def add_slash(directory):
    return directory if directory.endswith('/') else directory + '/'


def strand_complement(allele):
    """Return strand complement of a single-character allele string."""
    return ''.join(_STRAND_COMP.get(c.upper(), c) for c in allele)


def alleles_match(a, b):
    """True if a == b or a == strand complement of b (case-insensitive)."""
    a_u = a.upper()
    b_u = b.upper()
    return a_u == b_u or a_u == strand_complement(b_u)


def sign_correction_needed(mvp_a1, lookup_ref, lookup_alt):
    """
    Determine if mvp_r needs to be negated to match betanfer's allele coding.

    Betanfer uses ALT as the +1 effect allele.
    Plink --r square uses A1 (from BIM) as the +1 effect allele.

    If mvp_a1 matches lookup_alt → same effect allele → no correction
    If mvp_a1 matches lookup_ref → opposite effect allele → flip sign
    Returns True if sign flip is needed, False otherwise.
    """
    if alleles_match(mvp_a1, lookup_alt):
        return False
    if alleles_match(mvp_a1, lookup_ref):
        return True
    return False


def load_rsid_lookup(lookup_file):
    """
    Read the rsID position lookup JSON (built by build_rsid_lookup.py).
    Returns a dict: rsid -> (chr, pos, ref, alt)
    """
    with open(lookup_file, 'r') as fh:
        raw = json.load(fh)
    return {rsid: (int(v[0]), int(v[1]), v[2].upper(), v[3].upper())
            for rsid, v in raw.items()}


def load_betanfer_matrix(vcor1_path, vars_path):
    """Load betanfer matrix and variant list.
    Returns (n×n numpy float64 array, list of rsID strings).
    """
    with open(vars_path, 'r') as fh:
        variants = [line.strip() for line in fh if line.strip()]
    matrix = np.loadtxt(vcor1_path, delimiter='\t', dtype=np.float64)
    return matrix, variants


def load_mvp_matrix(vcor1_path, vars_path):
    """Load raw MVP plink matrix and variant list.
    Returns (n×n numpy float64 array, list of 'chr:pos:a1:a2' strings).
    """
    with open(vars_path, 'r') as fh:
        variants = [line.strip() for line in fh if line.strip()]
    # plink1.9 --r square output is space-delimited
    matrix = np.loadtxt(vcor1_path, dtype=np.float64)
    return matrix, variants


def find_shared_variants(betanfer_vars, mvp_vars):
    """
    Match betanfer variants (rsIDs) to MVP variants (chr:pos:a1:a2) by position.
    Uses the module-level _RSID_TO_INFO dict populated by the pool initializer.

    Returns a list of tuples:
      (b_idx, m_idx, rsid, lookup_ref, lookup_alt, mvp_a1, mvp_a2)
    for every variant that appears in both matrices.
    """
    mvp_pos_to_idx = {}
    for i, var_str in enumerate(mvp_vars):
        parts = var_str.split(':')
        if len(parts) < 4:
            continue
        try:
            pos = int(parts[1])
        except ValueError:
            continue
        a1 = parts[2].upper()
        a2 = parts[3].upper()
        if pos not in mvp_pos_to_idx:
            mvp_pos_to_idx[pos] = (i, a1, a2)

    shared = []
    for b_idx, rsid in enumerate(betanfer_vars):
        info = _RSID_TO_INFO.get(rsid)
        if info is None:
            continue
        _chrom, pos, lookup_ref, lookup_alt = info
        if pos not in mvp_pos_to_idx:
            continue
        m_idx, mvp_a1, mvp_a2 = mvp_pos_to_idx[pos]
        shared.append((b_idx, m_idx, rsid, lookup_ref, lookup_alt, mvp_a1, mvp_a2))

    return shared


###############################################################################
########################### POOL INITIALIZER + WORKER ########################
###############################################################################

def _init_worker(rsid_to_info):
    """Initialise each worker process with the shared lookup dict."""
    global _RSID_TO_INFO
    _RSID_TO_INFO = rsid_to_info


def process_locus_pop(work_item):
    """
    Process one (locus, population) pair: extract shared-variant submatrices
    and write 4 output files.

    Returns (success, n_shared, status_message).
    """
    locus_label  = work_item['locus_label']
    pop          = work_item['pop']
    b_vcor1      = work_item['b_vcor1']
    b_vars       = work_item['b_vars']
    m_vcor1      = work_item['m_vcor1']
    m_vars       = work_item['m_vars']
    out_b_vcor1  = work_item['out_b_vcor1']
    out_b_vars   = work_item['out_b_vars']
    out_m_vcor1  = work_item['out_m_vcor1']
    out_m_vars   = work_item['out_m_vars']

    tag = f"{locus_label}.{pop}"

    # Check all input files exist
    for path, label in [(b_vcor1, "betanfer vcor1"),
                        (b_vars,  "betanfer vars"),
                        (m_vcor1, "MVP vcor1"),
                        (m_vars,  "MVP vars")]:
        if not os.path.exists(path):
            return (False, 0, f"WARNING: Missing {label} for {tag} — skipping")

    # Load matrices
    b_matrix, b_variants = load_betanfer_matrix(b_vcor1, b_vars)
    m_matrix, m_variants = load_mvp_matrix(m_vcor1, m_vars)

    # Find shared variants using module-level lookup
    shared = find_shared_variants(b_variants, m_variants)

    if len(shared) < 2:
        return (False, len(shared),
                f"NOTE: Fewer than 2 shared variants for {tag} — skipping")

    # Extract submatrices
    b_idx_list = [s[0] for s in shared]
    m_idx_list = [s[1] for s in shared]
    rsids      = [s[2] for s in shared]

    sub_b = b_matrix[np.ix_(b_idx_list, b_idx_list)]
    sub_m = m_matrix[np.ix_(m_idx_list, m_idx_list)].copy()

    # Apply per-variant sign corrections to sub_m
    for k, (_, _, _, lookup_ref, lookup_alt, mvp_a1, _) in enumerate(shared):
        if sign_correction_needed(mvp_a1, lookup_ref, lookup_alt):
            sub_m[k, :] *= -1
            sub_m[:, k] *= -1

    # Write submatrices (tab-delimited, 6 significant figures)
    np.savetxt(out_b_vcor1, sub_b, delimiter='\t', fmt='%.6g')
    np.savetxt(out_m_vcor1, sub_m, delimiter='\t', fmt='%.6g')

    # Write shared variant lists
    with open(out_b_vars, 'w') as fh:
        fh.write('\n'.join(rsids) + '\n')
    with open(out_m_vars, 'w') as fh:
        fh.write('\n'.join(rsids) + '\n')

    n_shared = len(shared)
    return (True, n_shared,
            f"NOTE: {tag}: {n_shared} shared variants written")


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(betanfer_dir, mvp_dir, locus_csv, lookup_file,
           output_dir, populations, n_workers):

    betanfer_dir = add_slash(betanfer_dir)
    mvp_dir      = add_slash(mvp_dir)

    # ── Read loci ──────────────────────────────────────────────────────────
    print(f"NOTE: Reading locus file: {locus_csv}")
    loci = pd.read_csv(locus_csv)
    print(f"NOTE: {len(loci)} loci to process")

    # ── Load rsID lookup ───────────────────────────────────────────────────
    print(f"NOTE: Loading rsID lookup: {lookup_file}")
    rsid_to_info = load_rsid_lookup(lookup_file)
    print(f"NOTE: {len(rsid_to_info)} rsIDs in lookup")

    if len(rsid_to_info) == 0:
        print("ERROR: rsID lookup is empty — check the lookup file")
        sys.exit(1)

    # ── Create output population subdirectories ────────────────────────────
    for pop in populations:
        os.makedirs(os.path.join(output_dir, pop), exist_ok=True)

    # ── Build work item list (one per locus × population) ─────────────────
    work_items = []
    for _, locus in loci.iterrows():
        chr_num = int(locus['CHR'])
        start   = int(locus['START'])
        end     = int(locus['END'])
        locus_label = f"chr{chr_num}_{start}_{end}"

        for pop in populations:
            pop_out = os.path.join(output_dir, pop)
            work_items.append({
                'locus_label': locus_label,
                'pop':         pop,
                'b_vcor1':     betanfer_dir + pop + f"/{locus_label}.vcor1",
                'b_vars':      betanfer_dir + pop + f"/{locus_label}.vcor1.vars",
                'm_vcor1':     mvp_dir      + pop + f"/{locus_label}.unphased.vcor1",
                'm_vars':      mvp_dir      + pop + f"/{locus_label}.unphased.vcor1.vars",
                'out_b_vcor1': os.path.join(pop_out, f"{locus_label}.betanfer.vcor1"),
                'out_b_vars':  os.path.join(pop_out, f"{locus_label}.betanfer.vcor1.vars"),
                'out_m_vcor1': os.path.join(pop_out, f"{locus_label}.mvp.vcor1"),
                'out_m_vars':  os.path.join(pop_out, f"{locus_label}.mvp.vcor1.vars"),
            })

    print(f"NOTE: {len(work_items)} (locus, population) pairs to process "
          f"using {n_workers} workers")

    # ── Run parallel pool ──────────────────────────────────────────────────
    found   = 0
    skipped = 0

    with Pool(processes=n_workers,
              initializer=_init_worker,
              initargs=(rsid_to_info,)) as pool:

        for success, n_shared, msg in pool.imap_unordered(
                process_locus_pop, work_items):
            print(msg)
            if success:
                found += 1
            else:
                skipped += 1

    print(f"\nSUMMARY: {found}/{len(work_items)} pairs written to {output_dir}")
    if skipped > 0:
        print(f"WARNING: {skipped} pair(s) skipped — check warnings above")
    else:
        print("SUCCESS: All pairs processed without skipping")


###############################################################################
################################ CALL MAIN ###################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv)
