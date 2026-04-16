'''

compare_matrices.py

The purpose of this script is to compare Betanfer LD matrices (inferred from
summary statistics on Aurora) against raw plink-based MVP LD matrices (computed
from genotype data on Frontier) for 100 sampled PLATLAS loci.

For each (locus, population) pair:
  1. Load the betanfer .vcor1 matrix and .vcor1.vars file (rsIDs, one per line)
  2. Load the raw MVP .unphased.vcor1 matrix and .unphased.vcor1.vars file
     (chr:pos:a1:a2, one per line)
  3. Build a position-based lookup from the plink .pvar files to map
     rsID -> (CHR, POS, REF, ALT) for the betanfer variants
  4. Match betanfer and MVP variants by genomic position
  5. Apply sign correction: if the plink A1 allele differs from the betanfer
     ALT allele (i.e. the effect allele is coded oppositely), multiply the
     corresponding row and column of the MVP submatrix by -1
  6. Write all variant pairs (upper triangle, i < j) from the shared set to
     a single combined gzip TSV output file

Output columns:
  locus        - chr{CHR}_{START}_{END}
  population   - EUR or AFR
  variant_1    - betanfer rsID for the first variant
  variant_2    - betanfer rsID for the second variant
  betanfer_r   - correlation from betanfer matrix
  mvp_r        - correlation from raw MVP matrix (sign-corrected)
  delta_r      - betanfer_r - mvp_r

NOTE: For 750 kb loci with thousands of common variants, the number of output
rows can be very large (O(n^2/2) per locus).  Ensure sufficient disk space for
the output file before running.

Usage:
  python compare_matrices.py -b <betanfer_dir> -m <mvp_dir> -l <locus_csv>
                             -r <registry_file> -p <plink_dir> -o <output_file>
                             [-a <populations>]

***Run on Frontier (Python/3.9)***

'''

import sys
import getopt
import os
import csv
import gzip
import numpy as np
import pandas as pd

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_POPULATIONS = "AFR,EUR"
DEFAULT_PLINK_DIR   = "/gpfs/arx2/bif154_mde/proj-shared/RITM0309091/pgen"

_STRAND_COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def help(exit_num=1):
    print("""-----------------------------------------------------------------
ARGUMENTS
-b => Betanfer matrices base directory (with EUR/ and AFR/ subdirs) REQUIRED
-m => Raw MVP matrices base directory (with EUR/ and AFR/ subdirs) REQUIRED
-l => Sampled locus CSV with columns CHR,START,END REQUIRED
-r => Registry CSV (hg19.chunk_register.csv) REQUIRED
-p => Plink PGEN files directory
      (default: /gpfs/arx2/bif154_mde/proj-shared/RITM0309091/pgen) OPTIONAL
-o => Output gzip TSV file path REQUIRED
-a => Populations comma-separated (default: AFR,EUR) OPTIONAL
-h => Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, args = getopt.getopt(sys.argv[1:], "b:m:l:r:p:o:a:h")
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
        registry     = options_dict['-r']
        output_file  = options_dict['-o']
    except KeyError:
        print("ERROR: One or more required arguments are missing.")
        help()

    # Optional arguments
    plink_dir   = options_dict.get('-p', DEFAULT_PLINK_DIR)
    populations = [p.strip().upper()
                   for p in options_dict.get('-a', DEFAULT_POPULATIONS).split(',')]

    for path, label in [
        (betanfer_dir, "betanfer matrices directory"),
        (mvp_dir,      "raw MVP matrices directory"),
        (locus_csv,    "locus CSV"),
        (registry,     "registry file"),
        (plink_dir,    "plink directory"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    print("Acceptable Inputs Given")

    driver(betanfer_dir, mvp_dir, locus_csv, registry, plink_dir,
           output_file, populations)


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


def sign_correction_needed(mvp_a1, pvar_ref, pvar_alt):
    """
    Determine if mvp_r needs to be negated to match betanfer's allele coding.

    Betanfer uses ALT as the +1 effect allele.
    Plink --r square uses A1 (from BIM) as the +1 effect allele.

    If mvp_a1 matches pvar_alt  → same effect allele → no correction
    If mvp_a1 matches pvar_ref  → opposite effect allele → flip sign
    Returns True if sign flip is needed, False otherwise.
    """
    if alleles_match(mvp_a1, pvar_alt):
        return False
    if alleles_match(mvp_a1, pvar_ref):
        return True
    # Ambiguous allele pair (e.g. A/T or C/G strand-ambiguous) — default no flip
    return False


def load_pvar_region(pvar_path, chr_num, start, end):
    """
    Read a plink2 .pvar file line by line and return two dicts for variants
    that fall within [start, end] on chr_num:
      rsid_to_info : rsid -> (pos, ref, alt)
      pos_to_info  : pos  -> (rsid, ref, alt)

    Handles both compressed (.pvar.zst, .pvar.gz) and plain (.pvar) files.
    Lines starting with '#' are treated as header/metadata and skipped.
    """
    rsid_to_info = {}
    pos_to_info  = {}

    # Support compressed pvar files
    if pvar_path.endswith('.gz'):
        import gzip as _gz
        opener = lambda p: _gz.open(p, 'rt')
    elif pvar_path.endswith('.zst'):
        import zstandard as _zst
        def opener(p):
            ctx = _zst.ZstdDecompressor()
            return ctx.stream_reader(open(p, 'rb'), read_across_frames=True)
    else:
        opener = open

    try:
        with opener(pvar_path) as fh:
            for raw_line in fh:
                line = raw_line if isinstance(raw_line, str) else raw_line.decode('utf-8')
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 5:
                    continue
                try:
                    chrom = int(parts[0])
                    pos   = int(parts[1])
                except ValueError:
                    continue
                if chrom != chr_num:
                    continue
                if pos < start or pos > end:
                    continue
                var_id = parts[2]
                ref    = parts[3].upper()
                alt    = parts[4].split(',')[0].upper()   # first ALT only
                rsid_to_info[var_id] = (pos, ref, alt)
                # pos_to_info keeps the first entry at each position
                if pos not in pos_to_info:
                    pos_to_info[pos] = (var_id, ref, alt)
    except FileNotFoundError:
        print(f"WARNING: pvar file not found: {pvar_path}")

    return rsid_to_info, pos_to_info


def get_pvar_path(plink_dir, chunk_name):
    """
    Return the path to the .pvar file for a given chunk.
    Tries plain, .gz, and .zst extensions in order.
    chunk_name already includes '.dose' suffix.
    """
    base = os.path.join(plink_dir, chunk_name + ".pvar")
    for candidate in [base, base + ".gz", base + ".zst"]:
        if os.path.exists(candidate):
            return candidate
    return base   # return plain path even if missing (caller handles FileNotFoundError)


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


def find_shared_variants(betanfer_vars, mvp_vars, rsid_to_info):
    """
    Match betanfer variants (rsIDs) to MVP variants (chr:pos:a1:a2) by position.

    Returns a list of tuples:
      (b_idx, m_idx, rsid, pvar_ref, pvar_alt, mvp_a1, mvp_a2)
    for every variant that appears in both matrices.
    """
    # Build MVP pos → (index, a1, a2) lookup
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
        info = rsid_to_info.get(rsid)
        if info is None:
            # rsID not in pvar for this region — skip
            continue
        pos, pvar_ref, pvar_alt = info
        if pos not in mvp_pos_to_idx:
            continue
        m_idx, mvp_a1, mvp_a2 = mvp_pos_to_idx[pos]
        shared.append((b_idx, m_idx, rsid, pvar_ref, pvar_alt, mvp_a1, mvp_a2))

    return shared


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(betanfer_dir, mvp_dir, locus_csv, registry, plink_dir,
           output_file, populations):

    betanfer_dir = add_slash(betanfer_dir)
    mvp_dir      = add_slash(mvp_dir)
    plink_dir    = add_slash(plink_dir)

    # ── Read loci and registry ─────────────────────────────────────────────
    print(f"NOTE: Reading locus file: {locus_csv}")
    loci = pd.read_csv(locus_csv)
    print(f"NOTE: {len(loci)} loci to process")

    print(f"NOTE: Reading registry: {registry}")
    reg = pd.read_csv(registry, sep=',')
    reg.columns = reg.columns.str.strip()
    if 'pgen' in reg.columns:
        reg['Chunk'] = reg['pgen'].astype(str) + '.dose'
    elif 'Chunk' in reg.columns:
        reg['Chunk'] = reg['Chunk'].astype(str) + '.dose'

    # ── Open output file ───────────────────────────────────────────────────
    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # ── pvar cache: chunk_name → (rsid_to_info, pos_to_info) ──────────────
    pvar_cache = {}

    def get_pvar_for_chunk(chunk_name):
        if chunk_name not in pvar_cache:
            pvar_path = get_pvar_path(plink_dir.rstrip('/'), chunk_name)
            # We load lazily; the actual region filter happens inside load_pvar_region
            pvar_cache[chunk_name] = pvar_path
        return pvar_cache[chunk_name]

    total_pairs  = 0
    skipped_loci = 0

    with gzip.open(output_file, 'wt', newline='') as gz_fh:
        writer = csv.writer(gz_fh, delimiter='\t')
        writer.writerow(['locus', 'population', 'variant_1', 'variant_2',
                         'betanfer_r', 'mvp_r', 'delta_r'])

        for _, locus in loci.iterrows():
            chr_num = int(locus['CHR'])
            start   = int(locus['START'])
            end     = int(locus['END'])
            locus_label = f"chr{chr_num}_{start}_{end}"

            # Find plink chunks covering this locus
            chr_chunks = reg[reg['Chr'] == chr_num]
            if chr_chunks.empty:
                print(f"WARNING: No plink chunks for chr{chr_num} — skipping {locus_label}")
                skipped_loci += 1
                continue

            lower_cands = chr_chunks[chr_chunks['Start'] <= start]
            lower_row   = lower_cands.iloc[-1] if not lower_cands.empty else chr_chunks.iloc[0]
            upper_cands = chr_chunks[chr_chunks['End'] >= end]
            upper_row   = upper_cands.iloc[0] if not upper_cands.empty else chr_chunks.iloc[-1]

            lower_idx = lower_row.name
            upper_idx = upper_row.name
            chunk_names = [reg.loc[i, 'Chunk']
                           for i in range(lower_idx, upper_idx + 1)]

            # Build merged rsid/pos lookup from all relevant pvar chunks
            rsid_to_info = {}
            for chunk_name in chunk_names:
                pvar_path = get_pvar_for_chunk(chunk_name)
                r2i, _ = load_pvar_region(pvar_path, chr_num, start, end)
                rsid_to_info.update(r2i)

            if not rsid_to_info:
                print(f"WARNING: No variants found in pvar for {locus_label} — skipping")
                skipped_loci += 1
                continue

            for pop in populations:
                b_vcor1 = betanfer_dir + pop + f"/{locus_label}.vcor1"
                b_vars  = betanfer_dir + pop + f"/{locus_label}.vcor1.vars"
                m_vcor1 = mvp_dir      + pop + f"/{locus_label}.unphased.vcor1"
                m_vars  = mvp_dir      + pop + f"/{locus_label}.unphased.vcor1.vars"

                # Check that both matrix sets exist
                for path, label in [(b_vcor1, "betanfer vcor1"),
                                    (b_vars,  "betanfer vars"),
                                    (m_vcor1, "MVP vcor1"),
                                    (m_vars,  "MVP vars")]:
                    if not os.path.exists(path):
                        print(f"WARNING: Missing {label} for {locus_label}.{pop} — skipping")
                        break
                else:
                    # All files present — proceed
                    b_matrix, b_variants = load_betanfer_matrix(b_vcor1, b_vars)
                    m_matrix, m_variants = load_mvp_matrix(m_vcor1, m_vars)

                    shared = find_shared_variants(b_variants, m_variants, rsid_to_info)

                    if len(shared) < 2:
                        print(f"NOTE: Fewer than 2 shared variants for "
                              f"{locus_label}.{pop} — skipping")
                        continue

                    print(f"NOTE: {locus_label}.{pop}: "
                          f"{len(b_variants)} betanfer vars, "
                          f"{len(m_variants)} MVP vars, "
                          f"{len(shared)} shared")

                    # Extract submatrices for shared variants
                    b_idx_list = [s[0] for s in shared]
                    m_idx_list = [s[1] for s in shared]

                    sub_b = b_matrix[np.ix_(b_idx_list, b_idx_list)]
                    sub_m = m_matrix[np.ix_(m_idx_list, m_idx_list)].copy()

                    # Apply per-variant sign corrections to sub_m
                    for k, (_, _, _, pvar_ref, pvar_alt, mvp_a1, _) in enumerate(shared):
                        if sign_correction_needed(mvp_a1, pvar_ref, pvar_alt):
                            sub_m[k, :] *= -1
                            sub_m[:, k] *= -1

                    # Compute difference matrix
                    delta = sub_b - sub_m

                    # Write upper-triangle rows (i < j)
                    rsids = [s[2] for s in shared]
                    n = len(shared)
                    i_idx, j_idx = np.triu_indices(n, k=1)

                    for row_k in range(len(i_idx)):
                        i = i_idx[row_k]
                        j = j_idx[row_k]
                        writer.writerow([
                            locus_label,
                            pop,
                            rsids[i],
                            rsids[j],
                            f"{sub_b[i, j]:.6g}",
                            f"{sub_m[i, j]:.6g}",
                            f"{delta[i, j]:.6g}",
                        ])
                        total_pairs += 1

    print(f"\nSUMMARY: {total_pairs} variant pairs written to {output_file}")
    if skipped_loci > 0:
        print(f"WARNING: {skipped_loci} locus/population pair(s) skipped — "
              f"check warnings above")
    else:
        print("SUCCESS: All loci processed without skipping")


###############################################################################
################################ CALL MAIN ###################################
###############################################################################

if __name__ == "__main__":
    main(sys.argv)
