'''
merge_ld_matrices.py

Efficient sample-size-weighted merge of multiple plink2 --r-unphased square LD
matrices (.vcor1 format) into a single combined matrix.

Algorithm (small loci, estimated peak memory <= LARGE_LOCUS_GB):
  1. Load each non-empty input .vcor1 file (tab-delimited float32 square matrix)
     and its corresponding .vcor1.vars file (one variant ID per line).
  2. Build the union of all variant lists, preserving first-occurrence order.
  3. Expand each matrix to the union variant space using numpy fancy indexing;
     cells where a cohort has no data are filled with NaN.
  4. Compute a vectorized sample-size-weighted average:
       combined[i,j] = sum(N_k * M_k[i,j] for k where M_k[i,j] is not NaN)
                       / sum(N_k for k where M_k[i,j] is not NaN)
     Cells with no coverage across any cohort are set to 0.
  5. Write combined matrix as tab-delimited float32 and combined .vcor1.vars.

Algorithm (large loci, estimated peak memory > LARGE_LOCUS_GB):
  Uses streaming row-by-row I/O and np.memmap for the weighted_sum and
  total_weight accumulators to avoid loading full matrices into RAM.  Output
  is written in row-chunks of CHUNK_ROWS to avoid materialising the full
  combined matrix.  Intermediate memmap files are removed on completion.

Input .vcor1 files that are 0 bytes are silently skipped (betanfer found no
variants passing the p-value threshold for that cohort).

If only 1 non-empty file remains after filtering, that matrix is written as-is
(no arithmetic merge needed, but output is still produced so the fine-mapper
sees a uniform set of inputs).

If 0 non-empty files remain, the script exits 0 without writing any output.

Usage:
  python merge_ld_matrices.py \\
      --ld-files f1.vcor1 f2.vcor1 ... \\
      --sample-sizes N1 N2 ... \\
      --output-prefix /path/to/output_prefix

  Output files:
    {output_prefix}.merged.vcor1
    {output_prefix}.merged.vcor1.vars

***Run with ALCF frameworks module: module load frameworks***
'''

import sys
import getopt
import os
import shutil
import numpy as np

# Estimated peak memory (GB) above which the streaming/memmap path is used.
# The estimate is conservative (assumes no variant overlap between cohorts),
# so this threshold should be set well below the per-worker memory limit.
LARGE_LOCUS_GB = 50

# Number of output rows to materialise at once when writing the merged matrix
# in the large-locus path.  500 rows x 300k cols x 4 bytes ~ 600 MB.
CHUNK_ROWS = 500


def help(exit_num=0):
    print("""Usage: python merge_ld_matrices.py [OPTIONS]
Options:
  --ld-files        Space-separated list of .vcor1 input files (required)
  --sample-sizes    Space-separated list of sample sizes, one per LD file (required)
  --output-prefix   Output path prefix; .merged.vcor1 / .merged.vcor1.vars appended (required)
  -h, --help        Show this help message""")
    sys.exit(exit_num)


###############################################################################
################  I/O HELPERS  ################################################
###############################################################################

def load_vars(vars_file):
    """Return an ordered list of variant IDs from a .vcor1.vars file."""
    with open(vars_file, 'r') as fh:
        return [line.rstrip('\n').split('\t')[0] for line in fh if line.strip()]


def load_matrix(vcor1_file):
    """Load a tab-delimited .vcor1 square matrix as float32 numpy array."""
    return np.loadtxt(vcor1_file, dtype=np.float32)


def stream_matrix_rows(vcor1_file):
    """Yield one row at a time as a float32 numpy array.

    Avoids loading the full matrix into memory; suitable for large .vcor1 files
    where np.loadtxt would build multi-GB Python object intermediates.
    """
    with open(vcor1_file, 'r') as f:
        for line in f:
            yield np.fromstring(line, dtype=np.float32, sep='\t')


###############################################################################
################  CORE MERGE LOGIC  ###########################################
###############################################################################

def build_variant_union(all_variant_lists):
    """Return ordered union of variant IDs (first-occurrence order)."""
    seen   = set()
    union  = []
    for variants in all_variant_lists:
        for v in variants:
            if v not in seen:
                seen.add(v)
                union.append(v)
    return union


def expand_to_union(matrix, local_variants, union_index):
    """Expand matrix from local_variants space to the union variant space.

    Cells where a variant pair is not present in local_variants are set to NaN.

    Args:
        matrix:          numpy float32 array, shape (n_local, n_local)
        local_variants:  list of variant IDs in matrix row/column order
        union_index:     dict mapping variant_id -> position in union space

    Returns:
        numpy float32 array, shape (n_union, n_union), with NaN for missing pairs
    """
    n_union = len(union_index)
    expanded = np.full((n_union, n_union), np.nan, dtype=np.float32)

    # Map local positions to union positions; skip variants not in union (shouldn't happen)
    local_to_union = [union_index[v] for v in local_variants if v in union_index]
    # Also track which local indices are valid (all should be, but be safe)
    valid_local = [i for i, v in enumerate(local_variants) if v in union_index]

    if not valid_local:
        return expanded

    ui = np.array(local_to_union, dtype=np.intp)
    li = np.array(valid_local,    dtype=np.intp)

    # Use fancy indexing: expanded[ui[:, None], ui[None, :]] = matrix[li[:, None], li[None, :]]
    expanded[np.ix_(ui, ui)] = matrix[np.ix_(li, li)]
    return expanded


###############################################################################
################  META LOADING AND MEMORY ESTIMATION  #########################
###############################################################################

def _load_valid_meta(ld_files, sample_sizes):
    """Load variant lists for all non-empty files; return list of (vcor1_file, variants, n).

    Applies per-file validation (missing/empty .vars) and skips 0-byte files.
    Does NOT load the matrices themselves.
    """
    valid = [(f, n) for f, n in zip(ld_files, sample_sizes) if os.path.getsize(f) > 0]
    result = []
    for vcor1_file, n in valid:
        vars_file = vcor1_file + '.vars'
        if not os.path.exists(vars_file):
            print(f"WARNING: .vars file missing for {vcor1_file}, skipping this cohort")
            continue
        variants = load_vars(vars_file)
        if not variants:
            print(f"WARNING: Empty .vars file for {vcor1_file}, skipping this cohort")
            continue
        result.append((vcor1_file, variants, n))
    return result


def estimate_peak_memory_gb(valid_meta):
    """Conservative upper bound on peak RAM for the in-memory merge path (GB).

    Accounts for: loading the largest matrix (np.loadtxt) + four n_union-square
    float32 arrays (expanded, weighted_sum, total_weight, combined).  Uses the
    sum of all variant counts as an upper bound on n_union (no overlap assumed).
    """
    n_union_upper = sum(len(v) for _, v, _ in valid_meta)
    n_max_local   = max(len(v) for _, v, _ in valid_meta)
    return (n_max_local ** 2 + 4 * n_union_upper ** 2) * 4 / 1e9


###############################################################################
################  LARGE-LOCUS STREAMING / MEMMAP MERGE  #######################
###############################################################################

def merge_large(valid_meta, matrix_tmp, vars_tmp):
    """Memory-efficient merge via streaming I/O and np.memmap accumulators.

    Writes the merged matrix to matrix_tmp and variant list to vars_tmp.
    Intermediate memmap files ({matrix_tmp}.wsum, {matrix_tmp}.wtot) are
    created alongside the output and removed when done.

    Returns union_variants.
    """
    union_variants = build_variant_union([v for _, v, _ in valid_meta])
    union_index    = {v: i for i, v in enumerate(union_variants)}
    n_union        = len(union_variants)

    # Single-cohort: copy directly, no arithmetic needed
    if len(valid_meta) == 1:
        vcor1_file, _, _ = valid_meta[0]
        shutil.copyfile(vcor1_file, matrix_tmp)
        with open(vars_tmp, 'w') as fh:
            fh.write('\n'.join(union_variants) + '\n')
        return union_variants

    # Multi-cohort: stream each input into np.memmap accumulators
    wsum_path = matrix_tmp + '.wsum'
    wtot_path = matrix_tmp + '.wtot'

    for stale in (wsum_path, wtot_path):
        try:
            os.remove(stale)
        except FileNotFoundError:
            pass

    wsum = np.memmap(wsum_path, dtype=np.float32, mode='w+', shape=(n_union, n_union))
    wtot = np.memmap(wtot_path, dtype=np.float32, mode='w+', shape=(n_union, n_union))

    try:
        for vcor1_file, variants, sample_size in valid_meta:
            col_indices = np.array([union_index[v] for v in variants], dtype=np.intp)
            n_local     = len(variants)
            w           = np.float32(sample_size)

            for local_row, row_data in enumerate(stream_matrix_rows(vcor1_file)):
                if row_data.shape[0] != n_local:
                    print(f"WARNING: {vcor1_file}: row {local_row} has "
                          f"{row_data.shape[0]} values != {n_local} expected -- "
                          f"skipping remaining rows of this cohort")
                    break
                ui = union_index[variants[local_row]]
                wsum[ui, col_indices] += row_data * w
                wtot[ui, col_indices] += w

        # Write combined matrix in chunks to avoid materialising the full array
        with open(matrix_tmp, 'w') as out:
            for r0 in range(0, n_union, CHUNK_ROWS):
                r1 = min(r0 + CHUNK_ROWS, n_union)
                with np.errstate(invalid='ignore', divide='ignore'):
                    chunk = np.where(
                        wtot[r0:r1] > 0,
                        wsum[r0:r1] / wtot[r0:r1],
                        0.0,
                    ).astype(np.float32)
                np.savetxt(out, chunk, fmt='%.6f', delimiter='\t')

        with open(vars_tmp, 'w') as fh:
            fh.write('\n'.join(union_variants) + '\n')

        return union_variants

    finally:
        del wsum, wtot
        for p in (wsum_path, wtot_path):
            try:
                os.remove(p)
            except FileNotFoundError:
                pass


###############################################################################
################  SMALL-LOCUS IN-MEMORY MERGE  ################################
###############################################################################

def merge_matrices_small(valid_meta):
    """In-memory merge for small loci. Accepts pre-loaded valid_meta.

    Returns (combined_matrix, union_variants) or (None, None) if no valid input.
    """
    if not valid_meta:
        return None, None

    # Load matrices and validate shapes
    loaded = []
    for vcor1_file, variants, n in valid_meta:
        matrix = load_matrix(vcor1_file)
        if matrix.ndim == 0 or matrix.size == 0:
            print(f"WARNING: Could not load matrix from {vcor1_file}, skipping")
            continue
        # Ensure 2-D (single-variant locus loads as scalar or 1-D)
        if matrix.ndim == 1:
            matrix = matrix.reshape(1, -1)
        if matrix.shape[0] != len(variants):
            print(f"WARNING: {vcor1_file}: {matrix.shape[0]} matrix rows != "
                  f"{len(variants)} vars-file entries -- skipping this cohort")
            continue
        if matrix.shape[0] != matrix.shape[1]:
            print(f"WARNING: {vcor1_file}: non-square matrix {matrix.shape} -- skipping")
            continue
        loaded.append((matrix, variants, n))

    if not loaded:
        return None, None

    # Single-cohort pass-through: no arithmetic needed
    if len(loaded) == 1:
        matrix, variants, _ = loaded[0]
        return matrix, variants

    # Multi-cohort: build union and compute weighted average
    union_variants = build_variant_union([v for _, v, _ in loaded])
    union_index    = {v: i for i, v in enumerate(union_variants)}
    n_union        = len(union_variants)

    weighted_sum = np.zeros((n_union, n_union), dtype=np.float32)
    total_weight = np.zeros((n_union, n_union), dtype=np.float32)

    for matrix, variants, sample_size in loaded:
        expanded = expand_to_union(matrix, variants, union_index)
        mask     = ~np.isnan(expanded)
        w        = np.where(mask, float(sample_size), 0.0).astype(np.float32)
        weighted_sum += np.where(mask, expanded * float(sample_size), 0.0).astype(np.float32)
        total_weight += w

    with np.errstate(invalid='ignore', divide='ignore'):
        combined = np.where(total_weight > 0, weighted_sum / total_weight, 0.0).astype(np.float32)
    return combined, union_variants


###############################################################################
###############  COLLECT / CHECK ARGUMENTS AND RUN  ###########################
###############################################################################

def main():

    # -- Parse arguments -------------------------------------------------------
    # getopt doesn't natively handle multi-value options; collect them as raw tokens.
    # We split on '--ld-files' and '--sample-sizes' sentinels manually after getopt.
    try:
        opts, remaining = getopt.getopt(
            sys.argv[1:],
            "h",
            ["help", "ld-files=", "sample-sizes=", "output-prefix="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    ld_files_raw     = None
    sample_sizes_raw = None
    output_prefix    = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--ld-files":
            ld_files_raw = a
        elif o == "--sample-sizes":
            sample_sizes_raw = a
        elif o == "--output-prefix":
            output_prefix = a
        else:
            assert False, f"WARNING: Unhandled option: {o}"

    # Validate required arguments
    missing = [name for name, val in [
        ("--ld-files",      ld_files_raw),
        ("--sample-sizes",  sample_sizes_raw),
        ("--output-prefix", output_prefix),
    ] if val is None]
    if missing:
        print(f"ERROR: Missing required option(s): {', '.join(missing)}")
        sys.exit(2)

    ld_files     = ld_files_raw.split()
    sample_sizes = [int(x) for x in sample_sizes_raw.split()]

    if len(ld_files) != len(sample_sizes):
        print(f"ERROR: Number of --ld-files ({len(ld_files)}) must match "
              f"number of --sample-sizes ({len(sample_sizes)})")
        sys.exit(2)

    # Validate input files exist
    for f in ld_files:
        if not os.path.exists(f):
            print(f"ERROR: LD file not found: {f}")
            sys.exit(1)

    # -- Setup output paths ----------------------------------------------------
    out_dir = os.path.dirname(output_prefix)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    matrix_out = output_prefix + ".merged.vcor1"
    vars_out   = matrix_out + ".vars"
    matrix_tmp = matrix_out + ".tmp"
    vars_tmp   = vars_out + ".tmp"

    # Remove any stale .tmp files left by a prior killed run
    for stale in (matrix_tmp, vars_tmp):
        try:
            os.remove(stale)
        except FileNotFoundError:
            pass

    # -- Load variant lists (cheap) and choose merge path ----------------------
    valid_meta = _load_valid_meta(ld_files, sample_sizes)
    if not valid_meta:
        print("NOTE: No non-empty input matrices found; no output written.")
        sys.exit(0)

    est_gb = estimate_peak_memory_gb(valid_meta)

    try:
        if est_gb > LARGE_LOCUS_GB:
            sizes_str = ' + '.join(str(len(v)) for _, v, _ in valid_meta)
            print(f"NOTE: Large locus ({sizes_str} variants, ~{est_gb:.0f} GB estimated) "
                  f"-- using streaming/memmap merge")
            union_variants = merge_large(valid_meta, matrix_tmp, vars_tmp)
        else:
            combined, union_variants = merge_matrices_small(valid_meta)
            if combined is None:
                print("NOTE: No valid input matrices; no output written.")
                sys.exit(0)
            np.savetxt(matrix_tmp, combined, fmt='%.6f', delimiter='\t')
            with open(vars_tmp, 'w') as fh:
                fh.write('\n'.join(union_variants) + '\n')

        os.rename(matrix_tmp, matrix_out)
        os.rename(vars_tmp, vars_out)

    except Exception:
        for stale in (matrix_tmp, vars_tmp):
            try:
                os.remove(stale)
            except FileNotFoundError:
                pass
        raise

    print(f"SUCCESS: Wrote {len(union_variants)}-variant merged matrix to {matrix_out}")


if __name__ == "__main__":
    main()
