'''
merge_intersection_ld.py

Build the trait-aware, strict-intersection, sample-size-weighted LD matrix used
by RSparsePro for one GWAS meta-analysis locus.

The merged LD matrices produced by step 4 (merge_ld_matrices.py) are a UNION
across cohorts: a variant kept by any one cohort is kept in the output, and
cells with no coverage are filled in (0, after averaging only the cohorts that
have data). That heterogeneous coverage creates an LD mismatch problem for
RSparsePro, which assumes every entry of its LD matrix is estimated from the
same set of individuals.

This script instead keeps only variants present in EVERY contributing per-
cohort LD matrix, AND in every cohort's own trait-specific summary statistics
with a finite BETA and SE (and, optionally, a minimum MAF), AND in the meta-
analysis summary statistics with a finite z-score. The merged correlation on
that common set is the sample-size-weighted average

    R_meta[i,j] = sum_k(N_k * R_k[i,j]) / sum_k(N_k)

Output variant order follows the first --ld-files entry. Small loci use an
in-memory vectorized merge; loci whose estimated peak memory would exceed
--max-memory-gb use a row-streaming merge with a float64 memmap accumulator,
so peak memory scales with the (smaller) intersection rather than the full
per-cohort inputs.

If any selected LD cell is non-finite, only the SNPs involved are dropped from
the intersection and the merge is retried (rather than failing the locus).  A
locus whose intersection ends up empty publishes an explicit zero-byte LD
matrix/variant-list pair plus a `<output>.empty.json` sidecar recording why,
so empty results are distinguishable from a job that never ran.

Inputs:
  Per-cohort LD:  plink2 --r-unphased square matrices (.vcor1 + .vcor1.vars),
                  one pair per --ld-files entry, in the same order as
                  --sample-sizes (and --sumstats-files, if given).
  Trait sumstats: tabix-indexed per-cohort summary statistics with ID/BETA/SE
                  (and AAF, if --maf-threshold > 0) columns.
  Meta sumstats:  tabix-indexed meta-analysis summary statistics with
                  rs_number/z columns (GWAMA output).

Output:
  {output}       tab-delimited float64 square matrix, OR a zero-byte file for
                 an explicitly empty locus.
  {output}.vars  one variant ID per line in matrix row/column order, OR
                 zero-byte for an empty locus.
  {output}.empty.json  present only for an explicitly empty locus; records
                 {"status": "empty", "variants": 0, "reason": ..., chrom/
                 start/end}.

Usage:
  python merge_intersection_ld.py \\
      --ld-files f1.vcor1 f2.vcor1 ... \\
      --sumstats-files s1.out.gz s2.out.gz ... \\
      --sample-sizes N1 N2 ... \\
      --meta-sumstats meta.out.gz --chrom 1 --start 100000 --end 200000 \\
      --output /path/to/locus.trait-intersection-maf005.vcor1

***Run in finemap conda environment***
'''

import sys
import os
import getopt
import json
import math
import shutil
import subprocess
import time

# ASCII ----------------------------------------------------------------------
# 外层按 locus 并行；单个 locus 内部不再让 BLAS/numpy 开多线程，避免节点过量订阅。
for _name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np

TABIX_BIN = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/tabix"

# Use the streaming/memmap merge path when the in-memory estimate exceeds this.
DEFAULT_MAX_MEMORY_GB = 8.0
DEFAULT_CHUNK_ROWS    = 256


def help(exit_num=0):
    print("""Usage: python merge_intersection_ld.py [OPTIONS]
Options:
  --ld-files        Space-separated list of per-cohort .vcor1 files (required)
  --sample-sizes    Space-separated list of sample sizes, one per LD file (required)
  --sumstats-files  Space-separated list of per-cohort trait sumstats (.out.gz),
                     same order as --ld-files (optional; enables trait-coverage filtering)
  --meta-sumstats   Meta-analysis sumstats (.out.gz) used by RSparsePro (optional;
                     enables meta-coverage filtering)
  --tabix-bin       Path to tabix (required with --sumstats-files/--meta-sumstats;
                     default: {tabix})
  --chrom           Chromosome (required with --sumstats-files/--meta-sumstats)
  --start           Locus start bp   (required with --sumstats-files/--meta-sumstats)
  --end             Locus end bp     (required with --sumstats-files/--meta-sumstats)
  --maf-threshold   Minimum per-cohort MAF; 0 disables the filter (default: 0.0)
  --output          Output .vcor1 path (required)
  --max-memory-gb   Switch to streaming merge above this estimate (default: {mem})
  --chunk-rows      Rows written per chunk (default: {chunk})
  --force           Replace an existing, possibly incomplete output pair
  -h, --help        Show this help message""".format(
        tabix=TABIX_BIN, mem=DEFAULT_MAX_MEMORY_GB, chunk=DEFAULT_CHUNK_ROWS))
    sys.exit(exit_num)


class NonFiniteIntersectionError(ValueError):
    """Raised when selected LD cells for the current intersection are non-finite."""

    def __init__(self, ld_file, bad_variants):
        self.ld_file = ld_file
        self.bad_variants = list(bad_variants)
        super().__init__(
            f"{len(self.bad_variants)} intersection SNPs have non-finite LD entries in {ld_file}"
        )


###############################################################################
##########  VARIANT LIST / SUMSTATS I/O  ######################################
###############################################################################

def read_variants(ld_file):
    """Read the ordered, duplicate-free variant list paired with an LD matrix."""
    vars_file = ld_file + ".vars"
    if not os.path.exists(vars_file):
        raise FileNotFoundError(f"Variant list not found: {vars_file}")
    with open(vars_file) as fh:
        variants = [line.split('\t', 1)[0].strip() for line in fh if line.strip()]
    if not variants:
        raise ValueError(f"Variant list is empty: {vars_file}")
    if len(set(variants)) != len(variants):
        raise ValueError(f"Duplicate variant IDs found: {vars_file}")
    return variants


def ordered_intersection(variant_lists, required_sets):
    """Return the all-input SNP intersection, in first-input order."""
    common = set(variant_lists[0])
    for variants in variant_lists[1:]:
        common.intersection_update(variants)
    for required in required_sets:
        common.intersection_update(required)
    return [v for v in variant_lists[0] if v in common]


def tabix_lines(tabix_bin, sumstats_file, chrom, start, end):
    region = f"{chrom}:{start}-{end}"
    result = subprocess.run(
        [tabix_bin, "-h", sumstats_file, region],
        check=True, capture_output=True, text=True,
    )
    return [line for line in result.stdout.splitlines() if line]


def extract_sumstats_variants(sumstats_file, tabix_bin, chrom, start, end, maf_threshold):
    """Trait SNPs with finite BETA/SE (and optional minimum MAF) in this region."""
    lines = tabix_lines(tabix_bin, sumstats_file, chrom, start, end)
    if not lines:
        raise ValueError(f"No summary-stat data returned: {sumstats_file} {chrom}:{start}-{end}")

    header = lines[0].lstrip('#').split('\t')
    required = {"ID", "BETA", "SE"}
    missing = required.difference(header)
    if missing:
        raise ValueError(f"Missing summary-stat columns {sorted(missing)}: {sumstats_file}")
    if maf_threshold > 0 and "AAF" not in header:
        raise ValueError(f"Missing AAF column for MAF filter: {sumstats_file}")
    index = {col: i for i, col in enumerate(header)}

    variants = set()
    for line in lines[1:]:
        fields = line.split('\t')
        if len(fields) != len(header):
            continue
        try:
            beta = float(fields[index["BETA"]])
            se   = float(fields[index["SE"]])
        except ValueError:
            continue
        if not math.isfinite(beta) or not math.isfinite(se) or se <= 0:
            continue
        if maf_threshold > 0:
            try:
                aaf = float(fields[index["AAF"]])
            except ValueError:
                continue
            maf = min(aaf, 1.0 - aaf)
            if not math.isfinite(maf) or maf < maf_threshold:
                continue
        variants.add(fields[index["ID"]])

    if not variants:
        raise ValueError(f"No valid trait variants in {sumstats_file} for region {chrom}:{start}-{end}")
    return variants


def extract_meta_variants(sumstats_file, tabix_bin, chrom, start, end):
    """Meta-analysis SNPs with a finite z-score in this region."""
    lines = tabix_lines(tabix_bin, sumstats_file, chrom, start, end)
    if not lines:
        raise ValueError(f"No meta summary-stat data returned: {sumstats_file} {chrom}:{start}-{end}")

    header = lines[0].lstrip('#').split('\t')
    required = {"rs_number", "z"}
    missing = required.difference(header)
    if missing:
        raise ValueError(f"Missing meta summary-stat columns {sorted(missing)}: {sumstats_file}")
    index = {col: i for i, col in enumerate(header)}

    variants = set()
    for line in lines[1:]:
        fields = line.split('\t')
        if len(fields) != len(header):
            continue
        try:
            z = float(fields[index["z"]])
        except ValueError:
            continue
        if math.isfinite(z):
            variants.add(fields[index["rs_number"]])

    if not variants:
        raise ValueError(f"No valid meta variants in {sumstats_file} for region {chrom}:{start}-{end}")
    return variants


###############################################################################
##########  MERGE CORE  ########################################################
###############################################################################

def load_square_matrix(ld_file, n_variants):
    """Load a tab-delimited .vcor1 matrix and enforce its expected square shape."""
    matrix = np.loadtxt(ld_file, dtype=np.float32)
    if matrix.ndim == 0:
        matrix = matrix.reshape(1, 1)
    elif matrix.ndim == 1:
        matrix = matrix.reshape(1, -1)
    expected = (n_variants, n_variants)
    if matrix.shape != expected:
        raise ValueError(f"Matrix shape mismatch for {ld_file}: {matrix.shape} != {expected}")
    return matrix


def require_finite_intersection(matrix, common_variants, ld_file):
    """Reject only the SNPs implicated in non-finite selected LD cells."""
    finite = np.isfinite(matrix)
    if finite.all():
        return
    bad_rows, bad_cols = np.where(~finite)
    bad_positions = sorted(set(bad_rows.tolist()) | set(bad_cols.tolist()))
    raise NonFiniteIntersectionError(ld_file, [common_variants[i] for i in bad_positions])


def estimate_in_memory_gb(variant_lists, n_common):
    """Estimate the in-memory path's peak allocation (float32 inputs + float64 accumulator)."""
    largest_input = max(len(v) for v in variant_lists)
    # float32 input matrix + float32 selected submatrix + float64 accumulator, +25% margin.
    estimated_bytes = largest_input ** 2 * 4 + n_common ** 2 * 4 + n_common ** 2 * 8
    return estimated_bytes * 1.25 / 1e9


def merge_in_memory(ld_files, variant_lists, sample_sizes, common_variants):
    """Merge by loading one full per-cohort matrix at a time."""
    n = len(common_variants)
    weighted_sum = np.zeros((n, n), dtype=np.float64)
    for ld_file, variants, n_k in zip(ld_files, variant_lists, sample_sizes):
        local_index = {v: i for i, v in enumerate(variants)}
        keep = np.fromiter((local_index[v] for v in common_variants), dtype=np.intp, count=n)
        matrix = load_square_matrix(ld_file, len(variants))
        selected = matrix[np.ix_(keep, keep)]
        del matrix
        require_finite_intersection(selected, common_variants, ld_file)
        weighted_sum += selected * n_k
        del selected
    weighted_sum /= float(sum(sample_sizes))
    return weighted_sum


def accumulate_streaming(accumulator, ld_file, variants, common_variants, sample_size):
    """Stream one per-cohort matrix row-by-row into an intersection-sized accumulator."""
    local_index = {v: i for i, v in enumerate(variants)}
    columns = np.fromiter(
        (local_index[v] for v in common_variants), dtype=np.intp, count=len(common_variants)
    )
    row_to_output = {local_index[v]: i for i, v in enumerate(common_variants)}

    rows_seen = 0
    selected_seen = 0
    bad_positions = set()
    with open(ld_file) as fh:
        for local_row, line in enumerate(fh):
            rows_seen += 1
            output_row = row_to_output.get(local_row)
            if output_row is None:
                continue
            values = np.fromstring(line, dtype=np.float32, sep='\t')
            if values.size != len(variants):
                raise ValueError(
                    f"Row {local_row} in {ld_file} has {values.size} values; expected {len(variants)}"
                )
            selected = values[columns]
            finite = np.isfinite(selected)
            if not finite.all():
                bad_positions.add(output_row)
                bad_positions.update(np.flatnonzero(~finite).tolist())
            else:
                accumulator[output_row, :] += selected * sample_size
            selected_seen += 1

    if rows_seen != len(variants):
        raise ValueError(f"Matrix row count mismatch for {ld_file}: {rows_seen} != {len(variants)}")
    if selected_seen != len(common_variants):
        raise ValueError(
            f"Only {selected_seen} intersection rows were read from {ld_file}; "
            f"expected {len(common_variants)}"
        )
    if bad_positions:
        raise NonFiniteIntersectionError(ld_file, [common_variants[i] for i in sorted(bad_positions)])


def write_text_matrix(matrix, output_tmp, denominator=1.0, chunk_rows=DEFAULT_CHUNK_ROWS):
    """Write a dense matrix in row chunks, dividing by denominator if given."""
    with open(output_tmp, 'w') as fh:
        for row_start in range(0, matrix.shape[0], chunk_rows):
            row_end = min(row_start + chunk_rows, matrix.shape[0])
            chunk = np.asarray(matrix[row_start:row_end], dtype=np.float64)
            if denominator != 1.0:
                chunk = chunk / denominator
            np.savetxt(fh, chunk, fmt='%.6f', delimiter='\t')


def merge_streaming(ld_files, variant_lists, sample_sizes, common_variants, output_tmp, chunk_rows):
    """Merge large matrices with a float64 memmap accumulator sized to the intersection."""
    n = len(common_variants)
    accumulator_path = output_tmp + ".accumulator"
    accumulator = np.memmap(accumulator_path, dtype=np.float64, mode='w+', shape=(n, n))
    accumulator[:] = 0.0
    try:
        for ld_file, variants, n_k in zip(ld_files, variant_lists, sample_sizes):
            accumulate_streaming(accumulator, ld_file, variants, common_variants, n_k)
        accumulator.flush()
        write_text_matrix(accumulator, output_tmp, denominator=float(sum(sample_sizes)), chunk_rows=chunk_rows)
    finally:
        del accumulator
        try:
            os.remove(accumulator_path)
        except FileNotFoundError:
            pass


###############################################################################
##########  VALIDATION / ATOMIC PUBLISH  ######################################
###############################################################################

def validate_outputs(matrix_file, vars_file, n_common):
    if not os.path.exists(matrix_file) or os.path.getsize(matrix_file) == 0:
        raise ValueError(f"Output matrix is missing or empty: {matrix_file}")
    with open(vars_file) as fh:
        output_variants = [line.strip() for line in fh if line.strip()]
    if len(output_variants) != n_common:
        raise ValueError(f"Output variant count mismatch: {len(output_variants)} != {n_common}")


def validate_empty_outputs(matrix_file, vars_file, marker_file):
    if not os.path.exists(matrix_file) or os.path.getsize(matrix_file) != 0:
        raise ValueError(f"Empty LD matrix must be a zero-byte file: {matrix_file}")
    if not os.path.exists(vars_file) or os.path.getsize(vars_file) != 0:
        raise ValueError(f"Empty LD variant list must be zero-byte: {vars_file}")
    with open(marker_file) as fh:
        metadata = json.load(fh)
    if metadata.get("status") != "empty" or metadata.get("variants") != 0:
        raise ValueError(f"Malformed empty-LD metadata: {marker_file}")
    return metadata


def publish_empty_outputs(matrix_tmp, vars_tmp, marker_tmp, output_file, output_vars, output_marker, reason, chrom, start, end):
    with open(matrix_tmp, 'w'):
        pass
    with open(vars_tmp, 'w'):
        pass
    with open(marker_tmp, 'w') as fh:
        json.dump(
            {"status": "empty", "variants": 0, "reason": reason, "chrom": chrom, "start": start, "end": end},
            fh, sort_keys=True,
        )
        fh.write('\n')
    validate_empty_outputs(matrix_tmp, vars_tmp, marker_tmp)
    os.replace(matrix_tmp, output_file)
    os.replace(vars_tmp, output_vars)
    os.replace(marker_tmp, output_marker)


###############################################################################
##########  MAIN  #############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "ld-files=", "sample-sizes=", "sumstats-files=", "meta-sumstats=",
             "tabix-bin=", "chrom=", "start=", "end=", "maf-threshold=", "output=",
             "max-memory-gb=", "chunk-rows=", "force"]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    ld_files_raw = sample_sizes_raw = sumstats_files_raw = None
    meta_sumstats = tabix_bin = chrom = start = end = None
    maf_threshold = 0.0
    output = None
    max_memory_gb = DEFAULT_MAX_MEMORY_GB
    chunk_rows = DEFAULT_CHUNK_ROWS
    force = False

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--ld-files":
            ld_files_raw = a
        elif o == "--sample-sizes":
            sample_sizes_raw = a
        elif o == "--sumstats-files":
            sumstats_files_raw = a
        elif o == "--meta-sumstats":
            meta_sumstats = a
        elif o == "--tabix-bin":
            tabix_bin = a
        elif o == "--chrom":
            chrom = a
        elif o == "--start":
            start = int(a)
        elif o == "--end":
            end = int(a)
        elif o == "--maf-threshold":
            maf_threshold = float(a)
        elif o == "--output":
            output = a
        elif o == "--max-memory-gb":
            max_memory_gb = float(a)
        elif o == "--chunk-rows":
            chunk_rows = int(a)
        elif o == "--force":
            force = True
        else:
            assert False, f"Unhandled option: {o}"

    if tabix_bin is None:
        tabix_bin = TABIX_BIN

    missing = [name for name, val in [("--ld-files", ld_files_raw), ("--sample-sizes", sample_sizes_raw),
                                       ("--output", output)] if val is None]
    if missing:
        print(f"ERROR: Missing required arguments: {', '.join(missing)}")
        sys.exit(2)

    ld_files = ld_files_raw.split()
    sample_sizes = [float(x) for x in sample_sizes_raw.split()]
    sumstats_files = sumstats_files_raw.split() if sumstats_files_raw else []

    if len(ld_files) != len(sample_sizes):
        print(f"ERROR: {len(ld_files)} LD files but {len(sample_sizes)} sample sizes")
        sys.exit(2)
    if not ld_files:
        print("ERROR: At least one LD file is required")
        sys.exit(2)
    if any(n <= 0 for n in sample_sizes):
        print("ERROR: All sample sizes must be positive")
        sys.exit(2)
    if not 0 <= maf_threshold < 0.5:
        print("ERROR: --maf-threshold must be in [0, 0.5)")
        sys.exit(2)
    if chunk_rows < 1:
        print("ERROR: --chunk-rows must be positive")
        sys.exit(2)

    for ld_file in ld_files:
        if not os.path.isfile(ld_file):
            print(f"ERROR: LD matrix not found: {ld_file}")
            sys.exit(1)
        if os.path.getsize(ld_file) == 0:
            print(f"ERROR: Contributing LD matrix is empty: {ld_file}")
            sys.exit(1)

    coverage_sets = []
    if sumstats_files:
        if len(sumstats_files) != len(ld_files):
            print(f"ERROR: {len(ld_files)} LD files but {len(sumstats_files)} summary-stat files")
            sys.exit(2)
        if not os.path.isfile(tabix_bin):
            print("ERROR: --tabix-bin is required and must exist")
            sys.exit(1)
        if chrom is None or start is None or end is None:
            print("ERROR: --chrom, --start, and --end are required with --sumstats-files")
            sys.exit(2)
        if start > end:
            print("ERROR: --start must not exceed --end")
            sys.exit(2)
        for sumstats_file in sumstats_files:
            if not os.path.isfile(sumstats_file):
                print(f"ERROR: Summary-stat file not found: {sumstats_file}")
                sys.exit(1)
            coverage_sets.append(
                extract_sumstats_variants(sumstats_file, tabix_bin, chrom, start, end, maf_threshold)
            )

    meta_coverage = None
    if meta_sumstats:
        if not os.path.isfile(meta_sumstats):
            print(f"ERROR: Meta summary statistics not found: {meta_sumstats}")
            sys.exit(1)
        if not os.path.isfile(tabix_bin):
            print("ERROR: --tabix-bin is required and must exist")
            sys.exit(1)
        if chrom is None or start is None or end is None:
            print("ERROR: --chrom, --start, and --end are required with --meta-sumstats")
            sys.exit(2)
        meta_coverage = extract_meta_variants(meta_sumstats, tabix_bin, chrom, start, end)

    started = time.time()
    output_vars = output + ".vars"
    output_marker = output + ".empty.json"

    existing = (os.path.exists(output), os.path.exists(output_vars), os.path.exists(output_marker))
    if any(existing) and not force:
        if existing[0] and existing[1]:
            if existing[2]:
                validate_empty_outputs(output, output_vars, output_marker)
                print(f"SKIP: validated empty output already exists: {output}")
                sys.exit(0)
            if os.path.getsize(output) > 0 and os.path.getsize(output_vars) > 0:
                print(f"SKIP: complete output already exists: {output}")
                sys.exit(0)
        print(f"ERROR: Incomplete output pair exists; use --force after inspection: {output}")
        sys.exit(1)
    if force and os.path.exists(output_marker):
        os.remove(output_marker)

    variant_lists = [read_variants(f) for f in ld_files]
    required_sets = list(coverage_sets)
    if meta_coverage is not None:
        required_sets.append(meta_coverage)
    common_variants = ordered_intersection(variant_lists, required_sets)
    empty_reason = "no_common_variants" if not common_variants else None

    out_dir = os.path.dirname(os.path.abspath(output))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    token = str(os.getpid())
    matrix_tmp = f"{output}.tmp.{token}"
    vars_tmp = f"{output_vars}.tmp.{token}"
    marker_tmp = f"{output_marker}.tmp.{token}"
    for stale in (matrix_tmp, vars_tmp, marker_tmp):
        try:
            os.remove(stale)
        except FileNotFoundError:
            pass

    input_counts = ','.join(str(len(v)) for v in variant_lists)
    print(f"INPUT: strata={len(ld_files)} variants=[{input_counts}] intersection={len(common_variants)}")
    if coverage_sets:
        coverage_counts = ','.join(str(len(v)) for v in coverage_sets)
        print(f"COVERAGE: trait_variants=[{coverage_counts}] maf_threshold={maf_threshold:.6g}")
    if meta_coverage is not None:
        print(f"META_COVERAGE: variants={len(meta_coverage)}")

    dropped_nonfinite = set()
    try:
        while empty_reason is None:
            try:
                os.remove(matrix_tmp)
            except FileNotFoundError:
                pass
            with open(vars_tmp, 'w') as fh:
                fh.write('\n'.join(common_variants) + '\n')
            estimated_gb = estimate_in_memory_gb(variant_lists, len(common_variants))
            can_copy = len(ld_files) == 1 and common_variants == variant_lists[0]
            mode = "copy" if can_copy else ("memory" if estimated_gb <= max_memory_gb else "streaming")
            print(f"MODE: {mode} estimated_peak_gb={estimated_gb:.3f} limit_gb={max_memory_gb:.3f}")
            try:
                if mode == "copy":
                    shutil.copyfile(ld_files[0], matrix_tmp)
                elif mode == "memory":
                    merged = merge_in_memory(ld_files, variant_lists, sample_sizes, common_variants)
                    write_text_matrix(merged, matrix_tmp, chunk_rows=chunk_rows)
                    del merged
                else:
                    merge_streaming(ld_files, variant_lists, sample_sizes, common_variants, matrix_tmp, chunk_rows)
                break
            except NonFiniteIntersectionError as error:
                bad = set(error.bad_variants).intersection(common_variants)
                if not bad:
                    raise
                dropped_nonfinite.update(bad)
                common_variants = [v for v in common_variants if v not in bad]
                print(
                    f"FINITE_QC: source={error.ld_file} dropped={len(bad)} "
                    f"total_dropped={len(dropped_nonfinite)} remaining={len(common_variants)}"
                )
                if not common_variants:
                    empty_reason = "nonfinite_qc_removed_all"

        if empty_reason is not None:
            publish_empty_outputs(matrix_tmp, vars_tmp, marker_tmp, output, output_vars, output_marker,
                                   empty_reason, chrom, start, end)
            print(f"FINITE_QC: reason={empty_reason} variants=0")
        else:
            validate_outputs(matrix_tmp, vars_tmp, len(common_variants))
            os.replace(matrix_tmp, output)
            os.replace(vars_tmp, output_vars)
    except Exception:
        for stale in (matrix_tmp, vars_tmp, marker_tmp):
            try:
                os.remove(stale)
            except FileNotFoundError:
                pass
        raise

    elapsed = time.time() - started
    keep_fraction = len(common_variants) / max(len(v) for v in variant_lists)
    print(
        f"SUCCESS: output={output} status={'empty' if empty_reason is not None else 'success'} "
        f"variants={len(common_variants)} keep_fraction={keep_fraction:.4f} elapsed_seconds={elapsed:.2f}"
    )


if __name__ == "__main__":
    main()
