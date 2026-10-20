'''
run_rsparsepro_intersection.py

Per-locus RSparsePro fine-mapping wrapper that runs on the trait-aware
intersection LD matrix from merge_intersection_ld.py, with strict validation
and an atomic ".done" success sentinel.

Unlike run_rsparsepro.py (step 7, which runs on the UNION-merged LD from step
4 and always writes *some* output, even an empty one, on any failure), this
wrapper treats an explicitly empty LD locus (an
"<ld-file>.empty.json" marker from merge_intersection_ld.py) as the only valid
empty result; every other failure is a hard error so a broken locus cannot be
silently mistaken for one with no fine-mappable signal.

Because merge_intersection_ld.py already restricts the LD matrix to the exact
variant set validated against the meta-analysis sumstats (via --meta-sumstats),
the meta z-scores usually line up with the LD variant order directly ("fast
path": no LD rewrite needed). If they do not, only the needed rows/columns of
the LD matrix are streamed out to a subset file instead of reloading the
entire matrix.

RSparsePro itself (rsparsepro_ld.py) is invoked directly with --z/--ld/--save/
--K; its output columns RSID/Z/PIP/z_estimated/cs are renamed (RSID ->
VARIANT_ID) to match the fine-mapping result schema used by step 7.

Inputs:
  Meta sumstats: tabix-indexed meta-analysis summary statistics (rs_number, z)
  LD matrix:     {ld-file} / {ld-file}.vars from merge_intersection_ld.py
                 (or its zero-byte pair plus ".empty.json" marker)

Output columns:
  VARIANT_ID, Z, PIP, z_estimated, cs
Output files:
  {output-file}        tab-delimited result (or header-only for an empty locus)
  {output-file}.done   JSON success sentinel: {status, variants, elapsed_seconds}
  {output-dir}/timing/{output-file basename}  rsparsepro_time_seconds

Usage:
  python run_rsparsepro_intersection.py \\
      --sumstats meta.out.gz --ld-file locus.vcor1 --ld-vars locus.vcor1.vars \\
      --chrom 1 --start 100000 --end 200000 --output-file locus.rsparsepro.txt

***Run in finemap conda environment***
'''

import sys
import os
import getopt
import csv
import json
import math
import subprocess
import time
from io import StringIO

import numpy as np
import pandas as pd

# ASCII ----------------------------------------------------------------------
# 外层按 locus 并行；每个 RSparsePro 子进程固定单线程（巨型 locus 可用 --blas-threads 放宽）。
THREAD_VARIABLES = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")
REQUIRED_COLUMNS = {"VARIANT_ID", "Z", "PIP", "z_estimated", "cs"}

TABIX_BIN         = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/tabix"
FINEMAP_PYTHON    = "/flare/GeomicVar/mconery/tools/miniconda3/envs/finemap/bin/python"
RSPARSEPRO_SCRIPT = "/flare/GeomicVar/mconery/tools/RSparsePro_LD/src/rsparsepro_ld.py"


def help(exit_num=0):
    print("""Usage: python run_rsparsepro_intersection.py [OPTIONS]
Required:
  --sumstats    Path to bgzip+tabix meta-analysis sumstats (.out.gz)
  --ld-file     Path to intersection LD matrix (merge_intersection_ld.py output)
  --ld-vars     Path to intersection LD variant list ({ld-file}.vars)
  --chrom       Chromosome
  --start       Locus start position (bp)
  --end         Locus end position (bp)
  --output-file Full path for final output .txt file
Optional:
  --K               Max causal signals for RSparsePro (default: 10)
  --blas-threads    BLAS threads for the RSparsePro subprocess (default: 1;
                     use more for huge loci)
  --tabix-bin       Path to tabix (default: {tabix})
  --rsparsepro-script  Path to rsparsepro_ld.py (default: {script})
  --python-bin      Python interpreter for rsparsepro_ld.py (default: {python})
  --force           Re-run and overwrite a previously validated output
  -h, --help        Show this help message""".format(
        tabix=TABIX_BIN, script=RSPARSEPRO_SCRIPT, python=FINEMAP_PYTHON))
    sys.exit(exit_num)


###############################################################################
##########  LD / RESULT VALIDATION  ###########################################
###############################################################################

def read_ordered_ld_variants(path, allow_empty=False):
    with open(path) as fh:
        variants = [line.strip() for line in fh if line.strip()]
    if not variants and not allow_empty:
        raise ValueError(f"LD variant list is empty: {path}")
    if len(variants) != len(set(variants)):
        raise ValueError(f"Duplicate LD variants found: {path}")
    return variants


def validate_result(path, ld_variants, allow_empty=False):
    """Validate schema, row uniqueness, numerics, and LD membership."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        raise ValueError(f"RSparsePro output is missing or empty: {path}")

    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        columns = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS.difference(columns)
        if missing:
            raise ValueError(f"Missing output columns {sorted(missing)}: {path}")

        seen = set()
        for row_number, row in enumerate(reader, start=2):
            variant = row["VARIANT_ID"]
            if not variant or variant in seen:
                raise ValueError(f"Missing or duplicate variant at row {row_number}")
            if variant not in ld_variants:
                raise ValueError(f"Output variant absent from LD at row {row_number}: {variant}")
            seen.add(variant)
            for column in ("Z", "PIP", "z_estimated"):
                try:
                    value = float(row[column])
                except ValueError as error:
                    raise ValueError(f"Invalid {column} at row {row_number}: {row[column]}") from error
                if not math.isfinite(value):
                    raise ValueError(f"Non-finite {column} at row {row_number}: {value}")

    if not seen and not allow_empty:
        raise ValueError(f"RSparsePro output has no data rows: {path}")
    return len(seen)


def read_empty_ld_marker(ld_file, ld_vars):
    """Return validated empty-LD metadata, or None if this locus is not marked empty."""
    marker = ld_file + ".empty.json"
    if not os.path.exists(marker):
        return None
    with open(marker) as fh:
        metadata = json.load(fh)
    if metadata.get("status") != "empty" or metadata.get("variants") != 0:
        raise ValueError(f"Malformed empty-LD marker: {marker}")
    if os.path.getsize(ld_file) != 0 or os.path.getsize(ld_vars) != 0:
        raise ValueError("Marked empty LD must have zero-byte matrix and vars files")
    return metadata


def write_empty_result(output, timing):
    with open(output, 'w') as fh:
        fh.write("VARIANT_ID\tZ\tPIP\tz_estimated\tcs\n")
    os.makedirs(os.path.dirname(timing), exist_ok=True)
    with open(timing, 'w') as fh:
        fh.write("rsparsepro_time_seconds\n0.00\n")


###############################################################################
##########  META SUMSTATS / LD SUBSETTING  ####################################
###############################################################################

def extract_meta_sumstats(tabix_bin, sumstats, chrom, start, end):
    result = subprocess.run(
        [tabix_bin, "-h", sumstats, f"{chrom}:{start}-{end}"],
        check=True, capture_output=True, text=True,
    )
    lines = result.stdout.splitlines()
    if len(lines) < 2:
        raise ValueError(f"No meta summary statistics for {chrom}:{start}-{end}")
    content = lines[0].lstrip('#') + "\n" + "\n".join(lines[1:])
    return pd.read_csv(StringIO(content), sep='\t', low_memory=False)


def write_ld_subset(ld_file, n_variants, row_indices, output):
    """Stream only the needed LD rows/columns out in RSparsePro text format."""
    columns = np.asarray(row_indices, dtype=np.intp)
    row_to_output = {row: i for i, row in enumerate(row_indices)}
    rows_seen = 0
    selected_seen = 0
    buffered = []
    with open(ld_file) as source, open(output, 'w') as target:
        for local_row, line in enumerate(source):
            rows_seen += 1
            if local_row not in row_to_output:
                continue
            values = np.fromstring(line, dtype=np.float64, sep='\t')
            if values.size != n_variants:
                raise ValueError(f"Row {local_row} in {ld_file} has {values.size} values; expected {n_variants}")
            buffered.append(values[columns])
            selected_seen += 1
            if len(buffered) == 256:
                np.savetxt(target, np.vstack(buffered), fmt='%.6f', delimiter=' ')
                buffered.clear()
        if buffered:
            np.savetxt(target, np.vstack(buffered), fmt='%.6f', delimiter=' ')
    if rows_seen != n_variants or selected_seen != len(row_indices):
        raise ValueError(
            f"LD subset shape mismatch: rows={rows_seen}/{n_variants}, selected={selected_seen}/{len(row_indices)}"
        )


def run_rsparsepro(args, output):
    """Extract meta z-scores, hand off LD to RSparsePro, and write the final result."""
    started = time.time()
    frame = extract_meta_sumstats(args['tabix_bin'], args['sumstats'], args['chrom'], args['start'], args['end'])
    required = {"rs_number", "z"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing meta columns {sorted(missing)}: {args['sumstats']}")

    ld_vars = read_ordered_ld_variants(args['ld_vars'])
    ld_set = set(ld_vars)
    frame = frame[frame["rs_number"].astype(str).isin(ld_set)].copy()
    frame = frame.drop_duplicates(subset="rs_number", keep="first")
    frame["z"] = pd.to_numeric(frame["z"], errors="coerce")
    frame = frame.dropna(subset=["z"])
    if frame.empty:
        raise ValueError("No valid meta z-scores overlap the intersection LD")

    sumstats_set = set(frame["rs_number"].astype(str))
    ordered_vars = [v for v in ld_vars if v in sumstats_set]
    if not ordered_vars:
        raise ValueError("No ordered variants remain for RSparsePro")
    z_map = dict(zip(frame["rs_number"].astype(str), frame["z"]))

    output_prefix = os.path.splitext(str(output))[0]
    temporary_z = f"{output_prefix}.tmp.z.txt"
    temporary_ld = f"{output_prefix}.tmp.ld.txt"
    save_prefix = f"{output_prefix}.tmp"
    rs_result = f"{save_prefix}.rsparsepro.txt"
    rs_log = f"{save_prefix}.rsparsepro.log"
    timing = os.path.join(os.path.dirname(output), "timing", os.path.basename(output))
    os.makedirs(os.path.dirname(timing), exist_ok=True)

    z_frame = pd.DataFrame({"RSID": ordered_vars, "Z": [z_map[v] for v in ordered_vars]})
    z_frame.to_csv(temporary_z, sep='\t', index=False)
    direct_ld = ordered_vars == ld_vars
    rs_ld = args['ld_file']
    try:
        if direct_ld:
            print(f"LD_FAST_PATH: direct variants={len(ordered_vars)}")
        else:
            index = {v: i for i, v in enumerate(ld_vars)}
            row_indices = [index[v] for v in ordered_vars]
            write_ld_subset(args['ld_file'], len(ld_vars), row_indices, temporary_ld)
            rs_ld = temporary_ld
            print(f"LD_FAST_PATH: subset input={len(ld_vars)} output={len(ordered_vars)}")

        completed = subprocess.run(
            [args['python_bin'], args['rsparsepro_script'],
             "--z", temporary_z, "--ld", rs_ld, "--save", save_prefix, "--K", str(args['K'])],
            capture_output=True, text=True,
        )
        if completed.returncode != 0:
            raise RuntimeError(f"RSparsePro failed with exit code {completed.returncode}: {completed.stderr.strip()}")

        result = pd.read_csv(rs_result, sep='\t')
        result = result.rename(columns={"RSID": "VARIANT_ID"})
        result.to_csv(output, sep='\t', index=False, na_rep="NA")
        with open(timing, 'w') as fh:
            fh.write(f"rsparsepro_time_seconds\n{time.time() - started:.2f}\n")
    finally:
        for f in (temporary_z, temporary_ld, rs_result, rs_log):
            try:
                os.remove(f)
            except FileNotFoundError:
                pass


###############################################################################
##########  MAIN  #############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "sumstats=", "ld-file=", "ld-vars=", "chrom=", "start=", "end=",
             "output-file=", "K=", "blas-threads=", "tabix-bin=", "rsparsepro-script=",
             "python-bin=", "force"]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    sumstats = ld_file = ld_vars = chrom = start = end = output_file = None
    K = 10
    blas_threads = 1
    tabix_bin = TABIX_BIN
    rsparsepro_script = RSPARSEPRO_SCRIPT
    python_bin = FINEMAP_PYTHON
    force = False

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--sumstats":
            sumstats = a
        elif o == "--ld-file":
            ld_file = a
        elif o == "--ld-vars":
            ld_vars = a
        elif o == "--chrom":
            chrom = a
        elif o == "--start":
            start = int(a)
        elif o == "--end":
            end = int(a)
        elif o == "--output-file":
            output_file = a
        elif o == "--K":
            K = int(a)
        elif o == "--blas-threads":
            blas_threads = int(a)
        elif o == "--tabix-bin":
            tabix_bin = a
        elif o == "--rsparsepro-script":
            rsparsepro_script = a
        elif o == "--python-bin":
            python_bin = a
        elif o == "--force":
            force = True
        else:
            assert False, f"Unhandled option: {o}"

    required = {"--sumstats": sumstats, "--ld-file": ld_file, "--ld-vars": ld_vars,
                "--chrom": chrom, "--start": start, "--end": end, "--output-file": output_file}
    missing = [k for k, v in required.items() if v is None]
    if missing:
        print(f"ERROR: Missing required arguments: {', '.join(missing)}")
        sys.exit(2)
    if start > end:
        print("ERROR: --start must not exceed --end")
        sys.exit(2)
    if K < 1:
        print("ERROR: --K must be positive")
        sys.exit(2)
    if blas_threads < 1:
        print("ERROR: --blas-threads must be positive")
        sys.exit(2)
    for path, label in [(sumstats, "sumstats"), (ld_file, "LD file"), (ld_vars, "LD vars file"),
                         (tabix_bin, "tabix binary"), (rsparsepro_script, "rsparsepro script"),
                         (python_bin, "python binary")]:
        if not os.path.isfile(path):
            print(f"ERROR: Required input not found ({label}): {path}")
            sys.exit(1)

    output = os.path.abspath(output_file)
    done = f"{output}.done"
    started = time.time()

    empty_metadata = read_empty_ld_marker(ld_file, ld_vars)
    allow_empty = empty_metadata is not None
    ld_variant_set = set(read_ordered_ld_variants(ld_vars, allow_empty=allow_empty))

    if os.path.exists(output) and not force:
        n_variants = validate_result(output, ld_variant_set, allow_empty=allow_empty)
        if not os.path.exists(done):
            status = "empty" if allow_empty else "success"
            with open(done, 'w') as fh:
                json.dump({"status": status, "variants": n_variants,
                           "reason": empty_metadata.get("reason") if empty_metadata else None}, fh, sort_keys=True)
                fh.write('\n')
        print(f"SKIP: validated output already exists: {output}")
        return

    os.makedirs(os.path.dirname(output), exist_ok=True)
    token = str(os.getpid())
    temporary = os.path.join(os.path.dirname(output), f".{os.path.basename(output)}.tmp.{token}")
    try:
        os.remove(temporary)
    except FileNotFoundError:
        pass
    temporary_timing = os.path.join(os.path.dirname(output), "timing", os.path.basename(temporary))
    final_timing = os.path.join(os.path.dirname(output), "timing", os.path.basename(output))

    for name in THREAD_VARIABLES:
        os.environ[name] = str(blas_threads)

    try:
        if allow_empty:
            write_empty_result(temporary, temporary_timing)
            print(f"FINITE_QC: reason={empty_metadata.get('reason')} status=empty variants=0")
        else:
            run_rsparsepro(
                {'sumstats': sumstats, 'ld_file': ld_file, 'ld_vars': ld_vars, 'chrom': chrom,
                 'start': start, 'end': end, 'K': K, 'tabix_bin': tabix_bin,
                 'rsparsepro_script': rsparsepro_script, 'python_bin': python_bin},
                temporary,
            )
        n_variants = validate_result(temporary, ld_variant_set, allow_empty=allow_empty)
        os.replace(temporary, output)
        if os.path.exists(temporary_timing):
            os.makedirs(os.path.dirname(final_timing), exist_ok=True)
            os.replace(temporary_timing, final_timing)
        elapsed = time.time() - started
        status = "empty" if allow_empty else "success"
        done_tmp = f"{done}.tmp.{token}"
        with open(done_tmp, 'w') as fh:
            json.dump(
                {"status": status, "variants": n_variants, "elapsed_seconds": round(elapsed, 2),
                 "reason": empty_metadata.get("reason") if empty_metadata else None},
                fh, sort_keys=True,
            )
            fh.write('\n')
        os.replace(done_tmp, done)
    except Exception:
        for f in (temporary, temporary_timing):
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
        raise

    print(f"SUCCESS: output={output} status={status} variants={n_variants} elapsed_seconds={elapsed:.2f}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {type(error).__name__}: {error}", file=sys.stderr)
        sys.exit(1)
