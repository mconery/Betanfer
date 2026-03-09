#!/usr/bin/env python3
"""
find_missing_results.py

Identify missing fine-mapping results files by comparing what is recorded
in summary collation files against the complete expected set.

For each analysis directory, each method, and each LD matrix combination,
we expect results for all 1000 loci x 2 trait types x 3 causal variant counts.
Any combination absent from the summary file is reported as a missing replicate
(this catches both jobs that never ran AND jobs whose output files were unparseable).

Usage:
    python find_missing_results.py \
        --results-dir /path/to/fine-mapping_results \
        --summary-dir /path/to/summary_results \
        [--output missing_files.txt]
"""

import argparse
import csv
import multiprocessing
import os
import re
import sys

# ── Constants ─────────────────────────────────────────────────────────────────

RESULTS_BASE = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/fine-mapping_results"
SUMMARY_BASE = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro/SAIGE_simulations/summary_results"

# Top-level directories in fine-mapping_results that are NOT analysis directories
SKIP_DIRS = {"cleanup_jobs_carma", "command_files", "flux_output"}

# All method subdirectory names (with _results suffix)
METHOD_DIRS = [
    "carma_results",
    "mesusie_results",
    "multisusie_results",
    "rsparsepro_results",
    "susie_results",
    "susiex_results",
]

# sub-sample analysis directories have non-obvious summary file names
SPECIAL_SUMMARY_MAP = {
    "fixed_sub-sample_LD_3e-4":  "sub-sample_LD.fixed_homo_3e-4.results_file.txt",
    "random_sub-sample_LD_3e-4": "sub-sample_LD.random_homo_3e-4.results_file.txt",
}

# Expected simulation parameters
PHENO_NUMBERS = range(1, 1001)
CAUSAL_COUNTS = [1, 2, 3]
# (summary_trait_value, filename_trait_string)
TRAIT_SPECS = [
    ("binary",     "binary_0.25"),
    ("continuous", "continuous"),
]

# Regex to extract the filename suffix from an existing Pheno*.txt file
# Captures: group(1) = trait portion, group(2) = suffix (between trait string and .txt)
FNAME_RE = re.compile(
    r'^Pheno\d+\.\d+_causal\.(binary_[\d.]+|continuous)\.(.+)\.txt$'
)

# Summary TSV column indices (0-based)
COL_PHENO  = 0
COL_TRAIT  = 1
COL_CAUSAL = 2
COL_METHOD = 8
COL_LD     = 9


# ── Directory helpers ──────────────────────────────────────────────────────────

def build_summary_filename(dir_name):
    """Map an analysis directory name to its corresponding summary filename."""
    if dir_name in SPECIAL_SUMMARY_MAP:
        return SPECIAL_SUMMARY_MAP[dir_name]
    return f"{dir_name}.results_file.txt"


def get_analysis_dirs(results_dir):
    """
    Return list of (dir_name, abs_path) for all analysis directories in results_dir.
    Skips known non-analysis directories. Validates that a directory contains at
    least one method_results subdirectory before including it.
    """
    dirs = []
    try:
        entries = sorted(os.scandir(results_dir), key=lambda e: e.name)
    except PermissionError as exc:
        print(f"[ERROR] Cannot scan {results_dir}: {exc}", file=sys.stderr)
        return dirs

    for entry in entries:
        if not entry.is_dir(follow_symlinks=False):
            continue
        if entry.name in SKIP_DIRS:
            continue
        # Validate: must contain at least one method_results subdir
        has_method = any(
            os.path.isdir(os.path.join(entry.path, m)) for m in METHOD_DIRS
        )
        if has_method:
            dirs.append((entry.name, entry.path))
    return dirs


def get_method_ld_combos(analysis_dir_path):
    """
    Return list of (method_name, ld_dir_name, ld_dir_abs_path) for all
    method x ld_matrix combinations within an analysis directory.
    Skips non-directory entries (e.g., mesusie_commands.txt).
    """
    combos = []
    for method_full in METHOD_DIRS:
        method_name = method_full.replace("_results", "")
        method_dir = os.path.join(analysis_dir_path, method_full)
        if not os.path.isdir(method_dir):
            continue
        try:
            entries = sorted(os.scandir(method_dir), key=lambda e: e.name)
        except PermissionError:
            continue
        for entry in entries:
            # Only real directories; follow_symlinks=False avoids circular symlinks
            if entry.is_dir(follow_symlinks=False):
                combos.append((method_name, entry.name, entry.path))
    return combos


def extract_suffix(ld_dir_path):
    """
    Sample the first matching Pheno*.txt file in ld_dir_path and extract the
    filename suffix (the part between the trait string and '.txt').
    Returns None if no matching file is found (empty directory).
    """
    try:
        for entry in os.scandir(ld_dir_path):
            if not entry.is_file(follow_symlinks=False):
                continue
            m = FNAME_RE.match(entry.name)
            if m:
                return m.group(2)
    except PermissionError:
        pass
    return None


def _infer_suffix_fallback(method, ld_matrix_dir_name):
    """
    Fallback suffix inference for empty ld_matrix directories.
    Applies method-specific and dir-name-based rules.
    """
    if method == "multisusie":
        return f"{ld_matrix_dir_name}.multisusie"
    elif method in ("mesusie", "susiex"):
        return ld_matrix_dir_name
    else:
        # susie, carma, rsparsepro
        if '.' in ld_matrix_dir_name:
            # e.g., "beta-inferred_100_fixed_all.Meta" → "Meta"
            return ld_matrix_dir_name.rsplit('.', 1)[-1]
        else:
            # e.g., "PC-corrected", "raw_in-sample", "sub-sample_1000" → "weighted"
            return "weighted"


# ── Summary parsing ────────────────────────────────────────────────────────────

def parse_summary(args):
    """
    Parse a summary TSV file and return (analysis_dir_name, present_set).

    present_set contains tuples of:
        (method, ld_matrix_dir_name, phenotype_number, trait_type, causal_count)

    Accepts a single tuple argument for multiprocessing.Pool compatibility.
    """
    summary_path, analysis_dir_name = args
    present = set()

    if not os.path.isfile(summary_path):
        print(f"[WARNING] Summary file not found: {summary_path}", file=sys.stderr)
        return analysis_dir_name, present

    print(f"[INFO] Parsing {os.path.basename(summary_path)} ...", file=sys.stderr)
    try:
        with open(summary_path, newline='') as fh:
            reader = csv.reader(fh, delimiter='\t')
            header = next(reader)  # skip header row
            # Verify expected column positions
            expected = {
                COL_PHENO:  'phenotype_number',
                COL_TRAIT:  'trait_type',
                COL_CAUSAL: 'causal_variant_count',
                COL_METHOD: 'method',
                COL_LD:     'ld_matrix',
            }
            for idx, name in expected.items():
                if idx >= len(header) or header[idx] != name:
                    # Fall back to dynamic index lookup
                    col_map = {col: i for i, col in enumerate(header)}
                    COL_PHENO_dyn  = col_map.get('phenotype_number', COL_PHENO)
                    COL_TRAIT_dyn  = col_map.get('trait_type', COL_TRAIT)
                    COL_CAUSAL_dyn = col_map.get('causal_variant_count', COL_CAUSAL)
                    COL_METHOD_dyn = col_map.get('method', COL_METHOD)
                    COL_LD_dyn     = col_map.get('ld_matrix', COL_LD)
                    for row in reader:
                        try:
                            present.add((
                                row[COL_METHOD_dyn],
                                row[COL_LD_dyn],
                                int(row[COL_PHENO_dyn]),
                                row[COL_TRAIT_dyn],
                                int(row[COL_CAUSAL_dyn]),
                            ))
                        except (IndexError, ValueError):
                            pass
                    return analysis_dir_name, present

            # Fast path: known column indices
            for row in reader:
                try:
                    present.add((
                        row[COL_METHOD],
                        row[COL_LD],
                        int(row[COL_PHENO]),
                        row[COL_TRAIT],
                        int(row[COL_CAUSAL]),
                    ))
                except (IndexError, ValueError):
                    pass
    except (OSError, IOError) as exc:
        print(f"[ERROR] Could not read {summary_path}: {exc}", file=sys.stderr)

    return analysis_dir_name, present


# ── Missing file detection ─────────────────────────────────────────────────────

def find_missing_for_analysis(analysis_dir_path, present_set):
    """
    Return list of expected file paths whose (method, ld_matrix, pheno, trait, causal)
    combination is absent from present_set.
    """
    missing = []
    combos = get_method_ld_combos(analysis_dir_path)

    for method, ld_matrix_dir, ld_dir_path in combos:
        suffix = extract_suffix(ld_dir_path)
        if suffix is None:
            suffix = _infer_suffix_fallback(method, ld_matrix_dir)
            print(
                f"[WARNING] Empty directory, using fallback suffix '{suffix}': {ld_dir_path}",
                file=sys.stderr,
            )

        for pheno in PHENO_NUMBERS:
            for causal in CAUSAL_COUNTS:
                for trait_summary, trait_filename in TRAIT_SPECS:
                    key = (method, ld_matrix_dir, pheno, trait_summary, causal)
                    if key not in present_set:
                        fname = (
                            f"Pheno{pheno}.{causal}_causal"
                            f".{trait_filename}.{suffix}.txt"
                        )
                        missing.append(os.path.join(ld_dir_path, fname))

    return missing


# ── Argument parsing ───────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Find missing fine-mapping results files by comparing summary "
            "collation files against the complete expected set."
        )
    )
    parser.add_argument(
        "--results-dir",
        default=RESULTS_BASE,
        help="Path to fine-mapping_results/ directory (default: %(default)s)",
    )
    parser.add_argument(
        "--summary-dir",
        default=SUMMARY_BASE,
        help="Path to summary_results/ directory (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output file for missing file paths (default: stdout)",
    )
    return parser.parse_args()


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    results_dir = args.results_dir
    summary_dir = args.summary_dir

    print("=" * 60, file=sys.stderr)
    print("FINDING MISSING FINE-MAPPING RESULTS", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"Results directory : {results_dir}", file=sys.stderr)
    print(f"Summary directory : {summary_dir}", file=sys.stderr)
    print(file=sys.stderr)

    # Phase 1: Enumerate analysis directories
    print("[Phase 1] Enumerating analysis directories...", file=sys.stderr)
    analysis_dirs = get_analysis_dirs(results_dir)
    if not analysis_dirs:
        print("[ERROR] No analysis directories found.", file=sys.stderr)
        sys.exit(1)
    print(f"  Found {len(analysis_dirs)} analysis directories.", file=sys.stderr)
    print(file=sys.stderr)

    # Phase 2: Parse summary files in parallel
    print("[Phase 2] Parsing summary files in parallel...", file=sys.stderr)
    parse_args_list = []
    for dir_name, _dir_path in analysis_dirs:
        summary_filename = build_summary_filename(dir_name)
        summary_path = os.path.join(summary_dir, summary_filename)
        parse_args_list.append((summary_path, dir_name))

    num_workers = min(len(parse_args_list), multiprocessing.cpu_count())
    with multiprocessing.Pool(processes=num_workers) as pool:
        results = pool.map(parse_summary, parse_args_list)

    # Build lookup: analysis_dir_name → present_set
    present_sets = {dir_name: pset for dir_name, pset in results}
    total_present = sum(len(s) for s in present_sets.values())
    print(f"  Parsed {len(results)} summary files; {total_present:,} total entries found.", file=sys.stderr)
    print(file=sys.stderr)

    # Phase 3: Find missing entries and construct paths
    print("[Phase 3] Finding missing replicates...", file=sys.stderr)
    all_missing = []
    stats = {}

    for dir_name, dir_path in sorted(analysis_dirs, key=lambda x: x[0]):
        present_set = present_sets.get(dir_name, set())
        missing = find_missing_for_analysis(dir_path, present_set)
        stats[dir_name] = len(missing)
        all_missing.extend(missing)
        print(f"  {dir_name}: {len(missing):,} missing", file=sys.stderr)

    # Output missing paths
    output_content = "\n".join(all_missing)
    if args.output:
        with open(args.output, 'w') as fh:
            fh.write(output_content)
            if output_content:
                fh.write("\n")
        print(f"\n[INFO] Missing file paths written to: {args.output}", file=sys.stderr)
    else:
        if all_missing:
            print(output_content)

    # Summary statistics
    print(file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print("SUMMARY", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"{'Analysis Directory':<42} {'Missing':>10}", file=sys.stderr)
    print("-" * 54, file=sys.stderr)
    for dir_name, count in sorted(stats.items()):
        print(f"{dir_name:<42} {count:>10,}", file=sys.stderr)
    print("-" * 54, file=sys.stderr)
    print(f"{'TOTAL':<42} {sum(stats.values()):>10,}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)


if __name__ == "__main__":
    main()
