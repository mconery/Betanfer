#!/usr/bin/env python3
"""
regenerate_missing_results.py

Given the list of missing fine-mapping result file paths produced by
find_missing_results.py, this script:
  1. Deletes each file if it already exists (to remove corrupt/incomplete copies).
  2. Reverse-maps each output path back to its input files (GWAS summary
     statistics and LD matrix) using glob-based discovery.
  3. Constructs the appropriate regeneration command for each missing file.
  4. Writes six per-method command files (one per fine-mapping method).

The six command files are intended to be executed sequentially by the
companion shell script, each in its own conda environment.

Usage:
    python regenerate_missing_results.py \\
        --missing-list /path/to/missing_results_files.txt \\
        --output-dir   /path/to/command_files_dir \\
        [--populations AFR,EUR] \\
        [--plink-conv  CHR#.plink.maf]
"""

import argparse
import glob
import os
import re
import sys

# ── Constants ─────────────────────────────────────────────────────────────────

_BASE = "/lustre/orion/bif154/proj-shared/mconery/infeRSparsePro"
_SAIGE = f"{_BASE}/SAIGE_simulations"
_SCRIPTS = f"{_BASE}/Betanfer/paper_analyses/4_At-Scale_SAIGE_Simulations/3_Fine-Mapping_Scripts"
_TOOLS = "/lustre/orion/bif154/proj-shared/mconery/tools"

RESULTS_BASE   = f"{_SAIGE}/fine-mapping_results"
META_DIR       = f"{_SAIGE}/meta_analysis_results"
STEP2_DIR      = f"{_SAIGE}/step2_results"
LD_BASE        = f"{_SAIGE}/LD_matrices"
PHENO_DIR      = f"{_SAIGE}/simulated_phenotypes"
GENO_DIR       = f"{_SAIGE}/simulated_genotypes"
CMD_FILES_DIR  = f"{RESULTS_BASE}/command_files"

# Conda environment paths
ENVS = {
    "carma":      f"{_TOOLS}/environments/carma",
    "susie":      f"{_TOOLS}/environments/susie",
    "rsparsepro": f"{_TOOLS}/environments/rsparsepro",
    "mesusie":    f"{_TOOLS}/environments/carma",   # MESuSiE uses the carma env
    "susiex":     None,                             # uses default miniforge3
    "multisusie": f"{_TOOLS}/environments/multisusie",
}

# Wrapper script paths
WRAPPERS = {
    "carma":      f"{_SCRIPTS}/CARMA/finemap_CARMA.R",
    "susie":      f"{_SCRIPTS}/SuSiE/finemap_SuSiE.R",
    "rsparsepro": f"{_SCRIPTS}/RSparsePro/finemap_RSparsePro.sh",
    "mesusie":    f"{_SCRIPTS}/MESuSiE/finemap_MESuSiE.R",
    "susiex":     f"{_SCRIPTS}/SuSiEx/finemap_SuSiEx.py",
    "multisusie": f"{_SCRIPTS}/MultiSuSiE/finemap_MultiSuSiE.py",
}

SINGLE_POP_METHODS = {"carma", "susie", "rsparsepro"}
MULTI_POP_METHODS  = {"mesusie", "susiex", "multisusie"}

# Regex to parse output filenames
FNAME_RE = re.compile(
    r'^(Pheno\d+\.\d+_causal)\.(binary_[\d.]+|continuous)\.(.+)\.txt$'
)

# Regex to extract heritability value from a filename string.
# Uses a proper floating-point pattern to avoid greedily matching the trailing
# dot separator (e.g., avoids extracting "0.0003." from "heritability_0.0003.binary_...").
H_RE = re.compile(r'heritability_(\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)')

# Regex to extract chromosome from an LD filename
CHR_RE = re.compile(r'\.(chr\d+)\.')


# ── Analysis directory parsing ─────────────────────────────────────────────────

def parse_analysis_dir(analysis_dir):
    """
    Parse an analysis directory name into its constituent parameters.

    Handles all four analysis types:
      {fixed|random}_homo_{h}
      {fixed|random}_sub-sample_LD_{h}
      {fixed|random}_matched_to_{POP}_{h}
      {fixed|random}_drop_causal_{POP}_{h}

    Returns a dict with:
      analysis_type    : 'fixed' or 'random'
      heritability_label: e.g., '3e-4'
      heritability_value: float, e.g., 0.0003
      filter_keyword   : str or None (e.g., 'matched_to_AFR', 'drop_causal_EUR')
      is_sub_sample    : bool
    """
    HERITABILITY_MAP = {"3e-4": 0.0003, "5e-5": 0.00005}

    analysis_type = "fixed" if analysis_dir.startswith("fixed") else "random"

    # Heritability label is always the last underscore-delimited token
    h_label = analysis_dir.rsplit("_", 1)[-1]
    h_value = HERITABILITY_MAP.get(h_label)

    is_sub_sample = "sub-sample" in analysis_dir

    filter_keyword = None
    if "matched_to_" in analysis_dir:
        # e.g., fixed_matched_to_AFR_3e-4 → filter_keyword = matched_to_AFR
        idx = analysis_dir.index("matched_to_")
        tail = analysis_dir[idx:]                      # "matched_to_AFR_3e-4"
        pop_part = tail[len("matched_to_"):]           # "AFR_3e-4"
        pop = pop_part.rsplit("_", 1)[0]               # "AFR"
        filter_keyword = f"matched_to_{pop}"
    elif "drop_causal_" in analysis_dir:
        idx = analysis_dir.index("drop_causal_")
        tail = analysis_dir[idx:]
        pop_part = tail[len("drop_causal_"):]
        pop = pop_part.rsplit("_", 1)[0]
        filter_keyword = f"drop_causal_{pop}"

    return {
        "analysis_type":      analysis_type,
        "heritability_label": h_label,
        "heritability_value": h_value,
        "filter_keyword":     filter_keyword,
        "is_sub_sample":      is_sub_sample,
    }


# ── Filename parsing ───────────────────────────────────────────────────────────

def parse_missing_file_path(path):
    """
    Parse a missing output file path into constituent components.

    Expected path structure:
      {results_base}/{analysis_dir}/{method}_results/{ld_dir_name}/{filename}

    Returns a dict with keys:
      analysis_dir, method, ld_dir_name, pheno_causal, trait_type_str, suffix
    """
    norm = os.path.normpath(path)
    parts = norm.split(os.sep)
    # We need at least 4 path components after results_base
    if len(parts) < 5:
        raise ValueError(f"Path too short to parse: {path}")

    filename        = parts[-1]
    ld_dir_name     = parts[-2]
    method_results  = parts[-3]   # e.g., "carma_results"
    analysis_dir    = parts[-4]

    method = method_results.replace("_results", "")
    if method not in WRAPPERS:
        raise ValueError(f"Unknown method directory '{method_results}' in path: {path}")

    m = FNAME_RE.match(filename)
    if not m:
        raise ValueError(f"Cannot parse filename '{filename}' in path: {path}")

    return {
        "analysis_dir":   analysis_dir,
        "method":         method,
        "ld_dir_name":    ld_dir_name,
        "pheno_causal":   m.group(1),
        "trait_type_str": m.group(2),
        "suffix":         m.group(3),
    }


# ── Heritability helpers ───────────────────────────────────────────────────────

def _heritability_matches(filename, target_value):
    """Return True if the heritability encoded in filename equals target_value."""
    m = H_RE.search(filename)
    if not m:
        return False
    try:
        return abs(float(m.group(1)) - target_value) < 1e-10
    except (ValueError, TypeError):
        return False


# ── Chromosome extraction ──────────────────────────────────────────────────────

def extract_chromosome(filepath):
    """
    Extract the chromosome string (e.g., 'chr1') from an LD matrix file.

    First tries to extract from the filename (pattern: Pheno#.#_causal.chr#.vcor1).
    If the filename lacks a chr token (e.g., beta-inferred / pc-corrected / raw_in-sample
    LD files), falls back to reading the first line of the companion .vars file, which
    stores variant IDs in the format 'CHROM:POS:REF:ALT'.
    """
    m = CHR_RE.search(os.path.basename(filepath))
    if m:
        return m.group(1)

    # Fallback: read chromosome from the .vars file
    vars_path = filepath + ".vars"
    try:
        with open(vars_path) as fh:
            first_line = fh.readline().strip()
        chrom_raw = first_line.split(":")[0]  # e.g., "3" from "3:57261894:C:T"
        if not chrom_raw:
            raise ValueError("Empty chromosome field in .vars file")
        # Normalise to "chrN" format
        if not chrom_raw.startswith("chr"):
            chrom_raw = "chr" + chrom_raw
        return chrom_raw
    except (OSError, IndexError, ValueError) as exc:
        raise ValueError(
            f"Cannot extract chromosome from LD filename: {os.path.basename(filepath)}"
            f" (also failed to read {vars_path}: {exc})"
        ) from exc


# ── LD matrix file discovery ───────────────────────────────────────────────────

def find_ld_file_single_pop(ld_dir_name, pheno_causal, ld_base, suffix):
    """
    Locate the single-population LD matrix (.vcor1) file.

    ld_dir_name encodes both the matrix type and subfolder:
      - If it contains a '.', the last component is the subfolder
        (e.g., 'beta-inferred_100_fixed_all.Meta' → matrix='beta-inferred_100_fixed_all', sub='Meta')
      - Otherwise the subfolder is inferred from the filename suffix
        (e.g., 'raw_in-sample' with suffix 'weighted' → sub='weighted')

    Returns path string.  Raises FileNotFoundError if not found.
    """
    if "." in ld_dir_name:
        last_dot = ld_dir_name.rfind(".")
        matrix_type = ld_dir_name[:last_dot]
        subfolder   = ld_dir_name[last_dot + 1:]
    else:
        matrix_type = ld_dir_name
        subfolder   = suffix  # e.g., "weighted" or "Meta"

    search_dir = os.path.join(ld_base, matrix_type, subfolder)
    pattern    = os.path.join(search_dir, f"{pheno_causal}*.vcor1")
    matches    = glob.glob(pattern)

    if not matches:
        raise FileNotFoundError(
            f"No LD file found for {pheno_causal} in {search_dir}\n"
            f"  Pattern: {pattern}"
        )
    return matches[0]


def find_ld_files_multi_pop(ld_dir_name, pheno_causal, ld_base, populations):
    """
    Locate LD matrix files for each population.

    For multi-population methods, LD files live at:
      {ld_base}/{ld_dir_name}/{POP}/{pheno_causal}*.vcor1

    Returns dict: {pop: path}.  Raises FileNotFoundError if any is missing.
    """
    result = {}
    for pop in populations:
        search_dir = os.path.join(ld_base, ld_dir_name, pop)
        pattern    = os.path.join(search_dir, f"{pheno_causal}*.vcor1")
        matches    = glob.glob(pattern)
        if not matches:
            raise FileNotFoundError(
                f"No LD file found for pop={pop} {pheno_causal}\n"
                f"  Pattern: {pattern}"
            )
        result[pop] = matches[0]
    return result


# ── GWAS file discovery ────────────────────────────────────────────────────────

def find_gwas_file_single_pop(pheno_causal, trait_type_str, analysis_type,
                               heritability_value, chromosome, meta_dir,
                               filter_keyword):
    """
    Locate the single-population meta-analysis GWAS file for this replicate.

    GWAS files live in per-phenotype subdirectories:
      {meta_dir}/{pheno_causal}/{pheno_causal}.heritability_*.{trait_type_str}.{chromosome}[.{filter_keyword}].{analysis_type}_meta.out

    For standard analyses, the filter_keyword does NOT appear in the filename.
    For matched_to/drop_causal analyses, the filter_keyword appears AFTER the
    chromosome in the filename (with full population suffix, e.g. 'drop_causal_EUR').

    Returns the path string.  Raises FileNotFoundError if nothing is found.
    """
    # Subdirectory is named: {pheno_causal}.heritability_{h_value}.{trait_type_str}
    # Actual filename: {pheno_causal}.heritability_X.{trait_type}.[keyword.]{chromosome}.{analysis_type}_meta.out
    # The keyword (if any) appears between the trait_type and chromosome in the filename.
    # Use wildcards for both heritability and any intervening tokens, then post-filter.
    subdir_pattern = f"{pheno_causal}.heritability_*.{trait_type_str}"
    pattern = os.path.join(
        meta_dir,
        subdir_pattern,
        f"{pheno_causal}.*.{trait_type_str}.*.{analysis_type}_meta.out",
    )
    matches = glob.glob(pattern)

    # Post-filter by chromosome (chromosome token appears somewhere after the trait_type)
    if chromosome:
        matches = [p for p in matches if f".{chromosome}." in os.path.basename(p)]

    if filter_keyword:
        # Include only files that contain the full filter_keyword in their basename
        matches = [p for p in matches if filter_keyword in os.path.basename(p)]
    else:
        # Standard analysis: exclude files with any keyword pattern
        matches = [
            p for p in matches
            if "matched_to_" not in os.path.basename(p)
            and "drop_causal_" not in os.path.basename(p)
        ]

    # Post-filter by heritability value
    if heritability_value is not None:
        matches = [p for p in matches
                   if _heritability_matches(os.path.basename(p), heritability_value)]

    if not matches:
        raise FileNotFoundError(
            f"No single-pop GWAS file found for {pheno_causal} ({trait_type_str}) "
            f"analysis_type={analysis_type} chromosome={chromosome} filter={filter_keyword}\n"
            f"  Searched: {pattern}"
        )
    if len(matches) > 1:
        print(
            f"[WARNING] Multiple GWAS matches for {pheno_causal} {chromosome}; using first: "
            f"{os.path.basename(matches[0])}",
            file=sys.stderr,
        )
    return matches[0]


def find_gwas_files_multi_pop(pheno_causal, trait_type_str, heritability_value,
                               step2_dir, populations, filter_keyword):
    """
    Locate step2 GWAS files for each population.

    Multi-pop LD files do not carry chromosome information in their filenames,
    so chromosome is NOT used as a constraint here.  Each simulated phenotype
    has a single causal locus (one chromosome), so the phenotype-level glob
    will naturally return at most one matching file per population.

    File naming conventions differ by analysis type AND by population:

    Standard (no filter_keyword):
      All pops: {pheno_causal}.*.{trait_type}.*.step2.txt
      (exclude files with 'drop_causal' or 'matched_to_' in basename)

    drop_causal_X (central_pop = X):
      X  (central): {pheno_causal}.*.{trait_type}.*.drop_causal.step2.txt
      Others:       {pheno_causal}.*.{trait_type}.*.step2.txt (standard, no keyword)

    matched_to_X (central_pop = X):
      X  (central): {pheno_causal}.*.{trait_type}.*.step2.txt (standard, no keyword)
      Others:       {pheno_causal}.*.{trait_type}.matched_to_X.*.step2.txt
                    (note: keyword appears BEFORE chromosome in step2 filenames)

    Returns dict: {pop: path}.  Raises FileNotFoundError if any population is missing.
    """
    # Determine the central population from the filter_keyword
    central_pop = None
    keyword_type = None
    if filter_keyword:
        if filter_keyword.startswith("drop_causal_"):
            keyword_type = "drop_causal"
            central_pop = filter_keyword[len("drop_causal_"):]
        elif filter_keyword.startswith("matched_to_"):
            keyword_type = "matched_to"
            central_pop = filter_keyword[len("matched_to_"):]

    result = {}
    for pop in populations:
        pop_dir = os.path.join(step2_dir, pop)

        if keyword_type == "drop_causal":
            if pop == central_pop:
                # Central pop has the dropped-causal file (keyword after chromosome, no pop suffix)
                pattern = os.path.join(
                    pop_dir,
                    f"{pheno_causal}.*.{trait_type_str}.*.drop_causal.step2.txt",
                )
                matches = glob.glob(pattern)
            else:
                # Non-central pops use standard step2 files
                pattern = os.path.join(
                    pop_dir,
                    f"{pheno_causal}.*.{trait_type_str}.*.step2.txt",
                )
                matches = glob.glob(pattern)
                matches = [
                    p for p in matches
                    if "matched_to_" not in os.path.basename(p)
                    and "drop_causal" not in os.path.basename(p)
                ]

        elif keyword_type == "matched_to":
            if pop == central_pop:
                # Central pop uses standard step2 files
                pattern = os.path.join(
                    pop_dir,
                    f"{pheno_causal}.*.{trait_type_str}.*.step2.txt",
                )
                matches = glob.glob(pattern)
                matches = [
                    p for p in matches
                    if "matched_to_" not in os.path.basename(p)
                    and "drop_causal" not in os.path.basename(p)
                ]
            else:
                # Non-central pops have matched_to_X BEFORE the chromosome
                pattern = os.path.join(
                    pop_dir,
                    f"{pheno_causal}.*.{trait_type_str}.matched_to_{central_pop}.*.step2.txt",
                )
                matches = glob.glob(pattern)

        else:
            # Standard analysis — no keyword in filename
            pattern = os.path.join(
                pop_dir,
                f"{pheno_causal}.*.{trait_type_str}.*.step2.txt",
            )
            matches = glob.glob(pattern)
            matches = [
                p for p in matches
                if "matched_to_" not in os.path.basename(p)
                and "drop_causal" not in os.path.basename(p)
            ]

        # Post-filter by heritability
        if heritability_value is not None:
            matches = [p for p in matches
                       if _heritability_matches(os.path.basename(p), heritability_value)]

        if not matches:
            raise FileNotFoundError(
                f"No step2 GWAS file found for pop={pop} {pheno_causal} "
                f"({trait_type_str}) filter={filter_keyword}\n"
                f"  Searched: {pattern}"
            )
        if len(matches) > 1:
            print(
                f"[WARNING] Multiple GWAS matches for {pop} {pheno_causal}; "
                f"using first: {os.path.basename(matches[0])}",
                file=sys.stderr,
            )
        result[pop] = matches[0]

    return result


# ── Auxiliary file discovery (MultiSuSiE, SuSiEx) ─────────────────────────────

def find_pheno_files(pheno_causal, trait_type_str, pheno_dir, populations):
    """
    Locate per-population phenotype files needed by MultiSuSiE.

    Returns dict: {pop: path}.  Raises FileNotFoundError if any is missing.
    """
    result = {}
    for pop in populations:
        pop_dir = os.path.join(pheno_dir, pop)
        pattern = os.path.join(pop_dir, f"{pheno_causal}*.{trait_type_str}.txt")
        matches = glob.glob(pattern)
        if not matches:
            raise FileNotFoundError(
                f"No pheno file found for pop={pop} {pheno_causal} ({trait_type_str})\n"
                f"  Pattern: {pattern}"
            )
        result[pop] = matches[0]
    return result


def get_plink_paths(chromosome, geno_dir, populations, plink_conv="CHR#.plink.maf"):
    """
    Construct PLINK reference panel paths for SuSiEx.

    Uses the chromosome already extracted from the LD matrix file to build:
      {geno_dir}/{POP}/{chromosome}.plink.maf
    for each population (in populations order).

    Returns list of path strings.
    """
    plink_paths = []
    for pop in populations:
        panel_name = plink_conv.replace("CHR#", chromosome)
        plink_paths.append(os.path.join(geno_dir, pop, panel_name))
    return plink_paths


# ── Command construction ───────────────────────────────────────────────────────

def construct_command(method, gwas, ld, output_path,
                      populations, pheno_files=None, plink_paths=None):
    """
    Build the fine-mapping command string for a single replicate.

    Parameters
    ----------
    method      : one of the six method names
    gwas        : str (single-pop) or dict {pop: path} (multi-pop)
    ld          : str (single-pop) or dict {pop: path} (multi-pop)
    output_path : destination .txt file
    populations : ordered list of population names (for multi-pop ordering)
    pheno_files : dict {pop: path}  (MultiSuSiE only)
    plink_paths : list of str       (SuSiEx only)
    """
    wrapper = WRAPPERS[method]

    if method == "carma":
        return f"Rscript {wrapper} {gwas} {ld} {output_path}"

    if method == "susie":
        return f"Rscript {wrapper} {gwas} {ld} {output_path}"

    if method == "rsparsepro":
        return f"{wrapper} {gwas} {ld} {output_path}"

    # Multi-population methods — order populations consistently
    pops = [p for p in populations if p in gwas]
    n_gwas_str = ",".join(["40000"] * len(pops))
    labels_str = ",".join(pops)

    if method == "mesusie":
        gwas_str = ",".join(gwas[p] for p in pops)
        ld_str   = ",".join(ld[p]   for p in pops)
        return (
            f"Rscript {wrapper} "
            f"--gwas_files {gwas_str} "
            f"--ld_files {ld_str} "
            f"--n_gwas {n_gwas_str} "
            f"--ancestry_labels {labels_str} "
            f"--output_file {output_path}"
        )

    if method == "susiex":
        gwas_str   = ",".join(gwas[p] for p in pops)
        ld_str     = ",".join(ld[p]   for p in pops)
        panels_str = ",".join(plink_paths)
        return (
            f"python {wrapper} "
            f"--gwas_files {gwas_str} "
            f"--ld_files {ld_str} "
            f"--n_gwas {n_gwas_str} "
            f"--ref_panels {panels_str} "
            f"--output_file {output_path}"
        )

    if method == "multisusie":
        gwas_str  = ",".join(gwas[p]         for p in pops)
        ld_str    = ",".join(ld[p]           for p in pops)
        pheno_str = ",".join(pheno_files[p]  for p in pops)
        return (
            f"python {wrapper} "
            f"--gwas_files {gwas_str} "
            f"--ld_files {ld_str} "
            f"--pheno_files {pheno_str} "
            f"--output_file {output_path} "
            f"--ancestry_labels {labels_str}"
        )

    raise ValueError(f"Unknown method: {method}")


# ── Main ───────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Given a list of missing fine-mapping result file paths, delete any "
            "existing corrupt copies, construct the regeneration command for each, "
            "and write six per-method command files."
        )
    )
    parser.add_argument(
        "--missing-list", "-i",
        required=True,
        help="File containing missing result paths, one per line",
    )
    parser.add_argument(
        "--output-dir", "-o",
        default=CMD_FILES_DIR,
        help=f"Directory to write the six command files (default: {CMD_FILES_DIR})",
    )
    parser.add_argument(
        "--populations", "-pop",
        default="AFR,EUR",
        help="Comma-separated population names for multi-pop methods (default: AFR,EUR)",
    )
    parser.add_argument(
        "--plink-conv",
        default="CHR#.plink.maf",
        help="PLINK reference panel naming convention, CHR# replaced by e.g. chr1 (default: CHR#.plink.maf)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    populations = [p.strip().upper() for p in args.populations.split(",") if p.strip()]
    plink_conv  = args.plink_conv
    output_dir  = args.output_dir

    # Read missing file list
    with open(args.missing_list) as fh:
        missing_paths = [line.strip() for line in fh if line.strip()]

    print(f"[INFO] {len(missing_paths):,} missing file paths loaded.", file=sys.stderr)

    os.makedirs(output_dir, exist_ok=True)

    # ── Phase 1: Delete existing files ──────────────────────────────────────
    print("[Phase 1] Deleting existing (corrupt/incomplete) files ...", file=sys.stderr)
    deleted = 0
    for path in missing_paths:
        if os.path.exists(path):
            try:
                os.remove(path)
                deleted += 1
            except OSError as exc:
                print(f"[WARNING] Could not delete {path}: {exc}", file=sys.stderr)
    print(f"  Deleted {deleted:,} file(s).", file=sys.stderr)

    # ── Phase 2: Construct commands ──────────────────────────────────────────
    print("[Phase 2] Constructing regeneration commands ...", file=sys.stderr)

    commands_by_method = {m: [] for m in WRAPPERS}
    errors = []

    for path in missing_paths:
        # Parse the output path
        try:
            info = parse_missing_file_path(path)
        except ValueError as exc:
            errors.append((path, str(exc)))
            continue

        analysis_dir   = info["analysis_dir"]
        method         = info["method"]
        ld_dir_name    = info["ld_dir_name"]
        pheno_causal   = info["pheno_causal"]
        trait_type_str = info["trait_type_str"]
        suffix         = info["suffix"]

        # Parse the analysis directory name
        try:
            ainfo = parse_analysis_dir(analysis_dir)
        except Exception as exc:
            errors.append((path, f"Cannot parse analysis_dir '{analysis_dir}': {exc}"))
            continue

        analysis_type  = ainfo["analysis_type"]
        h_value        = ainfo["heritability_value"]
        filter_keyword = ainfo["filter_keyword"]
        is_sub_sample  = ainfo["is_sub_sample"]

        # LD base depends on whether this is a sub-sample analysis
        ld_base = os.path.join(LD_BASE, "raw_in-sample") if is_sub_sample else LD_BASE

        # Build command
        try:
            if method in SINGLE_POP_METHODS:
                # Step A: find LD file first to extract chromosome
                ld = find_ld_file_single_pop(ld_dir_name, pheno_causal, ld_base, suffix)
                chromosome = extract_chromosome(ld)

                # Step B: find GWAS file using chromosome
                gwas = find_gwas_file_single_pop(
                    pheno_causal, trait_type_str, analysis_type,
                    h_value, chromosome, META_DIR, filter_keyword,
                )
                cmd = construct_command(method, gwas, ld, path, populations)

            else:  # multi-pop
                # Step A: find LD files (multi-pop LD filenames don't carry chromosome)
                ld = find_ld_files_multi_pop(ld_dir_name, pheno_causal, ld_base, populations)

                # Step B: find GWAS files with per-population keyword logic
                # (no chromosome constraint — each phenotype has one causal locus)
                gwas = find_gwas_files_multi_pop(
                    pheno_causal, trait_type_str, h_value,
                    STEP2_DIR, populations, filter_keyword,
                )

                pheno_files = None
                plink_paths = None

                if method == "multisusie":
                    pheno_files = find_pheno_files(
                        pheno_causal, trait_type_str, PHENO_DIR, populations
                    )
                elif method == "susiex":
                    # Extract chromosome from one of the found GWAS files for PLINK path
                    sample_gwas = next(iter(gwas.values()))
                    chr_m = re.search(r'\.(chr\d+)\.', os.path.basename(sample_gwas))
                    if not chr_m:
                        raise ValueError(
                            f"Cannot extract chromosome from GWAS filename: "
                            f"{os.path.basename(sample_gwas)}"
                        )
                    chromosome = chr_m.group(1)
                    plink_paths = get_plink_paths(
                        chromosome, GENO_DIR, populations, plink_conv
                    )

                cmd = construct_command(
                    method, gwas, ld, path, populations,
                    pheno_files=pheno_files, plink_paths=plink_paths,
                )

            commands_by_method[method].append(cmd)

        except (FileNotFoundError, ValueError) as exc:
            errors.append((path, str(exc)))

    # ── Phase 3: Write command files ─────────────────────────────────────────
    print("[Phase 3] Writing command files ...", file=sys.stderr)

    total_cmds = 0
    for method, cmds in commands_by_method.items():
        cmd_file = os.path.join(output_dir, f"{method}_regenerate_commands.txt")
        with open(cmd_file, "w") as fh:
            for cmd in cmds:
                fh.write(cmd + "\n")
        count = len(cmds)
        total_cmds += count
        print(f"  {method:<12}: {count:>6,} commands → {cmd_file}", file=sys.stderr)

    # ── Summary ──────────────────────────────────────────────────────────────
    print(file=sys.stderr)
    print(f"Total commands written : {total_cmds:,}", file=sys.stderr)
    print(f"Paths with errors      : {len(errors):,}", file=sys.stderr)
    if errors:
        err_file = os.path.join(output_dir, "regenerate_errors.txt")
        with open(err_file, "w") as fh:
            for path, msg in errors:
                fh.write(f"{path}\n  ERROR: {msg}\n\n")
        print(f"Error details written to: {err_file}", file=sys.stderr)


if __name__ == "__main__":
    main()
