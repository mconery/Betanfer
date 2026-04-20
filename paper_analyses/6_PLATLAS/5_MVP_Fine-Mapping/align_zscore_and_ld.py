'''
align_zscore_and_ld.py

Preprocessing step for MVP RSparsePro fine-mapping.

rsparsepro_ld.py requires that the z-score vector and the LD correlation
matrix describe exactly the same variants in exactly the same order.  Two
mismatches arise upstream:

  1. Ordering mismatch: betanfer builds LD matrices with variants in
     lexicographic RSID order (sorted()), while prep_z-score_files.py
     outputs variants in SE-ascending order.  The orders differ, so
     z[i] is passed to the wrong LD row/column.

  2. Dimension mismatch: variants listed in the .vars file but absent
     from the sumstats are silently dropped by prep_z-score_files.py,
     leaving the z-score file with M rows while the LD matrix has N > M
     rows/columns.  rsparsepro_ld.py crashes with a ValueError at
     np.dot(R, beta_all_k).

This script resolves both issues in a single pass:
  - Reads the .vars file (authoritative RSID order) and the z-score file.
  - Computes the intersection of the two RSID sets, preserving .vars order.
  - Rewrites the z-score file in .vars order, keeping only intersection RSIDs.
  - Loads the full LD matrix, extracts the submatrix for intersection RSIDs,
    and writes it to a new file.

The output z-score file and output LD file are then safe to pass directly
to rsparsepro_ld.py.

Usage:
  python align_zscore_and_ld.py \\
      --zscore  <z-score file from prep_z-score_files.py> \\
      --vars    <.vcor1.vars file>  \\
      --ld      <.vcor1 LD matrix file> \\
      --out_zscore <aligned z-score output path> \\
      --out_ld     <aligned LD matrix output path>

***Run in finemap conda environment***
'''

import sys
import getopt
import os
import numpy as np
import pandas as pd


def _read_ld_submatrix(ld_file, keep_indices):
    """Stream only the needed rows/columns from a tab-delimited LD matrix."""
    keep_set = set(keep_indices)
    max_row  = keep_indices[-1]   # keep_indices is sorted ascending

    def _needed_rows():
        with open(ld_file) as fh:
            for i, line in enumerate(fh):
                if i > max_row:
                    return
                if i in keep_set:
                    yield line

    return np.loadtxt(
        _needed_rows(),
        delimiter='\t',
        usecols=keep_indices,
        dtype=np.float32,
    )


def help(exit_num=0):
    print("""Usage: python align_zscore_and_ld.py [OPTIONS]
Options:
  --zscore     Input z-score file (TSV with RSID and Z columns) (required)
  --vars       Input .vcor1.vars file (one RSID per line, authoritative order) (required)
  --ld         Input .vcor1 LD correlation matrix file (required)
  --out_zscore Output path for aligned z-score file (required)
  --out_ld     Output path for aligned LD matrix file (required)
  -h, --help   Show this help message""")
    sys.exit(exit_num)


def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "zscore=", "vars=", "ld=", "out_zscore=", "out_ld="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    zscore_file  = None
    vars_file    = None
    ld_file      = None
    out_zscore   = None
    out_ld       = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--zscore":
            zscore_file = a
        elif o == "--vars":
            vars_file = a
        elif o == "--ld":
            ld_file = a
        elif o == "--out_zscore":
            out_zscore = a
        elif o == "--out_ld":
            out_ld = a
        else:
            assert False, f"Unhandled option: {o}"

    missing = [name for name, val in [
        ("--zscore",     zscore_file),
        ("--vars",       vars_file),
        ("--ld",         ld_file),
        ("--out_zscore", out_zscore),
        ("--out_ld",     out_ld),
    ] if val is None]
    if missing:
        print(f"ERROR: Missing required option(s): {', '.join(missing)}")
        sys.exit(2)

    for path, label in [
        (zscore_file, "z-score file"),
        (vars_file,   ".vars file"),
        (ld_file,     "LD matrix file"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Load .vars (authoritative RSID order) ─────────────────────────────────
    with open(vars_file) as fh:
        vars_order = [line.strip() for line in fh if line.strip()]
    vars_index = {rsid: i for i, rsid in enumerate(vars_order)}

    # ── Load z-score file ─────────────────────────────────────────────────────
    zdf = pd.read_csv(zscore_file, sep='\t')
    if 'RSID' not in zdf.columns or 'Z' not in zdf.columns:
        print(f"ERROR: z-score file must have RSID and Z columns; found: {list(zdf.columns)}")
        sys.exit(1)

    # ── Compute intersection in .vars order ───────────────────────────────────
    z_rsids = set(zdf['RSID'])
    keep = [(vars_index[rsid], rsid) for rsid in vars_order if rsid in z_rsids]

    if not keep:
        print("ERROR: No overlapping RSIDs between z-score file and .vars file")
        sys.exit(1)

    keep_indices = [idx for idx, _ in keep]
    keep_rsids   = [rsid for _, rsid in keep]

    n_vars   = len(vars_order)
    n_zscore = len(zdf)
    n_keep   = len(keep)
    print(f"NOTE: .vars variants={n_vars}  z-score variants={n_zscore}  "
          f"intersection={n_keep}")
    if n_keep < n_zscore:
        print(f"NOTE: {n_zscore - n_keep} z-score variant(s) not found in .vars — dropped")
    if n_keep < n_vars:
        print(f"NOTE: {n_vars - n_keep} .vars variant(s) absent from sumstats — "
              f"LD matrix will be subsetted to {n_keep}×{n_keep}")

    # ── Write aligned z-score file (intersection in .vars order) ──────────────
    zdf_indexed  = zdf.set_index('RSID')
    zdf_aligned  = zdf_indexed.loc[keep_rsids].reset_index()[['RSID', 'Z']]
    zdf_aligned.to_csv(out_zscore, sep='\t', index=False)
    print(f"NOTE: Aligned z-score written to {out_zscore}")

    # ── Load LD submatrix (streaming — avoids loading full N×N into memory) ──
    with open(ld_file) as fh:
        n_lines = sum(1 for _ in fh)
    if n_lines != n_vars:
        print(f"WARNING: LD matrix has {n_lines} rows but .vars has {n_vars} entries; "
              f"proceeding with available dimensions")
        keep_indices = [i for i in keep_indices if i < n_lines]

    ld_sub = _read_ld_submatrix(ld_file, keep_indices)
    np.savetxt(out_ld, ld_sub, delimiter=' ', fmt='%.6f')
    print(f"NOTE: Aligned LD matrix ({ld_sub.shape[0]}×{ld_sub.shape[1]}) "
          f"written to {out_ld}")


if __name__ == "__main__":
    main()
