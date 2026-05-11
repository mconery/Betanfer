'''
harmonize_1000g.py

Preprocessing step for 1000G-reference fine-mapping of MVP GWAS loci.

The 1000G LD matrices use CHR:POS:REF:ALT variant IDs (plink2 --set-all-var-ids
"@:#:$r:$a"), while MVP summary statistics have an ID column that also contains
CHR:POS:REF:ALT IDs -- but with allele orientations from the MVP study, which may
differ from the 1000G panel at some loci.

This script:
  1. Reads the 1000G .unphased.vcor1.vars file (authoritative variant order).
  2. Reads a tabix-extracted MVP sumstats file (has #ID, BETA, SE, N columns).
  3. Matches .vars entries to sumstats ID values:
       - Exact match (same CHR:POS:REF:ALT): z = BETA / SE
       - Allele-flip match (CHR:POS:REF:ALT vs CHR:POS:ALT:REF): z = -BETA / SE
       - No match: variant excluded from analysis
       - Multiallelic: keep the sumstats row with the lowest SE per ID
  4. Streams the LD matrix and extracts the submatrix for matched variants
     (same streaming approach as align_zscore_and_ld.py).
  5. Writes requested output files (all flags optional; at least one required):
       --out_zscore   RSID (=1000G ID) and Z [tab-delimited; for RSparsePro]
       --out_ld_rsp   Aligned LD submatrix [space-delimited; for rsparsepro_ld.py]
       --out_sumstats VARIANT_ID, BETA (sign-corrected), SE, N [tab; for SuSiE]
       --out_ld_susie Aligned LD submatrix [tab-delimited; for finemap_SuSiE.R]
     When --out_ld_susie is given, {out_ld_susie}.vars is also written so that
     finemap_SuSiE.R can locate it via paste0(ld_file, ".vars").

Usage:
  python harmonize_1000g.py \\
      --sumstats <tabix-extracted MVP sumstats> \\
      --vars     <.unphased.vcor1.vars file> \\
      --ld       <.unphased.vcor1 LD matrix file> \\
      [--out_zscore   <z-score TSV for RSparsePro>] \\
      [--out_ld_rsp   <space-delimited LD for RSparsePro>] \\
      [--out_sumstats <sumstats TSV for SuSiE>] \\
      [--out_ld_susie <tab-delimited LD for SuSiE>]

***Run in finemap conda environment***
'''

import sys
import getopt
import os
import numpy as np
import pandas as pd


def _read_ld_submatrix(ld_file, keep_indices):
    """Stream LD matrix rows and extract submatrix for matched variant indices.

    Reads only through the last needed row (keep_indices[-1]) to avoid loading
    the full N×N matrix into memory, which can be several gigabytes for large loci.
    Assumes tab-delimited LD matrix (plink2 --r-unphased square output).
    """
    keep_set = set(keep_indices)
    max_row  = keep_indices[-1]

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
    print("""Usage: python harmonize_1000g.py [OPTIONS]
Options:
  --sumstats     Tabix-extracted MVP sumstats file (required)
                 Must have #ID, BETA, SE, N columns; '#' stripped from column 1
  --vars         .unphased.vcor1.vars file (one CHR:POS:REF:ALT per line) (required)
  --ld           .unphased.vcor1 LD matrix file (tab-delimited) (required)
  --out_zscore   Output z-score file (RSID, Z) for RSparsePro (optional)
  --out_ld_rsp   Output LD submatrix (space-delimited) for RSparsePro (optional)
  --out_sumstats Output sumstats (VARIANT_ID, BETA, SE, N) for SuSiE (optional)
  --out_ld_susie Output LD submatrix (tab-delimited) for SuSiE (optional)
                 {out_ld_susie}.vars is written automatically alongside
  -h, --help     Show this help message""")
    sys.exit(exit_num)


def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "sumstats=", "vars=", "ld=",
             "out_zscore=", "out_ld_rsp=", "out_sumstats=", "out_ld_susie="]
        )
    except getopt.GetoptError as err:
        print(str(err))
        sys.exit(2)

    sumstats_file = None
    vars_file     = None
    ld_file       = None
    out_zscore    = None
    out_ld_rsp    = None
    out_sumstats  = None
    out_ld_susie  = None

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--sumstats":
            sumstats_file = a
        elif o == "--vars":
            vars_file = a
        elif o == "--ld":
            ld_file = a
        elif o == "--out_zscore":
            out_zscore = a
        elif o == "--out_ld_rsp":
            out_ld_rsp = a
        elif o == "--out_sumstats":
            out_sumstats = a
        elif o == "--out_ld_susie":
            out_ld_susie = a
        else:
            assert False, f"Unhandled option: {o}"

    # Validate required inputs
    missing_inputs = [name for name, val in [
        ("--sumstats", sumstats_file),
        ("--vars",     vars_file),
        ("--ld",       ld_file),
    ] if val is None]
    if missing_inputs:
        print(f"ERROR: Missing required option(s): {', '.join(missing_inputs)}")
        sys.exit(2)

    # Require at least one output
    outputs = [out_zscore, out_ld_rsp, out_sumstats, out_ld_susie]
    if all(o is None for o in outputs):
        print("ERROR: At least one --out_* flag must be specified")
        sys.exit(2)

    for path, label in [
        (sumstats_file, "sumstats file"),
        (vars_file,     ".vars file"),
        (ld_file,       "LD matrix file"),
    ]:
        if not os.path.exists(path):
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    # ── Load .vars (authoritative 1000G variant order) ────────────────────────
    with open(vars_file) as fh:
        vars_order = [line.strip() for line in fh if line.strip()]
    print(f"NOTE: {len(vars_order)} variants in .vars file")

    # ── Load MVP sumstats ──────────────────────────────────────────────────────
    # Auto-detect delimiter; strip leading '#' from first column header.
    df = pd.read_csv(sumstats_file, sep=None, engine='python')
    df.columns = [df.columns[0].lstrip('#')] + list(df.columns[1:])

    for col in ('ID', 'BETA', 'SE', 'N'):
        if col not in df.columns:
            print(f"ERROR: Required column '{col}' not found in sumstats. "
                  f"Columns present: {list(df.columns)}")
            sys.exit(1)

    # Coerce numeric types; handle multiallelic: keep lowest-SE row per ID.
    df['SE']   = pd.to_numeric(df['SE'],   errors='coerce')
    df['BETA'] = pd.to_numeric(df['BETA'], errors='coerce')
    df['N']    = pd.to_numeric(df['N'],    errors='coerce')
    df = df.dropna(subset=['ID', 'BETA', 'SE'])
    n_se_zero = (df['SE'] <= 0).sum()
    if n_se_zero:
        print(f"NOTE: Removed {n_se_zero} variants with SE <= 0")
    df = df[df['SE'] > 0]
    df = df.sort_values('SE', ascending=True)
    df = df.drop_duplicates(subset=['ID'], keep='first')
    sumstats_by_id = df.set_index('ID')[['BETA', 'SE', 'N']]
    print(f"NOTE: {len(sumstats_by_id)} unique variants in sumstats after dedup")

    # ── Match .vars entries to sumstats ───────────────────────────────────────
    matched     = []   # list of (vars_idx, vars_id, z, beta_corrected, se, n)
    n_exact     = 0
    n_flipped   = 0
    n_unmatched = 0

    for i, vid in enumerate(vars_order):
        if vid in sumstats_by_id.index:
            row  = sumstats_by_id.loc[vid]
            beta = float(row['BETA'])
            se   = float(row['SE'])
            n    = float(row['N'])
            matched.append((i, vid, beta / se, beta, se, n))
            n_exact += 1
        else:
            parts = vid.split(':')
            if len(parts) == 4:
                chrom, pos, ref, alt = parts
                flipped_id = f"{chrom}:{pos}:{alt}:{ref}"
                if flipped_id in sumstats_by_id.index:
                    row  = sumstats_by_id.loc[flipped_id]
                    beta = float(row['BETA'])
                    se   = float(row['SE'])
                    n    = float(row['N'])
                    # Negate: sumstats effect allele = 1000G REF; flip to 1000G ALT frame
                    matched.append((i, vid, -beta / se, -beta, se, n))
                    n_flipped += 1
                else:
                    n_unmatched += 1
            else:
                n_unmatched += 1

    print(f"NOTE: Exact={n_exact}  Allele-flip={n_flipped}  "
          f"Unmatched={n_unmatched}  Total matched={len(matched)}")

    if not matched:
        print("ERROR: No variants matched between sumstats and .vars file")
        sys.exit(1)

    keep_indices = [t[0] for t in matched]
    keep_ids     = [t[1] for t in matched]
    z_values     = [t[2] for t in matched]
    beta_values  = [t[3] for t in matched]
    se_values    = [t[4] for t in matched]
    n_values     = [t[5] for t in matched]

    # ── Stream LD submatrix (only if an LD output is requested) ───────────────
    ld_sub = None
    if out_ld_rsp is not None or out_ld_susie is not None:
        # Verify LD matrix row count matches .vars count
        with open(ld_file) as fh:
            n_ld_rows = sum(1 for _ in fh)
        if n_ld_rows != len(vars_order):
            print(f"WARNING: LD matrix has {n_ld_rows} rows but .vars has "
                  f"{len(vars_order)} entries; subsetting to available rows")
            keep_indices_valid = [i for i in keep_indices if i < n_ld_rows]
            matched = [(i, vid, z, b, s, n_v)
                       for (i, vid, z, b, s, n_v) in matched if i < n_ld_rows]
            keep_indices = [t[0] for t in matched]
            keep_ids     = [t[1] for t in matched]
            z_values     = [t[2] for t in matched]
            beta_values  = [t[3] for t in matched]
            se_values    = [t[4] for t in matched]
            n_values     = [t[5] for t in matched]

        print(f"NOTE: Streaming LD matrix for {len(keep_indices)} matched variants")
        ld_sub = _read_ld_submatrix(ld_file, keep_indices)
        print(f"NOTE: LD submatrix: {ld_sub.shape[0]}×{ld_sub.shape[1]}")

    # ── Write outputs ─────────────────────────────────────────────────────────

    if out_zscore is not None:
        pd.DataFrame({'RSID': keep_ids, 'Z': z_values}).to_csv(
            out_zscore, sep='\t', index=False)
        print(f"NOTE: z-score file → {out_zscore}")

    if out_ld_rsp is not None:
        np.savetxt(out_ld_rsp, ld_sub, delimiter=' ', fmt='%.6f')
        print(f"NOTE: RSparsePro LD matrix → {out_ld_rsp}")

    if out_sumstats is not None:
        pd.DataFrame({
            'VARIANT_ID': keep_ids,
            'BETA':       beta_values,
            'SE':         se_values,
            'N':          n_values,
        }).to_csv(out_sumstats, sep='\t', index=False)
        print(f"NOTE: SuSiE sumstats → {out_sumstats}")

    if out_ld_susie is not None:
        np.savetxt(out_ld_susie, ld_sub, delimiter='\t', fmt='%.6f')
        print(f"NOTE: SuSiE LD matrix → {out_ld_susie}")
        # Write .vars alongside so finemap_SuSiE.R finds it via paste0(ld_file, ".vars")
        with open(out_ld_susie + ".vars", 'w') as fh:
            for vid in keep_ids:
                fh.write(vid + '\n')
        print(f"NOTE: SuSiE .vars → {out_ld_susie}.vars")


if __name__ == "__main__":
    main()
