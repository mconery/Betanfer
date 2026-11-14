# Example: testing your Betanfer install

This folder contains a small synthetic dataset and a walkthrough of the full
Betanfer workflow: indexing summary statistics, inferring an LD matrix, and
checking the result against the true LD. It takes a few minutes on a laptop
and needs no real GWAS data.

## Quick start

From the repository root, with your Betanfer environment active:

```bash
bash example/run_example.sh
```

The last line should be:

```
SUCCESS: Betanfer is installed and working.
```

Everything is written to `betanfer_example/` (or the directory you pass as
the first argument). Delete it when you're done.

Useful options:

| Variable | Default | Use |
|---|---|---|
| `BETANFER` | `betanfer` | Set to `"python -m betanfer"` to test a source checkout without installing |
| `PYTHON` | `python` | Python used for the helper scripts |
| `JOBS` | number of CPUs | Number of `betanfer index` jobs run in parallel |
| `INFER_ARGS` | (none) | Extra arguments for `betanfer infer`, e.g. `"--cpu"` |

Example: `INFER_ARGS="--cpu" bash example/run_example.sh my_test`

## The dataset

`make_example_data.py` generates the data with a fixed random seed, so
everyone gets identical files. It uses only the Python standard library.

- **100 variants** on chromosome 22 (20,001,000–20,199,000), in three LD
  blocks. Adjacent variants within a block have r = 0.95 and LD decays with
  distance; different blocks are independent.
- **150 traits**, one summary-statistics file each, with columns
  `CHR POS SNP A1 A2 MAF BETA SE P`.
  - 140 traits have no real association. Their z-scores are drawn from a
    multivariate normal with the true LD as its correlation matrix. This is
    the property Betanfer uses: across many GWAS, sampling noise in the betas
    is correlated in the same way as LD.
  - 10 traits also carry a real association at one variant (p < 1e-15), so
    the p-value filter has something to remove.
- **8 rare variants** (MAF ≈ 0.3%) have missing statistics (`NA`) in about
  half of the traits, mimicking rare variants lost to QC in many GWAS.
  Common variants are missing in about 2% of traits.

The generator also writes `true_ld.tsv` (the true LD matrix),
`rare_variants.txt` and `signal_traits.txt`.

## Step by step

These are the commands `run_example.sh` runs. Here `OUT=betanfer_example`.

### 1. Generate the data

```bash
python example/make_example_data.py --out betanfer_example
```

### 2. Index each trait file

`betanfer index` sorts, BGZF-compresses and tabix-indexes one file. The
`CHR` and `POS` columns are detected automatically.

```bash
betanfer index --file betanfer_example/sumstats/trait_001.txt --output-dir betanfer_example/indexed
```

Do this for all 150 files. `run_example.sh` runs them in parallel with
`xargs`; on a cluster you would normally use a job array.

### 3. Infer the LD matrix

```bash
betanfer infer --dir betanfer_example/indexed --chrom 22 --start 20000000 --end 20200000 --pval-threshold 5e-8 --min-traits 100 --min-pair-traits 2 --plink2-format --output betanfer_example/ld_minpair2
```

- `--pval-threshold 5e-8` drops any trait with a genome-wide significant
  association in the region, since real effects distort the beta
  correlations. The log should show `Retained 140 traits; filtered out 10 traits`.
- `--min-traits 100` skips the locus unless at least 100 traits remain. The
  default of 1000 is meant for real data; this example has only 150 traits.
- `--min-pair-traits 2` reports every variant pair that is observed together
  in at least 2 traits, which here is all of them.
- `--plink2-format` writes `ld_minpair2.vcor1` (the matrix) and
  `ld_minpair2.vcor1.vars` (variant IDs in matrix order).

### 4. Filter poorly supported pairs

Each rare variant is observed in only about 70 of the 140 traits, so its
LD estimates are based on fewer traits and are noisier. `--min-pair-traits`
sets pairs observed together in fewer than N traits to NaN:

```bash
betanfer infer --dir betanfer_example/indexed --chrom 22 --start 20000000 --end 20200000 --pval-threshold 5e-8 --min-traits 100 --min-pair-traits 100 --plink2-format --output betanfer_example/ld_minpair100
```

The log should report `764 / 4950` pairs set to NaN: every pair that involves
one of the 8 rare variants (8 × 92 + 28 = 764). Pairs of common variants are
observed together in at least 127 traits, so none of them are affected.

`--min-pair-traits` must be between 2 and `--min-traits`. Its default is
1000.

### 5. Check the result

```bash
python example/check_example.py --inferred betanfer_example/ld_minpair100 --data betanfer_example --min-pair-traits 100
```

Expected output (values may differ slightly in the last digit):

```
  All 100 variants present                                           PASS  (100 variants)
  Inferred vs true LD: r >= 0.90                                     PASS  (r = 0.966)
  Mean absolute error <= 0.10                                        PASS  (MAE = 0.063)
  NaN exactly for pairs with a rare variant (--min-pair-traits 100)  PASS  (764 / 4950 NaN, 764 expected)
RESULT: PASS
```

The inferred LD is close to but not exactly the true LD, because 140 traits
is a small sample. The estimates get more precise as the number of traits
grows, which is why real analyses use thousands of GWAS.

## Troubleshooting

- **`betanfer: command not found`**: activate the environment Betanfer is
  installed in, or run with `BETANFER="python -m betanfer"` from the
  repository root.
- **Indexing failed**: see `betanfer_example/index.log`. A
  `ModuleNotFoundError` (e.g. `pysam` or `torch`) means the environment is
  incomplete; see the installation steps in the main README.
- **Device shows `cpu` when you expected a GPU**: the run still passes, but
  check your PyTorch install for your GPU type.
- **A check fails**: send the full output and the logs in
  `betanfer_example/` to the maintainers.
