![Alt text](logo.png)

**Betanfer** infers linkage disequilibrium (LD) structure from GWAS summary statistics by computing pairwise correlations between variant effect sizes (betas) across many traits — no reference genotype panel required.

---

## How it works

1. **Index** your GWAS summary statistics files with tabix so they can be queried efficiently by genomic position.
2. **Infer** LD by extracting the effect estimate (beta) for each target variant from every indexed trait, then computing the pairwise Pearson correlation matrix across traits. High beta-correlation between two variants indicates they are in LD.

---

## Installation

Betanfer uses PyTorch for GPU-accelerated correlation computation. How you install it depends on your system:

### Step 1 — Install betanfer and non-pytorch dependencies

**With Conda (recommended approach)**

```bash
git clone https://github.com/your-org/betanfer.git
cd betanfer

conda env create -f environment.yml --prefix /full/path/to/env/betanfer
conda activate /full/path/to/env/betanfer
```

### Step 2 — Install pytorch

Many HPC managers provide instructions on how to use pytorch on their systems and we have tested a few systems as well. Here is some guidance for installing pytorch on systems with different types of GPUs:

| System | Instructions |
|---|---|
| Intel XPU (ALCF Aurora) | `pip install   torch==2.9.1+xpu   torchvision==0.24.1+xpu   torchaudio==2.9.1+xpu   --index-url https://download.pytorch.org/whl/xpu   --extra-index-url https://pypi.org/simple` |
| AMD ROCm (OLCF Frontier) | [PyTorch on OLCF Frontier](https://docs.olcf.ornl.gov/software/analytics/pytorch_frontier.html#installing-pytorch) - Make sure to reactivate environment after loading rocm module |
| NVIDIA CUDA (Most Systems) | `pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu###` - run with Cuda version appropriate to your system |

### Step 3 - Verify installation

```bash
betanfer --version
betanfer --help
```
If you intend to use GPUs for the correlation calculations we recommend you check the log file the first time you use the infer command to make sure the code hasn't defaulted to the CPU-based implementation.

---

## Usage

Betanfer has two subcommands: `index` and `infer`.

---

### `betanfer index`

Index a **single** GWAS summary statistics file.  Scale to many files with a
shell loop or an HPC job array:

```bash
betanfer index --help

# Auto-detect chromosome and position column names
betanfer index --file trait1.txt --output-dir indexed/

# Specify column names explicitly
betanfer index --file trait1.txt --chrom-col CHR --pos-col POS --output-dir indexed/

# Parse positions from variant IDs (e.g. rsid format "CHR:POS:REF:ALT")
betanfer index --file trait1.txt \
    --variant-col SNP --variant-id-format CHR:POS:REF:ALT \
    --output-dir indexed/

# Use a PLINK .bim file to map variant IDs to positions
betanfer index --file trait1.txt \
    --variant-col SNP --bim-file reference.bim \
    --output-dir indexed/

# Re-index (overwrite existing .gz / .tbi files)
betanfer index --file trait1.txt --chrom-col CHR --pos-col POS --force
```

**Column name auto-detection**

If you do not specify column names, Betanfer reads the header of your file and
looks for common column names (case-insensitive):

| Type       | Names checked                                              |
|------------|------------------------------------------------------------|
| Chromosome | `chr`, `chrom`, `chromosome`, `hm_chrom`, `contig`        |
| Position   | `pos`, `bp`, `position`, `base_pair_location`, `hm_pos`, `base_pair`, `bp_hg38`, `bp_hg19` |

If user-specified names are provided, they are validated against the file header
and an error is raised if they are not found.

**Key flags**

| Flag | Description |
|---|---|
| `--file SUMSTATS` | Input summary statistics file (required) |
| `--output-dir DIR` | Where to write `.gz` and `.tbi` files (default: same dir as `--file`) |
| `--chrom-col COL` | Chromosome column name |
| `--pos-col COL` | Position column name |
| `--variant-col COL` | Variant ID column (used with `--variant-id-format` or `--bim-file`) |
| `--variant-id-format FMT` | Format string, e.g. `CHR:POS:REF:ALT` |
| `--bim-file BIM [BIM...]` | One or more PLINK `.bim` files for position mapping |
| `--force` | Overwrite existing index files |

---

### `betanfer infer`

Compute a pairwise LD matrix from indexed summary statistics.

```bash
betanfer infer --help

# Variant file (contiguous region) — GPU used automatically if available
betanfer infer \
    --dir indexed/ \
    --variant-file variants.txt \
    --output ld_matrix.tsv

# Genomic range
betanfer infer \
    --dir indexed/ \
    --chrom 1 --start 1000000 --end 2000000 \
    --output ld_matrix.tsv

# From a text file listing indexed paths
betanfer infer \
    --file-list indexed_files.txt \
    --variant-file variants.txt \
    --output ld_matrix.tsv

# P-value filtering: exclude traits with significant associations in the region
betanfer infer \
    --dir indexed/ \
    --variant-file variants.txt \
    --output ld_matrix.tsv \
    --pval-threshold 5e-8 \
    --buffer 50000 \
    --stats-output filter_stats.tsv

# Save binary output (.vcor1 + .vcor1.vars) instead of TSV
betanfer infer \
    --dir indexed/ \
    --variant-file variants.txt \
    --output ld_prefix \
    --save-binary

# Non-contiguous variants; force CPU; 8 parallel processes for filtering
betanfer infer \
    --dir indexed/ \
    --variant-file scattered_variants.txt \
    --output ld_matrix.tsv \
    --non-contiguous --cpu --num-processes 8
```

**Variant specification**

Variants can be specified in two mutually exclusive ways:

| Method | Flags | Notes |
|---|---|---|
| Variant file | `--variant-file FILE` | One variant ID per line. Assumed contiguous unless `--non-contiguous` is set. |
| Genomic range | `--chrom CHR --start POS --end POS` | All variants found in the range are used. |

**Column name auto-detection**

Betanfer reads the header of the first indexed file and looks for:

| Type       | Names checked                                                     |
|------------|-------------------------------------------------------------------|
| Variant ID | `rs_number`, `rsid`, `snp`, `variant_id`, `variant`, `id`, `markerid` |
| Beta       | `beta`, `effect`, `effect_size`, `b`, `beta_hat`, `log_or`       |
| P-value    | `p_value`, `pvalue`, `p.value`, `pval`, `p`, `log10p`            |

**P-value filtering**

When `--pval-threshold` or `--pval-col` is set, traits with any variant at
`p < threshold` in the target region (optionally extended by `--buffer` bp)
are excluded before computing the LD matrix.  This avoids inflating the 
correlation estimates.

The filter is parallelised across CPU cores (see `--num-processes`).

**Key flags**

| Flag | Default | Description |
|---|---|---|
| `--dir DIR` | — | Directory of indexed files (mutually exclusive with `--file-list`) |
| `--file-list FILE` | — | Text file listing indexed paths |
| `--variant-file FILE` | — | Variant IDs, one per line |
| `--chrom / --start / --end` | — | Genomic range instead of variant file |
| `--output FILE` | — | Output path (required) |
| `--variant-col COL` | auto | Variant ID column name |
| `--variant-id-format FMT` | `CHR:POS:REF:ALT` | Variant ID format for CHR/POS parsing |
| `--beta-col COL` | auto | Beta column name |
| `--cpu` | off | Force CPU (default: use GPU if available) |
| `--non-contiguous` | off | Disable batch region fetch (use per-variant queries) |
| `--pval-col COL` | auto | P-value column for trait filtering |
| `--pval-threshold FLOAT` | 5e-8 (when active) | Significance threshold for filtering |
| `--buffer BP` | 0 | Extend p-value check region by this many bp |
| `--stats-output FILE` | — | Per-trait filtering statistics (TSV) |
| `--num-processes N` | all CPUs | Parallel processes for p-value filtering |
| `--save-binary` | off | Save `.vcor1` + `.vcor1.vars` instead of TSV |

---

## Output formats

| Format | Description |
|---|---|
| TSV (default) | Tab-separated matrix with variant IDs as row and column labels |
| Plink2 (`.vcor1`) | Raw float32 matrix without row or column labels consistent with PLINK2 `.vcor1` format |
| `.vcor1.vars` | Variant ID list matching the matrix order corresponding to PLINK2 format |

---

## Notes

- **GPU acceleration**: Correlation computation uses PyTorch. GPUs are used automatically when available. Pass `--cpu` to force CPU execution.
- **Input format**: Summary statistics files should be sorted by chromosome and position (required by tabix). The indexing step sorts automatically.
- **Header format**: Files with a `#`-prefixed header line (VCF-style) are handled correctly.
- **Missing data**: Variants absent from a particular file are masked for correlation computation.

---

## Citation

If you use Betanfer in your research, please cite: *(citation forthcoming)*
