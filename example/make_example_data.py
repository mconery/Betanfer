#!/usr/bin/env python3
"""
Generate a small synthetic GWAS dataset for testing a Betanfer install.

Uses only the Python standard library and a fixed random seed, so every
user gets byte-identical files.

What is simulated
-----------------
* 100 variants on chromosome 22 (positions 20,001,000-20,199,000, 2 kb apart)
  arranged in three LD blocks. Within a block, adjacent variants have
  r = 0.95 and LD decays geometrically with distance; variants in different
  blocks are uncorrelated. The true LD matrix is written to true_ld.tsv.
* 150 traits. For 140 "null" traits the z-scores are drawn from N(0, R),
  which is the property Betanfer relies on: across many traits with no
  real association, sampling noise in the betas is correlated like LD.
* 10 "signal" traits additionally carry one causal variant (|z| ~ 8-11),
  so `betanfer infer --pval-threshold 5e-8` should remove them.
* 8 rare variants (MAF ~0.3%) have missing summary statistics (NA) in about
  half of the traits, mimicking rare variants lost to QC in many GWAS.
  Common variants are missing in about 2% of traits.

Output layout
-------------
  OUT/sumstats/trait_001.txt ... trait_150.txt   tab-separated, one per trait
      columns: CHR POS SNP A1 A2 MAF BETA SE P
  OUT/true_ld.tsv          true LD (Pearson r) matrix, labelled
  OUT/rare_variants.txt    IDs of the 8 rare variants
  OUT/signal_traits.txt    trait files carrying a real association

Usage:  python make_example_data.py --out DIR
"""

import argparse
import math
import os
import random

SEED = 2026
CHROM = "22"
FIRST_POS = 20_001_000
SPACING = 2_000
BLOCK_SIZES = (40, 35, 25)          # 100 variants
RHO = 0.95                          # LD between adjacent variants in a block
N_NULL_TRAITS = 140
N_SIGNAL_TRAITS = 10
RARE_INDICES = (5, 17, 33, 48, 61, 70, 84, 95)
RARE_MISSING = 0.50                 # fraction of traits missing each rare variant
COMMON_MISSING = 0.02               # fraction of traits missing each common variant


def build_ld(block_sizes, rho):
    """Return (adjacent-rho list, full LD matrix) for a block AR(1) structure."""
    adj = []
    for size in block_sizes:
        adj.extend([0.0] + [rho] * (size - 1))   # 0 at the start of each block
    n = len(adj)
    R = [[1.0] * n for _ in range(n)]
    for i in range(n):
        r = 1.0
        for k in range(i + 1, n):
            r *= adj[k]
            R[i][k] = R[k][i] = r
    return adj, R


def draw_correlated_z(rng, adj):
    """Draw one z vector ~ N(0, R) using the Markov (AR(1)) structure."""
    z = []
    prev = 0.0
    for a in adj:
        val = a * prev + math.sqrt(1.0 - a * a) * rng.gauss(0.0, 1.0)
        z.append(val)
        prev = val
    return z


def pvalue(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--out", required=True, help="Output directory")
    args = ap.parse_args()

    rng = random.Random(SEED)
    adj, R = build_ld(BLOCK_SIZES, RHO)
    n_var = len(adj)

    # Variant annotation
    bases = "ACGT"
    variants = []
    for i in range(n_var):
        pos = FIRST_POS + i * SPACING
        a1 = rng.choice(bases)
        a2 = rng.choice([b for b in bases if b != a1])
        maf = rng.uniform(0.002, 0.004) if i in RARE_INDICES else rng.uniform(0.05, 0.5)
        variants.append({"pos": pos, "a1": a1, "a2": a2, "maf": maf,
                         "id": f"{CHROM}:{pos}:{a1}:{a2}"})

    # Signal traits are interleaved at random positions among all traits
    n_traits = N_NULL_TRAITS + N_SIGNAL_TRAITS
    signal_slots = set(rng.sample(range(n_traits), N_SIGNAL_TRAITS))
    common_idx = [i for i in range(n_var) if i not in RARE_INDICES]

    sumstats_dir = os.path.join(args.out, "sumstats")
    os.makedirs(sumstats_dir, exist_ok=True)

    signal_files = []
    for t in range(n_traits):
        name = f"trait_{t + 1:03d}.txt"
        n_samples = rng.uniform(50_000, 200_000)
        z = draw_correlated_z(rng, adj)
        if t in signal_slots:
            causal = rng.choice(common_idx)
            effect = rng.uniform(8.0, 11.0) * rng.choice((-1.0, 1.0))
            z = [z[i] + effect * R[i][causal] for i in range(n_var)]
            signal_files.append(name)

        lines = ["CHR\tPOS\tSNP\tA1\tA2\tMAF\tBETA\tSE\tP"]
        for i, v in enumerate(variants):
            miss = RARE_MISSING if i in RARE_INDICES else COMMON_MISSING
            prefix = f"{CHROM}\t{v['pos']}\t{v['id']}\t{v['a1']}\t{v['a2']}\t{v['maf']:.4f}"
            if rng.random() < miss:
                lines.append(f"{prefix}\tNA\tNA\tNA")
                continue
            se = 1.0 / math.sqrt(2.0 * v["maf"] * (1.0 - v["maf"]) * n_samples)
            lines.append(f"{prefix}\t{z[i] * se:.6g}\t{se:.6g}\t{pvalue(z[i]):.4g}")
        with open(os.path.join(sumstats_dir, name), "w") as fh:
            fh.write("\n".join(lines) + "\n")

    ids = [v["id"] for v in variants]
    with open(os.path.join(args.out, "true_ld.tsv"), "w") as fh:
        fh.write("SNP\t" + "\t".join(ids) + "\n")
        for vid, row in zip(ids, R):
            fh.write(vid + "\t" + "\t".join(f"{x:.6g}" for x in row) + "\n")
    with open(os.path.join(args.out, "rare_variants.txt"), "w") as fh:
        fh.write("\n".join(ids[i] for i in RARE_INDICES) + "\n")
    with open(os.path.join(args.out, "signal_traits.txt"), "w") as fh:
        fh.write("\n".join(sorted(signal_files)) + "\n")

    print(f"Wrote {n_traits} trait files ({n_var} variants each) to {sumstats_dir}")
    print(f"True LD matrix: {os.path.join(args.out, 'true_ld.tsv')}")


if __name__ == "__main__":
    main()
