'''
select_causal_variants.py

Pick user-defined sets of causal variants from chromosome-specific PLINK BIM
files while enforcing distance constraints.  Each set is written to a
tab-delimited file suitable for downstream GWAS simulation.

***Run in Conda finemap environment***
'''

import argparse
import glob
import os
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import numpy as np


# -----------------------------------------------------------------------------#
# Helper functions
# -----------------------------------------------------------------------------#
def load_bim_files(input_dir: Path, pattern: str) -> dict[int, pd.DataFrame]:
    """
    Read BIM files matching the pattern and return a dict keyed by chromosome.
    """
    files = glob.glob(str(input_dir / pattern))
    if not files:
        sys.exit(f"No BIM files matching pattern '{pattern}' found in {input_dir}")

    data: dict[int, pd.DataFrame] = {}
    for f in files:
        # Extract chromosome from filename by replacing '*' placeholder
        # User may provide pattern with '*' => replace with (\d+) regex group
        chrom = infer_chromosome_from_filename(f)
        df = pd.read_csv(
            f,
            sep=r"\s+",
            header=None,
            names=[
                "CHR",
                "SNP",
                "CM",
                "BP",
                "A1",
                "A2",
            ],
            dtype={"CHR": int, "SNP": str, "BP": int},
        )
        data[chrom] = df.sort_values("BP").reset_index(drop=True)
    return data


def infer_chromosome_from_filename(path: str) -> int:
    """
    Extract integer chromosome from filename (assumes 'chr<number>' in name).
    """
    import re

    m = re.search(r"chr(\d+)", os.path.basename(path))
    if not m:
        sys.exit(f"Cannot infer chromosome number from file name '{path}'")
    return int(m.group(1))


def pick_variant_set(
    df: pd.DataFrame,
    k: int,
    max_span: int,
    min_spacing: int,
    used_intervals: List[Tuple[int, int]],
    rng: random.Random,
    max_trials: int = 1000,
) -> Tuple[List[int], List[str]]:
    """
    Attempt to choose k variants from df that satisfy intra-set constraints and
    inter-set separation (no overlap with used_intervals). Returns positions and IDs.
    Raises ValueError if impossible after max_trials.
    """
    positions = df["BP"].to_numpy()
    ids = df["SNP"].to_numpy()
    n = len(df)

    for _ in range(max_trials):
        idx_anchor = rng.randrange(n)
        anchor_pos = positions[idx_anchor]
        window_end = anchor_pos + max_span

        # Candidates within window
        mask = (positions >= anchor_pos) & (positions <= window_end)
        cand_pos = positions[mask]
        cand_ids = ids[mask]
        if len(cand_pos) < k:
            continue

        # Greedy random spacing
        chosen_pos = [anchor_pos]
        chosen_ids = [ids[idx_anchor]]

        for _ in range(k - 1):
            valid = [
                (p, s)
                for p, s in zip(cand_pos, cand_ids)
                if all(abs(p - cp) >= min_spacing for cp in chosen_pos)
            ]
            if not valid:
                break
            p, s = rng.choice(valid)
            chosen_pos.append(p)
            chosen_ids.append(s)

        if len(chosen_pos) != k:
            continue

        # Check inter-set separation
        interval_start = min(chosen_pos)
        interval_end = max(chosen_pos)
        if any(
            not (interval_end + max_span < us or interval_start - max_span > ue)
            for us, ue in used_intervals
        ):
            continue

        used_intervals.append((interval_start, interval_end))
        return chosen_pos, chosen_ids

    raise ValueError("Failed to select variant set given the constraints.")


def write_set_file(
    out_dir: Path,
    pheno_idx: int,
    k: int,
    effect: float,
    max_span: int,
    min_spacing: int,
    variant_ids: List[str],
) -> None:
    """
    Write one file containing causal variants and effect sizes.
    """
    filename = (
        f"Pheno{pheno_idx}.{k}_causal."
        f"effect_size_{effect:g}.{max_span}_space{min_spacing}.txt"
    )
    out_path = out_dir / filename
    with open(out_path, "w") as fout:
        for vid in variant_ids:
            fout.write(f"{vid}\t{effect}\n")


# -----------------------------------------------------------------------------#
# Main
# -----------------------------------------------------------------------------#
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Select causal variant sets from PLINK BIM files."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        type=Path,
        help="Directory containing BIM files.",
    )
    parser.add_argument(
        "--naming-convention",
        default="chr*.plink.maf_0.05.shared.bim",
        help="Filename pattern with '*' as chromosome placeholder (glob syntax).",
    )
    parser.add_argument(
        "--num-sets",
        type=int,
        default=1000,
        help="Number of phenotype (causal) sets to create.",
    )
    parser.add_argument(
        "--variants-per-set",
        type=int,
        default=1,
        help="Number of causal variants in each set.",
    )
    parser.add_argument(
        "--effect-size",
        type=float,
        default=0.1,
        help="Effect size assigned to every causal variant.",
    )
    parser.add_argument(
        "--max-distance",
        type=int,
        default=500_000,
        help="Maximum genomic span (in bp) of variants within a set.",
    )
    parser.add_argument(
        "--min-spacing",
        type=int,
        default=100_000,
        help="Minimum distance (in bp) between any two variants within a set.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory to write output files.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducibility.",
    )

    args = parser.parse_args()

    # ---------------------------------------------------------------------- #
    # Validation
    # ---------------------------------------------------------------------- #
    if args.min_spacing > args.max_distance:
        sys.exit("min-spacing cannot exceed max-distance.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    # ---------------------------------------------------------------------- #
    # Load BIM data
    # ---------------------------------------------------------------------- #
    bim_data = load_bim_files(args.input_dir, args.naming_convention)

    # Flatten variant index per chromosome for weighted sampling
    chrom_weights = {c: len(df) for c, df in bim_data.items()}
    chroms = list(chrom_weights.keys())
    weights = np.array([chrom_weights[c] for c in chroms], dtype=float)
    weights /= weights.sum()

    used_intervals_by_chr: dict[int, List[Tuple[int, int]]] = defaultdict(list)

    # ---------------------------------------------------------------------- #
    # Main selection loop
    # ---------------------------------------------------------------------- #
    for pheno_idx in range(1, args.num_sets + 1):
        # Pick chromosome proportional to variant count
        chrom = rng.choices(chroms, weights=weights, k=1)[0]

        try:
            positions, vid_list = pick_variant_set(
                bim_data[chrom],
                args.variants_per_set,
                args.max_distance,
                args.min_spacing,
                used_intervals_by_chr[chrom],
                rng,
            )
        except ValueError as e:
            sys.exit(f"Phenotype {pheno_idx}: {e}")

        write_set_file(
            args.output_dir,
            pheno_idx,
            args.variants_per_set,
            args.effect_size,
            args.max_distance,
            args.min_spacing,
            vid_list,
        )

        if pheno_idx % 50 == 0 or pheno_idx == args.num_sets:
            print(f"Completed {pheno_idx}/{args.num_sets} sets", file=sys.stderr)


if __name__ == "__main__":
    main()
