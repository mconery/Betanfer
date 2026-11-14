#!/usr/bin/env python3
"""
Check a Betanfer LD matrix from the example dataset against the true LD.

Usage:
  python check_example.py --inferred PREFIX --data DIR [--min-pair-traits N]

  PREFIX  output prefix given to `betanfer infer --plink2-format --output`
          (reads PREFIX.vcor1 and PREFIX.vcor1.vars)
  DIR     directory written by make_example_data.py
  N       the --min-pair-traits value used for the run (default 2)

Checks
  1. All 100 simulated variants are in the matrix.
  2. Inferred LD tracks the true LD: Pearson r >= 0.90 and mean absolute
     error <= 0.10 over all reported (non-NaN) variant pairs.
  3. NaN pattern matches --min-pair-traits:
       N = 2   -> no NaN entries at all
       N > 2   -> every pair involving a rare variant is NaN, all other
                  pairs are reported (holds for this dataset with N = 100)

Exits with status 0 if all checks pass, 1 otherwise. Standard library only.
"""

import argparse
import math
import os
import sys

MIN_R = 0.90
MAX_MAE = 0.10


def read_matrix(prefix):
    with open(prefix + ".vcor1.vars") as fh:
        ids = [line.strip() for line in fh if line.strip()]
    rows = []
    with open(prefix + ".vcor1") as fh:
        for line in fh:
            if line.strip():
                rows.append([float(x) for x in line.split()])
    return ids, rows


def read_truth(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        truth = {}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            truth[parts[0]] = dict(zip(header, (float(x) for x in parts[1:])))
    return header, truth


def pearson(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    return sxy / math.sqrt(sxx * syy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inferred", required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--min-pair-traits", type=int, default=2)
    args = ap.parse_args()

    ids, M = read_matrix(args.inferred)
    true_ids, truth = read_truth(os.path.join(args.data, "true_ld.tsv"))
    with open(os.path.join(args.data, "rare_variants.txt")) as fh:
        rare = {line.strip() for line in fh if line.strip()}

    results = []

    # 1. Variant set
    missing = set(true_ids) - set(ids)
    results.append(("All %d variants present" % len(true_ids), not missing,
                    "missing %d" % len(missing) if missing else "%d variants" % len(ids)))

    # Collect off-diagonal pairs
    inferred, expected = [], []
    nan_pairs, nan_with_rare, rare_pairs, rare_pairs_nan = 0, 0, 0, 0
    n_pairs = 0
    for i in range(len(ids)):
        for k in range(i + 1, len(ids)):
            n_pairs += 1
            has_rare = ids[i] in rare or ids[k] in rare
            val = M[i][k]
            if has_rare:
                rare_pairs += 1
            if math.isnan(val):
                nan_pairs += 1
                nan_with_rare += has_rare
                rare_pairs_nan += has_rare
                continue
            inferred.append(val)
            expected.append(truth[ids[i]][ids[k]])

    # 2. Accuracy
    if len(inferred) >= 2:
        r = pearson(inferred, expected)
        mae = sum(abs(a - b) for a, b in zip(inferred, expected)) / len(inferred)
    else:
        r, mae = float("nan"), float("nan")
    results.append(("Inferred vs true LD: r >= %.2f" % MIN_R, r >= MIN_R, "r = %.3f" % r))
    results.append(("Mean absolute error <= %.2f" % MAX_MAE, mae <= MAX_MAE, "MAE = %.3f" % mae))

    # 3. NaN pattern
    if args.min_pair_traits <= 2:
        results.append(("No NaN pairs (--min-pair-traits %d)" % args.min_pair_traits,
                        nan_pairs == 0, "%d / %d NaN" % (nan_pairs, n_pairs)))
    else:
        ok = nan_pairs == nan_with_rare and rare_pairs_nan == rare_pairs
        results.append(("NaN exactly for pairs with a rare variant (--min-pair-traits %d)"
                        % args.min_pair_traits, ok,
                        "%d / %d NaN, %d expected" % (nan_pairs, n_pairs, rare_pairs)))

    width = max(len(name) for name, _, _ in results)
    all_ok = True
    for name, ok, detail in results:
        all_ok &= bool(ok)
        print("  %-*s  %s  (%s)" % (width, name, "PASS" if ok else "FAIL", detail))
    print("RESULT: %s" % ("PASS" if all_ok else "FAIL"))
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
