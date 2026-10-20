#!/usr/bin/env python3
'''
match_credible_sets_to_abf.py

Scores how well each new-method (SuSiEx/MultiSuSiE/RSparsePro) credible set
overlaps the historical Wakefield-ABF (GWAMA/MR-MEGA) credible set it is
meant to replace, using a PIP-weighted Jaccard overlap, and flags each
new-method CS as "shared" (confirmed against the ABF baseline) or not.

Five comparison scopes (pop, tool):
  (EUR, rsparsepro), (AFR, rsparsepro), (ALL, rsparsepro),
  (ALL, susiex), (ALL, multisusie)

ABF credible sets come straight from the raw per-trait GWAMA (EUR/AFR) or
MR-MEGA (ALL) output: the row's PIP_1 is the weight and CS95_1 flags
membership (the ABF arm only ever has one credible set, 'cs' == '1').

New-method credible sets and weights come from this step's own collation
output (--credible-sets-dir, same directory collate_credible_sets.py
writes to):
  - rsparsepro / multisusie: membership from gwas_cs_members.tsv, weight is
    that tool's PIP from gwas_finemapped_variants.tsv.
  - susiex (ALL only): membership and weight (pip_marginal, the marginal
    inclusion probability across all of a locus's components) both come
    from susiex_native_weights.tsv.gz -- SuSiEx's native 95%-coverage
    credible sets, not the single-assignment wrapper `cs` column.

Trait scope per pop: a trait is compared only if it has >=1 ABF credible
set (CS95_1 == TRUE somewhere) AND >=1 new-method credible set (any tool,
from gwas_cs_members.tsv) for that pop. This reproduces the historical
"fixed shared trait scope" without needing a separate static trait list:
it is exactly legacy-ABF-positive AND new-collation-positive, both
computed here from the two inputs directly.

Overlap score (PIP-weighted Jaccard), for an ABF CS with weights w_a and a
new-method CS with weights w_b over the union of their members:
  score = sum(min(w_a, w_b)) / sum(max(w_a, w_b))
Each new-method CS is scored against its single best-overlapping ABF CS
(within the same trait); score > 0.10 is reported as "shared".

Output (written to --credible-sets-dir):
  cs_match_status_native.tsv.gz  one row per (pop, tool, trait, locus, cs):
                                  best_pip_weighted_jaccard, shared

Usage:
  python match_credible_sets_to_abf.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import os
import getopt
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

# Same directory collate_credible_sets.py writes gwas_cs_members.tsv,
# gwas_finemapped_variants.tsv and susiex_native_weights.tsv.gz to.
CREDIBLE_SETS_DIR_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/credible_sets"
# Raw ABF (Wakefield) fine-mapping output: GWAMA per-ancestry, MR-MEGA cross-ancestry.
ABF_GWAMA_DIR_DEFAULT = "/lus/flare/projects/GeomicVar/jianhuawang/betainfer/data/raw/cs_gwama"
ABF_MRMEGA_DIR_DEFAULT = "/lus/flare/projects/GeomicVar/jianhuawang/betainfer/data/raw/cs_mrmega"

PAIRS = (("EUR", "rsparsepro"), ("AFR", "rsparsepro"), ("ALL", "rsparsepro"),
         ("ALL", "susiex"), ("ALL", "multisusie"))
SHARED_THRESHOLD = 0.10


def help(exit_num=0):
    print("""Usage: python match_credible_sets_to_abf.py [OPTIONS]
Optional:
  --credible-sets-dir  This step's own collation output dir (default: {csd})
  --abf-gwama-dir       Raw GWAMA per-ancestry CS directory (default: {gw})
  --abf-mrmega-dir      Raw MR-MEGA cross-ancestry CS directory (default: {mm})
  --workers             Parallel workers for ABF file scanning (default: min(16, cpu_count))
  -h, --help            Show this help message""".format(
        csd=CREDIBLE_SETS_DIR_DEFAULT, gw=ABF_GWAMA_DIR_DEFAULT, mm=ABF_MRMEGA_DIR_DEFAULT))
    sys.exit(exit_num)


###############################################################################
##########  ABF (WAKEFIELD) CREDIBLE SETS  ####################################
###############################################################################

def raw_abf(path):
    """Read one raw GWAMA/MR-MEGA CS file and return its 95% ABF credible
    set (there is always exactly one, 'cs' == '1'): trait, pop, locus, cs,
    variant_id, pip."""
    cols = ["#PHENO", "POP", "LOCUS", "ID", "PIP_1", "CS95_1"]
    d = pd.read_csv(path, sep="\t", usecols=cols)
    d = d[d.CS95_1.astype(str).str.upper().eq("TRUE")].copy()
    d = d.rename(columns={"#PHENO": "trait", "POP": "pop", "LOCUS": "locus",
                           "ID": "variant_id", "PIP_1": "pip"})
    d["cs"] = "1"
    return d[["trait", "pop", "locus", "cs", "variant_id", "pip"]]


def _abf_path(abf_gwama_dir, abf_mrmega_dir, pop, trait):
    if pop == "ALL":
        return abf_mrmega_dir / f"{trait}.ALL.cs.txt.gz"
    return abf_gwama_dir / f"{trait}.{pop}.gwama.cs.txt.gz"


###############################################################################
##########  NEW-METHOD (SUSIEX/MULTISUSIE/RSPARSEPRO) CREDIBLE SETS  #########
###############################################################################

def _new_union_traits(credible_sets_dir, pops):
    """Traits with >=1 credible-set member row for each pop, across any
    tool, from this step's own gwas_cs_members.tsv (SuSiEx membership in
    there is already native -- see collate_credible_sets.py)."""
    out = {p: set() for p in pops}
    path = credible_sets_dir / "gwas_cs_members.tsv"
    for chunk in pd.read_csv(path, sep="\t", usecols=["trait", "pop"], chunksize=1_000_000):
        for pop, traits in chunk.groupby("pop")["trait"]:
            if pop in out:
                out[pop].update(traits.unique())
    return out


def _load_members_with_pip(credible_sets_dir, tools):
    """rsparsepro/multisusie credible-set members with that tool's PIP as
    the matching weight (gwas_cs_members.tsv joined to
    gwas_finemapped_variants.tsv)."""
    keys = ["trait", "locus", "pop", "variant_id"]
    members = pd.read_csv(credible_sets_dir / "gwas_cs_members.tsv", sep="\t", dtype={"cs": str},
                           usecols=["trait", "locus", "pop", "tool", "cs", "variant_id"])
    members = members[members.tool.isin(tools)]
    pip_cols = [f"{t}_pip" for t in tools]
    pips = pd.read_csv(credible_sets_dir / "gwas_finemapped_variants.tsv", sep="\t",
                        usecols=keys + pip_cols)
    frames = []
    for tool in tools:
        subset = members[members.tool.eq(tool)]
        weights = pips[keys + [f"{tool}_pip"]].rename(columns={f"{tool}_pip": "pip"})
        d = subset.merge(weights, on=keys, how="left", validate="many_to_one")
        frames.append(d)
    return pd.concat(frames, ignore_index=True)


def _load_susiex_native(credible_sets_dir):
    """SuSiEx native (95%-coverage) credible-set members with the marginal
    inclusion probability (1 - prod(1 - alpha), across all of a locus's
    components) as the matching weight. ALL pop only."""
    d = pd.read_csv(credible_sets_dir / "susiex_native_weights.tsv.gz", sep="\t", dtype={"cs": str})
    d = d.rename(columns={"pip_marginal": "pip"})
    d["tool"] = "susiex"
    return d


###############################################################################
##########  PIP-WEIGHTED JACCARD MATCHING  ####################################
###############################################################################

def validate_weights(d):
    """Matching weights must be real numbers in [0, 1] and never missing --
    a missing weight silently treated as zero would bias every score."""
    assert d.pip.notna().all(), "Missing PIP must not be replaced with zero"
    assert np.isfinite(d.pip).all() and d.pip.between(0, 1).all()


def families(d):
    """One {variant_id: pip} weight dict per (locus, cs) credible set."""
    output = {}
    for key, rows in d.groupby(["locus", "cs"], sort=False):
        weights = dict(zip(rows.variant_id, rows.pip))
        assert sum(weights.values()) > 0
        output[key] = weights
    return output


def weighted_edges(a, b):
    """PIP-weighted Jaccard overlap between every credible set in a and
    every credible set in b that shares >=1 variant: (key_a, key_b, score)."""
    inverted = defaultdict(list)
    masses_a = {k: sum(v.values()) for k, v in a.items()}
    for key, weights in a.items():
        for variant, weight in weights.items():
            inverted[variant].append((key, weight))
    edges = []
    for kb, weights in b.items():
        shared = defaultdict(float)
        for variant, weight in weights.items():
            for ka, wa in inverted.get(variant, ()):
                shared[ka] += min(wa, weight)
        for ka, mass in shared.items():
            union_mass = masses_a[ka] + sum(weights.values()) - mass
            score = mass / union_mass
            assert -1e-8 <= score <= 1 + 1e-8
            edges.append((ka, kb, float(np.clip(score, 0, 1))))
    return edges


def _self_check():
    """Lightweight inline check of weighted_edges against hand-worked and
    randomized PIP-weighted-Jaccard cases. No file output."""
    key = ("locus", "1")
    score = lambda a, b: weighted_edges({key: a}, {key: b})[0][2]
    assert np.isclose(score({"a": .9, "b": .1}, {"a": .1, "b": .9}), 1 / 9)
    assert np.isclose(score({"a": .9, "b": .1}, {"a": .9, "b": .1}), 1)
    assert np.isclose(score({"a": 1, "b": 1}, {"a": 1}), .5)
    assert weighted_edges({key: {"a": 1}}, {key: {"b": 1}}) == []
    rng = np.random.default_rng(911)
    for _ in range(50):
        x, y = rng.random((2, 10))
        assert np.isclose(score(dict(enumerate(x)), dict(enumerate(y))),
                           np.minimum(x, y).sum() / np.maximum(x, y).sum())


###############################################################################
##########  MAIN  ##############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "credible-sets-dir=", "abf-gwama-dir=", "abf-mrmega-dir=", "workers="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    credible_sets_dir = CREDIBLE_SETS_DIR_DEFAULT
    abf_gwama_dir = ABF_GWAMA_DIR_DEFAULT
    abf_mrmega_dir = ABF_MRMEGA_DIR_DEFAULT
    workers = min(16, os.cpu_count() or 4)

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--credible-sets-dir":
            credible_sets_dir = a
        elif o == "--abf-gwama-dir":
            abf_gwama_dir = a
        elif o == "--abf-mrmega-dir":
            abf_mrmega_dir = a
        elif o == "--workers":
            workers = int(a)
        else:
            assert False, f"Unhandled option: {o}"

    credible_sets_dir = Path(credible_sets_dir)
    abf_gwama_dir = Path(abf_gwama_dir)
    abf_mrmega_dir = Path(abf_mrmega_dir)

    _self_check()
    print("NOTE: weighted_edges self-check passed")

    pops = sorted({p for p, _ in PAIRS})
    new_union = _new_union_traits(credible_sets_dir, pops)

    keep, abf_frames = {}, {}
    for pop in pops:
        candidates = sorted(new_union.get(pop, set()))
        paths = [p for t in candidates
                 if (p := _abf_path(abf_gwama_dir, abf_mrmega_dir, pop, t)).exists()]
        print(f"NOTE: {pop}: {len(candidates)} new-collation traits, {len(paths)} ABF files to scan")
        positive_frames = []
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for frame in ex.map(raw_abf, paths, chunksize=8):
                if len(frame):
                    positive_frames.append(frame)
        keep[pop] = {t for f in positive_frames for t in f.trait.unique()}
        abf_frames[pop] = pd.concat(positive_frames, ignore_index=True) if positive_frames else \
            pd.DataFrame(columns=["trait", "pop", "locus", "cs", "variant_id", "pip"])
        validate_weights(abf_frames[pop])
        print(f"NOTE: {pop}: {len(keep[pop])} traits kept (ABF CS95-positive and new-collation-positive)")

    non_susiex_tools = sorted({t for _, t in PAIRS if t != "susiex"})
    members = _load_members_with_pip(credible_sets_dir, non_susiex_tools)
    susiex = _load_susiex_native(credible_sets_dir)
    new_all = pd.concat([members, susiex], ignore_index=True)
    validate_weights(new_all)

    status = []
    for pop, tool in PAIRS:
        abf = abf_frames[pop]
        aa = {tr: families(d) for tr, d in abf.groupby("trait")}
        new_pop_tool = new_all[new_all["pop"].eq(pop) & new_all.tool.eq(tool)]
        bb = {tr: families(d) for tr, d in new_pop_tool.groupby("trait")}
        for trait in sorted(keep.get(pop, ())):
            a, b = aa.get(trait, {}), bb.get(trait, {})
            edges = weighted_edges(a, b)
            best = {k: 0. for k in b}
            for ka, kb, score in edges:
                best[kb] = max(best[kb], score)
            for key, score in best.items():
                status.append(dict(pop=pop, tool=tool, trait=trait, locus=key[0], cs=key[1],
                                    best_pip_weighted_jaccard=score, shared=score > SHARED_THRESHOLD))
        print(f"NOTE: {pop}/{tool}: {len(keep.get(pop, ()))} traits scored")

    out_path = credible_sets_dir / "cs_match_status_native.tsv.gz"
    pd.DataFrame(status, columns=["pop", "tool", "trait", "locus", "cs",
                                   "best_pip_weighted_jaccard", "shared"]).to_csv(
        out_path, sep="\t", index=False)
    print(f"NOTE: {len(status)} status rows -> {out_path}")


if __name__ == "__main__":
    main()
