#!/usr/bin/env python3
'''
collate_credible_sets.py

Collates PLATLAS GWAS-arm fine-mapping credible sets, across all three
tools (SuSiEx, MultiSuSiE, RSparsePro) and all three ancestry groups
(AFR, EUR, ALL), into the tidy tables the downstream gene-nomination and
summary-table/figure steps read.

Per-locus fine-mapping output format (one file per trait/locus/pop/tool,
written by run_susiex.py / finemap_SuSiEx.py / run_rsparsepro.py):
  Columns: VARIANT_ID (or RSID), Z, PIP, z_estimated, cs[, cs_all, LogBF_CS{k}_*]
  - MultiSuSiE and RSparsePro: cs is a single integer credible-set id
    (0 = not in any credible set).
  - SuSiEx: credible-set membership follows SuSiEx's own 95%-coverage rule,
    not the single cs column the upstream wrapper writes (that column assigns
    each variant to the first model component whose exported alpha exceeds
    0.01, which is neither 95%-coverage nor SuSiEx's own rule). If the file
    has cs_all (an updated upstream writer encoding native multi-CS
    membership directly as a comma-separated list, '0' = none), that is used.
    Otherwise each component k is rebuilt from its LogBF_CS{k}_<cohort>
    columns: alpha = softmax(sum of cohort LogBF); sort variants by alpha
    descending; keep adding until the cumulative alpha first exceeds 0.95.
    A variant can be a member of more than one component's credible set.
    PIP stays the exported value (sum of alpha across components); only
    credible-set membership is rebuilt.

Unit of work = one (trait, locus, pop) triple: each unit reads three small
per-tool files, so units are processed independently and in parallel.

Output (written to --out-dir):
  gwas_cs_members.tsv          one row per (trait, locus, pop, tool, cs,
                                variant_id) -- a SuSiEx variant belonging to
                                more than one native credible set gets one
                                row per credible set
  gwas_finemapped_variants.tsv one row per (trait, locus, pop, variant_id)
                                that passed the filter (in >=1 credible set,
                                or max PIP across tools >= 0.01), with each
                                tool's PIP, the cross-method median PIP, the
                                max PIP, and whether it is in any credible set
  susiex_native_weights.tsv.gz one row per (trait, locus, pop, cs, variant_id)
                                SuSiEx native credible-set member: alpha (that
                                component's posterior weight) and pip_marginal
                                (1 minus the product over components of
                                (1 - alpha), i.e. the per-variant marginal
                                inclusion probability). Needed by
                                match_credible_sets_to_abf.py because SuSiEx's
                                own exported PIP (sum of alpha) is not a valid
                                probability -- it can exceed 1.

Usage:
  python collate_credible_sets.py [OPTIONS]

***Run in finemap conda environment***
'''

import sys
import os
import csv
import getopt
import gzip
import re
import statistics
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

# SuSiEx and MultiSuSiE results live under the main GWAS fine-mapping results tree.
ROOT_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/fine-mapping_results/GWAS"
# RSparsePro results used here come from a separate rerun on a strict-
# intersection LD panel (the union LD panel gave excessive LD mismatch after
# meta-analysis across cohorts and barely split any credible sets).
RSPARSEPRO_ROOT_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/intersection_finemap"
OUT_DIR_DEFAULT = "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/application_results/credible_sets"

TOOLS = ("susiex", "multisusie", "rsparsepro")
PIP_KEEP = 0.01   # variant-retention threshold, matched to the gene-nomination PIP>=0.01 cutoff
LEVEL = 0.95      # SuSiEx native credible-set coverage target

_CS_LOGBF_RE = re.compile(r"^LogBF_CS(\d+)_(.+)$")


def help(exit_num=0):
    print("""Usage: python collate_credible_sets.py [OPTIONS]
Optional:
  --root             SuSiEx/MultiSuSiE results root
                      (default: {root})
  --rsparsepro-root  RSparsePro results root (default: strict-intersection-LD rerun)
                      (default: {rsp})
  --out-dir          Output directory (default: {out})
  --pops             Comma-separated ancestry groups (default: AFR,EUR,ALL)
  --workers          Parallel workers (default: min(16, cpu_count))
  --limit            Debug: only process the first N (trait,locus,pop) units
  -h, --help         Show this help message""".format(
        root=ROOT_DEFAULT, rsp=RSPARSEPRO_ROOT_DEFAULT, out=OUT_DIR_DEFAULT))
    sys.exit(exit_num)


###############################################################################
##########  PER-LOCUS FILE READING  ###########################################
###############################################################################

def _parse_cs(row, cs_i):
    """Parse the single integer credible-set id for one output row
    (MultiSuSiE / RSparsePro). Returns [] if the variant is in no CS."""
    try:
        v = int(float(row[cs_i]))
    except (ValueError, TypeError, IndexError):
        v = 0
    return [v] if v != 0 else []


def _read_tool(path):
    """Read one MultiSuSiE/RSparsePro output file. Returns
    {variant_id: (pip, cs_ids)}; cs_ids is a list of nonzero credible-set
    ids (see _parse_cs). Missing file / empty file / missing required
    columns all return {}."""
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        try:
            header = next(r)
        except StopIteration:
            return {}
        idx = {c: i for i, c in enumerate(header)}
        id_i = idx.get("VARIANT_ID", idx.get("RSID"))
        pip_i = idx.get("PIP")
        cs_i = idx.get("cs")
        if id_i is None or pip_i is None or cs_i is None:
            return {}
        hi = max(id_i, pip_i, cs_i)
        for row in r:
            if len(row) <= hi:
                continue
            try:
                pip = float(row[pip_i])
            except ValueError:
                pip = float("nan")
            out[row[id_i]] = (pip, _parse_cs(row, cs_i))
    return out


def _softmax(logbf):
    """Softmax over one component's per-variant log Bayes factors (already
    summed across cohorts). Matches collate_native.py exactly."""
    w = np.exp(logbf - logbf.max())
    return w / w.sum()


def _native_cs_indices(alpha):
    """Indices (into alpha) of the smallest alpha-descending prefix whose
    cumulative sum first exceeds LEVEL -- SuSiEx's own model.cpp 95%-coverage
    rule. side='right' is required: it is what makes the cutoff the first
    index where the cumulative sum EXCEEDS (not just reaches) LEVEL."""
    order = np.argsort(-alpha, kind="stable")
    acc = np.cumsum(alpha[order])
    n = min(int(np.searchsorted(acc, LEVEL, side="right")) + 1, len(alpha))
    return order[:n]


def _to_float(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return float("nan")


def _read_susiex_native(path):
    """Read one SuSiEx output file and rebuild native (95%-coverage)
    credible-set membership.

    Returns (pip_cs, weight_rows):
      pip_cs      {variant_id: (pip, cs_ids)} -- exported PIP, native cs ids
      weight_rows [(cs, variant_id, alpha, pip_marginal), ...] -- one row per
                  native credible-set member

    Component membership: if the file has cs_all (an updated writer that
    encodes native multi-CS membership directly as a comma-separated list,
    '0' = none), that is trusted as-is. Otherwise each component k's
    membership is rebuilt from its LogBF_CS{k}_<cohort> columns: alpha =
    softmax(sum of cohort LogBF); native members = _native_cs_indices(alpha).
    This LogBF-rebuild path reproduces collate_native.py exactly.

    alpha / pip_marginal come from the LogBF columns when present (whether
    or not cs_all is also present): alpha is that component's posterior
    weight for the variant, and pip_marginal = 1 - prod_k(1 - alpha_k) is
    the marginal inclusion probability across ALL of the locus's
    components -- not just the ones a variant is a native member of, which
    is why it is computed from the full per-component alpha matrix, exactly
    as collate_native.py does. If a file has cs_all but no LogBF columns,
    membership is still native but alpha/pip_marginal cannot be recovered
    and are left as NaN.
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}, []
    with open(path, newline="") as fh:
        r = csv.reader(fh, delimiter="\t")
        try:
            header = next(r)
        except StopIteration:
            return {}, []
        idx = {c: i for i, c in enumerate(header)}
        id_i = idx.get("VARIANT_ID", idx.get("RSID"))
        pip_i = idx.get("PIP")
        if id_i is None or pip_i is None:
            return {}, []
        rows = list(r)

    n = len(rows)
    vids = [row[id_i] if len(row) > id_i else "" for row in rows]
    pips = [_to_float(row[pip_i]) if len(row) > pip_i else float("nan") for row in rows]

    comp_cols = defaultdict(list)  # k -> [column index, ...]
    for col, i in idx.items():
        m = _CS_LOGBF_RE.match(col)
        if m:
            comp_cols[int(m.group(1))].append(i)

    comp_ks = sorted(comp_cols)
    alpha_by_comp = {}
    for k in comp_ks:
        cols = comp_cols[k]
        logbf = np.array([
            sum(_to_float(row[c]) if len(row) > c else float("nan") for c in cols)
            for row in rows
        ])
        alpha_by_comp[k] = _softmax(logbf)

    marginal = None
    if alpha_by_comp:
        alpha_mat = np.stack([alpha_by_comp[k] for k in comp_ks], axis=1)
        marginal = 1.0 - np.prod(1.0 - alpha_mat, axis=1)

    cs_all_i = idx.get("cs_all")
    member_cs = defaultdict(list)  # row index -> [cs, ...]
    if cs_all_i is not None:
        for i, row in enumerate(rows):
            val = row[cs_all_i] if len(row) > cs_all_i else ""
            for tok in val.split(","):
                tok = tok.strip()
                if not tok:
                    continue
                try:
                    v = int(float(tok))
                except ValueError:
                    continue
                if v != 0:
                    member_cs[i].append(v)
    elif alpha_by_comp:
        for k in comp_ks:
            for i in _native_cs_indices(alpha_by_comp[k]):
                member_cs[i].append(k)

    pip_cs = {vid: (pip, member_cs.get(i, [])) for i, (vid, pip) in enumerate(zip(vids, pips))}
    weight_rows = []
    for i, cs_ids in member_cs.items():
        vid = vids[i]
        pm = float(marginal[i]) if marginal is not None else float("nan")
        for k in cs_ids:
            a = float(alpha_by_comp[k][i]) if k in alpha_by_comp else float("nan")
            weight_rows.append((k, vid, a, pm))
    return pip_cs, weight_rows


###############################################################################
##########  PER-UNIT PROCESSING  ##############################################
###############################################################################

def _process_unit(arg):
    """Process one (trait, locus, pop): read all three tools, build CS-member
    rows, variant rows and SuSiEx native-weight rows. roots = {tool: root}."""
    trait, locus, pop, roots = arg
    per_tool = {}
    weight_rows = []
    for t in TOOLS:
        path = os.path.join(roots[t], t, pop, f"{trait}.{locus}.{t}.txt")
        if t == "susiex":
            pip_cs, w_rows = _read_susiex_native(path)
            per_tool[t] = pip_cs
            weight_rows.extend((trait, locus, pop, cs, vid, a, pm) for cs, vid, a, pm in w_rows)
        else:
            per_tool[t] = _read_tool(path)

    # ---- credible-set members: group each tool's variants by cs id ----
    member_rows = []
    for tool, d in per_tool.items():
        groups = {}
        for vid, (pip, cs_ids) in d.items():
            for cs in cs_ids:
                groups.setdefault(cs, []).append((vid, pip))
        for cs, members in groups.items():
            members.sort(key=lambda x: (x[1] if x[1] == x[1] else -1.0), reverse=True)
            sz = len(members)
            for mv, _mp in members:  # one row per CS member, for CS-size-bucketed analyses
                member_rows.append((trait, locus, pop, tool, cs, sz, mv))

    # ---- variant level: union variants across tools, cross-method median PIP, filter ----
    var_rows = []
    all_vids = set().union(*(set(d) for d in per_tool.values())) if per_tool else set()
    for vid in all_vids:
        pips, in_cs_any, tp = [], False, {}
        for tool in TOOLS:
            hit = per_tool[tool].get(vid)
            if hit is None:
                continue
            p, cs_ids = hit
            tp[tool] = p
            if p == p:
                pips.append(p)
            if cs_ids:
                in_cs_any = True
        if not pips:
            continue
        max_pip = max(pips)
        if max_pip < PIP_KEEP and not in_cs_any:
            continue
        var_rows.append((trait, locus, pop, vid,
                          tp.get("susiex", ""), tp.get("multisusie", ""), tp.get("rsparsepro", ""),
                          statistics.median(pips), max_pip, int(in_cs_any)))
    return var_rows, member_rows, weight_rows


def enumerate_units(roots, pops):
    """List (trait, locus, pop) triples as the union of per-tool output file names."""
    units = set()
    for pop in pops:
        for tool in TOOLS:
            d = Path(roots[tool]) / tool / pop
            if not d.is_dir():
                continue
            suf = f".{tool}.txt"
            for f in d.iterdir():
                if not f.name.endswith(suf):
                    continue
                stem = f.name[: -len(suf)]
                dot = stem.find(".")
                if dot > 0:
                    units.add((stem[:dot], stem[dot + 1:], pop))
    return sorted(units)


###############################################################################
##########  MAIN  ##############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "root=", "rsparsepro-root=", "out-dir=", "pops=", "workers=", "limit="]
        )
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    root = ROOT_DEFAULT
    rsparsepro_root = RSPARSEPRO_ROOT_DEFAULT
    out_dir = OUT_DIR_DEFAULT
    pops = ("AFR", "EUR", "ALL")
    workers = min(16, os.cpu_count() or 4)
    limit = 0

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--root":
            root = a
        elif o == "--rsparsepro-root":
            rsparsepro_root = a
        elif o == "--out-dir":
            out_dir = a
        elif o == "--pops":
            pops = tuple(p.strip() for p in a.split(",") if p.strip())
        elif o == "--workers":
            workers = int(a)
        elif o == "--limit":
            limit = int(a)
        else:
            assert False, f"Unhandled option: {o}"

    roots = {t: root for t in TOOLS}
    roots["rsparsepro"] = rsparsepro_root
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    units = enumerate_units(roots, pops)
    if limit:
        units = units[:limit]
    print(f"NOTE: {len(units)} (trait,locus,pop) units over pops={pops}, workers={workers}")
    print(f"NOTE: roots: susiex/multisusie={root}")
    print(f"NOTE:        rsparsepro={rsparsepro_root}")

    var_all, mem_all, weight_all = [], [], []
    tasks = [(t, l, p, roots) for (t, l, p) in units]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, (var_rows, member_rows, weight_rows) in enumerate(ex.map(_process_unit, tasks, chunksize=64), 1):
            var_all.extend(var_rows)
            mem_all.extend(member_rows)
            weight_all.extend(weight_rows)
            if i % 5000 == 0:
                print(f"NOTE:  {i}/{len(units)} units ...")

    mem_path = out_path / "gwas_cs_members.tsv"
    with open(mem_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["trait", "locus", "pop", "tool", "cs", "cs_size", "variant_id"])
        w.writerows(mem_all)

    var_path = out_path / "gwas_finemapped_variants.tsv"
    with open(var_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["trait", "locus", "pop", "variant_id",
                    "susiex_pip", "multisusie_pip", "rsparsepro_pip",
                    "median_pip", "max_pip", "in_cs_any"])
        w.writerows(var_all)

    weight_path = out_path / "susiex_native_weights.tsv.gz"
    with gzip.open(weight_path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["trait", "locus", "pop", "cs", "variant_id", "alpha", "pip_marginal"])
        w.writerows(weight_all)

    uniq_all = {r[3] for r in var_all}
    uniq_med01 = {r[3] for r in var_all if r[7] >= PIP_KEEP}
    uniq_incs = {r[3] for r in var_all if r[9] == 1}
    print(f"NOTE: finemapped-variant rows (kept): {len(var_all)} -> {var_path}")
    print(f"NOTE: UNIQUE variants (any kept)         : {len(uniq_all)}")
    print(f"NOTE: UNIQUE variants median_pip>=0.01   : {len(uniq_med01)}  <- gene-nomination/VEP target set")
    print(f"NOTE: UNIQUE variants in >=1 credible set: {len(uniq_incs)}")
    print(f"NOTE: SuSiEx native CS-member weight rows: {len(weight_all)} -> {weight_path}")
    print(f"NOTE: CS member rows (all tools)         : {len(mem_all)} -> {mem_path}")


if __name__ == "__main__":
    main()
