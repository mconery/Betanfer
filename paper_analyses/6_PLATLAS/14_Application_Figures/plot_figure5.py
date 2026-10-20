#!/usr/bin/env python3
'''
plot_figure5.py

Renders Figure 5: credible-set (CS) matching against the ABF baseline (a),
fine-mapping resolution by CS size (b), coding-variant composition by CS
size (c), nominated genes by CS-matching class (d), the per-gene pleiotropy
distribution (e), and the top pleiotropic genes' phecode-domain matrix with
their nomination-evidence breakdown (f).

Input:
  FIG5_SOURCE_DATA  directory with the tables written by build_figure5_tables.py:
    figure5a_cs_matching.tsv, cs_resolution_four_bins.tsv,
    annotation_composition.tsv, gene_ranking.tsv, gene_domain_matrix.tsv,
    f_top_bar_counts.tsv, d_genes_by_cs_class.tsv
  FONTS_DIR         directory with Arial.ttf (falls back to the default
                    sans-serif font if not found; see figure_style.register_arial)

Output:
  OUT_DIR/Figure5.png, OUT_DIR/Figure5.pdf

Usage:
  python plot_figure5.py [--source-data DIR] [--fonts-dir DIR] [--out-dir DIR]

***Run in finemap conda environment***
'''
import getopt
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

import figure_style as S

###############################################################################
##########  DEFAULT PATHS  ####################################################
###############################################################################

ROOT = Path("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping")
FIG5_SOURCE_DATA = ROOT / "application_results" / "figures" / "source_data" / "fig5"
FONTS_DIR = Path.home() / ".local/share/fonts"
OUT_DIR = ROOT / "application_results" / "figures"

###############################################################################
##########  PANEL A: credible-set matching against the ABF baseline  #########
###############################################################################

def draw_cs_matching(ax, frame):
    frame = frame.set_index(["pop", "tool"]).loc[S.PAIRS].reset_index()
    xx = np.arange(len(frame))
    for x, row in zip(xx, frame.itertuples()):
        shared_pct = 100 * row.new_matched / row.new_cs
        ax.bar(x, shared_pct, .75, color=S.SHARED_NEUTRAL)
        ax.bar(x, 100 - shared_pct, .75, bottom=shared_pct, color=S.TOOL_COLORS[row.tool])
        ax.text(x, shared_pct / 2, f"{row.new_matched:,}", ha="center", va="center", fontsize=6.2)
        ax.text(x, (100 + shared_pct) / 2, f"{row.new_unmatched:,}", ha="center", va="center", fontsize=6.2, color="white")
    ax.set_xticks(xx, [f"{row.pop}\n{S.TOOL_LABELS[row.tool]}" for row in frame.itertuples()], fontsize=6)
    ax.set_ylim(0, 100); ax.set_ylabel("Credible sets (%)")
    handles = [Patch(color=S.SHARED_NEUTRAL, label="Shared with ABF")] + \
        [Patch(color=S.TOOL_COLORS[t], label=f"Additional ({S.TOOL_LABELS[t]})") for t in S.TOOL_ORDER]
    ax.legend(handles=handles, ncol=2, frameon=False, fontsize=6.3, loc="upper center", bbox_to_anchor=(.5, -.23))
    S.clean_axis(ax); ax.tick_params(axis="x", labelsize=6)


###############################################################################
##########  PANEL B: fine-mapping resolution  ##################################
###############################################################################

def draw_resolution(ax, resolution):
    left = np.zeros(len(S.TOOL_ORDER))
    for size in S.RESOLUTION_BINS:
        vals = resolution[resolution.size_bin.eq(size)].set_index("tool").loc[S.TOOL_ORDER, "percent"].to_numpy()
        ax.barh(np.arange(len(S.TOOL_ORDER)), vals, .56, left=left, color=S.RESOLUTION_COLORS[size], label=size)
        for y, (value, lft) in enumerate(zip(vals, left)):
            ax.text(lft + value / 2, y, f"{value:.0f}", ha="center", va="center", fontsize=7)
        left += vals
    ax.set_yticks(range(len(S.TOOL_ORDER)), [S.TOOL_LABELS[t] for t in S.TOOL_ORDER])
    ax.invert_yaxis(); ax.set_xlim(0, 100)
    ax.set_xlabel("Credible sets (%)")
    ax.legend(ncol=4, title="Variants per CS", title_fontsize=7, fontsize=7,
              loc="upper center", bbox_to_anchor=(.5, -.22), handlelength=1)
    S.clean_axis(ax); ax.tick_params(labelsize=7, length=3, width=.6)


###############################################################################
##########  PANEL C: coding-variant composition by CS size  ###################
###############################################################################

def draw_coding_composition(ax, annotation):
    groups = ["11+", "6-10", "2-5", "1"]
    categories = [("Synonymous", "Synonymous"), ("Missense", "Missense"), ("Splice/Start/Stop Gain/Loss", "Splice/start/stop")]
    bottom = np.zeros(len(groups))
    for cat, label in categories:
        values = annotation[annotation.category.eq(cat)].set_index("group").loc[groups, "percent"].to_numpy()
        ax.bar(np.arange(len(groups)), values, .60, bottom=bottom, color=S.CODING_COLORS[cat], label=label)
        bottom += values
    for x, value in enumerate(bottom):
        ax.text(x, value + .12, f"{value:.2f}", ha="center", fontsize=7)
    all_coding = annotation[annotation.group.eq("All") & annotation.category.isin(S.CODING_CATS)].percent.sum()
    ax.axhline(all_coding, color="#777777", lw=.8, ls="--", label=f"All variants ({all_coding:.2f}%)")
    ax.set_xticks(range(len(groups)), ["11+", "6–10", "2–5", "1"])
    ax.set_xlabel("CS size"); ax.set_ylabel("Coding variants (%)")
    ax.set_ylim(0, 5.15)
    ax.legend(fontsize=6.8, loc="upper left", frameon=False)
    S.clean_axis(ax); ax.tick_params(labelsize=7, length=3, width=.6)


###############################################################################
##########  PANEL D: nominated genes by CS-matching class  #####################
###############################################################################

def draw_genes_by_cs_class(ax, counts):
    bottom = np.zeros(len(counts))
    for cls, color in zip(S.PANEL_D_CLASSES, S.PANEL_D_COLORS):
        ax.bar(range(len(counts)), counts[cls], .56, bottom=bottom, color=color, label=cls)
        for i, (v, b) in enumerate(zip(counts[cls], bottom)):
            ax.text(i, b + v / 2, f"{v:,}", ha="center", va="center", fontsize=7,
                    color="white" if cls == "Additional CS only" else "#202020")
        bottom += counts[cls].to_numpy()
    ax.set_xticks(range(len(counts)), counts.index)
    ax.set_xlabel("Nomination evidence"); ax.set_ylabel("Nominated genes")
    ax.set_xlim(-.6, len(counts) - .4)
    S.clean_axis(ax); ax.tick_params(labelsize=7, length=3, width=.6)
    ax.legend(frameon=False, fontsize=6.5, loc="upper right", bbox_to_anchor=(1.0, 1.0), handlelength=1.4)


###############################################################################
##########  PANEL E: per-gene pleiotropy distribution  #########################
###############################################################################

def draw_pleiotropy(ax, genes):
    dist = genes.groupby("phecode_groups").size().rename("genes").reset_index()
    ax.bar(dist.phecode_groups, dist.genes, .90, color=S.PLEIOTROPY_COLOR, linewidth=0)
    ax.set_yscale("log"); ax.set_ylim(.8, 6500); ax.set_xlim(0, 40)
    ax.set_xticks([1, 10, 20, 30, 40])
    ax.set_xlabel("Phecode groups per gene"); ax.set_ylabel("Genes")
    S.clean_axis(ax); ax.tick_params(labelsize=7, length=3, width=.6)


###############################################################################
##########  PANEL F: top pleiotropic genes x phecode-domain matrix  ############
###############################################################################

def draw_gene_domain_matrix(fig, slot, genes, cells, counts):
    inner = slot.subgridspec(2, 1, height_ratios=[.5, 1.65], hspace=.07)
    bar, ax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])
    names = genes.loc[genes.main_selected, "gene"].tolist()
    bottom = np.zeros(len(names))
    for kind, label in [("coding", "Coding"), ("abc", "ABC"), ("both", "Both")]:
        values = counts[counts.evidence.eq(kind)].set_index("gene").loc[names, "phecode_groups"].to_numpy()
        bar.bar(range(len(names)), values, .73, bottom=bottom, color=S.EVID_COLORS[kind], label=label)
        bottom += values
    bar.set_ylabel("Phecode\ngroups", fontsize=7); bar.set_xticks([])
    bar.set_yticks([0, 20, 40]); bar.set_ylim(0, 43)
    S.clean_axis(bar); bar.tick_params(labelsize=7, length=3, width=.6)
    bar.legend(ncol=3, loc="lower center", bbox_to_anchor=(.5, 1.02), frameon=False, fontsize=6)

    subset = cells[cells.gene.isin(names)]
    domains = sorted(set(subset.domain) - {"other", "unmapped", "unclassified"})
    domains += [d for d in ["other", "unmapped", "unclassified"] if d in set(subset.domain)]
    for y, domain in enumerate(domains):
        color = S.DOMAIN_SCATTER_COLORS[y % len(S.DOMAIN_SCATTER_COLORS)]
        ax.axhspan(y - .5, y + .5, color=color, alpha=.10, linewidth=0)
        part = subset[subset.domain.eq(domain)]
        ax.scatter([names.index(g) for g in part.gene], [y] * len(part), s=9 * part.phecode_groups, color=color, linewidth=0)
    ax.set_xticks(range(len(names)), names, rotation=65, ha="right", fontsize=6.4, fontstyle="italic")
    ax.set_yticks(range(len(domains)), [S.domain_label(d) for d in domains], fontsize=6.3)
    ax.set_xlim(-.6, len(names) - .4); bar.set_xlim(ax.get_xlim())
    ax.set_ylim(len(domains) - .5, -.5); ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    handles = [ax.scatter([], [], s=9 * n, color="#888888") for n in [1, 5, 10]]
    ax.legend(handles, ["1", "5", "10"], title="Phecode groups", title_fontsize=6.5, ncol=3, fontsize=6,
              frameon=False, loc="upper center", bbox_to_anchor=(.5, -.25))
    return bar, ax


###############################################################################
##########  ASSEMBLY  ##########################################################
###############################################################################

def build_figure(source_data):
    matching = pd.read_csv(source_data / "figure5a_cs_matching.tsv", sep="\t")
    resolution = pd.read_csv(source_data / "cs_resolution_four_bins.tsv", sep="\t")
    annotation = pd.read_csv(source_data / "annotation_composition.tsv", sep="\t")
    class_counts = pd.read_csv(source_data / "d_genes_by_cs_class.tsv", sep="\t", index_col=0)
    genes = pd.read_csv(source_data / "gene_ranking.tsv", sep="\t")
    domain_matrix = pd.read_csv(source_data / "gene_domain_matrix.tsv", sep="\t")
    bar_counts = pd.read_csv(source_data / "f_top_bar_counts.tsv", sep="\t")

    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8})
    fig = plt.figure(figsize=(9.6, 9.2))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.18], left=.085, right=.975,
                           top=.965, bottom=.12, wspace=.34, hspace=.54)

    a = fig.add_subplot(gs[0, 0]); draw_cs_matching(a, matching)
    b = fig.add_subplot(gs[0, 1]); draw_resolution(b, resolution)
    c = fig.add_subplot(gs[1, 0]); draw_coding_composition(c, annotation)
    d = fig.add_subplot(gs[1, 1]); draw_genes_by_cs_class(d, class_counts[S.PANEL_D_CLASSES])
    e = fig.add_subplot(gs[2, 0]); draw_pleiotropy(e, genes)
    fbar, f = draw_gene_domain_matrix(fig, gs[2, 1], genes, domain_matrix, bar_counts)

    right = gs[0, 1].get_position(fig).x0 - .073
    for ax, letter in [(a, "a"), (b, "b"), (c, "c"), (d, "d"), (e, "e"), (fbar, "f")]:
        fig.text(.03 if letter in "ace" else right, ax.get_position().y1 + .015, letter.upper(),
                  fontsize=12, fontweight="bold")
    return fig


def usage():
    print(__doc__)


def main(argv):
    source_data, fonts_dir, out_dir = FIG5_SOURCE_DATA, FONTS_DIR, OUT_DIR
    try:
        opts, _ = getopt.getopt(argv, "h", ["help", "source-data=", "fonts-dir=", "out-dir="])
    except getopt.GetoptError as e:
        print(f"ERROR: {e}"); usage(); sys.exit(2)
    for opt, val in opts:
        if opt in ("-h", "--help"):
            usage(); sys.exit(0)
        elif opt == "--source-data":
            source_data = Path(val)
        elif opt == "--fonts-dir":
            fonts_dir = Path(val)
        elif opt == "--out-dir":
            out_dir = Path(val)

    if not source_data.exists():
        print(f"ERROR: source-data directory not found: {source_data}"); sys.exit(1)
    out_dir.mkdir(parents=True, exist_ok=True)

    S.register_arial(fonts_dir)
    S.apply_rc()
    fig = build_figure(source_data)
    for ext in ["png", "pdf"]:
        fig.savefig(out_dir / f"Figure5.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"NOTE: wrote {out_dir / 'Figure5.png'} and Figure5.pdf")


if __name__ == "__main__":
    main(sys.argv[1:])
