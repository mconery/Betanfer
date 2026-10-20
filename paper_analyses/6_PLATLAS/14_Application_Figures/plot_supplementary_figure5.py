#!/usr/bin/env python3
'''
plot_supplementary_figure5.py

Renders Supplementary Figure 5: credible-set matching sensitivity across
PIP-weighted-Jaccard thresholds (a), VEP annotation composition by CS size
(b), gene-nomination evidence summary (c), and the GO/KEGG pathway
enrichment dot plot for nominated genes, one selected term per phecode
domain (d).

Panel d reads the already-selected, non-redundant enrichment terms written
by build_figureS5_tables.py; if its two input tables are not found, panel d
is left blank with a short placeholder note instead of failing the whole
figure.

Input:
  FIGS5_SOURCE_DATA  directory with the tables written by build_figureS5_tables.py:
    a_matching_sensitivity.tsv, b_annotation_exact.tsv, c_gene_evidence.tsv,
    and (optional) d_enrichment_display_order.tsv, d_enrichment_visible_cells.tsv
  FONTS_DIR          directory with Arial.ttf (see figure_style.register_arial)

Output:
  OUT_DIR/SupplementaryFigure5.png, OUT_DIR/SupplementaryFigure5.pdf

Usage:
  python plot_supplementary_figure5.py [--source-data DIR] [--fonts-dir DIR] [--out-dir DIR]

***Run in finemap conda environment***
'''
import getopt
import re
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.text import Text
import numpy as np
import pandas as pd

import figure_style as S

# Arial has no subscript-digit glyphs, which the panel-d colorbar label
# ("-log10 FDR") needs; rewrite unicode subscript digits as mathtext
# subscripts right before saving so the glyph renders instead of falling
# back to a different font.
SUBSCRIPT_DIGITS = re.compile("[₀-₉]+")
DIGIT_TRANSLATION = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def fix_subscript_glyphs(fig):
    for text in fig.findobj(Text):
        label = text.get_text()
        if SUBSCRIPT_DIGITS.search(label):
            text.set_text(SUBSCRIPT_DIGITS.sub(lambda m: "$_{%s}$" % m.group().translate(DIGIT_TRANSLATION), label))

###############################################################################
##########  DEFAULT PATHS  ####################################################
###############################################################################

ROOT = Path("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping")
FIGS5_SOURCE_DATA = ROOT / "application_results" / "figures" / "source_data" / "figs5"
FONTS_DIR = Path.home() / ".local/share/fonts"
OUT_DIR = ROOT / "application_results" / "figures"

###############################################################################
##########  PANEL A: credible-set matching sensitivity  #######################
###############################################################################

def draw_sensitivity(ax, sensitivity):
    markers, styles = ["o", "s", "^", "o", "o"], ["-", "--", ":", "-", "-"]
    for (pop, tool), marker, style in zip(S.PAIRS, markers, styles):
        d = sensitivity[sensitivity["pop"].eq(pop) & sensitivity.tool.eq(tool)].sort_values("threshold")
        pct = 100 * d.new_matched / d.new_cs
        ax.plot(d.threshold, pct, marker=marker, ms=4, ls=style, lw=1.2,
                color=S.TOOL_COLORS[tool], label=f"{pop} {S.TOOL_LABELS[tool]}")
    ax.axvline(S.PRIMARY_THRESHOLD, color="#aaaaaa", lw=.7, ls="--", zorder=0)
    ax.set_xlabel("PIP-weighted Jaccard threshold"); ax.set_ylabel("CS shared with ABF (%)")
    ax.set_xticks([.05, .10, .20, .30, .50], ["0.05", "0.10", "0.20", "0.30", "0.50"])
    ax.set_xlim(.035, .515); ax.set_ylim(0, 65); ax.set_yticks([0, 20, 40, 60])
    ax.legend(ncol=2, frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(.5, -.24),
              columnspacing=1, handlelength=1.8)
    S.clean_axis(ax); ax.tick_params(labelsize=7, length=3, width=.6)


###############################################################################
##########  PANEL B: VEP annotation composition by CS size  ###################
###############################################################################

def draw_annotation(ax, annotation):
    """Collapses the 10 raw VEP categories into figure_style.VEP_GROUPS at draw
    time, straight from the exact per-category table (no precomputed display
    table); a display category with zero variants in every group is skipped."""
    groups = ["All", "11+", "6-10", "2-5", "1"]
    bottom = np.zeros(len(groups))
    for label, categories in S.VEP_GROUPS:
        values = annotation[annotation.category.isin(categories)].groupby("group")["percent"].sum().loc[groups].to_numpy()
        if not np.any(values):
            continue
        display_label = label if label == "UTR / ncRNA exon" else label.capitalize()
        ax.bar(range(5), values, .67, bottom=bottom, color=S.VEP_GROUP_COLORS[label], label=display_label)
        bottom += values
    ax.set_xticks(range(5), ["All", "11+", "6–10", "2–5", "1"])
    ax.set_ylim(0, 100); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("CS size"); ax.set_ylabel("Annotated variants (%)")
    ax.legend(ncol=2, frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(.5, -.24),
              columnspacing=1, handlelength=1)
    S.clean_axis(ax); ax.tick_params(labelsize=7, length=3, width=.6)


###############################################################################
##########  PANEL C: gene-nomination evidence summary  ########################
###############################################################################

def draw_gene_evidence(ax, evidence):
    order = ["coding", "abc", "both"]
    evidence = evidence.set_index("evidence").loc[order].reset_index()
    ax.barh(range(3), evidence.percent, height=.5, color=[S.EVID_COLORS[e] for e in order])
    for y, row in evidence.iterrows():
        ax.text(row.percent + 1.5, y, f"{row.genes:,} ({row.percent:.1f}%)", va="center", fontsize=8)
    ax.set_yticks(range(3), evidence.label)
    ax.set_ylim(2.6, -.6); ax.set_xlim(0, 100); ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Nominated genes (%)")
    S.clean_axis(ax); ax.tick_params(labelsize=7, length=3, width=.6)


###############################################################################
##########  PANEL D: GO/KEGG pathway enrichment dot plot  #####################
###############################################################################

def draw_enrichment(fig, slot, selected, visible):
    """Two-column layout inside `slot`: a text axis with wrapped term labels
    on the left, the domain x term dot-plot scatter on the right. Returns the
    (label_axis, scatter_axis) pair so the caller can place the panel letter."""
    if selected is None or visible is None:
        ax = fig.add_subplot(slot)
        ax.text(.5, .5, "GO/KEGG enrichment selection\nnot available", ha="center", va="center",
                fontsize=7.5, color=S.INK_SECONDARY, transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        return ax, ax

    selected = selected.sort_values("display_row").reset_index(drop=True)
    domains = sorted(visible.domain.unique())
    domain_text = dict(zip(visible.domain, visible.domain_label))
    keys = list(zip(selected.library, selected.term))
    inner = slot.subgridspec(1, 2, width_ratios=[1.25, 1], wspace=.04)
    labels, ax = fig.add_subplot(inner[0]), fig.add_subplot(inner[1])

    for y, row in selected.iterrows():
        label = row.term_label
        if row.library.startswith("GO_"):
            label = label[0].upper() + label[1:].lower()
        labels.text(.98, y, textwrap.fill(label, 33, break_long_words=False),
                    ha="right", va="center", fontsize=6.3)
        if y % 2 == 0:
            for axis in [labels, ax]:
                axis.axhspan(y - .5, y + .5, color="#f3f4f5", lw=0, zorder=0)

    points = ax.scatter([domains.index(v) for v in visible.domain],
                         [keys.index((r.library, r.term)) for r in visible.itertuples()],
                         s=5 * visible.n_overlap, c=-np.log10(visible.q_global), cmap=S.enrichment_cmap(),
                         norm=Normalize(-np.log10(.05), 3.1), edgecolors="white", linewidths=.35)
    ax.set_xticks(range(len(domains)), [domain_text[v] for v in domains], rotation=60, ha="right", fontsize=6.2)
    ax.set_yticks([]); ax.set_xlim(-.55, len(domains) - .45)
    labels.set_xlim(0, 1); labels.axis("off")
    for axis in [labels, ax]:
        axis.set_ylim(len(selected) - .5, -.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)

    cax = ax.inset_axes([1.07, .52, .035, .46])
    colorbar = fig.colorbar(points, cax=cax)
    colorbar.set_label("−log₁₀ FDR", fontsize=6.5)
    colorbar.ax.tick_params(labelsize=6, length=2)

    handles = [ax.scatter([], [], s=5 * n, color="#77828d") for n in [3, 15, 30]]
    ax.legend(handles, ["3", "15", "30"], title="Genes", title_fontsize=6.5, fontsize=6,
              frameon=False, loc="upper left", bbox_to_anchor=(1.035, .43),
              borderpad=0, labelspacing=1.3, handletextpad=.5)
    return labels, ax


###############################################################################
##########  ASSEMBLY  ##########################################################
###############################################################################

def build_figure(source_data):
    sensitivity = pd.read_csv(source_data / "a_matching_sensitivity.tsv", sep="\t")
    annotation = pd.read_csv(source_data / "b_annotation_exact.tsv", sep="\t")
    evidence = pd.read_csv(source_data / "c_gene_evidence.tsv", sep="\t")
    display_order_path = source_data / "d_enrichment_display_order.tsv"
    visible_cells_path = source_data / "d_enrichment_visible_cells.tsv"
    display_order = pd.read_csv(display_order_path, sep="\t") if display_order_path.exists() else None
    visible_cells = pd.read_csv(visible_cells_path, sep="\t") if visible_cells_path.exists() else None
    if display_order is None or visible_cells is None:
        print("NOTE: GO/KEGG enrichment tables not found; panel d will be left blank")

    plt.rcParams.update({"font.size": 8, "axes.labelsize": 8})
    fig = plt.figure(figsize=(10.8, 8.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.18], left=.095, right=.945,
                           top=.96, bottom=.20, wspace=.30, hspace=.64)

    a = fig.add_subplot(gs[0, 0]); draw_sensitivity(a, sensitivity)
    b = fig.add_subplot(gs[0, 1]); draw_annotation(b, annotation)
    for container in b.containers:
        for group, patch in zip(["All", "11+", "6-10", "2-5", "1"], container.patches):
            value = patch.get_height()
            if value < 3:
                continue
            b.text(patch.get_x() + patch.get_width() / 2, patch.get_y() + value / 2,
                   f"{value:.1f}%", ha="center", va="center",
                   fontsize=7 if value >= 10 else 5.5, color="#202020")
    c = fig.add_subplot(gs[1, 0]); draw_gene_evidence(c, evidence)
    d, _ = draw_enrichment(fig, gs[1, 1], display_order, visible_cells)

    right = gs[0, 1].get_position(fig).x0 - .035
    for ax, letter, x in [(a, "a", .025), (b, "b", right), (c, "c", .025), (d, "d", right)]:
        fig.text(x, ax.get_position().y1 + .02, letter.upper(), fontsize=12, fontweight="bold")
    return fig


def usage():
    print(__doc__)


def main(argv):
    source_data, fonts_dir, out_dir = FIGS5_SOURCE_DATA, FONTS_DIR, OUT_DIR
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
    fix_subscript_glyphs(fig)
    for ext in ["png", "pdf"]:
        fig.savefig(out_dir / f"SupplementaryFigure5.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"NOTE: wrote {out_dir / 'SupplementaryFigure5.png'} and SupplementaryFigure5.pdf")


if __name__ == "__main__":
    main(sys.argv[1:])
