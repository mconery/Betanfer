#!/usr/bin/env python3

'''

plot_mvp_teff_vs_t.py

The purpose of this script is to plot the effective number of independent
traits (T_eff, from collate_effective_trait_counts.py) as a function of the
number of traits actually included, for the beta, z, and normalized-z
Pearson correlation estimators. This is the Supplementary Note figure
showing how far below the T_eff = T line each estimator falls.

Input:  <in_dir>/mvp_teff_vs_t.tsv, <in_dir>/mvp_teff_summary.tsv
        (collate_effective_trait_counts.py).
Output: <in_dir>/mvp_teff_vs_t.{pdf,png}

Usage:
  python plot_mvp_teff_vs_t.py [-i <in_dir>]

***Run in the betainfer conda environment on Aurora (needs matplotlib)***

'''

import sys
import getopt
import os

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import pandas as pd

# ── Defaults ──────────────────────────────────────────────────────────────────
DEFAULT_IN_DIR = ("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/"
                   "mvp_results/effective_traits")
SHOW = {"beta": ("β", "#0072B2"), "z": ("z", "#D55E00"), "z_norm": ("normalized z", "#009E73")}

# ── Font: Arial, with Liberation Sans standing in for the faces not present
# on this system (metric-compatible with Arial) ─────────────────────────────
FONT_DIR = Path.home() / ".local/share/fonts"
FONT_FACES = [("Arial.ttf", "normal", 400), ("arialbd.ttf", "normal", 700),
              ("ariali.ttf", "italic", 400), ("arialbi.ttf", "italic", 700)]
FONT_FALLBACK = {"arialbd.ttf": "LiberationSans-Bold.ttf", "ariali.ttf": "LiberationSans-Italic.ttf",
                  "arialbi.ttf": "LiberationSans-BoldItalic.ttf"}
FONT_RC = {"font.family": "Arial", "pdf.fonttype": 42, "mathtext.fontset": "custom",
           "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"}

# ── Style shared with the application figures (Figures 5 and 6) ───────────────
INK, INK_SECONDARY = "#1A1A1A", "#4D4D4D"
STYLE_RC = {"font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7,
            "ytick.labelsize": 7, "legend.fontsize": 7, "text.color": INK, "axes.labelcolor": INK,
            "xtick.color": INK, "ytick.color": INK, "axes.edgecolor": INK_SECONDARY,
            "axes.linewidth": 0.6, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
            "xtick.major.size": 2.8, "ytick.major.size": 2.8}


def _register_arial_fonts():
    """Register Arial under its own family name so normal font lookups resolve
    to it, falling back to Liberation Sans for any face whose Arial file is
    not installed, then apply the font and style rcParams."""
    fm.fontManager.ttflist = [f for f in fm.fontManager.ttflist if f.name != "Arial"]
    for name, style, weight in FONT_FACES:
        path = FONT_DIR / name if (FONT_DIR / name).exists() else FONT_DIR / FONT_FALLBACK[name]
        if path.exists():
            fm.fontManager.ttflist.append(fm.FontEntry(fname=str(path), name="Arial", style=style, weight=weight))
    plt.rcParams.update(FONT_RC)
    plt.rcParams.update(STYLE_RC)


_register_arial_fonts()


def help(exit_num=0):
    print("""-----------------------------------------------------------------
ARGUMENTS
    -i  Directory with mvp_teff_vs_t.tsv and mvp_teff_summary.tsv, and where
        mvp_teff_vs_t.{pdf,png} are written
        (default: mvp_results/effective_traits/)
    -h  Print this help message
-----------------------------------------------------------------""")
    sys.exit(exit_num)


def clean(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#e6e6e6", lw=0.6)
    ax.set_axisbelow(True)


###############################################################################
################################# DRIVER #####################################
###############################################################################

def driver(in_dir):
    cur = pd.read_csv(os.path.join(in_dir, "mvp_teff_vs_t.tsv"), sep="\t")
    cur["arm"] = cur.estimator.str.replace("_pearson", "")
    summ = pd.read_csv(os.path.join(in_dir, "mvp_teff_summary.tsv"), sep="\t")
    full = summ[summ.variants.eq("all_common")].set_index("estimator").corrected_mean
    T = cur["T"].max()

    fig, ax = plt.subplots(figsize=(3.6, 3.0), constrained_layout=True)
    ax.plot([0, 500], [0, 500], color="#9e9e9e", ls="--", lw=0.8)
    ax.text(330, 440, "$T_{\\mathrm{eff}} = T$", color="#616161", ha="right", va="center")
    for arm, (lab, col) in SHOW.items():
        g = cur[cur.arm.eq(arm)].groupby("T").corrected.agg(["mean", "min", "max"])
        ax.fill_between(g.index, g["min"], g["max"], color=col, alpha=0.18, lw=0)
        ax.plot(g.index, g["mean"], "-o", color=col, lw=1.5, ms=3.5, mec="white", mew=0.6, label=lab)
        ax.annotate(f"{full[arm + '_pearson']:.0f}", (T, full[arm + "_pearson"]), xytext=(5, 0),
                    textcoords="offset points", va="center", color="#212121")
    ax.set_xlim(0, T * 1.1)
    ax.set_ylim(0, 500)
    ax.set_xlabel("Number of traits, $T$")
    ax.set_ylabel("Effective number of independent traits, $T_{\\mathrm{eff}}$")
    clean(ax)
    ax.legend(frameon=False, loc="lower right")

    for ext, kw in (("pdf", {}), ("png", {"dpi": 300})):
        out_path = os.path.join(in_dir, f"mvp_teff_vs_t.{ext}")
        fig.savefig(out_path, **kw)
        print(f"NOTE: Written to {out_path}")
    plt.close(fig)


###############################################################################
########################### COLLECT AND CHECK ARGUMENTS ######################
###############################################################################

def main(argv):
    try:
        opts, _ = getopt.getopt(argv[1:], "i:h")
    except getopt.GetoptError as e:
        print(f"ERROR: {e}")
        help(1)

    options_dict = dict(opts)
    if '-h' in options_dict:
        help()

    in_dir = options_dict.get('-i', DEFAULT_IN_DIR)

    for fname in ("mvp_teff_vs_t.tsv", "mvp_teff_summary.tsv"):
        if not os.path.exists(os.path.join(in_dir, fname)):
            print(f"ERROR: required input not found: {os.path.join(in_dir, fname)}")
            sys.exit(1)

    print("Acceptable Inputs Given")
    driver(in_dir)


if __name__ == "__main__":
    main(sys.argv)
