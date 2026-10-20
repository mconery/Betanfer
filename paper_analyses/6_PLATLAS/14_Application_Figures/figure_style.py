#!/usr/bin/env python3
'''
figure_style.py

Shared constants and small helpers for Figure 5 and Supplementary Figure 5
(credible-set matching, fine-mapping resolution, gene nomination and
pleiotropy panels). Not a script: imported by the build_* and plot_* files
in this directory. One colour means one thing in every panel of every
figure that uses this module.

Usage: import figure_style as S   (from a script in this same directory)
'''

import re

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

###############################################################################
##########  TOOLS, POPULATIONS, EVIDENCE  #####################################
###############################################################################

TOOL_ORDER = ["rsparsepro", "susiex", "multisusie"]
TOOL_COLORS = {"rsparsepro": "#B96F45", "susiex": "#547A96", "multisusie": "#428679"}
TOOL_LABELS = {"rsparsepro": "RSparsePro", "susiex": "SuSiEx", "multisusie": "MultiSuSiE"}
SHARED_NEUTRAL = "#BBC7CA"  # panel-a "no ABF counterpart" bar segment

# (population, tool) pairs the new fine-mapping is compared against the original
# ABF credible sets for. SuSiEx and MultiSuSiE were only run on the cross-ancestry
# "ALL" meta-analysis; RSparsePro was also run per-ancestry.
PAIRS = [("EUR", "rsparsepro"), ("AFR", "rsparsepro"), ("ALL", "rsparsepro"),
         ("ALL", "susiex"), ("ALL", "multisusie")]

# Credible-set matching: PIP-weighted Jaccard threshold above which a new-method CS
# counts as "shared with" (recovering) an ABF CS. 0.10 is the threshold used for the
# headline figures; the other values are the supplementary sensitivity curve.
MATCH_THRESHOLDS = [0.05, 0.10, 0.20, 0.30, 0.50]
PRIMARY_THRESHOLD = 0.10

EVID_COLORS = {"coding": "#B96F45", "abc": "#428679", "both": "#89708F"}
EVID_LABELS = {"coding": "Coding", "abc": "ABC", "both": "Both"}

# Fig5 new panel d: nominated genes classified by whether the RSparsePro credible
# set(s) linking them to a locus have an ABF counterpart (PIP-weighted Jaccard
# matching, see MATCH_THRESHOLDS below).
PANEL_D_CLASSES = ["Shared CS only", "Shared and additional CS", "Additional CS only"]
PANEL_D_COLORS = ["#BCC6CA", "#DDB7A2", "#B96F45"]

# Fig5 panel e (pleiotropy histogram) bar colour.
PLEIOTROPY_COLOR = "#C3A064"

# Fig5 panel f gene x phecode-domain scatter colour cycle (one colour per domain,
# assigned in the order domains first appear).
DOMAIN_SCATTER_COLORS = ["#547A96", "#428679", "#C3A064", "#A65F57", "#89708F", "#81989F"]

# S5 panel d (GO/KEGG enrichment dot plot) colour ramp, low -> high -log10(FDR).
ENRICHMENT_CMAP_COLORS = ["#E7DFB8", "#72A79B", "#345C79"]


def enrichment_cmap():
    return LinearSegmentedColormap.from_list("enrichment", ENRICHMENT_CMAP_COLORS)

###############################################################################
##########  VEP CONSEQUENCE CATEGORIES  #######################################
###############################################################################

# Variant -> consequence category (compute_vep_enrichment.py in 12_Gene_Nomination
# defines the same mapping; duplicated here because scripts in this directory must
# stay self-contained and not import across step directories).
CODING_CATS = ("Splice/Start/Stop Gain/Loss", "Missense", "Synonymous")
# Fig5 panel c stacked-bar colours, one per coding category (stacking order
# Synonymous / Missense / Splice-start-stop, bottom to top).
CODING_COLORS = {"Synonymous": "#D8BB7C", "Missense": "#B96F45", "Splice/Start/Stop Gain/Loss": "#A65F57"}

MVP_CATS = [
    ("Regulatory Region", "#EF5FA0"), ("TF Binding Site", "#D24FBF"),
    ("Non-Gene Unknown", "#9B84D9"), ("Intron", "#33A0E6"),
    ("Non-Coding Exon", "#1FC3C3"), ("3/5 prime UTR", "#00A878"),
    ("Mature miRNA", "#3CB043"), ("Synonymous", "#9E9327"),
    ("Missense", "#C77B2E"), ("Splice/Start/Stop Gain/Loss", "#F1746A"),
]
CAT_NAMES = [c for c, _ in MVP_CATS]
CAT_COLOR = dict(MVP_CATS)

SO2CAT = {}
for _t in ("splice_acceptor_variant", "splice_donor_variant", "stop_gained", "stop_lost",
           "start_lost", "frameshift_variant", "transcript_ablation", "transcript_amplification"):
    SO2CAT[_t] = "Splice/Start/Stop Gain/Loss"
for _t in ("missense_variant", "protein_altering_variant", "inframe_insertion", "inframe_deletion"):
    SO2CAT[_t] = "Missense"
for _t in ("synonymous_variant", "stop_retained_variant", "start_retained_variant",
           "incomplete_terminal_codon_variant"):
    SO2CAT[_t] = "Synonymous"
SO2CAT["mature_miRNA_variant"] = "Mature miRNA"
for _t in ("5_prime_UTR_variant", "3_prime_UTR_variant"):
    SO2CAT[_t] = "3/5 prime UTR"
for _t in ("non_coding_transcript_exon_variant", "non_coding_transcript_variant", "NMD_transcript_variant"):
    SO2CAT[_t] = "Non-Coding Exon"
for _t in ("intron_variant", "splice_region_variant", "splice_donor_region_variant",
           "splice_polypyrimidine_tract_variant", "splice_donor_5th_base_variant"):
    SO2CAT[_t] = "Intron"
for _t in ("intergenic_variant", "upstream_gene_variant", "downstream_gene_variant"):
    SO2CAT[_t] = "Non-Gene Unknown"
for _t in ("TF_binding_site_variant", "TFBS_ablation", "TFBS_amplification"):
    SO2CAT[_t] = "TF Binding Site"
for _t in ("regulatory_region_variant", "regulatory_region_ablation", "regulatory_region_amplification"):
    SO2CAT[_t] = "Regulatory Region"


def size_bin(sz):
    """Credible-set-size bucket used by the VEP annotation-composition panels (5c/S5b)."""
    if sz == 1:
        return "1"
    if sz <= 5:
        return "2-5"
    if sz <= 10:
        return "6-10"
    return "11+"


def resolution_bin(n_variants):
    """Credible-set-size bucket used by the resolution panel (5b): wider bins, en-dash labels."""
    if n_variants <= 1:
        return "1"
    if n_variants <= 5:
        return "2–5"
    if n_variants <= 20:
        return "6–20"
    return ">20"


RESOLUTION_BINS = ["1", "2–5", "6–20", ">20"]
RESOLUTION_COLORS = {"1": "#428679", "2–5": "#93B8AC", "6–20": "#E2D1A5", ">20": "#C9947D"}

# 10 raw VEP categories collapsed to 7 display groups (S5b); "Regulatory Region"
# and "TF Binding Site" are pooled into one "other non-coding" group.
VEP_GROUPS = [
    ("intergenic", ["Non-Gene Unknown"]),
    ("intron", ["Intron"]),
    ("UTR / ncRNA exon", ["Non-Coding Exon", "3/5 prime UTR", "Mature miRNA"]),
    ("other non-coding", ["Regulatory Region", "TF Binding Site"]),
    ("synonymous", ["Synonymous"]),
    ("missense", ["Missense"]),
    ("splice/start/stop", ["Splice/Start/Stop Gain/Loss"]),
]
VEP_GROUP_COLORS = {
    "intergenic": "#BBC7CA", "intron": "#81989F", "UTR / ncRNA exon": "#428679",
    "other non-coding": "#B39DDB",
    "synonymous": "#D8BB7C", "missense": "#B96F45", "splice/start/stop": "#A65F57",
}
CODING_GROUPS = ("synonymous", "missense", "splice/start/stop")

###############################################################################
##########  PHECODES / GENE DOMAINS  ##########################################
###############################################################################

# A phecode subtype ("Phe_250_2") rolls up to its parent phecode ("Phe_250") for every
# pleiotropy / trait-count statistic. Anchored at both ends: a bare "_\d+$" pattern
# would also eat the phecode's own trailing digits (e.g. "Phe_038" -> "Phe"), silently
# collapsing every subtype-free phecode into one trait and understating pleiotropy.
_SUBTYPE = re.compile(r"^(Phe_\d+)_\d+$")


def base_trait(trait):
    return _SUBTYPE.sub(r"\1", trait)


# Raw phecode "group" category (phecode_info.tsv) -> short domain slug. Only
# categories with an unwieldy raw name are aliased; every other phecode group
# (e.g. "endocrine/metabolic", "sense organs", "genitourinary") passes through
# unchanged as its own slug. "unclassified" phecodes become domain "other".
DOMAIN_NAMES = {
    "circulatory system": "circulatory", "hematopoietic": "blood",
    "infectious diseases": "infectious", "injuries & poisonings": "injury",
    "mental disorders": "mental", "pregnancy complications": "pregnancy",
    "congenital anomalies": "congenital", "unclassified": "other",
}


def domain_label(domain):
    """Display label for a domain slug (panel f / enrichment-panel axis ticks)."""
    overrides = {"injury": "Injury/poisoning", "other": "Other", "unmapped": "Unmapped",
                 "unclassified": "Unclassified"}
    return overrides.get(domain, domain.capitalize())

###############################################################################
##########  MATPLOTLIB STYLE  ##################################################
###############################################################################

INK = "#1A1A1A"
INK_SECONDARY = "#4D4D4D"
RULE = "#C8C8C8"


def register_arial(fonts_dir):
    """Register Arial for matplotlib; bold/italic/bold-italic fall back to the metric-
    compatible Liberation Sans faces (registered under the name "Arial") when the real
    Arial face is not installed. If even the regular face is missing, leaves the default
    sans-serif font in place."""
    faces = [("Arial.ttf", "normal", 400, None), ("arialbd.ttf", "normal", 700, "LiberationSans-Bold.ttf"),
             ("ariali.ttf", "italic", 400, "LiberationSans-Italic.ttf"),
             ("arialbi.ttf", "italic", 700, "LiberationSans-BoldItalic.ttf")]
    if not (fonts_dir / "Arial.ttf").exists():
        print(f"NOTE: {fonts_dir / 'Arial.ttf'} not found; figure will use the default sans-serif font")
        return
    fm.fontManager.ttflist = [f for f in fm.fontManager.ttflist if f.name != "Arial"]
    for name, style_, weight, fallback in faces:
        path = fonts_dir / name if (fonts_dir / name).exists() else fonts_dir / fallback
        if not path.exists():
            print(f"NOTE: no face available for Arial {style_}/{weight}; skipping")
            continue
        fm.fontManager.ttflist.append(fm.FontEntry(fname=str(path), name="Arial", style=style_, weight=weight))
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
                          "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                          "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})


def apply_rc():
    """Journal-grade rcParams shared by Figure 5 and Supplementary Figure 5."""
    plt.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK, "ytick.color": INK,
        "axes.edgecolor": INK_SECONDARY, "axes.linewidth": 0.7,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "xtick.major.size": 2.8, "ytick.major.size": 2.8,
        "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


def clean_axis(ax, grid_axis=None):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.set_axisbelow(True)
        ax.grid(axis=grid_axis, color=RULE, linewidth=0.5, alpha=0.7)
    ax.tick_params(axis="both", labelsize=7, length=3, width=0.6)


def panel_letter(fig, ax, letter, x=None, upper=False):
    """Bold panel letter above the top-left of `ax`; `x` overrides the default left margin."""
    text = letter.upper() if upper else letter
    xpos = 0.03 if x is None else x
    fig.text(xpos, ax.get_position().y1 + 0.015, text, fontsize=12, fontweight="bold")
