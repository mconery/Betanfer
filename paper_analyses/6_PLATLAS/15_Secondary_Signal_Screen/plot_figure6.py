#!/usr/bin/env python3
'''
plot_figure6.py

Draws Figure 6: three loci where a secondary signal is resolved only when LD between
variants is modelled (RSparsePro, strict-intersection LD panel), compared against the
published single-causal ABF credible set that cannot separate it.

The three panels (gene, trait, locus window, which variant is the ABF lead, which is
the secondary lead, and any extra labelled variant) are a fixed, hand-picked selection
and are read from --panels rather than hard-coded; see figure6_panels.tsv in this
directory for the default set (BMP2/glaucoma, BAIAP2L2/hearing loss, TGFB2/abdominal
aortic aneurysm).

Each panel has three tracks:
  1. -log10(P) association plot, points coloured by LD with the ABF lead (blue) and
     with the secondary lead (magenta); grey = r2 < 0.2 with both, or not fine-mapped.
  2. Credible-set membership: "Without LD" = the original ABF 95% CS; "With LD" = the
     RSparsePro 95% CS, split into the CS matching the ABF signal (blue) and the
     secondary CS the ABF fit cannot see (magenta). A CS from a multi-causal fit
     matches an ABF CS when their PIP-weighted Jaccard similarity exceeds MATCH_THRESHOLD
     (sum(min(pip_a, pip_b)) / (mass_a + mass_b - sum(min)), same definition used to
     screen for secondary signals).
  3. Protein-coding gene track (GENCODE, Ensembl-canonical transcripts).

INPUT:
  --panels            Hand-picked panel definitions (default: figure6_panels.tsv, this directory)
  --cs-members        Credible-set membership TSV (output of 11_Credible_Set_Collation)
  --rsparsepro-dir     Per-locus RSparsePro output, files <trait>.<locus>.rsparsepro.txt
  --abf-dir           Directory of per-trait original ABF CS files, <trait>.<pop>.cs.txt.gz
  --meta-analysis-dir Per-population meta-analysis summary statistics, <pop>/<trait>.out.gz (tabix)
  --ld-dir            Per-trait/population LD matrices, <trait>.<pop>/<locus>*.vcor1(.vars)
  --gtf               GENCODE GTF (tabix-indexed), used for the gene track
  --fonts-dir         Directory holding Arial (regular) and, optionally, bold/italic/bold-italic
                      faces; missing bold/italic/bold-italic faces fall back to the metric-
                      compatible Liberation Sans faces registered under the name "Arial".
                      If Arial itself is absent, falls back to the default sans-serif font.
  --pop               Population (default: ALL)

OUTPUT:
  --output-prefix  Figure path prefix; writes <prefix>.png and <prefix>.pdf
  --output-points  Per-panel plotted CS points (both tracks), gzipped TSV

Usage:
  python plot_figure6.py [OPTIONS]

***Run with an environment that has pandas, numpy, matplotlib and pysam
   (e.g. $HOME/miniforge3/envs/betainfer)***
'''

import sys
import os
import getopt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import pysam

HERE = Path(__file__).resolve().parent
UPSTREAM = Path("/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping")
PANELS_DEFAULT = HERE / "figure6_panels.tsv"
CS_MEMBERS_DEFAULT = UPSTREAM / "application_results/credible_sets/gwas_cs_members.tsv"
RSPARSEPRO_DIR_DEFAULT = UPSTREAM / "application_results/intersection_finemap/rsparsepro/ALL"
ABF_DIR_DEFAULT = Path("/lus/flare/projects/GeomicVar/jianhuawang/betainfer/data/raw/cs_mrmega")
META_ANALYSIS_DIR_DEFAULT = UPSTREAM / "meta-analysis_results"
LD_DIR_DEFAULT = UPSTREAM / "LD_matrices"
GTF_DEFAULT = Path("/lus/flare/projects/GeomicVar/alphagenome_dataset/gencode_v48/gencode.v48.annotation.sorted.gtf.gz")
FONTS_DIR_DEFAULT = Path.home() / ".local/share/fonts"
OUTPUT_PREFIX_DEFAULT = UPSTREAM / "application_results/figures/Figure6"
OUTPUT_POINTS_DEFAULT = UPSTREAM / "application_results/figures/figure6_points.tsv.gz"

MATCH_THRESHOLD = 0.10
TOOL = "rsparsepro"   # the only multi-causal method drawn ("With LD" track)

INK, INK_SECONDARY, INK_MUTED = "#1A1A1A", "#4D4D4D", "#8A8A8A"
PRIMARY = "#1f5fa8"      # CS matching the ABF CS
SECONDARY = "#d81b60"    # CS the ABF fit cannot see: the figure's accent colour
LIGHT = {PRIMARY: "#8fb4dc", SECONDARY: "#f0a3c0"}
GREY = "#9a9a9a"
FS, FS_TICK = 8, 7
TRACK_WIDTH_IN = 2.95     # width of the association track, for gene-label packing
FIXED_ROW_RATIOS = [1.7, .42, .26, .05]   # assoc / x-axis gap / CS track / spacer
GENE_ROW_RATIO = .2                        # per gene-track row; panel height scales with gene rows


def help(exit_num=0):
    print("""Usage: python plot_figure6.py [OPTIONS]
Options:
  --panels             Panel definitions TSV (default: %s)
  --cs-members         Credible-set membership TSV (default: %s)
  --rsparsepro-dir     Per-locus RSparsePro output directory (default: %s)
  --abf-dir            Per-trait original ABF CS directory (default: %s)
  --meta-analysis-dir  Per-population meta-analysis sumstats directory (default: %s)
  --ld-dir             Per-trait/population LD matrix directory (default: %s)
  --gtf                GENCODE GTF, tabix-indexed (default: %s)
  --fonts-dir          Directory with Arial / Liberation Sans faces (default: %s)
  --pop                Population (default: ALL)
  --output-prefix      Figure output path prefix (default: %s)
  --output-points      Plotted CS points output path (default: %s)
  -h, --help           Show this help message""" % (
        PANELS_DEFAULT, CS_MEMBERS_DEFAULT, RSPARSEPRO_DIR_DEFAULT, ABF_DIR_DEFAULT,
        META_ANALYSIS_DIR_DEFAULT, LD_DIR_DEFAULT, GTF_DEFAULT, FONTS_DIR_DEFAULT,
        OUTPUT_PREFIX_DEFAULT, OUTPUT_POINTS_DEFAULT))
    sys.exit(exit_num)


###############################################################################
##########  PANEL DEFINITIONS  ################################################
###############################################################################

def load_panels(path):
    """Read the hand-picked panel table into the dict shape the rest of this script uses."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    panels = []
    for _, row in df.iterrows():
        ex = dict(panel=row["panel"], gene=row["gene"], trait=row["trait"], title=row["title"],
                  locus=row["locus"], abf_locus=row["abf_locus"],
                  window=(int(row["window_start"]), int(row["window_end"])),
                  abf_lead=(row["abf_lead_variant"], row["abf_lead_rsid"], row["abf_lead_label"]),
                  secondary=(row["secondary_variant"], row["secondary_rsid"], row["secondary_label"]),
                  gene_rows=int(row["gene_rows"]) if row["gene_rows"] else 1)
        if row["show_genes"]:
            ex["show_genes"] = tuple(row["show_genes"].split(","))
        if row["extra"]:
            ex["extra"] = [(v, l, (float(dx), float(dy)))
                            for v, l, dx, dy in (item.split("|") for item in row["extra"].split(";"))]
        panels.append(ex)
    return panels


def ratios(ex):
    return FIXED_ROW_RATIOS + [GENE_ROW_RATIO * ex.get("gene_rows", 1) + .05]


###############################################################################
##########  GENERIC TABIX / GTF HELPERS  ######################################
###############################################################################

def tabix_frame(path, chrom, start, end):
    with pysam.TabixFile(str(path)) as tabix:
        header = tabix.header[-1].lstrip("#")
        lines = list(tabix.fetch(chrom, start - 1, end))
    if not lines:
        return pd.DataFrame()
    return pd.read_csv(pd.io.common.StringIO(header + "\n" + "\n".join(lines)), sep="\t")


def parse_attributes(blob, key):
    token = f'{key} "'
    if token not in blob:
        return None
    return blob.split(token, 1)[1].split('"', 1)[0]


def gene_features(gtf_path, chrom, start, end):
    """Ensembl-canonical transcript of each named protein-coding gene in the window."""
    spans, exons = {}, {}
    with pysam.TabixFile(str(gtf_path)) as gtf:
        for line in gtf.fetch("chr" + chrom, start - 1, end):
            f = line.split("\t")
            if f[2] not in ("transcript", "exon") or 'gene_type "protein_coding"' not in f[8] \
                    or 'tag "Ensembl_canonical"' not in f[8]:
                continue
            name = parse_attributes(f[8], "gene_name")
            if not name or name.startswith("ENSG"):
                continue
            if f[2] == "transcript":
                spans[name] = (int(f[3]), int(f[4]))
            else:
                exons.setdefault(name, []).append((int(f[3]), int(f[4])))
    return spans, exons


def pack_genes(spans, start, end):
    """Greedy rows; a gene's footprint is its line plus its centred label, kept inside the window."""
    bp_per_char = (end - start) / TRACK_WIDTH_IN * (.62 * FS_TICK / 72)
    rows, placed = [], []
    for name, (a, b) in sorted(spans.items(), key=lambda kv: kv[1][0]):
        x1, x2 = max(a, start), min(b, end)
        half = len(name) * bp_per_char / 2
        mid = min(max((x1 + x2) / 2, start + half), end - half)
        lo, hi = min(x1, mid - half), max(x2, mid + half)
        row = next((i for i, right in enumerate(rows) if lo > right + bp_per_char), len(rows))
        rows[row:row + 1] = [hi]
        placed.append(dict(name=name, x1=x1, x2=x2, mid=mid, row=row))
    return placed


###############################################################################
##########  READ: ASSOCIATION, LD, CS MEMBERSHIP  #############################
###############################################################################

def read_assoc(meta_analysis_dir, trait, chrom, start, end, pop):
    path = Path(meta_analysis_dir) / pop / f"{trait}.out.gz"
    frame = tabix_frame(path, chrom, start, end).rename(
        columns={"rs_number": "variant_id", "POS": "pos", "p-value": "p_value", "_-log10_p-value": "neglog10p"})
    frame = frame[["variant_id", "pos", "p_value", "neglog10p"]].dropna()
    return frame.sort_values("p_value").drop_duplicates("variant_id")


def lead_ld(ld_dir, trait, locus, lead, assoc, pop):
    """r2 with `lead` from the sample-size-weighted multi-cohort ("merged") LD panel --
    the panel whose LD with each fine-mapped variant equals the strict-intersection LD
    RSparsePro used. Ties among merged panels are broken by variant coverage, then MVP."""
    ranked = []
    for matrix in (Path(ld_dir) / f"{trait}.{pop}").glob(f"{locus}*.vcor1"):
        if ".merged." not in matrix.name:
            continue
        names = Path(str(matrix) + ".vars").read_text().split()
        if lead in names:
            ranked.append((int(assoc.variant_id.isin(names).sum()), int(".MVP." in matrix.name), matrix, names))
    _, _, matrix, names = max(ranked, key=lambda r: (r[0], r[1]))
    row = names.index(lead)
    with matrix.open() as handle:
        for i, line in enumerate(handle):
            if i == row:
                r2 = np.array(line.split(), dtype=float) ** 2
                break
    return dict(zip(names, r2))


def abf_pips(ex, abf_dir, pop):
    a = pd.read_csv(Path(abf_dir) / f"{ex['trait']}.{pop}.cs.txt.gz", sep="\t",
                     usecols=["LOCUS", "ID", "PIP_1", "CS95_1"])
    a = a[a.LOCUS.eq(ex["abf_locus"])]
    return pd.DataFrame(dict(variant_id=a.ID, pip=a.PIP_1, cs=np.where(a.CS95_1, "abf", None)))


def rsparsepro_pips(rsparsepro_dir, ex, members):
    path = Path(rsparsepro_dir) / f"{ex['trait']}.{ex['locus']}.rsparsepro.txt"
    raw = pd.read_csv(path, sep="\t", usecols=["VARIANT_ID", "PIP"]).rename(
        columns={"VARIANT_ID": "variant_id", "PIP": "pip"})
    cs = members[members.tool.eq(TOOL)][["variant_id", "cs"]].astype({"cs": str})
    return raw.merge(cs, on="variant_id", how="left")


def classify(frames, secondary_lead):
    """Label each RSparsePro CS: matched to the ABF CS (Jaccard > MATCH_THRESHOLD), holding
    the secondary lead, or other; the ABF CS itself is simply "primary"."""
    abf = frames["abf"].query("cs == 'abf'").set_index("variant_id").pip
    out = {"abf": frames["abf"].assign(kind=lambda d: np.where(d.cs.eq("abf"), "primary", None))}
    d = frames[TOOL]
    kind = {}
    for cs, g in d.dropna(subset=["cs"]).groupby("cs"):
        w = g.set_index("variant_id").pip
        shared = np.minimum(w, abf.reindex(w.index).fillna(0)).sum()
        j = shared / (w.sum() + abf.sum() - shared)
        kind[cs] = "primary" if j > MATCH_THRESHOLD else ("secondary" if secondary_lead in set(g.variant_id) else "other")
    out[TOOL] = d.assign(kind=d.cs.map(kind))
    return out


def prepare(ex, members, paths):
    sec = ex["secondary"][0]
    chrom = ex["abf_lead"][0].split(":")[0]
    start, end = ex["window"]
    assoc = read_assoc(paths["meta_analysis_dir"], ex["trait"], chrom, start, end, paths["pop"])
    for key, lead in (("r2_abf", ex["abf_lead"][0]), ("r2_sec", sec)):
        r2 = lead_ld(paths["ld_dir"], ex["trait"], ex["locus"], lead, assoc, paths["pop"])
        assoc[key] = assoc.variant_id.map(r2).astype(float)
    mem = members[members.trait.eq(ex["trait"]) & members.locus.eq(ex["locus"])]
    frames = {"abf": abf_pips(ex, paths["abf_dir"], paths["pop"]),
              TOOL: rsparsepro_pips(paths["rsparsepro_dir"], ex, mem)}
    pips = {m: d.assign(pos=d.variant_id.str.split(":").str[1].astype(int)).query("@start <= pos <= @end")
            for m, d in classify(frames, sec).items()}
    spans, exons = gene_features(paths["gtf"], chrom, start, end)
    return dict(ex=ex, assoc=assoc, pips=pips,
                genes=dict(start=start, end=end, chrom=chrom, spans=spans, exons=exons, gene=ex["gene"]))


###############################################################################
##########  DRAWING  ###########################################################
###############################################################################

def style(ax):
    ax.tick_params(labelsize=FS_TICK, width=.7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(.7)


def two_lead_colors(a):
    color = np.full(len(a), GREY, dtype=object)
    for key, dark in (("r2_abf", PRIMARY), ("r2_sec", SECONDARY)):
        color[(a[key] >= .2).to_numpy()] = LIGHT[dark]
        color[(a[key] >= .6).to_numpy()] = dark
    return color


def draw_assoc(ax, item):
    a = item["assoc"]
    used = a.variant_id.isin(item["pips"][TOOL].variant_id)
    # Variants not used in fine-mapping are folded into the r2 < 0.2 class (not a separate legend entry).
    ax.scatter(a.pos[~used] / 1e6, a.neglog10p[~used], s=4, color=GREY, linewidths=0, rasterized=True, zorder=1)
    u = a[used]
    ax.scatter(u.pos / 1e6, u.neglog10p, s=7, c=two_lead_colors(u), linewidths=0, rasterized=True, zorder=2)
    x_of, y_of = dict(zip(a.variant_id, a.pos / 1e6)).get, dict(zip(a.variant_id, a.neglog10p)).get
    lo, hi = (w / 1e6 for w in item["ex"]["window"])
    # Both leads: filled diamond + a short leader line to the rsID.
    for (variant, rsid, _), color in ((item["ex"]["abf_lead"], PRIMARY), (item["ex"]["secondary"], SECONDARY)):
        x, y = x_of(variant), y_of(variant)
        ax.scatter([x], [y], s=42, marker="D", color=color, edgecolors="white", linewidths=.6, zorder=7)
        right = x > lo + .7 * (hi - lo)
        ax.annotate(rsid, (x, y), xytext=(-16 if right else 16, 9), textcoords="offset points",
                    fontsize=FS_TICK, color=color, ha="right" if right else "left", va="center",
                    arrowprops=dict(arrowstyle="-", color=color, lw=.8, shrinkA=1, shrinkB=3), zorder=8)
    # Other coding variants inside the secondary CS worth calling out: small diamond, same color.
    for variant, label, offset in item["ex"].get("extra", []):
        x, y = x_of(variant), y_of(variant)
        ax.scatter([x], [y], s=24, marker="D", color=SECONDARY, edgecolors="white", linewidths=.6, zorder=7)
        ax.annotate(label, (x, y), xytext=offset, textcoords="offset points", fontsize=FS_TICK, color=SECONDARY,
                    ha="right" if offset[0] < 0 else "left", va="center",
                    arrowprops=dict(arrowstyle="-", color=SECONDARY, lw=.8, shrinkA=1, shrinkB=3), zorder=8)
    ax.axhline(-np.log10(5e-8), color=INK_MUTED, lw=.7, ls=(0, (3, 2)), zorder=1)
    ax.set_ylim(-.3, a.neglog10p.max() * 1.3)
    ax.set_ylabel(r"$-\log_{10}P$", fontsize=FS)
    style(ax)


def draw_cs(ax, item):
    """CS membership track: one row per method; vertical ticks mark 95% CS members."""
    d = item["pips"][TOOL].dropna(subset=["kind"])
    n = {k: d.loc[d.kind.eq(k), "cs"].nunique() for k in ("primary", "secondary", "other")}
    if n != dict(primary=1, secondary=1, other=0):
        print(f"ERROR: expected exactly one primary and one secondary {TOOL} CS, got {n}")
        sys.exit(1)
    f = item["pips"]["abf"].dropna(subset=["kind"])
    ax.vlines(f.pos / 1e6, .62, .98, colors=INK, lw=.9)
    ax.vlines(d.pos / 1e6, .02, .38, colors=d.kind.map({"primary": PRIMARY, "secondary": SECONDARY}), lw=.9)
    ax.set_ylim(-.1, 1.1)
    ax.set_yticks([.8, .2], ["Without LD", "With LD"], fontsize=FS_TICK)
    ax.tick_params(axis="y", length=0)
    style(ax)
    ax.spines[["top", "right"]].set_visible(True)   # the CS track gets a full box
    for side in ax.spines.values():
        side.set_linewidth(.6)


def draw_genes(ax, g, rows):
    """Gene track: every gene in black except the panel's target gene, in the accent color."""
    for p in pack_genes(g["spans"], g["start"], g["end"]):
        hit = p["name"] == g["gene"]
        color, y = (SECONDARY if hit else "black"), -p["row"]
        ax.plot([p["x1"] / 1e6, p["x2"] / 1e6], [y, y], color=color, lw=.8, solid_capstyle="butt")
        for a, b in g["exons"].get(p["name"], []):
            a, b = max(a, g["start"]), min(b, g["end"])
            if a <= b:
                ax.add_patch(Rectangle((a / 1e6, y - .13), (b - a) / 1e6, .26, color=color, lw=0))
        ax.text(p["mid"] / 1e6, y + .22, p["name"], ha="center", va="bottom", fontsize=FS_TICK,
                style="italic", color=color, fontweight="bold" if hit else "normal")
    ax.set_ylim(-(rows - 1) - .5, .95)
    ax.set_yticks([])
    ax.spines[["left", "right", "top", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)


def draw_panel(sub, item, letter):
    ex = item["ex"]
    gs = sub.add_gridspec(5, 1, left=.09, right=.98, top=.95, bottom=.02, height_ratios=ratios(ex), hspace=0)
    sub.text(.005, .995, letter, fontsize=13, weight="bold", va="top")
    top = sub.add_subplot(gs[0])
    draw_assoc(top, item)
    top.text(.012, .97, ex["title"], transform=top.transAxes, fontsize=FS, va="top", ha="left", color="black")
    cs = sub.add_subplot(gs[2], sharex=top)
    draw_cs(cs, item)
    genes = sub.add_subplot(gs[4], sharex=top)
    g = item["genes"]
    keep = ex.get("show_genes")
    spans = {k: v for k, v in g["spans"].items() if keep is None or k in keep}
    draw_genes(genes, dict(g, spans=spans), ex.get("gene_rows", 1))
    top.set_xlim(ex["window"][0] / 1e6, ex["window"][1] / 1e6)
    top.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(7))
    # Mb position is only shown under the association track; the CS and gene tracks share its x-axis.
    top.set_xlabel(f"Chr {g['chrom']} position (Mb)", fontsize=FS)
    top.tick_params(axis="x", labelbottom=True, labelsize=FS_TICK)
    cs.tick_params(axis="x", bottom=False, labelbottom=False)


def legend(fig):
    dot = lambda c, l, **k: Line2D([], [], marker="o", ls="none", color=c, ms=5, label=l, **k)
    handles = [dot(PRIMARY, r"$r^2 \geq 0.6$ with ABF lead"), dot(LIGHT[PRIMARY], "0.2–0.6"),
               dot(SECONDARY, r"$r^2 \geq 0.6$ with secondary lead"), dot(LIGHT[SECONDARY], "0.2–0.6"),
               dot(GREY, r"$r^2$ < 0.2 with both")]
    kw = dict(frameon=False, fontsize=FS_TICK, handletextpad=.3, columnspacing=1, loc="lower center")
    fig.legend(handles=handles, ncol=5, bbox_to_anchor=(.5, -.25 / fig.get_figheight()), **kw)


###############################################################################
##########  FONT SETUP  ########################################################
###############################################################################

def register_arial(fonts_dir):
    """Register Arial for matplotlib; bold/italic/bold-italic fall back to the metric-
    compatible Liberation Sans faces (registered under the name "Arial") when the real
    Arial face is not installed. If even the regular face is missing, leaves the default
    sans-serif font in place."""
    fonts_dir = Path(fonts_dir)
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
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "mathtext.fontset": "custom",
                          "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold"})


###############################################################################
##########  MAIN  #############################################################
###############################################################################

def main():
    try:
        opts, _ = getopt.getopt(
            sys.argv[1:], "h",
            ["help", "panels=", "cs-members=", "rsparsepro-dir=", "abf-dir=", "meta-analysis-dir=",
             "ld-dir=", "gtf=", "fonts-dir=", "pop=", "output-prefix=", "output-points="])
    except getopt.GetoptError as err:
        print(f"ERROR: {err}")
        sys.exit(2)

    panels_path = PANELS_DEFAULT
    paths = dict(cs_members=CS_MEMBERS_DEFAULT, rsparsepro_dir=RSPARSEPRO_DIR_DEFAULT, abf_dir=ABF_DIR_DEFAULT,
                 meta_analysis_dir=META_ANALYSIS_DIR_DEFAULT, ld_dir=LD_DIR_DEFAULT, gtf=GTF_DEFAULT, pop="ALL")
    fonts_dir = FONTS_DIR_DEFAULT
    output_prefix = OUTPUT_PREFIX_DEFAULT
    output_points = OUTPUT_POINTS_DEFAULT

    for o, a in opts:
        if o in ("-h", "--help"):
            help()
        elif o == "--panels":
            panels_path = a
        elif o == "--cs-members":
            paths["cs_members"] = a
        elif o == "--rsparsepro-dir":
            paths["rsparsepro_dir"] = a
        elif o == "--abf-dir":
            paths["abf_dir"] = a
        elif o == "--meta-analysis-dir":
            paths["meta_analysis_dir"] = a
        elif o == "--ld-dir":
            paths["ld_dir"] = a
        elif o == "--gtf":
            paths["gtf"] = a
        elif o == "--fonts-dir":
            fonts_dir = a
        elif o == "--pop":
            paths["pop"] = a
        elif o == "--output-prefix":
            output_prefix = a
        elif o == "--output-points":
            output_points = a
        else:
            assert False, f"Unhandled option: {o}"

    if not os.path.exists(panels_path):
        print(f"ERROR: --panels not found: {panels_path}")
        sys.exit(1)

    examples = load_panels(panels_path)
    traits = {e["trait"] for e in examples}
    print(f"NOTE: loading credible-set membership for {len(traits)} trait(s)")
    chunks = pd.read_csv(paths["cs_members"], sep="\t", chunksize=2_000_000, dtype={"cs": str})
    members = pd.concat(c[c["pop"].eq(paths["pop"]) & c.trait.isin(traits)] for c in chunks)

    print("NOTE: preparing panels")
    items = [prepare(ex, members, paths) for ex in examples]

    register_arial(fonts_dir)
    plt.rcParams.update({"font.size": FS, "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK,
                          "ytick.color": INK, "axes.edgecolor": INK_SECONDARY, "pdf.fonttype": 42})
    heights = [sum(ratios(i["ex"])) for i in items]
    fig = plt.figure(figsize=(7.2, 1.09 * sum(heights)))
    subfigs = fig.subfigures(len(items), 1, hspace=.02, height_ratios=heights)
    for item, sub in zip(items, np.atleast_1d(subfigs)):
        draw_panel(sub, item, item["ex"]["panel"])
    legend(fig)

    os.makedirs(os.path.dirname(str(output_prefix)), exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{output_prefix}.{ext}", dpi=220, bbox_inches="tight", facecolor="white")
    print(f"NOTE: figure written to {output_prefix}.png / .pdf")

    rows = [d.assign(gene=i["ex"]["gene"], method=m) for i in items for m, d in i["pips"].items()]
    os.makedirs(os.path.dirname(str(output_points)), exist_ok=True)
    pd.concat(rows).to_csv(output_points, sep="\t", index=False)
    print(f"NOTE: plotted CS points written to {output_points}")


if __name__ == "__main__":
    main()
