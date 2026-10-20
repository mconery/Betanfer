#!/usr/bin/env Rscript
####################################################################################
#plot_mvp_ld_error_scatter.R
#purpose: scatter of betanfer beta correlation (y) against real MVP genotype
#         correlation (x), one point per sampled variant pair, annotated with
#         the exact Pearson r computed over every pair.
#environment: Conda r_plotting on Aurora (ggplot2, data.table)
####################################################################################
#
# Input:  <in_dir>/mvp_ld_error_sample.tsv.gz    (extract_mvp_ld_error_sample.py)
#         <in_dir>/mvp_ld_error_moments.tsv      (extract_mvp_ld_error_sample.py)
# Output: <in_dir>/mvp_ld_error_scatter.{png,tiff}
#
# Usage:  Rscript plot_mvp_ld_error_scatter.R [png|tiff]
#         png = quick preview at 300 dpi (default); tiff = final 600 dpi LZW

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

# ── Paths ──────────────────────────────────────────────────────────────────────
in_dir <- "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545"

args <- commandArgs(trailingOnly = TRUE)
fmt  <- if (length(args)) args[1] else "png"

smp <- fread(cmd = paste("zcat", file.path(in_dir, "mvp_ld_error_sample.tsv.gz")))
mom <- fread(file.path(in_dir, "mvp_ld_error_moments.tsv"))

# ── Plot style ───────────────────────────────────────────────────────────────────
NAVY  <- "#132B43"
GUIDE <- "grey80"
FONT  <- "Liberation Sans"
theme_mitch <- theme_classic(base_size = 12, base_family = FONT) +
  theme(axis.text = element_text(colour = "black", size = 11),
        axis.title = element_text(colour = "black", size = 12),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        plot.margin = margin(8, 10, 6, 6))

###############################################################################
## Step 1: Scatterplot of betanfer beta correlation vs MVP genotype correlation
###############################################################################

r_all <- mom[population == "ALL", pearson_r]
p <- ggplot(smp, aes(mvp_r, betanfer_r)) +
  geom_point(colour = NAVY, size = 0.35, stroke = 0) +
  geom_abline(slope = 1, intercept = 0, colour = GUIDE, linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = -1, y = 1, hjust = 0, vjust = 1, family = FONT, size = 11 / .pt,
           label = sprintf("r = %.2f", r_all)) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), labels = sprintf("%.1f", seq(-1, 1, 0.5))) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), labels = sprintf("%.1f", seq(-1, 1, 0.5))) +
  coord_fixed() +
  labs(x = "Genotype Correlation", y = "Beta Correlation") +
  theme_mitch

f <- file.path(in_dir, paste0("mvp_ld_error_scatter.", fmt))
if (fmt == "tiff") {
  ggsave(f, p, width = 3.4, height = 3.3, units = "in", dpi = 600, device = "tiff",
         compression = "lzw", type = "cairo")
} else {
  ggsave(f, p, width = 3.4, height = 3.3, units = "in", dpi = 300, type = "cairo")
}
message("wrote ", f)
