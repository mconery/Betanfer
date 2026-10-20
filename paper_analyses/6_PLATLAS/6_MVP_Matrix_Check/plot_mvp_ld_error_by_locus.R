#!/usr/bin/env Rscript
####################################################################################
#plot_mvp_ld_error_by_locus.R
#purpose: mean absolute error of betanfer LD against real MVP genotype LD, one
#         point per locus and population, coloured by the number of traits
#         betanfer used at that locus. Figure 1C layout: theme_minimal, boxed
#         legend inside the top-right corner.
#environment: Conda r_plotting on Aurora (ggplot2, data.table)
####################################################################################
#
# Input:  <in_dir>/mvp_ld_error_by_locus_bin.tsv  (extract_mvp_ld_error_by_locus_bin.py)
#         <in_dir>/mvp_traits_per_locus.tsv        (collate_mvp_traits_per_locus.py)
# Output: <in_dir>/mvp_ld_error_by_locus_bin_with_traits.tsv  (merged source data)
#         <in_dir>/mvp_ld_error_by_locus.{png,tiff}
#
# Usage:  Rscript plot_mvp_ld_error_by_locus.R [png|tiff]
#         png = quick preview at 300 dpi (default); tiff = final 600 dpi LZW

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

# ── Paths ──────────────────────────────────────────────────────────────────────
in_dir <- "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545"

args <- commandArgs(trailingOnly = TRUE)
fmt  <- if (length(args)) args[1] else "png"

###############################################################################
## Step 1: Merge the per-bin MAE table with the per-locus trait counts
###############################################################################

d <- merge(fread(file.path(in_dir, "mvp_ld_error_by_locus_bin.tsv")),
           fread(file.path(in_dir, "mvp_traits_per_locus.tsv")),
           by = c("locus", "population"))[order(bin, population, locus)]
fwrite(d, file.path(in_dir, "mvp_ld_error_by_locus_bin_with_traits.tsv"), sep = "\t")

###############################################################################
## Step 2: Plot — one point per locus x population, coloured by # traits
###############################################################################

set.seed(1)
d <- d[sample(.N)]   # shuffle plotting order so one population does not hide the other

p <- ggplot(d, aes(x = bin, y = mae, color = n_traits)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.30), breaks = seq(0, 0.30, 0.05), labels = sprintf("%.2f", seq(0, 0.30, 0.05))) +
  guides(color = guide_colorbar(barheight = unit(0.9, "cm"), barwidth = unit(0.35, "cm"))) +
  labs(x = "Genotype Correlation", y = "Mean Absolute Error (MAE)", color = "# Traits") +
  theme(axis.title = element_text(size = 11),
        strip.text = element_blank(),
        legend.position = c(1, 1),
        legend.title = element_text(size = 11, hjust = 0.75),
        legend.text = element_text(size = 8), legend.margin = margin(3, 4, 3, 4),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.25),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

f <- file.path(in_dir, paste0("mvp_ld_error_by_locus.", fmt))
if (fmt == "tiff") {
  ggsave(f, p, width = 2.25, height = 3, units = "in", dpi = 600, compression = "lzw")
} else {
  ggsave(f, p, width = 2.25, height = 3, units = "in", dpi = 300)
}
message("wrote ", f)
