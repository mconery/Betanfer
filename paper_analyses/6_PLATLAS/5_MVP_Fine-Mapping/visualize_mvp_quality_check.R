# MVP Fine-Mapping Quality Check Visualizations
# Author: Mitch Conery
# Date: May 2026
#
# Produces two figures comparing original MVP fine-mapping results against
# new SuSiE and RSparsePro analyses using updated LD matrices (AFR + EUR only):
#   1. Pie chart of credible-set signal recovery (including signals not recovered)
#   2. PIP correlation scatterplot, faceted by method

# 0) Libraries ====

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# 1) Paths ====

inp_dir <- "C:/Users/mitch/Documents/Argonne/Meta-Analysis_Fine-Mapping/Data/MVP_Quality_Check"
out_dir <- "C:/Users/mitch/Documents/Argonne/Meta-Analysis_Fine-Mapping/Figures/MVP_Check"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

signal_file_betanfer <- file.path(inp_dir, "mvp_collation.signal_recovery.tsv")
pip_file_betanfer    <- file.path(inp_dir, "mvp_collation.pip_comparison.tsv")
signal_file_1000g    <- file.path(inp_dir, "mvp_1000g_collation.signal_recovery.tsv")
pip_file_1000g       <- file.path(inp_dir, "mvp_1000g_collation.pip_comparison.tsv")

# 2) Read data ====

signal_df <- bind_rows(
  read_tsv(signal_file_betanfer, show_col_types = FALSE) %>% mutate(ld_source = "Betanfer"),
  read_tsv(signal_file_1000g,    show_col_types = FALSE) %>% mutate(ld_source = "1000G")
) %>% mutate(ld_source = factor(ld_source, levels = c("Betanfer", "1000G")))

pip_df <- bind_rows(
  read_tsv(pip_file_betanfer, show_col_types = FALSE, na = "NA") %>% mutate(ld_source = "Betanfer"),
  read_tsv(pip_file_1000g,    show_col_types = FALSE, na = "NA") %>% mutate(ld_source = "1000G")
) %>% mutate(ld_source = factor(ld_source, levels = c("Betanfer", "1000G")))

# 3) Figure 1 — Pie chart of credible-set signal recovery ====

pie_df <- signal_df %>%
  mutate(
    category = case_when(
      any_variant_in_susie_cs  & any_variant_in_rsparsepro_cs  ~ "Both",
      any_variant_in_susie_cs  & !any_variant_in_rsparsepro_cs ~ "SuSiE Only",
      !any_variant_in_susie_cs & any_variant_in_rsparsepro_cs  ~ "RSparsePro Only",
      TRUE ~ "Neither"
    ),
    category = factor(category, levels = c("Both", "SuSiE Only", "RSparsePro Only", "Neither"))
  ) %>%
  count(ld_source, category) %>%
  group_by(ld_source) %>%
  mutate(
    pct   = n / sum(n),
    label = paste0(category, "\n", n, " (", round(pct * 100, 1), "%)")
  ) %>%
  ungroup()

# x = 2 sets the pie ring at radius 1.5–2.5; labels at x = 3.5 sit outside it.
# xlim lower bound of -1 prevents label clipping on the left.
p_pie <- ggplot(pie_df, aes(x = 1, y = pct, fill = category)) +
  geom_col(width = 2, color = "white", linewidth = 0.6) +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(x = 3.2, label = label),
    position = position_stack(vjust = 0.5),
    size     = 4,
    color    = "black"
  ) +
  xlim(c(0, 5)) +
  facet_wrap(~ ld_source) +
  scale_fill_manual(
    values = c(
      "Both"            = "darkgreen",
      "SuSiE Only"      = "#7ad151",
      "RSparsePro Only" = "#22a884",
      "Neither"         = "lightgrey"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 14, face = "bold")
  )

ggsave(
  file.path(out_dir, "mvp_cs_overlap_pie.tiff"),
  plot   = p_pie,
  width  = 14,
  height = 7,
  dpi    = 300
)

# 3b) Figure 1b — Pie chart faceted by median credible-set MAF ====

pie_df_maf <- signal_df %>%
  filter(!is.na(median_maf)) %>%
  mutate(
    maf_group = if_else(
      median_maf < 0.05,
      "Low Frequency (Median MAF < 0.05)",
      "Common (Median MAF ≥ 0.05)"
    ),
    maf_group = factor(maf_group, levels = c(
      "Low Frequency (Median MAF < 0.05)",
      "Common (Median MAF ≥ 0.05)"
    )),
    category = case_when(
      any_variant_in_susie_cs  & any_variant_in_rsparsepro_cs  ~ "Both",
      any_variant_in_susie_cs  & !any_variant_in_rsparsepro_cs ~ "SuSiE Only",
      !any_variant_in_susie_cs & any_variant_in_rsparsepro_cs  ~ "RSparsePro Only",
      TRUE ~ "Neither"
    ),
    category = factor(
      category,
      levels = c("Both", "SuSiE Only", "RSparsePro Only", "Neither")
    )
  ) %>%
  count(ld_source, maf_group, category) %>%
  group_by(ld_source, maf_group) %>%
  mutate(
    pct   = n / sum(n),
    label = paste0(category, "\n", n, " (", round(pct * 100, 1), "%)")
  ) %>%
  ungroup()

p_pie_maf <- ggplot(pie_df_maf, aes(x = 1, y = pct, fill = category)) +
  geom_col(width = 2, color = "white", linewidth = 0.6) +
  coord_polar(theta = "y", start = 0) +
  geom_text(
    aes(x = 3.2, label = label),
    position = position_stack(vjust = 0.5),
    size     = 4,
    color    = "black"
  ) +
  xlim(c(0, 5)) +
  facet_grid(maf_group ~ ld_source) +
  scale_fill_manual(
    values = c(
      "Both"            = "darkgreen",
      "SuSiE Only"      = "#7ad151",
      "RSparsePro Only" = "#22a884",
      "Neither"         = "lightgrey"
    )
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 14, face = "bold")
  )

ggsave(
  file.path(out_dir, "mvp_cs_overlap_pie_by_maf.tiff"),
  plot   = p_pie_maf,
  width  = 28,
  height = 14,
  dpi    = 300
)

# 4) Figure 2 — PIP correlation scatterplot ====

pip_long <- pip_df %>%
  pivot_longer(
    cols      = c(susie_pip, rsparsepro_pip),
    names_to  = "method",
    values_to = "new_pip"
  ) %>%
  mutate(method = case_when(
    method == "susie_pip"      ~ "SuSiE",
    method == "rsparsepro_pip" ~ "RSparsePro"
  )) %>%
  filter(!is.na(original_pip), !is.na(new_pip))

cor_stats <- pip_long %>%
  group_by(ld_source, method) %>%
  summarize(
    r  = cor(original_pip, new_pip, use = "complete.obs"),
    r2 = cor(original_pip, new_pip, use = "complete.obs")^2,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", round(r, 3), "\nR² = ", round(r2, 3)))

set.seed(5)
pip_plot <- pip_long %>% slice_sample(n = min(20000, nrow(pip_long)))

pip_theme <- theme_minimal() +
  theme(
    axis.title       = element_text(size = 18),
    axis.text        = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black"),
    strip.text       = element_text(size = 17, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing.x  = unit(1.5, "cm"),
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5)
  )

pip_base <- list(
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black"),
  geom_text(
    data        = cor_stats,
    aes(x = 0.05, y = 0.92, label = label),
    hjust       = 0,
    size        = 5,
    inherit.aes = FALSE
  ),
  facet_grid(ld_source ~ method),
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)),
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)),
  labs(x = "Original PIP", y = "New PIP",
       title = "PIP Correlation: Original vs. New MVP Fine-Mapping")
)

p_pip <- ggplot(pip_plot, aes(x = original_pip, y = new_pip)) +
  geom_point(aes(color = method), alpha = 0.4, size = 0.8) +
  scale_color_manual(
    values = c("SuSiE" = "#7ad151", "RSparsePro" = "#22a884")
  ) +
  pip_base +
  pip_theme +
  theme(legend.position = "none")

ggsave(
  file.path(out_dir, "mvp_pip_correlation.tiff"),
  plot   = p_pip,
  width  = 20,
  height = 10,
  dpi    = 300
)

p_pip_maf <- ggplot(pip_plot, aes(x = original_pip, y = new_pip)) +
  geom_point(aes(color = maf), alpha = 0.4, size = 0.8) +
  scale_color_viridis_c(name = "MAF", option = "plasma") +
  pip_base +
  pip_theme +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 14)
  )

ggsave(
  file.path(out_dir, "mvp_pip_correlation_maf.tiff"),
  plot   = p_pip_maf,
  width  = 20,
  height = 10,
  dpi    = 300
)

# 5) SuSiE convergence failure summary ====

convergence_df <- pip_df %>%
  group_by(ld_source, Trait, Population, Locus) %>%
  summarize(
    mapped       = any(!is.na(rsparsepro_pip)),
    no_susie_pip = all(is.na(susie_pip)),
    has_cs_zero  = any(!is.na(susie_cs) & susie_cs == 0),
    .groups      = "drop"
  ) %>%
  filter(mapped) %>%
  mutate(failed = no_susie_pip & has_cs_zero)

convergence_df %>%
  group_by(ld_source) %>%
  summarize(
    n_total  = n(),
    n_failed = sum(failed),
    pct      = round(n_failed / n_total * 100, 1),
    .groups  = "drop"
  ) %>%
  rowwise() %>%
  group_walk(~ message(sprintf(
    "[%s] SuSiE convergence failures: %d / %d mapped loci (%.1f%%)",
    .y$ld_source, .x$n_failed, .x$n_total, .x$pct
  )))

message("Done. Figures written to: ", out_dir)
