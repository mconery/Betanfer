library(ggplot2)
library(dplyr)
library(ggpubr)
library(data.table)

# ── Paths ──────────────────────────────────────────────────────────────────────
input_file <- "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545/mvp_matrix_comparison.tsv"
out_dir    <- "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545"

# ── Read data ──────────────────────────────────────────────────────────────────
df <- fread(input_file)

set.seed(42)
df_scatter <- df[sample(.N, min(1e6L, .N))]

# ── Join min MAF per variant pair ──────────────────────────────────────────────
maf_file <- "/flare/GeomicVar/mconery/Meta-Analysis_Fine-Mapping/mvp_results/RITM0322545/maf_lookup.tsv"
maf_df   <- fread(maf_file, na.strings = "NA")

# Pivot to long form (population, maf) for a single join key
maf_long <- melt(maf_df,
                 id.vars       = "rsid",
                 measure.vars  = c("EUR_maf", "AFR_maf"),
                 variable.name = "pop_col",
                 value.name    = "maf")
maf_long[, population := sub("_maf$", "", pop_col)]
maf_long[, pop_col    := NULL]

df_scatter[maf_long, on = c(variant_1 = "rsid", "population"), maf_v1 := i.maf]
df_scatter[maf_long, on = c(variant_2 = "rsid", "population"), maf_v2 := i.maf]
df_scatter[, min_maf := pmin(maf_v1, maf_v2)]

r2_val <- cor(df_scatter$mvp_r, df_scatter$betanfer_r, use = "complete.obs")^2

maf_breaks  <- c(0, 0.001, 0.01, 0.1, Inf)
maf_labels  <- c("<0.001", "0.001-<0.01", "0.01-<0.1", ">=0.1")
maf_palette <- setNames(scales::viridis_pal(direction = -1)(4), maf_labels)
df_scatter[, maf_bucket := cut(min_maf, breaks = maf_breaks, labels = maf_labels,
                                include.lowest = TRUE, right = FALSE)]

###############################################################################
# Plot 1: Scatterplot — betanfer_r (y) vs mvp_r (x), colored by min MAF
###############################################################################

p1 <- ggplot(df_scatter, aes(x = mvp_r, y = betanfer_r, color = min_maf)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_viridis_c(name = "Min MAF", na.value = "grey70") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("R² = %.4f", r2_val),
           hjust = -0.1, vjust = 1.5, size = 4) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Genotype Correlation",
    y = "Beta Correlation"
  ) +
  theme(legend.position = "right")

ggsave(
  filename    = file.path(out_dir, "scatter.betanfer_vs_mvp.tiff"),
  plot        = p1,
  width       = 200,
  height      = 80,
  units       = "mm",
  dpi         = 300,
  compression = "lzw"
)

###############################################################################
# Plot 1b: Scatterplot — same as p1 but MAF color bucketed by order of magnitude
###############################################################################

p1b <- ggplot(df_scatter, aes(x = mvp_r, y = betanfer_r, color = maf_bucket)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(name = "Min MAF", values = maf_palette, na.value = "grey70") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("R² = %.4f", r2_val),
           hjust = -0.1, vjust = 1.5, size = 4) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Genotype Correlation",
    y = "Beta Correlation"
  ) +
  theme(legend.position = "right")

ggsave(
  filename    = file.path(out_dir, "scatter.betanfer_vs_mvp.maf_bucket.tiff"),
  plot        = p1b,
  width       = 200,
  height      = 80,
  units       = "mm",
  dpi         = 300,
  compression = "lzw"
)

###############################################################################
# Plot 1c: Scatterplot — min MAF >= 0.01 only, R² computed on filtered points
###############################################################################

df_maf01  <- df_scatter[!is.na(min_maf) & min_maf >= 0.01]
r2_maf01  <- cor(df_maf01$mvp_r, df_maf01$betanfer_r, use = "complete.obs")^2

p1c <- ggplot(df_maf01, aes(x = mvp_r, y = betanfer_r, color = maf_bucket)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(name = "Min MAF", values = maf_palette, na.value = "grey70") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("R² = %.4f", r2_maf01),
           hjust = -0.1, vjust = 1.5, size = 4) +
  theme_pubr(base_size = 12) +
  labs(
    x = "Genotype Correlation",
    y = "Beta Correlation"
  ) +
  theme(legend.position = "right")

ggsave(
  filename    = file.path(out_dir, "scatter.betanfer_vs_mvp.maf_filtered.tiff"),
  plot        = p1c,
  width       = 200,
  height      = 80,
  units       = "mm",
  dpi         = 300,
  compression = "lzw"
)

###############################################################################
# Plot 2: MAE by absolute genotype correlation bucket, colored by population
###############################################################################

bucket_breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
bucket_labels <- c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")

mae_df <- df %>%
  mutate(bucket_range = cut(abs(mvp_r),
                            breaks        = bucket_breaks,
                            labels        = bucket_labels,
                            include.lowest = TRUE)) %>%
  filter(!is.na(bucket_range)) %>%
  group_by(population, bucket_range) %>%
  summarize(mae = mean(abs(delta_r), na.rm = TRUE), .groups = "drop") %>%
  mutate(bucket_range = factor(bucket_range, levels = bucket_labels))

p2 <- ggplot(mae_df, aes(x = bucket_range, y = mae, color = population)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_pubr(base_size = 12) +
  labs(
    x     = "Genotype Correlation",
    y     = "Mean Absolute Error (MAE)",
    color = "Population"
  ) +
  theme(
    legend.position  = "bottom",
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.margin      = margin(t = 7, r = 5, b = 5, l = 5, unit = "mm")
  )

ggsave(
  filename    = file.path(out_dir, "MAE.by_genotype_bucket.tiff"),
  plot        = p2,
  width       = 133,
  height      = 100,
  units       = "mm",
  dpi         = 300,
  compression = "lzw"
)
