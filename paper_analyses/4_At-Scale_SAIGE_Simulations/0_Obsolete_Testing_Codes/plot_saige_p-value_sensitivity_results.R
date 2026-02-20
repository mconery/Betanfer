library(tidyverse)

sensitivity_file <- "C:/Users/mconery/Box/Mitchell_Work_Files/Meta-Analysis_Fine-Mapping/Data/SAIGE_Simulations/combined_causal_results.true_ld-based_selection.maf-balanced.variance-normalized.csv"
df <- read.table(sensitivity_file, sep=",", header = TRUE)

# Ensure columns are correctly parsed
df$num_causal_variants <- as.factor(df$num_causal_variants)
df$population <- as.factor(df$population)

# Create the faceted scatter plot
ggplot(df, aes(x = heritability, y = p.value, color = population)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~num_causal_variants, nrow = 1) +
  scale_y_log10() +
  scale_x_continuous(labels = scales::scientific_format(digits = 2)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Heritability",
    y = "P-value (log scale)",
    color = "Population",
    title = "Faceted Scatter Plot by Causal Variant Number"
  )