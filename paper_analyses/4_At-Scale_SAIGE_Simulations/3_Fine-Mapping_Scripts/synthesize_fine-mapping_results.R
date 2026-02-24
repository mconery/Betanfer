# Genomic Fine-mapping Performance Analysis Script - Summary Level Data with True AUPRC and Timing
# Author: Mitch Conery
# Date: September 26, 2025

# 0) Read required libraries ====

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(tidyr)
library(stringr)
library(viridis)
library(PRROC)  # For precision-recall curves

# 1) Read in files ====

# Set input directory
inp_dir = "C:/Users/mconery/Box/Mitchell_Work_Files/Meta-Analysis_Fine-Mapping/Data/SAIGE_Simulations/fine-mapping_summary"
# Set output directory
out_dir = "C:/Users/mconery/Box/Mitchell_Work_Files/Meta-Analysis_Fine-Mapping/Figures/SAIGE_Simulations"

# Read in files
method_df <- read.csv(paste(inp_dir, "summary_by_method_ld.txt", sep = "/"), sep = "\t")
causal_method_df <- read.csv(paste(inp_dir, "summary_by_causal_method_ld.txt", sep = "/"), sep = "\t")
trait_method_df <- read.csv(paste(inp_dir, "summary_by_trait_method_ld.txt", sep = "/"), sep = "\t")
trait_causal_method_df <- read.csv(paste(inp_dir, "summary_by_trait_causal_method_ld.txt", sep = "/"), sep = "\t")

# 2) Set basic functions ==== 

# Function to format method names with proper capitalization
format_method_names <- function(data) {
  data %>%
    mutate(method = case_when(
      method == "susie" ~ "SuSiE",
      method == "rsparsepro" ~ "RSparsePro",
      method == "carma" ~ "CARMA",
      method == "mesusie" ~ "MESuSiE",
      method == "multisusie" ~ "MultiSuSiE",
      method == "susiex" ~ "SuSiEx",
      TRUE ~ method
    ))
}

# Function to format LD matrix names with camel case and spaces
format_ld_matrix_names <- function(data) {
  data %>%
    mutate(ld_matrix = str_replace_all(ld_matrix, "_fixed", "")) %>%
    mutate(ld_matrix = str_replace_all(ld_matrix, "external", "(External)")) %>%
    mutate(ld_matrix = str_replace_all(ld_matrix, "_", " ")) %>%
    mutate(ld_matrix = str_to_title(ld_matrix)) %>%
    mutate(ld_matrix = str_replace_all(ld_matrix, "Mvp", "MVP")) %>%
    mutate(ld_matrix = str_replace_all(ld_matrix, "Pc", "PC"))
}

# Function to create coverage vs power scatter plot with filtering
plot_coverage_vs_power_filtered <- function(
    metrics_data,
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    methods = NULL,
    n_traits_filter = NULL,
    matrix_type_filter = NULL,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Filter for specific n_traits for beta-inferred matrices, but always keep raw_in-sample and PC-corrected
  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>%
      filter(
        # Keep raw_in-sample and PC-corrected
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        # Or keep beta-inferred with specified n_traits
        (is_beta_inferred == TRUE & n_traits == n_traits_filter)
      )
  }

  # Filter by matrix trait types (continuous, binary, all)
  # Always keep raw_in-sample and PC-corrected
  if (!is.null(matrix_type_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & trait_types == matrix_type_filter)
      )
  }

  # Create a matrix type label for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      # Factor for consistent ordering
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      ))
    )

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,                # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,        # filled square
    "Beta-Inferred (Weighted)" = 18,    # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8  # asterisk
  )

  # Create the plot
  p <- ggplot(plot_data, aes(x = cs_power, y = cs_coverage)) +
    geom_point(aes(color = method, shape = matrix_type), size = 10, alpha = 0.8) +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(nrow = 3, byrow = TRUE),
      color = guide_legend(nrow = 2)
    ) +
    labs(
      y = "Coverage",
      x = "Power"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 32),
      axis.text = element_text(size = 30),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

  return(p)
}

# Function to create a standalone legend for fine-mapping methods and LD matrix types
create_finemapping_legend <- function() {

  # Define all fine-mapping methods
  all_methods <- c("CARMA", "SuSiE", "RSparsePro", "MESuSiE", "MultiSuSiE", "SuSiEx")

  # Define all matrix types
  all_matrix_types <- c(
    "Raw In-Sample",
    "PC-Corrected In-Sample",
    "MVP (External)",
    "1000-Null Beta-Inferred (Meta)",
    "1000-Null Beta-Inferred (Weighted)",
    "1000-Null Beta-Inferred (Pop-Specific)",
    "Realistic Beta-Inferred (Weighted)",
    "Realistic Beta-Inferred (Pop-Specific)"
  )

  # Create dummy data for legend
  dummy_data <- expand.grid(
    method = factor(all_methods, levels = all_methods),
    matrix_type = factor(all_matrix_types, levels = all_matrix_types),
    x = 0.5,
    y = 0.5
  )

  # Define matrix shapes (consistent with other plots)
  matrix_shapes <- c(
    "Raw In-Sample" = 16,                          # filled circle
    "PC-Corrected In-Sample" = 17,                 # filled triangle
    "MVP (External)" = 4,                          # X shape
    "1000-Null Beta-Inferred (Meta)" = 15,                   # filled square
    "1000-Null Beta-Inferred (Weighted)" = 18,               # filled diamond
    "1000-Null Beta-Inferred (Pop-Specific)" = 8,            # asterisk
    "Realistic Beta-Inferred (Weighted)" = 0,      # open square
    "Realistic Beta-Inferred (Pop-Specific)" = 2   # open triangle
  )

  # Create a dummy plot to extract legend
  p <- ggplot(dummy_data, aes(x = x, y = y, color = method, shape = matrix_type)) +
    geom_point(size = 5, alpha = 0.8) +
    scale_color_viridis_d(name = "Fine-mapping Method", drop = FALSE) +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(ncol = 1, byrow = TRUE, order = 1, override.aes = list(size = 4)),
      color = guide_legend(ncol = 1, order = 2)
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 19, face = "bold"),
      legend.title.align = 0,
      legend.box.just = "left"
    )

  # Extract legend using ggplot functionality
  # Get the grob
  g <- ggplotGrob(p)
  legend_grob <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]

  return(legend_grob)
}

# Function to create faceted coverage vs power plot by trait type and causal variant count
plot_coverage_vs_power_faceted <- function(
    metrics_data,
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    methods = NULL,
    n_traits_filter = NULL,
    matrix_type_filter = NULL,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Filter for specific n_traits for beta-inferred matrices, but always keep raw_in-sample and PC-corrected
  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>%
      filter(
        # Keep raw_in-sample and PC-corrected
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        # Or keep beta-inferred with specified n_traits
        (is_beta_inferred == TRUE & n_traits == n_traits_filter)
      )
  }

  # Filter by matrix trait types (continuous, binary, all)
  # Always keep raw_in-sample and PC-corrected
  if (!is.null(matrix_type_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & trait_types == matrix_type_filter)
      )
  }

  # Create a matrix type label for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      # Factor for consistent ordering
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      )),
      # Format trait type labels
      trait_type_label = case_when(
        trait_type == "binary" ~ "Binary",
        trait_type == "continuous" ~ "Continuous",
        TRUE ~ trait_type
      ),
      # Format causal variant count labels
      causal_label = paste0(causal_variant_count, " Causal")
    )

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,                # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,        # filled square
    "Beta-Inferred (Weighted)" = 18,    # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8  # asterisk
  )

  # Create the faceted plot
  p <- ggplot(plot_data, aes(x = cs_power, y = cs_coverage)) +
    geom_point(aes(color = method, shape = matrix_type), size = 5, alpha = 0.8) +
    facet_grid(trait_type_label ~ causal_label, scales = "free") +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = TRUE
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE),
      color = guide_legend(nrow = 1)
    ) +
    labs(
      y = "Coverage",
      x = "Power"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 17, face = "bold"),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

  return(p)
}

# Function to create faceted coverage vs power plot by heritability (horizontal faceting)
plot_coverage_vs_power_by_heritability <- function(
    metrics_data,
    simulation_type = NULL,
    gwas_meta_type = NULL,
    methods = NULL,
    n_traits_filter = NULL,
    matrix_type_filter = NULL,
    causal_variant_count_filter = NULL,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data

  # Exclude drop_causal simulation results
  plot_data <- plot_data %>% filter(simulation_type != "drop_causal")

  # Apply filters
  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  if (!is.null(causal_variant_count_filter)) {
    plot_data <- plot_data %>% filter(causal_variant_count == causal_variant_count_filter)
  }

  # Filter for specific n_traits for beta-inferred matrices, but always keep raw_in-sample and PC-corrected
  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>%
      filter(
        # Keep raw_in-sample and PC-corrected
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        # Or keep beta-inferred with specified n_traits
        (is_beta_inferred == TRUE & n_traits == n_traits_filter)
      )
  }

  # Filter by matrix trait types (continuous, binary, all)
  # Always keep raw_in-sample and PC-corrected
  if (!is.null(matrix_type_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & trait_types == matrix_type_filter)
      )
  }

  # Create a matrix type label for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      # Factor for consistent ordering
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      )),
      # Format heritability labels
      heritability_label = paste0("h² = ", snp_heritability),
      # Format simulation type labels
      simulation_type_label = case_when(
        simulation_type == "homo" ~ "Homogenous Heritability",
        simulation_type == "matched_to" ~ "Homogenous Effect Sizes",
        TRUE ~ simulation_type
      ),
      # Factor for consistent ordering (homo above matched_to)
      simulation_type_label = factor(simulation_type_label, levels = c(
        "Homogenous Heritability",
        "Homogenous Effect Sizes"
      ))
    )

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,                # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,        # filled square
    "Beta-Inferred (Weighted)" = 18,    # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8  # asterisk
  )

  # Create the faceted plot (vertical faceting by simulation type, horizontal by heritability)
  p <- ggplot(plot_data, aes(x = cs_power, y = cs_coverage)) +
    geom_point(aes(color = method, shape = matrix_type), size = 5, alpha = 0.8) +
    facet_grid(simulation_type_label ~ heritability_label) +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE),
      color = guide_legend(nrow = 1)
    ) +
    labs(
      y = "Coverage",
      x = "Power"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 17, face = "bold"),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

  return(p)
}

# Function to create faceted coverage vs power plot by GWAS meta-analysis type (horizontal faceting)
# Only for meta-analysis mapping methods: CARMA, SuSiE, RSparsePro (by default)
plot_coverage_vs_power_by_meta_type <- function(
    metrics_data,
    simulation_type = NULL,
    snp_heritability = NULL,
    methods = c("CARMA", "SuSiE", "RSparsePro"),
    n_traits_filter = NULL,
    matrix_type_filter = NULL,
    causal_variant_count_filter = NULL,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data

  # Get all unique methods in the dataset to preserve color scheme
  all_methods <- unique(metrics_data$method)

  # Convert method to factor with all levels before filtering
  # This ensures colors stay consistent even after filtering
  plot_data <- plot_data %>%
    mutate(method = factor(method, levels = all_methods))

  # Exclude drop_causal simulation results
  plot_data <- plot_data %>% filter(simulation_type != "drop_causal")

  # Filter for specified fine-mapping methods
  if (!is.null(methods)) {
    plot_data <- plot_data %>%
      filter(method %in% methods)
  }

  # Apply filters
  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(causal_variant_count_filter)) {
    plot_data <- plot_data %>% filter(causal_variant_count == causal_variant_count_filter)
  }

  # Filter for specific n_traits for beta-inferred matrices, but always keep raw_in-sample and PC-corrected
  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>%
      filter(
        # Keep raw_in-sample and PC-corrected
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        # Or keep beta-inferred with specified n_traits
        (is_beta_inferred == TRUE & n_traits == n_traits_filter)
      )
  }

  # Filter by matrix trait types (continuous, binary, all)
  # Always keep raw_in-sample and PC-corrected
  if (!is.null(matrix_type_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & trait_types == matrix_type_filter)
      )
  }

  # Create a matrix type label for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      # Factor for consistent ordering
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      )),
      # Format meta-analysis type labels
      meta_type_label = case_when(
        meta_analysis_type == "fixed" ~ "Fixed-Effects",
        meta_analysis_type == "random" ~ "Random-Effects",
        TRUE ~ meta_analysis_type
      ),
      # Format simulation type labels
      simulation_type_label = case_when(
        simulation_type == "homo" ~ "Homogenous Heritability",
        simulation_type == "matched_to" ~ "Homogenous Effect Sizes",
        TRUE ~ simulation_type
      ),
      # Factor for consistent ordering (homo above matched_to)
      simulation_type_label = factor(simulation_type_label, levels = c(
        "Homogenous Heritability",
        "Homogenous Effect Sizes"
      ))
    )

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,                # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,        # filled square
    "Beta-Inferred (Weighted)" = 18,    # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8  # asterisk
  )

  # Get methods and matrix types present in filtered data for legend
  methods_present <- unique(plot_data$method)
  methods_present <- methods_present[!is.na(methods_present)]

  matrix_types_present <- unique(plot_data$matrix_type)
  matrix_types_present <- matrix_types_present[!is.na(matrix_types_present)]

  # Create the faceted plot (vertical faceting by simulation type, horizontal by meta-analysis type)
  p <- ggplot(plot_data, aes(x = cs_power, y = cs_coverage)) +
    geom_point(aes(color = method, shape = matrix_type), size = 5, alpha = 0.8) +
    facet_grid(simulation_type_label ~ meta_type_label) +
    scale_color_viridis_d(name = "Fine-mapping Method", drop = FALSE, breaks = methods_present) +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE,
      breaks = matrix_types_present
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE),
      color = guide_legend(nrow = 1)
    ) +
    labs(
      y = "Coverage",
      x = "Power"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 17, face = "bold"),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

  return(p)
}

# Function to plot performance metrics for realistic LD matrices vs p-value threshold
plot_performance_realistic_ld <- function(
    metrics_data,
    gwas_meta_type = NULL,
    snp_heritability = NULL,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data %>% filter(simulation_type == "homo")

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  # Filter for realistic LD matrices with 0kb padding, 1000-trait matrices, and raw_in-sample
  plot_data <- plot_data %>%
    filter(
      (str_starts(ld_matrix, "beta-inferred_realistic_") & str_detect(ld_matrix, "_0kb")) |
      (is_beta_inferred == TRUE & n_traits == 1000 & trait_types=="all") |
      ld_matrix == "raw_in-sample"
    )

  # Filter out "Meta" construction for meta-analysis methods (keep only weighted)
  plot_data <- plot_data %>%
    filter(matrix_construction != "meta" | is.na(matrix_construction))

  # Separate in-sample, realistic, and 1000-trait data
  insample_data <- plot_data %>% filter(ld_matrix == "raw_in-sample")
  realistic_data <- plot_data %>% filter(str_starts(ld_matrix, "beta-inferred_realistic_"))
  trait_1000_data <- plot_data %>% filter(is_beta_inferred == TRUE & n_traits == 1000)

  # Extract p-value threshold from realistic matrices
  realistic_data <- realistic_data %>%
    mutate(
      # Extract p-value from the matrix name (format: beta-inferred_realistic_<pval>_0kb[_weighted])
      p_value_str = str_extract(ld_matrix, "(?<=beta-inferred_realistic_)[^_]+"),
      p_value = as.numeric(p_value_str),
      # Create categorical x-axis label
      x_category = format(p_value, scientific = TRUE),
      # Determine matrix shape based on method
      matrix_shape = case_when(
        method %in% c("MESuSiE", "MultiSuSiE", "SuSiEx") ~ "Realistic Beta-Inferred (Pop-Specific)",
        method %in% c("CARMA", "RSparsePro", "SuSiE") ~ "Realistic Beta-Inferred (Weighted)",
        TRUE ~ "Other"
      )
    )

  # Prepare 1000-trait data
  trait_1000_data <- trait_1000_data %>%
    mutate(
      x_category = "1000-trait",
      matrix_shape = case_when(
        method %in% c("MESuSiE", "MultiSuSiE", "SuSiEx") ~ "1000-Null Beta-Inferred (Pop-Specific)",
        method %in% c("CARMA", "RSparsePro", "SuSiE") ~ "1000-Null Beta-Inferred (Weighted)",
        TRUE ~ "Other"
      )
    )

  # Define the desired order for x-axis categories
  x_order <- c("0e+00", "5e-08", "5e-06", "5e-04", "1000-trait", "In-Sample")

  # Prepare realistic data
  realistic_long <- realistic_data %>%
    select(method, matrix_shape, x_category, cs_coverage, cs_power, variant_auprc) %>%
    pivot_longer(
      cols = c(cs_coverage, cs_power, variant_auprc),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("cs_coverage", "cs_power", "variant_auprc")),
      metric_label = case_when(
        metric == "cs_coverage" ~ "Coverage",
        metric == "cs_power" ~ "Power",
        metric == "variant_auprc" ~ "AUPRC",
        TRUE ~ metric
      ),
      x_category = factor(x_category, levels = x_order)
    )

  # Prepare 1000-trait data
  trait_1000_long <- trait_1000_data %>%
    select(method, matrix_shape, x_category, cs_coverage, cs_power, variant_auprc) %>%
    pivot_longer(
      cols = c(cs_coverage, cs_power, variant_auprc),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("cs_coverage", "cs_power", "variant_auprc")),
      metric_label = case_when(
        metric == "cs_coverage" ~ "Coverage",
        metric == "cs_power" ~ "Power",
        metric == "variant_auprc" ~ "AUPRC",
        TRUE ~ metric
      ),
      x_category = factor(x_category, levels = x_order)
    )

  # Prepare in-sample data
  insample_long <- insample_data %>%
    mutate(
      x_category = "In-Sample",
      matrix_shape = "Raw In-Sample"
    ) %>%
    select(method, matrix_shape, x_category, cs_coverage, cs_power, variant_auprc) %>%
    pivot_longer(
      cols = c(cs_coverage, cs_power, variant_auprc),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("cs_coverage", "cs_power", "variant_auprc")),
      metric_label = case_when(
        metric == "cs_coverage" ~ "Coverage",
        metric == "cs_power" ~ "Power",
        metric == "variant_auprc" ~ "AUPRC",
        TRUE ~ metric
      ),
      x_category = factor(x_category, levels = x_order)
    )

  # Create x-axis labels (convert scientific notation to cleaner format)
  x_labels <- c("0", "5e-8", "5e-6", "5e-4", "1000-Null", "In-Sample")

  # Define matrix shapes (using new shapes distinct from previous: 16, 17, 15, 18, 8)
  matrix_shapes <- c(
    "Raw In-Sample" = 16,                          # filled circle (for reference)
    "1000-Null Beta-Inferred (Weighted)" = 18,               # filled diamond
    "1000-Null Beta-Inferred (Pop-Specific)" = 8,            # asterisk
    "Realistic Beta-Inferred (Weighted)" = 0,      # open square
    "Realistic Beta-Inferred (Pop-Specific)" = 2   # open triangle
  )

  # Create the plot with faceting by metric
  p <- ggplot() +
    # Realistic data with lines
    geom_line(data = realistic_long, aes(x = x_category, y = value, color = method, group = method),
              alpha = 0.6, linewidth = 3) +
    geom_point(data = realistic_long, aes(x = x_category, y = value, color = method, shape = matrix_shape),
               size = 8, alpha = 0.8) +
    # 1000-trait data as disconnected points
    geom_point(data = trait_1000_long, aes(x = x_category, y = value, color = method, shape = matrix_shape),
               size = 8, alpha = 0.8) +
    # In-sample data as disconnected points
    geom_point(data = insample_long, aes(x = x_category, y = value, color = method, shape = matrix_shape),
               size = 8.5, alpha = 0.8) +
    facet_grid(. ~ metric_label, scales = "free_y") +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(nrow = 1, byrow = TRUE),
      color = guide_legend(order = 2)
    ) +
    scale_x_discrete(
      labels = x_labels
    ) +
    labs(
      x = "P-value Threshold for Trait Filtering",
      y = "Performance Metric"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 28),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 28, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    )

  return(p)
}

# Function to plot mean CS size and failure rate for realistic LD matrices vs p-value threshold
plot_cs_size_failure_realistic_ld <- function(
    metrics_data,
    gwas_meta_type = NULL,
    snp_heritability = NULL,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data %>% filter(simulation_type == "homo")

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  # Filter for realistic LD matrices with 0kb padding, 1000-trait matrices, and raw_in-sample
  plot_data <- plot_data %>%
    filter(
      (str_starts(ld_matrix, "beta-inferred_realistic_") & str_detect(ld_matrix, "_0kb")) |
      (is_beta_inferred == TRUE & n_traits == 1000 & trait_types=="") |
      ld_matrix == "raw_in-sample"
    )

  # Filter out "Meta" construction for meta-analysis methods (keep only weighted)
  plot_data <- plot_data %>%
    filter(matrix_construction != "meta" | is.na(matrix_construction))

  # Calculate failure rate
  plot_data <- plot_data %>%
    mutate(failure_rate = (n_failed_loci / n_replicates) * 100)

  # Separate in-sample, realistic, and 1000-trait data
  insample_data <- plot_data %>% filter(ld_matrix == "raw_in-sample")
  realistic_data <- plot_data %>% filter(str_starts(ld_matrix, "beta-inferred_realistic_"))
  trait_1000_data <- plot_data %>% filter(is_beta_inferred == TRUE & n_traits == 1000)

  # Extract p-value threshold from realistic matrices
  realistic_data <- realistic_data %>%
    mutate(
      # Extract p-value from the matrix name (format: beta-inferred_realistic_<pval>_0kb[_weighted])
      p_value_str = str_extract(ld_matrix, "(?<=beta-inferred_realistic_)[^_]+"),
      p_value = as.numeric(p_value_str),
      # Create categorical x-axis label
      x_category = format(p_value, scientific = TRUE),
      # Determine matrix shape based on method
      matrix_shape = case_when(
        method %in% c("MESuSiE", "MultiSuSiE", "SuSiEx") ~ "Realistic Beta-Inferred (Pop-Specific)",
        method %in% c("CARMA", "RSparsePro", "SuSiE") ~ "Realistic Beta-Inferred (Weighted)",
        TRUE ~ "Other"
      )
    )

  # Prepare 1000-trait data
  trait_1000_data <- trait_1000_data %>%
    mutate(
      x_category = "1000-trait",
      matrix_shape = case_when(
        method %in% c("MESuSiE", "MultiSuSiE", "SuSiEx") ~ "1000-Null Beta-Inferred (Pop-Specific)",
        method %in% c("CARMA", "RSparsePro", "SuSiE") ~ "1000-Null Beta-Inferred (Weighted)",
        TRUE ~ "Other"
      )
    )

  # Define the desired order for x-axis categories
  x_order <- c("0e+00", "5e-08", "5e-06", "5e-04", "1000-trait", "In-Sample")

  # Prepare realistic data
  realistic_long <- realistic_data %>%
    select(method, matrix_shape, x_category, mean_cs_size, failure_rate) %>%
    pivot_longer(
      cols = c(mean_cs_size, failure_rate),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("mean_cs_size", "failure_rate")),
      metric_label = case_when(
        metric == "mean_cs_size" ~ "Mean CS Size",
        metric == "failure_rate" ~ "Failure Rate (%)",
        TRUE ~ metric
      ),
      x_category = factor(x_category, levels = x_order)
    )

  # Prepare 1000-trait data
  trait_1000_long <- trait_1000_data %>%
    select(method, matrix_shape, x_category, mean_cs_size, failure_rate) %>%
    pivot_longer(
      cols = c(mean_cs_size, failure_rate),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("mean_cs_size", "failure_rate")),
      metric_label = case_when(
        metric == "mean_cs_size" ~ "Mean CS Size",
        metric == "failure_rate" ~ "Failure Rate (%)",
        TRUE ~ metric
      ),
      x_category = factor(x_category, levels = x_order)
    )

  # Prepare in-sample data
  insample_long <- insample_data %>%
    mutate(
      x_category = "In-Sample",
      matrix_shape = "Raw In-Sample"
    ) %>%
    select(method, matrix_shape, x_category, mean_cs_size, failure_rate) %>%
    pivot_longer(
      cols = c(mean_cs_size, failure_rate),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("mean_cs_size", "failure_rate")),
      metric_label = case_when(
        metric == "mean_cs_size" ~ "Mean CS Size",
        metric == "failure_rate" ~ "Failure Rate (%)",
        TRUE ~ metric
      ),
      x_category = factor(x_category, levels = x_order)
    )

  # Create x-axis labels (convert scientific notation to cleaner format)
  x_labels <- c("0", "5e-8", "5e-6", "5e-4", "1000-Null", "In-Sample")

  # Define matrix shapes (using new shapes distinct from previous: 16, 17, 15, 18, 8)
  matrix_shapes <- c(
    "Raw In-Sample" = 16,                          # filled circle (for reference)
    "1000-Null Beta-Inferred (Weighted)" = 18,               # filled diamond
    "1000-Null Beta-Inferred (Pop-Specific)" = 8,            # asterisk
    "Realistic Beta-Inferred (Weighted)" = 0,      # open square
    "Realistic Beta-Inferred (Pop-Specific)" = 2   # open triangle
  )

  # Create the plot with faceting by metric
  p <- ggplot() +
    # Realistic data with lines
    geom_line(data = realistic_long, aes(x = x_category, y = value, color = method, group = method),
              alpha = 0.6, linewidth = 1.8) +
    geom_point(data = realistic_long, aes(x = x_category, y = value, color = method, shape = matrix_shape),
               size = 5, alpha = 0.8) +
    # 1000-trait data as disconnected points
    geom_point(data = trait_1000_long, aes(x = x_category, y = value, color = method, shape = matrix_shape),
               size = 5, alpha = 0.8) +
    # In-sample data as disconnected points
    geom_point(data = insample_long, aes(x = x_category, y = value, color = method, shape = matrix_shape),
               size = 5.5, alpha = 0.8) +
    facet_grid(. ~ metric_label, scales = "free_y") +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE, order = 1),
      color = guide_legend(nrow = 1, order = 2)
    ) +
    scale_x_discrete(
      labels = x_labels
    ) +
    labs(
      x = "P-value Threshold for Trait Filtering",
      y = "Metric Value"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 19, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 19, face = "bold"),
      legend.title.align = 0,
      legend.box.just = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    )

  return(p)
}

# Function to plot performance vs number of traits with raw_in-sample reference
plot_performance_vs_ntraits <- function(
    metrics_data,
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    methods = NULL,
    trait_types = NULL,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(trait_types)) {
    plot_data <- plot_data %>% filter(trait_types == !!trait_types)
  }
  
  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Filter for beta-inferred matrices, raw_in-sample, PC-corrected, and MVP (External)
  plot_data <- plot_data %>%
    filter(is_beta_inferred == TRUE | ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external"))

  # Separate in-sample and beta-inferred data
  insample_data <- plot_data %>% filter(ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external"))
  beta_data <- plot_data %>% filter(is_beta_inferred == TRUE)

  # Get max n_traits to position in-sample on the far right
  max_traits <- max(beta_data$n_traits, na.rm = TRUE)
  raw_position <- max_traits * 1.10  # Position raw_in-sample 10% beyond max
  pc_position <- max_traits * 1.20   # Position PC-corrected 20% beyond max
  mvp_position <- max_traits * 1.30  # Position MVP (External) 30% beyond max

  # Prepare beta-inferred data with matrix construction labels
  # Population-specific methods (SuSiEx, MESuSiE, MultiSuSiE) should appear in both Weighted and Meta facets
  beta_long <- beta_data %>%
    mutate(
      construction_facet = case_when(
        matrix_construction == "weighted" ~ "Weighted",
        matrix_construction == "meta" ~ "Meta",
        matrix_construction == "population-specific" ~ "population-specific",  # Will be duplicated
        TRUE ~ matrix_construction
      ),
      method_type = case_when(
        method %in% c("MESuSiE", "MultiSuSiE", "SuSiEx") ~ "Population-Specific",
        method %in% c("CARMA", "RSparsePro", "SuSiE") ~ "Weighted Average",
        TRUE ~ "Other"
      )
    ) %>%
    select(method, method_type, n_traits, construction_facet, cs_coverage, cs_power, variant_auprc)

  # Duplicate population-specific methods to both Weighted and Meta
  pop_specific <- beta_long %>% filter(construction_facet == "population-specific")

  # Create two copies: one for Weighted, one for Meta
  pop_specific_weighted <- pop_specific %>% mutate(construction_facet = "Weighted")
  pop_specific_meta <- pop_specific %>% mutate(construction_facet = "Meta")

  # Combine: replace population-specific with both copies
  beta_long <- beta_long %>%
    filter(construction_facet != "population-specific") %>%
    bind_rows(pop_specific_weighted, pop_specific_meta) %>%
    mutate(
      construction_facet = factor(construction_facet, levels = c("Weighted", "Meta")),
      construction_facet_label = case_when(
        construction_facet == "Weighted" ~ "Weighted Mean of Population-Specific Matrices",
        construction_facet == "Meta" ~ "Matrices Constructed from Meta-Analyses",
        TRUE ~ as.character(construction_facet)
      ),
      # Shape based on matrix construction type within each facet
      matrix_shape = case_when(
        method_type == "Population-Specific" ~ "Beta-Inferred (Pop-Specific)",
        construction_facet == "Weighted" ~ "Beta-Inferred (Weighted)",
        construction_facet == "Meta" ~ "Beta-Inferred (Meta)",
        TRUE ~ "Other"
      )
    ) %>%
    pivot_longer(
      cols = c(cs_coverage, cs_power, variant_auprc),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("cs_coverage", "cs_power", "variant_auprc")),
      metric_label = case_when(
        metric == "cs_coverage" ~ "Coverage",
        metric == "cs_power" ~ "Power",
        metric == "variant_auprc" ~ "AUPRC",
        TRUE ~ metric
      )
    )

  # Prepare in-sample data (raw_in-sample and PC-corrected) with adjusted positions
  # Separate positions for each type, duplicate to both Weighted and Meta facets
  construction_facets <- c("Weighted", "Meta")
  insample_long <- insample_data %>%
    mutate(
      n_traits = case_when(
        ld_matrix == "raw_in-sample" ~ raw_position,
        ld_matrix == "PC-corrected" ~ pc_position,
        ld_matrix == "mvp_external" ~ mvp_position,
        TRUE ~ NA_real_
      ),
      matrix_shape = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ "Other"
      )
    ) %>%
    select(method, n_traits, ld_matrix, matrix_shape, cs_coverage, cs_power, variant_auprc) %>%
    crossing(construction_facet = construction_facets) %>%
    mutate(
      construction_facet = factor(construction_facet, levels = c("Weighted", "Meta")),
      construction_facet_label = case_when(
        construction_facet == "Weighted" ~ "Weighted Mean of Population-Specific Matrices",
        construction_facet == "Meta" ~ "Matrices Constructed from Meta-Analyses",
        TRUE ~ as.character(construction_facet)
      )
    ) %>%
    pivot_longer(
      cols = c(cs_coverage, cs_power, variant_auprc),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("cs_coverage", "cs_power", "variant_auprc")),
      metric_label = case_when(
        metric == "cs_coverage" ~ "Coverage",
        metric == "cs_power" ~ "Power",
        metric == "variant_auprc" ~ "AUPRC",
        TRUE ~ metric
      )
    )

  # Create x-axis breaks at every 500 traits
  beta_breaks <- seq(0, max_traits, by = 500)
  beta_breaks <- beta_breaks[beta_breaks > 0]  # Remove 0 if present
  all_breaks <- c(beta_breaks, raw_position, pc_position, mvp_position)
  all_labels <- c(as.character(beta_breaks), "Raw\nIn-Sample", "PC-\nCorrected", "MVP\n(External)")

  # Define matrix shapes (consistent with other plots)
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,                # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,        # filled square
    "Beta-Inferred (Weighted)" = 18,    # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8  # asterisk
  )

  # Create the plot with faceting by metric and construction method
  p <- ggplot() +
    # Beta-inferred data with lines
    geom_line(data = beta_long, aes(x = n_traits, y = value, color = method, group = method),
              alpha = 0.6, linewidth = 1.8) +
    geom_point(data = beta_long, aes(x = n_traits, y = value, color = method, shape = matrix_shape),
               size = 5, alpha = 0.8) +
    # In-sample data as disconnected points
    geom_point(data = insample_long, aes(x = n_traits, y = value, color = method, shape = matrix_shape),
               size = 5.5, alpha = 0.8) +
    facet_grid(metric_label ~ construction_facet_label, scales = "free_y") +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(order = 1),
      color = guide_legend(order = 2)
    ) +
    scale_x_continuous(
      breaks = all_breaks,
      labels = all_labels
    ) +
    labs(
      x = "Number of Traits in LD Matrix",
      y = "Performance Metric"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 19, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 19, face = "bold"),
      legend.title.align = 0,
      legend.box.just = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    )

  return(p)
}

# Function to plot performance vs number of traits faceted by trait types
plot_performance_vs_ntraits_by_trait_type <- function(
    metrics_data,
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    methods = NULL,
    average_across_populations = FALSE,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Average across populations if requested
  if (average_across_populations) {
    plot_data <- plot_data %>%
      group_by(
        meta_analysis_type, simulation_type, snp_heritability, method,
        ld_matrix, n_traits, trait_types, matrix_construction, is_beta_inferred
      ) %>%
      summarise(
        cs_power = mean(cs_power, na.rm = TRUE),
        cs_coverage = mean(cs_coverage, na.rm = TRUE),
        variant_auprc = mean(variant_auprc, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Filter for beta-inferred matrices and raw_in-sample
  plot_data <- plot_data %>%
    filter(is_beta_inferred == TRUE | ld_matrix == "raw_in-sample")

  # Separate in-sample and beta-inferred data
  insample_data <- plot_data %>% filter(ld_matrix == "raw_in-sample")
  beta_data <- plot_data %>% filter(is_beta_inferred == TRUE)

  # Filter for weighted matrices only
  beta_data <- beta_data %>% filter(matrix_construction == "weighted" | matrix_construction == "population-specific")

  # Remove NA trait_types
  beta_data <- beta_data %>% filter(!is.na(trait_types))

  # Cap all data at 2000 traits
  beta_data <- beta_data %>% filter(n_traits <= 2000)

  # Get max n_traits to position in-sample on the far right
  max_traits <- 2000
  raw_position <- max_traits * 1.10  # Position raw_in-sample 10% beyond max

  # Prepare beta-inferred data with trait type labels
  beta_long <- beta_data %>%
    mutate(
      trait_type_facet = case_when(
        trait_types == "all" ~ "All Traits",
        trait_types == "binary" ~ "Binary Traits",
        trait_types == "continuous" ~ "Continuous Traits",
        TRUE ~ trait_types
      ),
      trait_type_facet = factor(trait_type_facet,
                                levels = c("All Traits", "Binary Traits", "Continuous Traits")),
      # Shape based on matrix construction type (all are weighted in this function)
      matrix_shape = case_when(
        method %in% c("MESuSiE", "MultiSuSiE", "SuSiEx") ~ "Beta-Inferred (Pop-Specific)",
        method %in% c("CARMA", "RSparsePro", "SuSiE") ~ "Beta-Inferred (Weighted)",
        TRUE ~ "Other"
      )
    ) %>%
    select(method, matrix_shape, n_traits, trait_type_facet, cs_coverage, cs_power, variant_auprc) %>%
    pivot_longer(
      cols = c(cs_coverage, cs_power, variant_auprc),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("cs_coverage", "cs_power", "variant_auprc")),
      metric_label = case_when(
        metric == "cs_coverage" ~ "Coverage",
        metric == "cs_power" ~ "Power",
        metric == "variant_auprc" ~ "AUPRC",
        TRUE ~ metric
      )
    )

  # Prepare in-sample data (raw_in-sample) with adjusted position
  # Duplicate to all three trait type facets
  trait_type_facets <- c("All Traits", "Binary Traits", "Continuous Traits")
  insample_long <- insample_data %>%
    mutate(
      n_traits = raw_position,
      matrix_shape = "Raw In-Sample"
    ) %>%
    select(method, matrix_shape, n_traits, ld_matrix, cs_coverage, cs_power, variant_auprc) %>%
    crossing(trait_type_facet = trait_type_facets) %>%
    mutate(
      trait_type_facet = factor(trait_type_facet,
                                levels = c("All Traits", "Binary Traits", "Continuous Traits"))
    ) %>%
    pivot_longer(
      cols = c(cs_coverage, cs_power, variant_auprc),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("cs_coverage", "cs_power", "variant_auprc")),
      metric_label = case_when(
        metric == "cs_coverage" ~ "Coverage",
        metric == "cs_power" ~ "Power",
        metric == "variant_auprc" ~ "AUPRC",
        TRUE ~ metric
      )
    )

  # Create x-axis breaks at every 500 traits
  beta_breaks <- seq(0, max_traits, by = 500)
  beta_breaks <- beta_breaks[beta_breaks > 0]  # Remove 0 if present
  all_breaks <- c(beta_breaks, raw_position)
  all_labels <- c(as.character(beta_breaks), "In-Sample")

  # Define matrix shapes (consistent with other plots)
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,                # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,        # filled square
    "Beta-Inferred (Weighted)" = 18,    # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8  # asterisk
  )

  # Create the plot with faceting by metric and trait type
  p <- ggplot() +
    # Beta-inferred data with lines
    geom_line(data = beta_long, aes(x = n_traits, y = value, color = method, group = method),
              alpha = 0.6, linewidth = 1.8) +
    geom_point(data = beta_long, aes(x = n_traits, y = value, color = method, shape = matrix_shape),
               size = 5, alpha = 0.8) +
    # In-sample data as disconnected points
    geom_point(data = insample_long, aes(x = n_traits, y = value, color = method, shape = matrix_shape),
               size = 5.5, alpha = 0.8) +
    facet_grid(metric_label ~ trait_type_facet, scales = "free_y") +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    scale_x_continuous(
      breaks = all_breaks,
      labels = all_labels
    ) +
    labs(
      x = "Number of Traits in LD Matrix",
      y = "Performance Metric"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 19, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 19, face = "bold"),
      legend.title.align = 0,
      legend.box.just = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    )

  return(p)
}

# Function to plot failure rate vs number of traits faceted by trait types
plot_failure_rate_vs_ntraits_by_trait_type <- function(
    metrics_data,
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    methods = NULL,
    average_across_populations = FALSE,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Average across populations if requested
  if (average_across_populations) {
    plot_data <- plot_data %>%
      group_by(
        meta_analysis_type, simulation_type, snp_heritability, method,
        ld_matrix, n_traits, trait_types, matrix_construction, is_beta_inferred
      ) %>%
      summarise(
        n_failed_loci = sum(n_failed_loci, na.rm = TRUE),
        n_replicates = sum(n_replicates, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Calculate failure rate as a percentage
  plot_data <- plot_data %>%
    mutate(failure_rate = (n_failed_loci / n_replicates) * 100)

  # Filter for beta-inferred matrices and raw_in-sample
  plot_data <- plot_data %>%
    filter(is_beta_inferred == TRUE | ld_matrix == "raw_in-sample")

  # Separate in-sample and beta-inferred data
  insample_data <- plot_data %>% filter(ld_matrix == "raw_in-sample")
  beta_data <- plot_data %>% filter(is_beta_inferred == TRUE)

  # Filter for weighted matrices only
  beta_data <- beta_data %>% filter(matrix_construction == "weighted" | matrix_construction == "population-specific")

  # Remove NA trait_types
  beta_data <- beta_data %>% filter(!is.na(trait_types))

  # Cap all data at 2000 traits
  beta_data <- beta_data %>% filter(n_traits <= 2000)

  # Get max n_traits to position in-sample on the far right
  max_traits <- 2000
  raw_position <- max_traits * 1.10  # Position raw_in-sample 10% beyond max

  # Prepare beta-inferred data with trait type labels
  beta_long <- beta_data %>%
    mutate(
      trait_type_facet = case_when(
        trait_types == "all" ~ "All Traits",
        trait_types == "binary" ~ "Binary Traits",
        trait_types == "continuous" ~ "Continuous Traits",
        TRUE ~ trait_types
      ),
      trait_type_facet = factor(trait_type_facet,
                                levels = c("All Traits", "Binary Traits", "Continuous Traits")),
      # Shape based on matrix construction type
      matrix_shape = case_when(
        method %in% c("MESuSiE", "MultiSuSiE", "SuSiEx") ~ "Beta-Inferred (Pop-Specific)",
        method %in% c("CARMA", "RSparsePro", "SuSiE") ~ "Beta-Inferred (Weighted)",
        TRUE ~ "Other"
      )
    )

  # Prepare in-sample data (raw_in-sample) with adjusted position
  # Duplicate to all three trait type facets
  trait_type_facets <- c("All Traits", "Binary Traits", "Continuous Traits")
  insample_long <- insample_data %>%
    mutate(
      n_traits = raw_position,
      matrix_shape = "Raw In-Sample"
    ) %>%
    crossing(trait_type_facet = trait_type_facets) %>%
    mutate(
      trait_type_facet = factor(trait_type_facet,
                                levels = c("All Traits", "Binary Traits", "Continuous Traits"))
    )

  # Create x-axis breaks at every 500 traits
  beta_breaks <- seq(0, max_traits, by = 500)
  beta_breaks <- beta_breaks[beta_breaks > 0]  # Remove 0 if present
  all_breaks <- c(beta_breaks, raw_position)
  all_labels <- c(as.character(beta_breaks), "In-Sample")

  # Define matrix shapes (consistent with other plots)
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,     # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,       # filled square
    "Beta-Inferred (Weighted)" = 18,   # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8 # asterisk
  )

  # Create the plot with faceting by trait type
  p <- ggplot() +
    # Beta-inferred data with lines
    geom_line(data = beta_long, aes(x = n_traits, y = failure_rate, color = method, group = method),
              alpha = 0.6, linewidth = 1.8) +
    geom_point(data = beta_long, aes(x = n_traits, y = failure_rate, color = method, shape = matrix_shape),
               size = 5, alpha = 0.8) +
    # In-sample data as disconnected points
    geom_point(data = insample_long, aes(x = n_traits, y = failure_rate, color = method, shape = matrix_shape),
               size = 5.5, alpha = 0.8) +
    facet_grid(. ~ trait_type_facet) +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    scale_x_continuous(
      breaks = all_breaks,
      labels = all_labels
    ) +
    labs(
      x = "Number of Traits in LD Matrix",
      y = "Failure Rate (%)"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 19, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 19, face = "bold"),
      legend.title.align = 0,
      legend.box.just = "left",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    )

  return(p)
}

# Function to create bar plots for various metrics with extensive filtering
plot_metric_barplot <- function(
    metrics_data,
    metric = "variant_auprc",  # Options: "variant_auprc", "mean_execution_time_seconds", "mean_cs_size"
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    n_traits_filter = NULL,
    matrix_type_filter = NULL,  # Filter for trait types: "all", "continuous", "binary"
    ld_construction_filter = NULL,  # Filter: "raw_in-sample", "PC-corrected", "beta-inferred"
    methods = NULL,
    average_across_populations = FALSE,
    title_suffix = ""
) {

  # Validate metric parameter
  valid_metrics <- c("variant_auprc", "mean_execution_time_seconds", "mean_cs_size")
  if (!(metric %in% valid_metrics)) {
    stop(paste("metric must be one of:", paste(valid_metrics, collapse = ", ")))
  }

  # Start with the full dataset
  plot_data <- metrics_data

  # Get all unique methods in the dataset to preserve color scheme
  all_methods <- unique(metrics_data$method)

  # Convert method to factor with all levels before filtering
  # This ensures colors stay consistent even after filtering
  plot_data <- plot_data %>%
    mutate(method = factor(method, levels = all_methods))

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Filter for specific n_traits for beta-inferred matrices, but always keep raw_in-sample and PC-corrected
  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>%
      filter(
        # Keep raw_in-sample and PC-corrected
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        # Or keep beta-inferred with specified n_traits
        (is_beta_inferred == TRUE & n_traits == n_traits_filter)
      )
  }

  # Filter by matrix trait types (continuous, binary, all)
  # Always keep raw_in-sample and PC-corrected
  if (!is.null(matrix_type_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & trait_types == matrix_type_filter)
      )
  }

  # Filter by LD construction type
  if (!is.null(ld_construction_filter)) {
    if (ld_construction_filter == "raw_in-sample") {
      plot_data <- plot_data %>% filter(ld_matrix == "raw_in-sample")
    } else if (ld_construction_filter == "PC-corrected") {
      plot_data <- plot_data %>% filter(ld_matrix == "PC-corrected")
    } else if (ld_construction_filter == "beta-inferred") {
      plot_data <- plot_data %>% filter(is_beta_inferred == TRUE)
    }
  }

  # Average across populations if requested
  if (average_across_populations) {
    # Determine grouping columns based on what's in the data
    grouping_cols <- c("meta_analysis_type", "simulation_type", "snp_heritability",
                       "method", "ld_matrix", "matrix_construction", "is_beta_inferred")

    # Add optional grouping columns if they exist and aren't filtered to a single value
    if ("n_traits" %in% names(plot_data)) {
      grouping_cols <- c(grouping_cols, "n_traits")
    }
    if ("trait_types" %in% names(plot_data)) {
      grouping_cols <- c(grouping_cols, "trait_types")
    }

    plot_data <- plot_data %>%
      group_by(across(all_of(grouping_cols))) %>%
      summarise(
        variant_auprc = mean(variant_auprc, na.rm = TRUE),
        mean_execution_time_seconds = mean(mean_execution_time_seconds, na.rm = TRUE),
        mean_cs_size = mean(mean_cs_size, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Create a matrix type label for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      # Factor for consistent ordering
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      ))
    )

  # Get methods present in filtered data for legend
  methods_present <- unique(plot_data$method)
  methods_present <- methods_present[!is.na(methods_present)]

  # Set up metric-specific parameters
  if (metric == "variant_auprc") {
    y_column <- "variant_auprc"
    y_label <- "Variant-level AUPRC"
  } else if (metric == "mean_execution_time_seconds") {
    y_column <- "mean_execution_time_seconds"  # Actual column name in data
    y_label <- "Mean Execution Time (seconds)"
  } else if (metric == "mean_cs_size") {
    y_column <- "mean_cs_size"
    y_label <- "Mean Credible Set Size"
  }

  # Determine x-axis label format based on number of unique matrix types
  # If only one matrix type, show just method names; otherwise show method + matrix type
  unique_matrix_types <- unique(plot_data$matrix_type)
  unique_matrix_types <- unique_matrix_types[!is.na(unique_matrix_types)]

  if (length(unique_matrix_types) == 1) {
    # Only one matrix type - just show method names
    plot_data <- plot_data %>%
      mutate(x_label = as.character(method))
    x_axis_label <- "Fine-mapping Method"
  } else {
    # Multiple matrix types - show method + matrix type
    plot_data <- plot_data %>%
      mutate(x_label = paste(method, matrix_type, sep = "\n"))
    x_axis_label <- "Method + LD Matrix Type"
  }

  # Create the bar plot
  p <- ggplot(plot_data, aes(x = x_label, y = .data[[y_column]], fill = method)) +
    geom_col(alpha = 0.8, position = "dodge") +
    scale_fill_viridis_d(name = "Fine-mapping Method", drop = FALSE, breaks = methods_present) +
    labs(
      x = x_axis_label,
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 22, color = "black"),
      axis.text = element_text(size = 20, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.text = element_text(size = 19, color = "black"),
      legend.title = element_text(size = 20, face = "bold", color = "black"),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.y = element_line(color = "black"),
      axis.ticks.length.y = unit(0.15, "cm")
    )

  # Add pseudo-log scale for execution time if needed
  # Pseudo-log scale is linear near 0 and logarithmic for larger values
  if (metric == "mean_execution_time_seconds") {
    # Labels at: 0, 0.5, 1, 5, 10, 50 (extended for larger values)
    custom_breaks <- c(
      seq(0.1, 0.9, by = 0.1),   # 0.1, 0.2, ..., 0.9
      1:10,                       # 1, 2, 3, ..., 10
      seq(20, 100, by = 10),      # 20, 30, ..., 100
      seq(200, 1000, by = 100),   # 200, 300, ..., 1000 (extended)
      seq(2000, 10000, by = 1000) # 2000, 3000, ..., 10000 (extended)
    )

    # Tick marks at: 0.1-0.9, 1-10, and every 10 through 100
    custom_labels <- ifelse(custom_breaks %in% c(0, 0.1, 1, 10, 100, 1000), custom_breaks, "")

    p <- p + scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = custom_breaks,
      labels = custom_labels, limits = c(0,100)
    )
  }

  return(p)
}

# Function to create scatter plots for metrics by method and matrix type
plot_metric_scatterplot <- function(
    metrics_data,
    metric = "variant_auprc",  # Options: "variant_auprc", "mean_cs_size"
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    n_traits_filter = NULL,
    matrix_type_filter = NULL,
    ld_construction_filter = NULL,
    methods = NULL,
    average_across_populations = FALSE,
    title_suffix = ""
) {

  # Validate metric parameter
  valid_metrics <- c("variant_auprc", "mean_cs_size")
  if (!(metric %in% valid_metrics)) {
    stop(paste("metric must be one of:", paste(valid_metrics, collapse = ", ")))
  }

  # Start with the full dataset
  plot_data <- metrics_data

  # Get all unique methods in the dataset to preserve color scheme
  all_methods <- unique(metrics_data$method)

  # Convert method to factor with all levels before filtering
  # This ensures colors stay consistent even after filtering
  plot_data <- plot_data %>%
    mutate(method = factor(method, levels = all_methods))

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Filter for specific n_traits for beta-inferred matrices, but always keep raw_in-sample and PC-corrected
  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>%
      filter(
        # Keep raw_in-sample and PC-corrected
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        # Or keep beta-inferred with specified n_traits
        (is_beta_inferred == TRUE & n_traits == n_traits_filter)
      )
  }

  # Filter by matrix trait types (continuous, binary, all)
  # Always keep raw_in-sample and PC-corrected
  if (!is.null(matrix_type_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & trait_types == matrix_type_filter)
      )
  }

  # Filter by LD construction type
  if (!is.null(ld_construction_filter)) {
    if (ld_construction_filter == "raw_in-sample") {
      plot_data <- plot_data %>% filter(ld_matrix == "raw_in-sample")
    } else if (ld_construction_filter == "PC-corrected") {
      plot_data <- plot_data %>% filter(ld_matrix == "PC-corrected")
    } else if (ld_construction_filter == "beta-inferred") {
      plot_data <- plot_data %>% filter(is_beta_inferred == TRUE)
    }
  }

  # Average across populations if requested
  if (average_across_populations) {
    # Determine grouping columns based on what's in the data
    grouping_cols <- c("meta_analysis_type", "simulation_type", "snp_heritability",
                       "method", "ld_matrix", "matrix_construction", "is_beta_inferred")

    # Add optional grouping columns if they exist and aren't filtered to a single value
    if ("n_traits" %in% names(plot_data)) {
      grouping_cols <- c(grouping_cols, "n_traits")
    }
    if ("trait_types" %in% names(plot_data)) {
      grouping_cols <- c(grouping_cols, "trait_types")
    }

    plot_data <- plot_data %>%
      group_by(across(all_of(grouping_cols))) %>%
      summarise(
        variant_auprc = mean(variant_auprc, na.rm = TRUE),
        mean_cs_size = mean(mean_cs_size, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Create a matrix type label for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      # Factor for consistent ordering
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      ))
    )

  # Get methods and matrix types present in filtered data for legend
  methods_present <- unique(plot_data$method)
  methods_present <- methods_present[!is.na(methods_present)]

  matrix_types_present <- unique(plot_data$matrix_type)
  matrix_types_present <- matrix_types_present[!is.na(matrix_types_present)]

  # Set up metric-specific parameters
  if (metric == "variant_auprc") {
    y_column <- "variant_auprc"
    y_label <- "Variant-level AUPRC"
  } else if (metric == "mean_cs_size") {
    y_column <- "mean_cs_size"
    y_label <- "Mean Credible Set Size"
  }

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,     # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,        # filled square
    "Beta-Inferred (Weighted)" = 18,    # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8  # asterisk
  )

  # Create the scatter plot
  p <- ggplot(plot_data, aes(x = method, y = .data[[y_column]], color = method, shape = matrix_type)) +
    geom_point(size = 12, alpha = 0.8, position = position_dodge(width = 0.5)) +
    scale_color_viridis_d(name = "Fine-mapping Method", drop = FALSE, breaks = methods_present) +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE,
      breaks = matrix_types_present
    ) +
    guides(
      shape = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(size = 5))
    ) +
    labs(
      x = "",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 32, color = "black"),
      axis.text = element_text(size = 30, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      plot.margin = margin(t = 20, r = 5, b = 5, l = 5)
    )

  return(p)
}

# Function to create faceted scatter plots by heritability (horizontal faceting)
plot_metric_scatterplot_by_heritability <- function(
    metrics_data,
    metric = "variant_auprc",  # Options: "variant_auprc", "mean_cs_size"
    simulation_type = NULL,
    gwas_meta_type = NULL,
    methods = NULL,
    n_traits_filter = NULL,
    matrix_type_filter = NULL,
    causal_variant_count_filter = NULL,
    title_suffix = ""
) {

  # Validate metric parameter
  valid_metrics <- c("variant_auprc", "mean_cs_size")
  if (!(metric %in% valid_metrics)) {
    stop(paste("metric must be one of:", paste(valid_metrics, collapse = ", ")))
  }

  # Start with the full dataset
  plot_data <- metrics_data

  # Exclude drop_causal simulation results
  plot_data <- plot_data %>% filter(simulation_type != "drop_causal")

  # Get all unique methods in the dataset to preserve color scheme
  all_methods <- unique(metrics_data$method)

  # Convert method to factor with all levels before filtering
  plot_data <- plot_data %>%
    mutate(method = factor(method, levels = all_methods))

  # Apply filters
  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  if (!is.null(causal_variant_count_filter)) {
    plot_data <- plot_data %>% filter(causal_variant_count == causal_variant_count_filter)
  }

  # Filter for specific n_traits for beta-inferred matrices, but always keep raw_in-sample and PC-corrected
  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & n_traits == n_traits_filter)
      )
  }

  # Filter by matrix trait types
  if (!is.null(matrix_type_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & trait_types == matrix_type_filter)
      )
  }

  # Create labels for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      )),
      heritability_label = paste0("h² = ", snp_heritability),
      simulation_type_label = case_when(
        simulation_type == "homo" ~ "Homogenous Heritability",
        simulation_type == "matched_to" ~ "Homogenous Effect Sizes",
        TRUE ~ simulation_type
      ),
      simulation_type_label = factor(simulation_type_label, levels = c(
        "Homogenous Heritability",
        "Homogenous Effect Sizes"
      ))
    )

  # Get methods and matrix types present in filtered data
  methods_present <- unique(plot_data$method)
  methods_present <- methods_present[!is.na(methods_present)]

  matrix_types_present <- unique(plot_data$matrix_type)
  matrix_types_present <- matrix_types_present[!is.na(matrix_types_present)]

  # Set up metric-specific parameters
  if (metric == "variant_auprc") {
    y_column <- "variant_auprc"
    y_label <- "Variant-level AUPRC"
  } else if (metric == "mean_cs_size") {
    y_column <- "mean_cs_size"
    y_label <- "Mean Credible Set Size"
  }

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,
    "PC-Corrected In-Sample" = 17,
    "MVP (External)" = 4,
    "Beta-Inferred (Meta)" = 15,
    "Beta-Inferred (Weighted)" = 18,
    "Beta-Inferred (Pop-Specific)" = 8
  )

  # Create the faceted plot
  p <- ggplot(plot_data, aes(x = method, y = .data[[y_column]], color = method, shape = matrix_type)) +
    geom_point(size = 5, alpha = 0.8, position = position_dodge(width = 0.5)) +
    facet_grid(simulation_type_label ~ heritability_label) +
    scale_color_viridis_d(name = "Fine-mapping Method", drop = FALSE, breaks = methods_present) +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE,
      breaks = matrix_types_present
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE),
      color = guide_legend(nrow = 1)
    ) +
    labs(
      x = "Fine-mapping Method",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 17, face = "bold", color = "black"),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    )

  return(p)
}

# Function to create faceted scatter plots by GWAS meta-analysis type (horizontal faceting)
plot_metric_scatterplot_by_meta_type <- function(
    metrics_data,
    metric = "variant_auprc",  # Options: "variant_auprc", "mean_cs_size"
    simulation_type = NULL,
    snp_heritability = NULL,
    methods = c("CARMA", "SuSiE", "RSparsePro"),
    n_traits_filter = NULL,
    matrix_type_filter = NULL,
    causal_variant_count_filter = NULL,
    title_suffix = ""
) {

  # Validate metric parameter
  valid_metrics <- c("variant_auprc", "mean_cs_size")
  if (!(metric %in% valid_metrics)) {
    stop(paste("metric must be one of:", paste(valid_metrics, collapse = ", ")))
  }

  # Start with the full dataset
  plot_data <- metrics_data

  # Get all unique methods in the dataset to preserve color scheme
  all_methods <- unique(metrics_data$method)

  # Convert method to factor with all levels before filtering
  plot_data <- plot_data %>%
    mutate(method = factor(method, levels = all_methods))

  # Exclude drop_causal simulation results
  plot_data <- plot_data %>% filter(simulation_type != "drop_causal")

  # Filter for specified fine-mapping methods
  if (!is.null(methods)) {
    plot_data <- plot_data %>%
      filter(method %in% methods)
  }

  # Apply filters
  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(causal_variant_count_filter)) {
    plot_data <- plot_data %>% filter(causal_variant_count == causal_variant_count_filter)
  }

  # Filter for specific n_traits for beta-inferred matrices
  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & n_traits == n_traits_filter)
      )
  }

  # Filter by matrix trait types
  if (!is.null(matrix_type_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & trait_types == matrix_type_filter)
      )
  }

  # Create labels for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      )),
      meta_type_label = case_when(
        meta_analysis_type == "fixed" ~ "Fixed-Effects",
        meta_analysis_type == "random" ~ "Random-Effects",
        TRUE ~ meta_analysis_type
      ),
      simulation_type_label = case_when(
        simulation_type == "homo" ~ "Homogenous Heritability",
        simulation_type == "matched_to" ~ "Homogenous Effect Sizes",
        TRUE ~ simulation_type
      ),
      simulation_type_label = factor(simulation_type_label, levels = c(
        "Homogenous Heritability",
        "Homogenous Effect Sizes"
      ))
    )

  # Get methods and matrix types present in filtered data
  methods_present <- unique(plot_data$method)
  methods_present <- methods_present[!is.na(methods_present)]

  matrix_types_present <- unique(plot_data$matrix_type)
  matrix_types_present <- matrix_types_present[!is.na(matrix_types_present)]

  # Set up metric-specific parameters
  if (metric == "variant_auprc") {
    y_column <- "variant_auprc"
    y_label <- "Variant-level AUPRC"
  } else if (metric == "mean_cs_size") {
    y_column <- "mean_cs_size"
    y_label <- "Mean Credible Set Size"
  }

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,
    "PC-Corrected In-Sample" = 17,
    "MVP (External)" = 4,
    "Beta-Inferred (Meta)" = 15,
    "Beta-Inferred (Weighted)" = 18,
    "Beta-Inferred (Pop-Specific)" = 8
  )

  # Create the faceted plot
  p <- ggplot(plot_data, aes(x = method, y = .data[[y_column]], color = method, shape = matrix_type)) +
    geom_point(size = 5, alpha = 0.8, position = position_dodge(width = 0.5)) +
    facet_grid(simulation_type_label ~ meta_type_label) +
    scale_color_viridis_d(name = "Fine-mapping Method", drop = FALSE, breaks = methods_present) +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE,
      breaks = matrix_types_present
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE),
      color = guide_legend(nrow = 1)
    ) +
    labs(
      x = "Fine-mapping Method",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 17, face = "bold", color = "black"),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    )

  return(p)
}

# Function to create faceted metric scatterplot by trait type and causal variant count
plot_metric_scatterplot_faceted <- function(
    metrics_data,
    metric = "variant_auprc",  # Options: "variant_auprc", "mean_cs_size"
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    methods = NULL,
    n_traits_filter = NULL,
    matrix_type_filter = NULL,
    title_suffix = ""
) {

  # Validate metric parameter
  valid_metrics <- c("variant_auprc", "mean_cs_size")
  if (!(metric %in% valid_metrics)) {
    stop(paste("metric must be one of:", paste(valid_metrics, collapse = ", ")))
  }

  # Start with the full dataset
  plot_data <- metrics_data

  # Get all unique methods in the dataset to preserve color scheme
  all_methods <- unique(metrics_data$method)

  # Convert method to factor with all levels before filtering
  plot_data <- plot_data %>%
    mutate(method = factor(method, levels = all_methods))

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Filter for specific n_traits for beta-inferred matrices, but always keep raw_in-sample and PC-corrected
  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>%
      filter(
        # Keep raw_in-sample and PC-corrected
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        # Or keep beta-inferred with specified n_traits
        (is_beta_inferred == TRUE & n_traits == n_traits_filter)
      )
  }

  # Filter by matrix trait types (continuous, binary, all)
  # Always keep raw_in-sample and PC-corrected
  if (!is.null(matrix_type_filter)) {
    plot_data <- plot_data %>%
      filter(
        ld_matrix %in% c("raw_in-sample", "PC-corrected", "mvp_external") |
        (is_beta_inferred == TRUE & trait_types == matrix_type_filter)
      )
  }

  # Create labels for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      # Factor for consistent ordering
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      )),
      # Format trait type labels
      trait_type_label = case_when(
        trait_type == "binary" ~ "Binary",
        trait_type == "continuous" ~ "Continuous",
        TRUE ~ trait_type
      ),
      # Format causal variant count labels
      causal_label = paste0(causal_variant_count, " Causal")
    )

  # Get methods and matrix types present in filtered data
  methods_present <- unique(plot_data$method)
  methods_present <- methods_present[!is.na(methods_present)]

  matrix_types_present <- unique(plot_data$matrix_type)
  matrix_types_present <- matrix_types_present[!is.na(matrix_types_present)]

  # Set up metric-specific parameters
  if (metric == "variant_auprc") {
    y_column <- "variant_auprc"
    y_label <- "Variant-level AUPRC"
  } else if (metric == "mean_cs_size") {
    y_column <- "mean_cs_size"
    y_label <- "Mean Credible Set Size"
  }

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,     # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,       # filled square
    "Beta-Inferred (Weighted)" = 18,   # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8 # asterisk
  )

  # Create the faceted plot
  p <- ggplot(plot_data, aes(x = method, y = .data[[y_column]], color = method, shape = matrix_type)) +
    geom_point(size = 5, alpha = 0.8, position = position_dodge(width = 0.5)) +
    facet_grid(trait_type_label ~ causal_label, scales = "free") +
    scale_color_viridis_d(name = "Fine-mapping Method", drop = FALSE, breaks = methods_present) +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE,
      breaks = matrix_types_present
    ) +
    guides(
      shape = guide_legend(nrow = 2, byrow = TRUE),
      color = guide_legend(nrow = 1)
    ) +
    labs(
      x = "Fine-mapping Method",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 17, face = "bold", color = "black"),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    )

  return(p)
}

# Function to create faceted coverage vs power plot by GWAS trait type and LD matrix trait types
plot_coverage_vs_power_by_trait_types <- function(
    metrics_data,
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    n_traits_filter = NULL,
    methods = NULL,
    title_suffix = ""
) {

  # Start with the full dataset
  plot_data <- metrics_data

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>% filter(n_traits == n_traits_filter)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Create labels for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      # Factor for consistent ordering
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      )),
      # Format GWAS trait type labels
      gwas_trait_label = case_when(
        trait_type == "binary" ~ "Binary GWAS",
        trait_type == "continuous" ~ "Continuous GWAS",
        TRUE ~ trait_type
      ),
      # Format LD matrix trait types labels
      ld_trait_label = case_when(
        trait_types == "binary" ~ "Binary LD",
        trait_types == "continuous" ~ "Continuous LD",
        trait_types == "all" ~ "All-Trait LD",
        TRUE ~ trait_types
      )
    )

  # Get matrix types present in filtered data
  matrix_types_present <- unique(plot_data$matrix_type)
  matrix_types_present <- matrix_types_present[!is.na(matrix_types_present)]

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,     # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,       # filled square
    "Beta-Inferred (Weighted)" = 18,   # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8 # asterisk
  )

  # Create the faceted plot
  p <- ggplot(plot_data, aes(x = cs_power, y = cs_coverage)) +
    geom_point(aes(color = method, shape = matrix_type), size = 5, alpha = 0.8) +
    facet_grid(ld_trait_label ~ gwas_trait_label, scales = "free") +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE,
      breaks = matrix_types_present
    ) +
    guides(
      shape = guide_legend(nrow = 1, byrow = TRUE),
      color = guide_legend(nrow = 1)
    ) +
    labs(
      y = "Coverage",
      x = "Power"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 17, face = "bold"),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

  return(p)
}

# Function to create faceted metric scatterplot by GWAS trait type and LD matrix trait types
plot_metric_scatterplot_by_trait_types <- function(
    metrics_data,
    metric = "variant_auprc",  # Options: "variant_auprc", "mean_cs_size"
    gwas_meta_type = NULL,
    simulation_type = NULL,
    snp_heritability = NULL,
    n_traits_filter = NULL,
    methods = NULL,
    title_suffix = ""
) {

  # Validate metric parameter
  valid_metrics <- c("variant_auprc", "mean_cs_size")
  if (!(metric %in% valid_metrics)) {
    stop(paste("metric must be one of:", paste(valid_metrics, collapse = ", ")))
  }

  # Start with the full dataset
  plot_data <- metrics_data

  # Get all unique methods in the dataset to preserve color scheme
  all_methods <- unique(metrics_data$method)

  # Convert method to factor with all levels before filtering
  plot_data <- plot_data %>%
    mutate(method = factor(method, levels = all_methods))

  # Apply filters
  if (!is.null(gwas_meta_type)) {
    plot_data <- plot_data %>% filter(meta_analysis_type == gwas_meta_type)
  }

  if (!is.null(simulation_type)) {
    plot_data <- plot_data %>% filter(simulation_type == !!simulation_type)
  }

  if (!is.null(snp_heritability)) {
    plot_data <- plot_data %>% filter(snp_heritability == !!snp_heritability)
  }

  if (!is.null(n_traits_filter)) {
    plot_data <- plot_data %>% filter(n_traits == n_traits_filter)
  }

  if (!is.null(methods)) {
    plot_data <- plot_data %>% filter(method %in% methods)
  }

  # Create labels for plotting
  plot_data <- plot_data %>%
    mutate(
      matrix_type = case_when(
        ld_matrix == "raw_in-sample" ~ "Raw In-Sample",
        ld_matrix == "PC-corrected" ~ "PC-Corrected In-Sample",
        is_beta_inferred & matrix_construction == "meta" ~ "Beta-Inferred (Meta)",
        is_beta_inferred & matrix_construction == "weighted" ~ "Beta-Inferred (Weighted)",
        is_beta_inferred & matrix_construction == "population-specific" ~ "Beta-Inferred (Pop-Specific)",
        ld_matrix == "mvp_external" ~ "MVP (External)",
        TRUE ~ ld_matrix
      ),
      # Factor for consistent ordering
      matrix_type = factor(matrix_type, levels = c(
        "Raw In-Sample", "PC-Corrected In-Sample", "MVP (External)",
        "Beta-Inferred (Meta)", "Beta-Inferred (Weighted)", "Beta-Inferred (Pop-Specific)"
      )),
      # Format GWAS trait type labels
      gwas_trait_label = case_when(
        trait_type == "binary" ~ "Binary GWAS",
        trait_type == "continuous" ~ "Continuous GWAS",
        TRUE ~ trait_type
      ),
      # Format LD matrix trait types labels
      ld_trait_label = case_when(
        trait_types == "binary" ~ "Binary LD",
        trait_types == "continuous" ~ "Continuous LD",
        trait_types == "all" ~ "All-Trait LD",
        TRUE ~ trait_types
      )
    )

  # Get methods and matrix types present in filtered data
  methods_present <- unique(plot_data$method)
  methods_present <- methods_present[!is.na(methods_present)]

  matrix_types_present <- unique(plot_data$matrix_type)
  matrix_types_present <- matrix_types_present[!is.na(matrix_types_present)]

  # Set up metric-specific parameters
  if (metric == "variant_auprc") {
    y_column <- "variant_auprc"
    y_label <- "Variant-level AUPRC"
  } else if (metric == "mean_cs_size") {
    y_column <- "mean_cs_size"
    y_label <- "Mean Credible Set Size"
  }

  # Define shapes for different matrix types
  matrix_shapes <- c(
    "Raw In-Sample" = 16,              # filled circle
    "PC-Corrected In-Sample" = 17,     # filled triangle
    "MVP (External)" = 4,              # X shape
    "Beta-Inferred (Meta)" = 15,       # filled square
    "Beta-Inferred (Weighted)" = 18,   # filled diamond
    "Beta-Inferred (Pop-Specific)" = 8 # asterisk
  )

  # Create the faceted plot
  p <- ggplot(plot_data, aes(x = method, y = .data[[y_column]], color = method, shape = matrix_type)) +
    geom_point(size = 5, alpha = 0.8, position = position_dodge(width = 0.5)) +
    facet_grid(ld_trait_label ~ gwas_trait_label, scales = "free") +
    scale_color_viridis_d(name = "Fine-mapping Method", drop = FALSE, breaks = methods_present) +
    scale_shape_manual(
      name = "LD Matrix Type",
      values = matrix_shapes,
      drop = FALSE,
      breaks = matrix_types_present
    ) +
    guides(
      shape = guide_legend(nrow = 1, byrow = TRUE),
      color = guide_legend(nrow = 1)
    ) +
    labs(
      x = "Fine-mapping Method",
      y = y_label
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 18, color = "black"),
      axis.text = element_text(size = 16, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 17, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 16, color = "black"),
      legend.title = element_text(size = 17, face = "bold", color = "black"),
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = element_line(color = "black")
    )

  return(p)
}

# Function to process beta-inferred LD matrix names
# Convention: beta-inferred_<#traits>_<meta-analysis_type>_<trait_types>.<Meta or weighted>
# Only carma, susie, and rsparsepro have the .<Meta or weighted> suffix
process_beta_inferred_ld <- function(data) {
  data %>%
    mutate(
      # Check if it's a beta-inferred matrix
      is_beta_inferred = grepl("^beta-inferred", ld_matrix, ignore.case = TRUE),

      # Create a lowercase version for pattern matching
      ld_matrix_lower = tolower(ld_matrix),

      # Extract number of traits (second field)
      # Pattern: beta-inferred_<number>_...
      n_traits = ifelse(is_beta_inferred,
                        as.numeric(str_extract(ld_matrix_lower, "(?<=beta-inferred_)\\d+")),
                        NA),

      # Extract meta-analysis type used to create the LD matrix (third field: fixed or random)
      # This is different from the existing meta_analysis_type column which is for the GWAS
      # Pattern: _<number>_(fixed|random)_
      # Use str_match to capture the group, then extract the second column
      ld_matrix_meta_type = ifelse(is_beta_inferred,
                                    {
                                      matches <- str_match(ld_matrix_lower, "_\\d+_(fixed|random)_")
                                      matches[, 2]
                                    },
                                    NA_character_),

      # Extract trait types (fourth field: continuous, binary, or all)
      # Pattern: _(fixed|random)_(continuous|binary|all)
      trait_types = ifelse(is_beta_inferred,
                           {
                             matches <- str_match(ld_matrix_lower, "_(fixed|random)_(continuous|binary|all)")
                             matches[, 3]
                           },
                           NA_character_),

      # Extract matrix construction method (Meta or weighted, or population-specific for other methods)
      # Only carma, susie, and rsparsepro have this suffix
      matrix_construction = case_when(
        !is_beta_inferred ~ NA_character_,
        method %in% c("carma", "susie", "rsparsepro", "CARMA", "SuSiE", "RSparsePro") &
          grepl("\\.(meta|Meta)$", ld_matrix, ignore.case = TRUE) ~ "meta",
        method %in% c("carma", "susie", "rsparsepro", "CARMA", "SuSiE", "RSparsePro") &
          grepl("\\.(weighted|Weighted)$", ld_matrix, ignore.case = TRUE) ~ "weighted",
        is_beta_inferred ~ "population-specific",
        TRUE ~ NA_character_
      )
    ) %>%
    select(-ld_matrix_lower)  # Remove temporary column
}

# 3) Process dataframes ====

# Apply formatting and processing to all dataframes
cat("Processing dataframes...\n")
method_df <- method_df %>% format_method_names() %>% process_beta_inferred_ld()
causal_method_df <- causal_method_df %>% format_method_names() %>% process_beta_inferred_ld()
trait_method_df <- trait_method_df %>% format_method_names() %>% process_beta_inferred_ld()
trait_causal_method_df <- trait_causal_method_df %>% format_method_names() %>% process_beta_inferred_ld()

# 4) Make Coverage vs Power Plots for All-Trait Matrices ====

# Create output subdirectories
coverage_power_dir <- file.path(out_dir, "coverage_vs_power")
if (!dir.exists(coverage_power_dir)) {
  dir.create(coverage_power_dir, recursive = TRUE)
}

cat("Creating plots...\n")

# Coverage vs Power plots
cat("\n--- Coverage vs Power Plots ---\n")

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 2000-all-trait matrices
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 2000 all traits...\n")
p_cvp_homo_fixed_5e5_2000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 2000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=5e-5, 2000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.homo.fixed_5e-5.2000_all_traits.png"),
       p_cvp_homo_fixed_5e5_2000, width = 10, height = 8, dpi = 300)

# Homogenous effect sizes, Fixed-effects, 5e-5 heritability, 2000-all-trait matrices
cat("  - Homogenous Effect Sizes, Fixed-effects, 5e-5 heritability, 2000 all traits...\n")
p_cvp_matched_fixed_5e5_2000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 2000,
  simulation_type = "matched_to",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Effect Sizes, Fixed-Effects GWAS, h²=5e-5, 2000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.matched_effects.fixed_5e-5.2000_all_traits.png"),
       p_cvp_matched_fixed_5e5_2000, width = 10, height = 8, dpi = 300)

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 1000-all-trait matrices
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 1000 all traits...\n")
p_cvp_homo_fixed_5e5_1000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=5e-5, 1000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.homo.fixed_5e-5.1000_all_traits.png"),
       p_cvp_homo_fixed_5e5_1000, width = 10, height = 8, dpi = 300)

# Homogenous effect sizes, Fixed-effects, 5e-5 heritability, 1000-all-trait matrices
cat("  - Homogenous Effect Sizes, Fixed-effects, 5e-5 heritability, 1000 all traits...\n")
p_cvp_matched_fixed_5e5_1000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  simulation_type = "matched_to",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Effect Sizes, Fixed-Effects GWAS, h²=5e-5, 1000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.matched_effects.fixed_5e-5.1000_all_traits.png"),
       p_cvp_matched_fixed_5e5_1000, width = 10, height = 8, dpi = 300)

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 4000-all-trait matrices
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 4000 all traits...\n")
p_cvp_homo_fixed_5e5_4000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 4000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=5e-5, 4000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.homo.fixed_5e-5.4000_all_traits.png"),
       p_cvp_homo_fixed_5e5_4000, width = 10, height = 8, dpi = 300)

cat("\n--- Faceted Coverage vs Power Plots (by Trait Type and Causal Variant Count) ---\n")

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 2000-all-trait matrices, faceted
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 2000 all traits (faceted)...\n")
p_cvp_faceted_homo_fixed_5e5_2000 <- plot_coverage_vs_power_faceted(
  trait_causal_method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 2000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=5e-5, 2000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.faceted.homo.fixed_5e-5.2000_all_traits.png"),
       p_cvp_faceted_homo_fixed_5e5_2000, width = 14, height = 8, dpi = 300)

# Homogenous effect sizes, Fixed-effects, 5e-5 heritability, 2000-all-trait matrices, faceted
cat("  - Homogenous Effect Sizes, Fixed-effects, 5e-5 heritability, 2000 all traits (faceted)...\n")
p_cvp_faceted_matched_fixed_5e5_2000 <- plot_coverage_vs_power_faceted(
  trait_causal_method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 2000,
  simulation_type = "matched_to",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Effect Sizes, Fixed-Effects GWAS, h²=5e-5, 2000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.faceted.matched_effects.fixed_5e-5.2000_all_traits.png"),
       p_cvp_faceted_matched_fixed_5e5_2000, width = 14, height = 8, dpi = 300)

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 1000-all-trait matrices, faceted
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 1000 all traits (faceted)...\n")
p_cvp_faceted_homo_fixed_5e5_1000 <- plot_coverage_vs_power_faceted(
  trait_causal_method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=5e-5, 1000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.faceted.homo.fixed_5e-5.1000_all_traits.png"),
       p_cvp_faceted_homo_fixed_5e5_1000, width = 14, height = 8, dpi = 300)

# Homogenous effect sizes, Fixed-effects, 5e-5 heritability, 1000-all-trait matrices, faceted
cat("  - Homogenous Effect Sizes, Fixed-effects, 5e-5 heritability, 1000 all traits (faceted)...\n")
p_cvp_faceted_matched_fixed_5e5_1000 <- plot_coverage_vs_power_faceted(
  trait_causal_method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  simulation_type = "matched_to",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Effect Sizes, Fixed-Effects GWAS, h²=5e-5, 1000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.faceted.matched_effects.fixed_5e-5.1000_all_traits.png"),
       p_cvp_faceted_matched_fixed_5e5_1000, width = 14, height = 8, dpi = 300)

# Coverage vs power faceted by heritability and simulation type (1000 traits)
cat("  - Both simulation types, 1000 traits, faceted by simulation type and heritability, fixed-effects...\n")
p_cvp_herit_1000_fixed <- plot_coverage_vs_power_by_heritability(
  method_df,
  gwas_meta_type = "fixed",
  n_traits_filter = 1000,
  matrix_type_filter = "all"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.by_heritability.fixed.1000_all_traits.png"),
       p_cvp_herit_1000_fixed, width = 16, height = 10, dpi = 300)

# Coverage vs power faceted by GWAS meta-analysis type and simulation type (1000 traits, meta-analysis methods only)
cat("  - Both simulation types, 1000 traits, faceted by simulation type and meta-analysis type, h²=5e-5...\n")
p_cvp_meta_1000_5e5 <- plot_coverage_vs_power_by_meta_type(
  method_df,
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  matrix_type_filter = "all",
  methods = c("CARMA", "SuSiE", "RSparsePro")
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.by_meta_type.5e-5.1000_all_traits.png"),
       p_cvp_meta_1000_5e5, width = 16, height = 10, dpi = 300)

# Coverage vs power faceted by GWAS trait type and LD matrix trait types
cat("  Creating coverage vs power plots faceted by GWAS trait type and LD matrix trait types...\n")
cat("  - Homogenous Heritability, Fixed-effects, 5e-5, 1000 traits...\n")
p_cvp_by_trait_types_homo_fixed_5e5_1000 <- plot_coverage_vs_power_by_trait_types(
  trait_method_df,
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 5e-5,
  n_traits_filter = 1000
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.by_trait_types.homo.fixed_5e-5.1000.png"),
       p_cvp_by_trait_types_homo_fixed_5e5_1000, width = 14, height = 10, dpi = 300)

# Metric scatterplots faceted by GWAS trait type and LD matrix trait types
cat("  Creating metric scatterplots faceted by GWAS trait type and LD matrix trait types...\n")

# Performance metrics for realistic LD matrices vs p-value threshold
cat("  Creating performance plots for realistic LD matrices...\n")

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, realistic LD
cat("  - Homogenous Heritability, Fixed-effects, 5e-5, realistic LD matrices (0kb padding)...\n")
p_perf_realistic_homo_fixed_5e5 <- plot_performance_realistic_ld(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5
)
ggsave(file.path(coverage_power_dir, "performance_realistic_ld.homo.fixed_5e-5.png"),
       p_perf_realistic_homo_fixed_5e5, width = 18, height = 8, dpi = 300)

# Mean CS size and failure rate for realistic LD matrices
cat("  - Mean CS Size and Failure Rate, realistic LD matrices (0kb padding)...\n")
p_cs_failure_realistic_homo_fixed_5e5 <- plot_cs_size_failure_realistic_ld(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5
)
ggsave(file.path(coverage_power_dir, "cs_size_failure_realistic_ld.homo.fixed_5e-5.png"),
       p_cs_failure_realistic_homo_fixed_5e5, width = 14, height = 8, dpi = 300)

# Create and save standalone legend
cat("  Creating standalone legend for fine-mapping methods and LD matrix types...\n")
legend_grob <- create_finemapping_legend()
png(file.path(out_dir, "finemapping_legend.png"), width = 7, height = 5, units = "in", res = 300)
grid::grid.draw(legend_grob)
dev.off()

# 5) Make Performance vs. Number of Traits Plots ====

# Create output subdirectories
performance_ntraits_dir <- file.path(out_dir, "performance_vs_ntraits")
if (!dir.exists(performance_ntraits_dir)) {
  dir.create(performance_ntraits_dir, recursive = TRUE)
}

cat("\n--- Performance vs. Number of Traits Plots ---\n")

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, all trait types faceted
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, all trait types...\n")
p_perf_ntraits_homo_fixed_5e5 <- plot_performance_vs_ntraits(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  simulation_type = "homo",
  trait_types = "all"
)
ggsave(file.path(performance_ntraits_dir, "performance_vs_ntraits.homo.fixed_5e-5.png"),
       p_perf_ntraits_homo_fixed_5e5, width = 14, height = 12, dpi = 300)

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, faceted by trait types
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, faceted by trait types...\n")
p_perf_ntraits_by_trait_type_homo_fixed_5e5 <- plot_performance_vs_ntraits_by_trait_type(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  simulation_type = "homo"
)
ggsave(file.path(performance_ntraits_dir, "performance_vs_ntraits_by_trait_type.homo.fixed_5e-5.png"),
       p_perf_ntraits_by_trait_type_homo_fixed_5e5, width = 18, height = 12, dpi = 300)

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, dropping causal variant, faceted by trait types
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, dropping causal variant, faceted by trait types...\n")
p_perf_ntraits_by_trait_type_drop_fixed_5e5 <- plot_performance_vs_ntraits_by_trait_type(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  simulation_type = "drop_causal",
  average_across_populations = TRUE
)
ggsave(file.path(performance_ntraits_dir, "performance_vs_ntraits_by_trait_type.drop_causal.fixed_5e-5.png"),
       p_perf_ntraits_by_trait_type_drop_fixed_5e5, width = 18, height = 12, dpi = 300)

# Failure rate vs number of traits, faceted by trait types
cat("  Creating failure rate vs number of traits plots...\n")

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, faceted by trait types
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, faceted by trait types...\n")
p_failure_rate_ntraits_by_trait_type_homo_fixed_5e5 <- plot_failure_rate_vs_ntraits_by_trait_type(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  simulation_type = "homo"
)
ggsave(file.path(performance_ntraits_dir, "failure_rate_vs_ntraits_by_trait_type.homo.fixed_5e-5.png"),
       p_failure_rate_ntraits_by_trait_type_homo_fixed_5e5, width = 18, height = 8, dpi = 300)

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, dropping causal variant, faceted by trait types
cat("  - Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, dropping causal variant, faceted by trait types...\n")
p_failure_rate_ntraits_by_trait_type_drop_fixed_5e5 <- plot_failure_rate_vs_ntraits_by_trait_type(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  simulation_type = "drop_causal",
  average_across_populations = TRUE
)
ggsave(file.path(performance_ntraits_dir, "failure_rate_vs_ntraits_by_trait_type.drop_causal.fixed_5e-5.png"),
       p_failure_rate_ntraits_by_trait_type_drop_fixed_5e5, width = 18, height = 8, dpi = 300)

cat("\nAll plots complete!\n")

# 6) Make Bar Plots ====

# Create output subdirectories
metrics_barplot_dir <- file.path(out_dir, "metric_barplots")
if (!dir.exists(metrics_barplot_dir)) {
  dir.create(metrics_barplot_dir, recursive = TRUE)
}

# Plot variant AUPRC for homogenous heritability, fixed-effects, 1000 traits
p_execution_time <- plot_metric_barplot(
  method_df,
  metric = "mean_execution_time_seconds",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  matrix_type_filter = "all",
  ld_construction_filter = "raw_in-sample"
)

ggsave(file.path(metrics_barplot_dir, "execution_times.raw_in-sample.fixed_5e-5.png"),
       p_execution_time, width = 18, height = 12, dpi = 300)

# 7) Make Scatter Plots ====

# Create output subdirectories
metrics_scatterplot_dir <- file.path(out_dir, "metric_scatterplots")
if (!dir.exists(metrics_scatterplot_dir)) {
  dir.create(metrics_scatterplot_dir, recursive = TRUE)
}

cat("\n--- Metric Scatter Plots ---\n")

# Plot variant AUPRC for homogenous heritability, fixed-effects, 1000 traits
cat("  - Variant AUPRC, homogenous heritability, fixed-effects, 5e-5 heritability, 1000 all traits...\n")
p_auprc_scatter <- plot_metric_scatterplot(
  method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.homo.fixed_5e-5.1000_all_traits.png"),
       p_auprc_scatter, width = 12, height = 10.1, dpi = 300)

# Plot mean credible set size for homogenous heritability, fixed-effects, 1000 traits
cat("  - Mean credible set size, homogenous heritability, fixed-effects, 5e-5 heritability, 1000 all traits...\n")
p_cs_size_scatter <- plot_metric_scatterplot(
  method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.homo.fixed_5e-5.1000_all_traits.png"),
       p_cs_size_scatter, width = 12, height = 10, dpi = 300)

# Plot variant AUPRC for matched_to (homogenous effect sizes), averaged across populations
cat("  - Variant AUPRC, homogenous effect sizes, fixed-effects, 5e-5 heritability, 1000 all traits (averaged)...\n")
p_auprc_scatter_matched <- plot_metric_scatterplot(
  method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  simulation_type = "matched_to",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  matrix_type_filter = "all",
  average_across_populations = TRUE
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.matched_to.fixed_5e-5.1000_all_traits.png"),
       p_auprc_scatter_matched, width = 12, height = 10, dpi = 300)

# Plot mean credible set size for matched_to (homogenous effect sizes), averaged across populations
cat("  - Mean credible set size, homogenous effect sizes, fixed-effects, 5e-5 heritability, 1000 all traits (averaged)...\n")
p_cs_size_scatter_matched <- plot_metric_scatterplot(
  method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  simulation_type = "matched_to",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  matrix_type_filter = "all",
  average_across_populations = TRUE
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.matched_to.fixed_5e-5.1000_all_traits.png"),
       p_cs_size_scatter_matched, width = 12, height = 10, dpi = 300)

# Variant AUPRC faceted by heritability and simulation type (1000 traits)
cat("  - Variant AUPRC, both simulation types, 1000 traits, faceted by simulation type and heritability, fixed-effects...\n")
p_auprc_herit_1000_fixed <- plot_metric_scatterplot_by_heritability(
  method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  n_traits_filter = 1000,
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.by_heritability.fixed.1000_all_traits.png"),
       p_auprc_herit_1000_fixed, width = 16, height = 10, dpi = 300)

# Mean credible set size faceted by heritability and simulation type (1000 traits)
cat("  - Mean credible set size, both simulation types, 1000 traits, faceted by simulation type and heritability, fixed-effects...\n")
p_cs_size_herit_1000_fixed <- plot_metric_scatterplot_by_heritability(
  method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  n_traits_filter = 1000,
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.by_heritability.fixed.1000_all_traits.png"),
       p_cs_size_herit_1000_fixed, width = 16, height = 10, dpi = 300)

# Variant AUPRC faceted by GWAS meta-analysis type and simulation type (1000 traits, meta-analysis methods only)
cat("  - Variant AUPRC, both simulation types, 1000 traits, faceted by simulation type and meta-analysis type, h²=5e-5...\n")
p_auprc_meta_1000_5e5 <- plot_metric_scatterplot_by_meta_type(
  method_df,
  metric = "variant_auprc",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  matrix_type_filter = "all",
  methods = c("CARMA", "SuSiE", "RSparsePro")
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.by_meta_type.5e-5.1000_all_traits.png"),
       p_auprc_meta_1000_5e5, width = 16, height = 10, dpi = 300)

# Mean credible set size faceted by GWAS meta-analysis type and simulation type (1000 traits, meta-analysis methods only)
cat("  - Mean credible set size, both simulation types, 1000 traits, faceted by simulation type and meta-analysis type, h²=5e-5...\n")
p_cs_size_meta_1000_5e5 <- plot_metric_scatterplot_by_meta_type(
  method_df,
  metric = "mean_cs_size",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  matrix_type_filter = "all",
  methods = c("CARMA", "SuSiE", "RSparsePro")
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.by_meta_type.5e-5.1000_all_traits.png"),
       p_cs_size_meta_1000_5e5, width = 16, height = 10, dpi = 300)

# Faceted scatter plots by trait type and causal variant count
cat("  Creating faceted scatter plots by trait type and causal variant count...\n")

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 2000 all traits - variant AUPRC
cat("  - Homogenous Heritability, Fixed-effects, 5e-5, 2000 all traits (faceted) - variant AUPRC...\n")
p_auprc_faceted_homo_fixed_5e5_2000 <- plot_metric_scatterplot_faceted(
  trait_causal_method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 2000,
  simulation_type = "homo",
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.faceted.homo.fixed_5e-5.2000_all_traits.png"),
       p_auprc_faceted_homo_fixed_5e5_2000, width = 14, height = 8, dpi = 300)

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 2000 all traits - mean CS size
cat("  - Homogenous Heritability, Fixed-effects, 5e-5, 2000 all traits (faceted) - mean CS size...\n")
p_cs_size_faceted_homo_fixed_5e5_2000 <- plot_metric_scatterplot_faceted(
  trait_causal_method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 2000,
  simulation_type = "homo",
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.faceted.homo.fixed_5e-5.2000_all_traits.png"),
       p_cs_size_faceted_homo_fixed_5e5_2000, width = 14, height = 8, dpi = 300)

# Homogenous effect sizes, Fixed-effects, 5e-5 heritability, 2000 all traits - variant AUPRC
cat("  - Homogenous Effect Sizes, Fixed-effects, 5e-5, 2000 all traits (faceted) - variant AUPRC...\n")
p_auprc_faceted_matched_fixed_5e5_2000 <- plot_metric_scatterplot_faceted(
  trait_causal_method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 2000,
  simulation_type = "matched_to",
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.faceted.matched_effects.fixed_5e-5.2000_all_traits.png"),
       p_auprc_faceted_matched_fixed_5e5_2000, width = 14, height = 8, dpi = 300)

# Homogenous effect sizes, Fixed-effects, 5e-5 heritability, 2000 all traits - mean CS size
cat("  - Homogenous Effect Sizes, Fixed-effects, 5e-5, 2000 all traits (faceted) - mean CS size...\n")
p_cs_size_faceted_matched_fixed_5e5_2000 <- plot_metric_scatterplot_faceted(
  trait_causal_method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 2000,
  simulation_type = "matched_to",
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.faceted.matched_effects.fixed_5e-5.2000_all_traits.png"),
       p_cs_size_faceted_matched_fixed_5e5_2000, width = 14, height = 8, dpi = 300)

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 1000 all traits - variant AUPRC
cat("  - Homogenous Heritability, Fixed-effects, 5e-5, 1000 all traits (faceted) - variant AUPRC...\n")
p_auprc_faceted_homo_fixed_5e5_1000 <- plot_metric_scatterplot_faceted(
  trait_causal_method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  simulation_type = "homo",
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.faceted.homo.fixed_5e-5.1000_all_traits.png"),
       p_auprc_faceted_homo_fixed_5e5_1000, width = 14, height = 8, dpi = 300)

# Homogenous per-SNP heritability, Fixed-effects, 5e-5 heritability, 1000 all traits - mean CS size
cat("  - Homogenous Heritability, Fixed-effects, 5e-5, 1000 all traits (faceted) - mean CS size...\n")
p_cs_size_faceted_homo_fixed_5e5_1000 <- plot_metric_scatterplot_faceted(
  trait_causal_method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  simulation_type = "homo",
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.faceted.homo.fixed_5e-5.1000_all_traits.png"),
       p_cs_size_faceted_homo_fixed_5e5_1000, width = 14, height = 8, dpi = 300)

# Homogenous effect sizes, Fixed-effects, 5e-5 heritability, 1000 all traits - variant AUPRC
cat("  - Homogenous Effect Sizes, Fixed-effects, 5e-5, 1000 all traits (faceted) - variant AUPRC...\n")
p_auprc_faceted_matched_fixed_5e5_1000 <- plot_metric_scatterplot_faceted(
  trait_causal_method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  simulation_type = "matched_to",
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.faceted.matched_effects.fixed_5e-5.1000_all_traits.png"),
       p_auprc_faceted_matched_fixed_5e5_1000, width = 14, height = 8, dpi = 300)

# Homogenous effect sizes, Fixed-effects, 5e-5 heritability, 1000 all traits - mean CS size
cat("  - Homogenous Effect Sizes, Fixed-effects, 5e-5, 1000 all traits (faceted) - mean CS size...\n")
p_cs_size_faceted_matched_fixed_5e5_1000 <- plot_metric_scatterplot_faceted(
  trait_causal_method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  snp_heritability = 5e-5,
  n_traits_filter = 1000,
  simulation_type = "matched_to",
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.faceted.matched_effects.fixed_5e-5.1000_all_traits.png"),
       p_cs_size_faceted_matched_fixed_5e5_1000, width = 14, height = 8, dpi = 300)

# Variant AUPRC
cat("  - Variant AUPRC - Homogenous Heritability, Fixed-effects, 5e-5, 1000 traits...\n")
p_auprc_by_trait_types_homo_fixed_5e5_1000 <- plot_metric_scatterplot_by_trait_types(
  trait_method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 5e-5,
  n_traits_filter = 1000
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.by_trait_types.homo.fixed_5e-5.1000.png"),
       p_auprc_by_trait_types_homo_fixed_5e5_1000, width = 14, height = 10, dpi = 300)

# Mean credible set size
cat("  - Mean CS Size - Homogenous Heritability, Fixed-effects, 5e-5, 1000 traits...\n")
p_cs_size_by_trait_types_homo_fixed_5e5_1000 <- plot_metric_scatterplot_by_trait_types(
  trait_method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 5e-5,
  n_traits_filter = 1000
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.by_trait_types.homo.fixed_5e-5.1000.png"),
       p_cs_size_by_trait_types_homo_fixed_5e5_1000, width = 14, height = 10, dpi = 300)

cat("\nScatter plots complete!\n")

# ---------------------------------------------------------------------
# Also generate the same plots for snp_heritability = 3e-4 (where applicable)
# ---------------------------------------------------------------------
cat("\n--- Repeating selected plots for snp_heritability = 3e-4 ---\n")

# Coverage vs Power plots (repeat several key views for h2 = 3e-4)
cat("  - Coverage vs Power: Homogenous h2=3e-4, Fixed-effects, 2000 all traits...\n")
p_cvp_homo_fixed_3e4_2000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  n_traits_filter = 2000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=3e-4, 2000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.homo.fixed_3e-4.2000_all_traits.png"),
       p_cvp_homo_fixed_3e4_2000, width = 10, height = 8, dpi = 300)

cat("  - Coverage vs Power: Matched effects h2=3e-4, Fixed-effects, 2000 all traits...\n")
p_cvp_matched_fixed_3e4_2000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  n_traits_filter = 2000,
  simulation_type = "matched_to",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Effect Sizes, Fixed-Effects GWAS, h²=3e-4, 2000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.matched_effects.fixed_3e-4.2000_all_traits.png"),
       p_cvp_matched_fixed_3e4_2000, width = 10, height = 8, dpi = 300)

cat("  - Coverage vs Power: Homogenous h2=3e-4, Fixed-effects, 1000 all traits...\n")
p_cvp_homo_fixed_3e4_1000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  n_traits_filter = 1000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=3e-4, 1000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.homo.fixed_3e-4.1000_all_traits.png"),
       p_cvp_homo_fixed_3e4_1000, width = 10, height = 8, dpi = 300)

cat("  - Coverage vs Power: Matched effects h2=3e-4, Fixed-effects, 1000 all traits...\n")
p_cvp_matched_fixed_3e4_1000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  n_traits_filter = 1000,
  simulation_type = "matched_to",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Effect Sizes, Fixed-Effects GWAS, h²=3e-4, 1000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.matched_effects.fixed_3e-4.1000_all_traits.png"),
       p_cvp_matched_fixed_3e4_1000, width = 10, height = 8, dpi = 300)

cat("  - Coverage vs Power: Homogenous h2=3e-4, Fixed-effects, 4000 all traits...\n")
p_cvp_homo_fixed_3e4_4000 <- plot_coverage_vs_power_filtered(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  n_traits_filter = 4000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=3e-4, 4000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.homo.fixed_3e-4.4000_all_traits.png"),
       p_cvp_homo_fixed_3e4_4000, width = 10, height = 8, dpi = 300)

# Faceted coverage vs power (trait/causal) for h2 = 3e-4
cat("  - Faceted coverage vs power (h2=3e-4), 2000 and 1000 traits...\n")
p_cvp_faceted_homo_fixed_3e4_2000 <- plot_coverage_vs_power_faceted(
  trait_causal_method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  n_traits_filter = 2000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=3e-4, 2000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.faceted.homo.fixed_3e-4.2000_all_traits.png"),
       p_cvp_faceted_homo_fixed_3e4_2000, width = 14, height = 8, dpi = 300)

p_cvp_faceted_matched_fixed_3e4_2000 <- plot_coverage_vs_power_faceted(
  trait_causal_method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  n_traits_filter = 2000,
  simulation_type = "matched_to",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Effect Sizes, Fixed-Effects GWAS, h²=3e-4, 2000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.faceted.matched_effects.fixed_3e-4.2000_all_traits.png"),
       p_cvp_faceted_matched_fixed_3e4_2000, width = 14, height = 8, dpi = 300)

p_cvp_faceted_homo_fixed_3e4_1000 <- plot_coverage_vs_power_faceted(
  trait_causal_method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  n_traits_filter = 1000,
  simulation_type = "homo",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Heritability, Fixed-Effects GWAS, h²=3e-4, 1000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.faceted.homo.fixed_3e-4.1000_all_traits.png"),
       p_cvp_faceted_homo_fixed_3e4_1000, width = 14, height = 8, dpi = 300)

p_cvp_faceted_matched_fixed_3e4_1000 <- plot_coverage_vs_power_faceted(
  trait_causal_method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  n_traits_filter = 1000,
  simulation_type = "matched_to",
  matrix_type_filter = "all",
  title_suffix = "\nHomogenous Effect Sizes, Fixed-Effects GWAS, h²=3e-4, 1000 All Traits"
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.faceted.matched_effects.fixed_3e-4.1000_all_traits.png"),
       p_cvp_faceted_matched_fixed_3e4_1000, width = 14, height = 8, dpi = 300)

# Realistic LD performance and CS/failure plots for h2 = 3e-4
cat("  - Performance and CS/failure (realistic LD) for h2=3e-4...\n")
p_perf_realistic_homo_fixed_3e4 <- plot_performance_realistic_ld(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4
)
ggsave(file.path(coverage_power_dir, "performance_realistic_ld.homo.fixed_3e-4.png"),
       p_perf_realistic_homo_fixed_3e4, width = 18, height = 8, dpi = 300)

p_cs_failure_realistic_homo_fixed_3e4 <- plot_cs_size_failure_realistic_ld(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4
)
ggsave(file.path(coverage_power_dir, "cs_size_failure_realistic_ld.homo.fixed_3e-4.png"),
       p_cs_failure_realistic_homo_fixed_3e4, width = 14, height = 8, dpi = 300)

# Performance vs number of traits for h2 = 3e-4
cat("  - Performance vs n_traits for h2=3e-4...\n")
p_perf_ntraits_homo_fixed_3e4 <- plot_performance_vs_ntraits(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  simulation_type = "homo",
  trait_types = "all"
)
ggsave(file.path(performance_ntraits_dir, "performance_vs_ntraits.homo.fixed_3e-4.png"),
       p_perf_ntraits_homo_fixed_3e4, width = 14, height = 12, dpi = 300)

p_perf_ntraits_by_trait_type_homo_fixed_3e4 <- plot_performance_vs_ntraits_by_trait_type(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  simulation_type = "homo"
)
ggsave(file.path(performance_ntraits_dir, "performance_vs_ntraits_by_trait_type.homo.fixed_3e-4.png"),
       p_perf_ntraits_by_trait_type_homo_fixed_3e4, width = 18, height = 12, dpi = 300)

p_perf_ntraits_by_trait_type_drop_fixed_3e4 <- plot_performance_vs_ntraits_by_trait_type(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  simulation_type = "drop_causal",
  average_across_populations = TRUE
)
ggsave(file.path(performance_ntraits_dir, "performance_vs_ntraits_by_trait_type.drop_causal.fixed_3e-4.png"),
       p_perf_ntraits_by_trait_type_drop_fixed_3e4, width = 18, height = 12, dpi = 300)

# Failure rate vs n_traits for h2 = 3e-4
cat("  - Failure rate vs n_traits for h2=3e-4...\n")
p_failure_rate_ntraits_by_trait_type_homo_fixed_3e4 <- plot_failure_rate_vs_ntraits_by_trait_type(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  simulation_type = "homo"
)
ggsave(file.path(performance_ntraits_dir, "failure_rate_vs_ntraits_by_trait_type.homo.fixed_3e-4.png"),
       p_failure_rate_ntraits_by_trait_type_homo_fixed_3e4, width = 18, height = 8, dpi = 300)

p_failure_rate_ntraits_by_trait_type_drop_fixed_3e4 <- plot_failure_rate_vs_ntraits_by_trait_type(
  method_df,
  gwas_meta_type = "fixed",
  snp_heritability = 3e-4,
  simulation_type = "drop_causal",
  average_across_populations = TRUE
)
ggsave(file.path(performance_ntraits_dir, "failure_rate_vs_ntraits_by_trait_type.drop_causal.fixed_3e-4.png"),
       p_failure_rate_ntraits_by_trait_type_drop_fixed_3e4, width = 18, height = 8, dpi = 300)

# Metric barplot (execution time) for h2 = 3e-4 (raw in-sample)
cat("  - Execution time barplot for h2=3e-4...\n")
p_execution_time_3e4 <- plot_metric_barplot(
  method_df,
  metric = "mean_execution_time_seconds",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 3e-4,
  n_traits_filter = 1000,
  matrix_type_filter = "all",
  ld_construction_filter = "raw_in-sample"
)
ggsave(file.path(metrics_barplot_dir, "execution_times.raw_in-sample.fixed_3e-4.png"),
       p_execution_time_3e4, width = 18, height = 12, dpi = 300)

# Metric scatterplots for h2 = 3e-4
cat("  - Metric scatterplots for h2=3e-4...\n")
p_auprc_scatter_3e4 <- plot_metric_scatterplot(
  method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 3e-4,
  n_traits_filter = 1000,
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.homo.fixed_3e-4.1000_all_traits.png"),
       p_auprc_scatter_3e4, width = 12, height = 10, dpi = 300)

p_cs_size_scatter_3e4 <- plot_metric_scatterplot(
  method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 3e-4,
  n_traits_filter = 1000,
  matrix_type_filter = "all"
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.homo.fixed_3e-4.1000_all_traits.png"),
       p_cs_size_scatter_3e4, width = 12, height = 10, dpi = 300)

p_auprc_scatter_matched_3e4 <- plot_metric_scatterplot(
  method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  simulation_type = "matched_to",
  snp_heritability = 3e-4,
  n_traits_filter = 1000,
  matrix_type_filter = "all",
  average_across_populations = TRUE
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.matched_to.fixed_3e-4.1000_all_traits.png"),
       p_auprc_scatter_matched_3e4, width = 12, height = 10, dpi = 300)

p_cs_size_scatter_matched_3e4 <- plot_metric_scatterplot(
  method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  simulation_type = "matched_to",
  snp_heritability = 3e-4,
  n_traits_filter = 1000,
  matrix_type_filter = "all",
  average_across_populations = TRUE
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.matched_to.fixed_3e-4.1000_all_traits.png"),
       p_cs_size_scatter_matched_3e4, width = 12, height = 10, dpi = 300)

# Coverage vs power faceted by GWAS trait type and LD matrix trait types
cat("  Creating coverage vs power plots faceted by GWAS trait type and LD matrix trait types...\n")
cat("  - Homogenous Heritability, Fixed-effects, 3e-4, 1000 traits...\n")
p_cvp_by_trait_types_homo_fixed_3e4_1000 <- plot_coverage_vs_power_by_trait_types(
  trait_method_df,
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 3e-4,
  n_traits_filter = 1000
)
ggsave(file.path(coverage_power_dir, "coverage_vs_power.by_trait_types.homo.fixed_3e-4.1000.png"),
       p_cvp_by_trait_types_homo_fixed_3e4_1000, width = 14, height = 10, dpi = 300)

# Variant AUPRC
cat("  - Variant AUPRC - Homogenous Heritability, Fixed-effects, 3e-4, 1000 traits...\n")
p_auprc_by_trait_types_homo_fixed_3e4_1000 <- plot_metric_scatterplot_by_trait_types(
  trait_method_df,
  metric = "variant_auprc",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 3e-4,
  n_traits_filter = 1000
)
ggsave(file.path(metrics_scatterplot_dir, "variant_auprc.by_trait_types.homo.fixed_3e-4.1000.png"),
       p_auprc_by_trait_types_homo_fixed_3e4_1000, width = 14, height = 10, dpi = 300)

# Mean credible set size
cat("  - Mean CS Size - Homogenous Heritability, Fixed-effects, 3e-4, 1000 traits...\n")
p_cs_size_by_trait_types_homo_fixed_3e4_1000 <- plot_metric_scatterplot_by_trait_types(
  trait_method_df,
  metric = "mean_cs_size",
  gwas_meta_type = "fixed",
  simulation_type = "homo",
  snp_heritability = 3e-4,
  n_traits_filter = 1000
)
ggsave(file.path(metrics_scatterplot_dir, "mean_cs_size.by_trait_types.homo.fixed_3e-4.1000.png"),
       p_cs_size_by_trait_types_homo_fixed_3e4_1000, width = 14, height = 10, dpi = 300)

cat("\nCompleted plots for snp_heritability = 3e-4\n")
