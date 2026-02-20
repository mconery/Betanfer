# Genomic Fine-mapping Performance Analysis Script - Summary Level Data with True AUPRC and Timing
# Author: Mitch Conery
# Date: September 26, 2025

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(tidyr)
library(stringr)
library(viridis)
library(PRROC)  # For precision-recall curves

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
    mutate(ld_matrix = str_replace_all(ld_matrix, "_fixed","")) %>% 
    mutate(ld_size = as.numeric(str_replace_all(ld_matrix, "sub-sample_", ""))) %>% 
    mutate(ld_matrix = str_replace_all(ld_matrix, "_", " ") %>%
           str_to_title()) %>% mutate(ld_matrix = str_replace_all(ld_matrix, "Pc", "PC"))
}

# Function to read and prepare summary-level fine-mapping data
read_finemapping_summary <- function(file_path) {
  # Read the data (adjust based on your file format)
  if (grepl("\\.csv$", file_path)) {
    data <- read_csv(file_path)
  } else if (grepl("\\.(txt|tsv)$", file_path)) {
    data <- read_tsv(file_path)
  } else {
    stop("Unsupported file format. Please use CSV or TSV.")
  }
  # Check for required columns
  required_cols <- c("phenotype_number", "causal_variant_count", "causal_in_credible_sets",
                     "avg_variants_per_cs", "causal_variant_pips", "method", "ld_matrix", "total_variants")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  # Clean the data and handle failed convergence cases
  data <- data %>%
  mutate(
      # Mark convergence failures
      converged = !grepl("nan", causal_variant_pips, ignore.case = TRUE) & avg_variants_per_cs > 0,
      # Parse PIPs (handle comma-separated values and convert to numeric)
      causal_pips_list = map(causal_variant_pips, function(x) {
        if (grepl("nan", x, ignore.case = TRUE) || is.na(x)) {
          return(numeric(0))
        }
        as.numeric(str_split(as.character(x), ",")[[1]])
      })
    )
  # Apply formatting to method and LD matrix names
  data <- format_method_names(data)
  data <- format_ld_matrix_names(data)

  return(data)
}

# Function to calculate coverage
# Coverage = proportion of loci where at least one causal variant is in a credible set
calculate_coverage <- function(data) {
  coverage_stats <- data %>%
    # Only consider converged cases
    group_by(method, ld_size, causal_variant_count) %>%
    summarise(
      # Coverage: proportion of credible sets with at least one causal variant
      total_credible_sets_with_causal = sum(credible_sets_with_causal, na.rm = TRUE),
      total_credible_sets = sum(num_credible_sets, na.rm = TRUE),
      coverage = total_credible_sets_with_causal/total_credible_sets,
      n_loci = n(),
      .groups = "drop"
    )
  return(coverage_stats)
}

# Function to calculate power
# Power = proportion of all true causal variants that are captured in credible sets
calculate_power <- function(data) {
  power_stats <- data %>%
    # Only consider converged cases
    filter(converged == TRUE) %>%
    group_by(method, ld_size, causal_variant_count) %>%
    summarise(
      # Power: total causal variants in credible sets divided by total causal variants
      total_causal_variants = sum(causal_variant_count, na.rm = TRUE),
      total_causal_in_cs = sum(causal_in_credible_sets, na.rm = TRUE),
      power = total_causal_in_cs / total_causal_variants,
      .groups = "drop"
    )
  return(power_stats)
}

# Function to calculate credible set sizes
calculate_credible_set_sizes <- function(data) {
  size_stats <- data %>%
    # Only consider converged cases with credible sets
    filter(converged == TRUE, avg_variants_per_cs > 0) %>%
    group_by(method, ld_size, causal_variant_count) %>%
    summarise(
      total_variants_in_cs = sum(avg_variants_per_cs * num_credible_sets, na.rm = TRUE),
      total_cs = sum(num_credible_sets, na.rm = TRUE),
      mean_credible_set_size = total_variants_in_cs / total_cs,
      n_loci_with_cs = n(),
      .groups = "drop"
    )
  return(size_stats)
}

# Function to calculate true AUPRC using the new non-causal PIP count columns
calculate_true_auprc <- function(data) {
  auprc_stats <- data %>%
    mutate(causal_variant_pips = str_replace_all(causal_variant_pips, "nan", "0")) %>%
    filter(converged == TRUE) %>%
    group_by(method, ld_size, causal_variant_count) %>%
    do({
      df <- .
      if (nrow(df) == 0) {
        data.frame(auprc = NA, n_loci = 0)
      } else {
        # Initialize vectors for precision and recall calculations
        precisions <- numeric()
        recalls <- numeric()
        thresholds <- seq(0.01, 1.00, by = 0.01)
        # Get total numbers across all loci in this group
        total_causal_variants <- sum(df$causal_variant_count, na.rm = TRUE)
        for (thresh in thresholds) {
          # Count causal variants with PIP >= threshold
          causal_above_threshold <- 0
          for (i in 1:nrow(df)) {
            pips <- df$causal_pips_list[[i]]
            if (length(pips) > 0) {
              causal_above_threshold <- causal_above_threshold + sum(pips >= thresh, na.rm = TRUE)
            }
          }
          # Count non-causal variants with PIP >= threshold
          # Use the appropriate column name based on threshold
          col_name <- paste0("non_causal_pip_ge_", sprintf("%.2f", thresh))
          if (col_name %in% names(df)) {
            non_causal_above_threshold <- sum(df[[col_name]], na.rm = TRUE)
          } else {
            non_causal_above_threshold <- 0
          }
          # Calculate precision and recall
          tp <- causal_above_threshold
          fp <- non_causal_above_threshold
          precision <- if (tp + fp > 0) tp / (tp + fp) else 1
          recall <- if (total_causal_variants > 0) tp / total_causal_variants else 0
          precisions <- c(precisions, precision)
          recalls <- c(recalls, recall)
        }
        # Calculate AUPRC using trapezoidal rule
        # Sort by recall for proper integration
        order_idx <- order(recalls)
        recall_sorted <- recalls[order_idx]
        precision_sorted <- precisions[order_idx]
        # Calculate AUC using trapezoidal rule
        auprc <- 0
        if (length(recall_sorted) > 1) {
          for (k in 2:length(recall_sorted)) {
            auprc <- auprc + (recall_sorted[k] - recall_sorted[k-1]) *
              (precision_sorted[k] + precision_sorted[k-1]) / 2
          }
        }
        data.frame(auprc = auprc, n_loci = nrow(df))
      }
    }) %>%
    ungroup()
  return(auprc_stats)
}

# Function to calculate convergence rates
calculate_convergence_rates <- function(data) {
  convergence_stats <- data %>%
    group_by(method, ld_size, causal_variant_count) %>%
    summarise(
      total_loci = n(),
      converged_loci = sum(converged, na.rm = TRUE),
      convergence_rate = mean(converged, na.rm = TRUE),
      .groups = "drop"
    )
  return(convergence_stats)
}

# NEW FUNCTION: Calculate execution time statistics
calculate_execution_times <- function(data) {
  # Check if execution time data is available
  if (!"execution_time_seconds" %in% names(data)) {
    cat("Warning: execution_time_seconds column not found. Skipping timing analysis.\n")
    return(NULL)
  }
  timing_stats <- data %>%
    filter(!is.na(execution_time_seconds), converged == TRUE) %>%
    group_by(method, ld_size, causal_variant_count) %>%
    summarise(
      mean_execution_time = mean(execution_time_seconds, na.rm = TRUE),
      median_execution_time = median(execution_time_seconds, na.rm = TRUE),
      min_execution_time = min(execution_time_seconds, na.rm = TRUE),
      max_execution_time = max(execution_time_seconds, na.rm = TRUE),
      sd_execution_time = sd(execution_time_seconds, na.rm = TRUE),
      n_loci_with_timing = n(),
      .groups = "drop"
    )
  return(timing_stats)
}

# Function to create comprehensive performance metrics from summary data
calculate_all_metrics_summary <- function(data, causal_count_filter = NULL) {
  # Optional filtering by number of causal variants
  if (!is.null(causal_count_filter)) {
    data <- data %>%
      filter(causal_variant_count %in% causal_count_filter)
  }
  # Calculate all metrics
  coverage <- calculate_coverage(data)
  power <- calculate_power(data)
  sizes <- calculate_credible_set_sizes(data)
  auprc <- calculate_true_auprc(data)  # Use the new true AUPRC calculation
  convergence <- calculate_convergence_rates(data)
  timing <- calculate_execution_times(data)  # NEW: Calculate timing statistics
  # Combine all metrics
  all_metrics <- coverage %>%
    full_join(power, by = c("method", "ld_size", "causal_variant_count")) %>%
    full_join(sizes, by = c("method", "ld_size", "causal_variant_count")) %>%
    full_join(auprc, by = c("method", "ld_size", "causal_variant_count")) %>%
    full_join(convergence, by = c("method", "ld_size", "causal_variant_count"))
  # Add timing metrics if available
  if (!is.null(timing)) {
    all_metrics <- all_metrics %>%
      full_join(timing, by = c("method", "ld_size", "causal_variant_count"))
  }
  # Reorder columns for better readability
  base_cols <- c("method", "ld_size", "causal_variant_count", "coverage", "power", "auprc",
                 "mean_credible_set_size", "convergence_rate")
  # Add timing columns if they exist
  timing_cols <- c("mean_execution_time", "median_execution_time", "min_execution_time", 
                   "max_execution_time", "sd_execution_time", "n_loci_with_timing")
  existing_timing_cols <- intersect(timing_cols, names(all_metrics))
  final_cols <- c(base_cols, existing_timing_cols)
  other_cols <- setdiff(names(all_metrics), final_cols)
  all_metrics <- all_metrics %>%
    select(all_of(c(final_cols, other_cols))) %>%
    rename(n_loci = n_loci.x)  # Handle duplicate n_loci columns
  return(all_metrics)
}

# Function to create coverage scatter plot
plot_coverage_summary <- function(metrics_data, title_suffix = "") {
  p <- ggplot(metrics_data, aes(x = ld_size, y = coverage)) +
    geom_point(aes(color = method, size = factor(causal_variant_count)),
               alpha = 0.8) +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(name = "LD Type",
                       values = c(16, 17, 15, 18, 8, 13, 10, 12, 14, 11)) +
    scale_size_manual(name = "Causal Variants",
                      values = c("1" = 4, "2" = 6, "3" = 8)) +
    labs(
      y = "Coverage",
      x = "Sample Size of LD Matrix"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size=18),
      axis.text = element_text(size=18),
      legend.position = "bottom",
      legend.text = element_text(size=16),
      legend.title  = element_text(size=16),
      legend.box = "vertical",
      legend.direction = "horizontal",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(0, 1))
  return(p)
}

# Function to create power scatter plot
plot_power_summary <- function(metrics_data, title_suffix = "") {
  p <- ggplot(metrics_data, aes(x = ld_size, y = power)) +
    geom_point(aes(color = method, size = factor(causal_variant_count)),
               alpha = 0.8) +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(name = "LD Type",
                       values = c(16, 17, 15, 18, 8, 13, 10, 12, 14, 11)) +
    scale_size_manual(name = "Causal Variants",
                      values = c("1" = 4, "2" = 6, "3" = 8)) +
    labs(
      y = "Power",
      x = "Sample Size of LD Matrix"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size=18),
      axis.text = element_text(size=18),
      legend.position = "bottom",
      legend.text = element_text(size=16),
      legend.title  = element_text(size=16),
      legend.box = "vertical",
      legend.direction = "horizontal",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(0, 1))
  return(p)
}

# Function to create auprc scatter plot
plot_auprc_summary <- function(metrics_data, title_suffix = "") {
  p <- ggplot(metrics_data, aes(x = ld_size, y = auprc)) +
    geom_point(aes(color = method, size = factor(causal_variant_count)),
               alpha = 0.8) +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(name = "LD Type",
                       values = c(16, 17, 15, 18, 8, 13, 10, 12, 14, 11)) +
    scale_size_manual(name = "Causal Variants",
                      values = c("1" = 4, "2" = 6, "3" = 8)) +
    labs(
      y = "AUPRC",
      x = "Sample Size of LD Matrix"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size=18),
      axis.text = element_text(size=18),
      legend.position = "bottom",
      legend.text = element_text(size=16),
      legend.title  = element_text(size=16),
      legend.box = "vertical",
      legend.direction = "horizontal",
      panel.grid.minor = element_blank()
    ) +
    coord_cartesian(ylim = c(0, 1))
  return(p)
}

# Function to create credible set size plot
plot_credible_set_sizes_summary <- function(metrics_data, title_suffix = "") {
  p <- ggplot(metrics_data, aes(x = ld_size, y = mean_credible_set_size)) +
    geom_point(aes(color = method, size = factor(causal_variant_count)),
               alpha = 0.8) +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(name = "LD Type",
                       values = c(16, 17, 15, 18, 8, 13, 10, 12, 14, 11)) +
    scale_size_manual(name = "Causal Variants",
                      values = c("1" = 4, "2" = 6, "3" = 8)) +
    labs(
      y = "Mean Credible Set Size",
      x = "Sample Size of LD Matrix"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size=18),
      axis.text = element_text(size=18),
      legend.position = "bottom",
      legend.text = element_text(size=16),
      legend.title  = element_text(size=16),
      legend.box = "vertical",
      legend.direction = "horizontal",
      panel.grid.minor = element_blank()
    ) 
  return(p)
}

# Function to create credible set size plot
plot_execution_times <- function(metrics_data, title_suffix = "") {
  p <- ggplot(metrics_data, aes(x = ld_size, y = median_execution_time)) +
    geom_point(aes(color = method, size = factor(causal_variant_count)),
               alpha = 0.8) +
    scale_color_viridis_d(name = "Fine-mapping Method") +
    scale_shape_manual(name = "LD Type",
                       values = c(16, 17, 15, 18, 8, 13, 10, 12, 14, 11)) +
    scale_size_manual(name = "Causal Variants",
                      values = c("1" = 4, "2" = 6, "3" = 8)) +
    labs(
      y = "Median Execution Times",
      x = "Sample Size of LD Matrix"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size=18),
      axis.text = element_text(size=18),
      legend.position = "bottom",
      legend.text = element_text(size=16),
      legend.title  = element_text(size=16),
      legend.box = "vertical",
      legend.direction = "horizontal",
      panel.grid.minor = element_blank()
    ) 
  return(p)
}

# Main analysis function for summary data
run_finemapping_analysis_summary <- function(file_path, output_dir = "finemapping_results") {
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  # Read data
  cat("Reading fine-mapping summary data...\n")
  data <- read_finemapping_summary(file_path)
  # Print data summary
  cat("Data summary:\n")
  cat(" Total rows:", nrow(data), "\n")
  cat(" Unique phenotypes:", length(unique(data$phenotype_number)), "\n")
  cat(" Methods:", paste(unique(data$method), collapse = ", "), "\n")
  cat(" LD matrices:", paste(unique(data$ld_size), collapse = ", "), "\n")
  cat(" Causal variant counts:", paste(sort(unique(data$causal_variant_count)), collapse = ", "), "\n")
  cat(" Overall convergence rate:", round(mean(data$converged, na.rm = TRUE), 3), "\n")
  # Print execution time summary if available
  if ("execution_time_seconds" %in% names(data)) {
    exec_times <- data$execution_time_seconds[!is.na(data$execution_time_seconds)]
    if (length(exec_times) > 0) {
      cat(" Execution time summary:\n")
      cat("   Mean:", round(mean(exec_times), 2), "seconds\n")
      cat("   Median:", round(median(exec_times), 2), "seconds\n")
      cat("   Range:", round(min(exec_times), 2), "-", round(max(exec_times), 2), "seconds\n")
      cat("   Missing timing data:", sum(is.na(data$execution_time_seconds)), "files\n")
    }
  }
  cat("\n")
  # Calculate metrics for all data
  cat("Calculating performance metrics...\n")
  all_metrics <- calculate_all_metrics_summary(data)
  # Save metrics
  write_csv(all_metrics, file.path(output_dir, "performance_metrics_summary.csv"))
  # Create plots for all data
  cat("Creating plots for all data...\n")
  # Coverage vs Power scatter plot
  p1 <- plot_coverage_summary(all_metrics)
  ggsave(file.path(output_dir, "coverage_all.png"), p1,
         width = 7.5, height = 10, dpi = 300)
  p2 <- plot_power_summary(all_metrics)
  ggsave(file.path(output_dir, "power_all.png"), p2,
         width = 7.5, height = 10, dpi = 300)
  # AUPRC bar plot
  p3 <- plot_auprc_summary(all_metrics)
  ggsave(file.path(output_dir, "auprc_all.png"), p3,
         width = 14, height = 8, dpi = 300)
  # Credible set sizes bar plot
  p4 <- plot_credible_set_sizes_summary(all_metrics)
  ggsave(file.path(output_dir, "credible_set_sizes_all.png"), p4,
         width = 14, height = 8, dpi = 300)
  # NEW: Execution times plot
  p5 <- plot_execution_times(all_metrics)
  if (!is.null(p5)) {
    ggsave(file.path(output_dir, "execution_times_all.png"), p5,
           width = 14, height = 8, dpi = 300)
  }
  # Analysis by individual causal variant counts
  causal_counts <- sort(unique(data$causal_variant_count))
  for (causal_count in causal_counts) {
    cat(paste0("Creating analysis for ", causal_count, " causal variant(s)...\n"))
    # Calculate metrics for this stratum
    stratified_metrics <- calculate_all_metrics_summary(data, causal_count_filter = causal_count)
    if (nrow(stratified_metrics) > 0) {
      # Save stratified metrics
      write_csv(stratified_metrics,
                file.path(output_dir, paste0("performance_metrics_", causal_count, "_causal.csv")))
      # Create stratified plots
      suffix <- paste0(" (", causal_count, " causal variant",
                       if(causal_count > 1) "s" else "", ")")
      p1_strat <- plot_coverage_summary(stratified_metrics, suffix)
      ggsave(file.path(output_dir, paste0("coverage_", causal_count, "_causal.png")),
             p1_strat, width = 6, height = 10, dpi = 300)
      p2_strat <- plot_power_summary(stratified_metrics, suffix)
      ggsave(file.path(output_dir, paste0("power_", causal_count, "_causal.png")),
             p2_strat, width = 6, height = 10, dpi = 300)
      # For individual causal counts, create simpler bar plots without faceting
      plot_data_simple <- stratified_metrics %>%
        mutate(method_ld = paste(method, ld_size, sep = " + "))
      # AUPRC plot
      p2_simple <- ggplot(plot_data_simple, aes(x = method_ld, y = auprc, fill = method)) +
        geom_col(alpha = 0.8) +
        scale_fill_viridis_d(name = "Method") +
        labs(title = paste0("True AUPRC", suffix), x = "Method + LD Matrix", y = "AUPRC") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ylim(0, 1)
      ggsave(file.path(output_dir, paste0("auprc_", causal_count, "_causal.png")),
             p2_simple, width = 12, height = 6, dpi = 300)
      # Credible set size plot
      p3_simple <- ggplot(plot_data_simple, aes(x = method_ld, y = mean_credible_set_size, fill = method)) +
        geom_col(alpha = 0.8) +
        scale_fill_viridis_d(name = "Method") +
        labs(title = paste0("Mean Credible Set Size", suffix),
             x = "Method + LD Matrix", y = "Mean Size") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(output_dir, paste0("credible_set_sizes_", causal_count, "_causal.png")),
             p3_simple, width = 12, height = 6, dpi = 300)
      # NEW: Execution time plot
      p_timing_simple <- plot_execution_times(stratified_metrics, suffix)
      if (!is.null(p_timing_simple)) {
        ggsave(file.path(output_dir, paste0("execution_times_", causal_count, "_causal.png")),
               p_timing_simple, width = 12, height = 6, dpi = 300)
      }
    }
  }
  cat("Analysis complete! Results saved to:", output_dir, "\n")
  return(all_metrics)
}

# Run the full analysis for random effects
out_dir = "C:/Users/mconery/Box/Mitchell_Work_Files/Meta-Analysis_Fine-Mapping/Figures/SAIGE_Simulations/sub-sample_LD.random_homo_3e-4"
homogeneity_results_file = "C:/Users/mconery/Box/Mitchell_Work_Files/Meta-Analysis_Fine-Mapping/Data/SAIGE_Simulations/sub-sample/sub-sample_LD.random_homo_3e-4.results_file.txt"
results <- run_finemapping_analysis_summary(homogeneity_results_file, output_dir = out_dir)

# Run the full analysis for fixed effects
out_dir = "C:/Users/mconery/Box/Mitchell_Work_Files/Meta-Analysis_Fine-Mapping/Figures/SAIGE_Simulations/sub-sample_LD.fixed_homo_3e-4"
homogeneity_results_file = "C:/Users/mconery/Box/Mitchell_Work_Files/Meta-Analysis_Fine-Mapping/Data/SAIGE_Simulations/sub-sample/sub-sample_LD.fixed_homo_3e-4.results_file.txt"
results <- run_finemapping_analysis_summary(homogeneity_results_file, output_dir = out_dir)