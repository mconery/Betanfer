#!/usr/bin/env Rscript

#' MESuSiE Multi-Population Fine-Mapping Script
#'
#' R script designed to take input files and parameters and use them to run MESuSiE
#' Equivalent functionality to the Python SuSiEx script
#'
#' @author Generated based on SuSiEx Python script structure

library(argparse)
library(data.table)
library(dplyr)
library(MESuSiE)

# Parse command line arguments
parse_arguments <- function() {
  parser <- ArgumentParser(
    description = 'Multi-population fine-mapping using MESuSiE (auto chr/bp from LD matrix .vars)'
  )
  
  parser$add_argument('--gwas_files', '-g', required = TRUE,
                      help = 'Comma-delimited list of GWAS summary statistics files')
  parser$add_argument('--ld_files', '-l', required = TRUE,
                      help = 'Comma-delimited list of LD matrix file prefixes (each must have a corresponding .vars file)')
  parser$add_argument('--output_file', '-o', required = TRUE,
                      help = 'Full path for final output file')
  parser$add_argument('--n_gwas', '-n', required = TRUE,
                      help = 'Comma-delimited list of GWAS sample sizes, same order as GWAS files')
  parser$add_argument('--ancestry_labels', '-a', required = TRUE,
                      help = 'Comma-delimited list of ancestry labels (e.g., EUR,AFR), same order as GWAS files')
  parser$add_argument('--maf', type = 'double', default = 0.005,
                      help = 'Minor allele frequency threshold (default: 0.005)')
  parser$add_argument('--L', type = 'integer', default = 10,
                      help = 'Maximum number of causal effects (default: 10)')
  parser$add_argument('--coverage', type = 'double', default = 0.95,
                      help = 'Coverage for credible sets (default: 0.95)')
  parser$add_argument('--min_purity', type = 'double', default = 0.5,
                      help = 'Minimum purity for credible sets (default: 0.5)')
  parser$add_argument('--estimate_residual_variance', action = 'store_true',
                      help = 'Estimate residual variance (default: FALSE)')
  parser$add_argument('--cor_method', default = 'min_abs_corr',
                      help = 'Correlation method for credible sets (default: min_abs_corr)')
  parser$add_argument('--cor_threshold', type = 'double', default = 0.5,
                      help = 'Correlation threshold for credible sets (default: 0.5)')
  parser$add_argument('--prior_weights', default = NULL,
                      help = 'Prior weights for variants (default: uniform)')
  parser$add_argument('--ancestry_weight', default = NULL,
                      help = 'Ancestry weights as comma-separated values (default: c(1,1,3))')
  
  return(parser$parse_args())
}

# Read variant IDs from .vars file
read_vars_file <- function(vars_file) {
  cat(sprintf("Reading variants from %s\n", vars_file))
  
  if (!file.exists(vars_file)) {
    stop(sprintf("ERROR: .vars file not found: %s", vars_file))
  }
  
  variant_ids <- readLines(vars_file)
  variant_ids <- variant_ids[variant_ids != ""]
  
  # Extract chromosome and positions
  chrom <- NULL
  positions <- c()
  
  for (vid in variant_ids) {
    # Match CHR:POS:REF:ALT format
    match_result <- regmatches(vid, regexpr("^(\\w+):(\\d+):[^:]+:[^:]+$", vid, perl = TRUE))
    if (length(match_result) > 0) {
      parts <- strsplit(match_result, ":")[[1]]
      c <- parts[1]
      pos <- as.numeric(parts[2])
      
      if (is.null(chrom)) chrom <- c
      positions <- c(positions, pos)
    } else {
      cat(sprintf("WARNING: Variant ID does not match CHR:POS:REF:ALT: %s\n", vid))
    }
  }
  
  if (length(positions) == 0 || is.null(chrom)) {
    stop(sprintf("ERROR: Could not parse chromosome and positions from %s", vars_file))
  }
  
  return(list(
    chrom = chrom,
    min_bp = min(positions),
    max_bp = max(positions),
    variant_ids = variant_ids
  ))
}

# Find intersection of variants across all populations
find_common_variants <- function(ld_files, gwas_files, ancestry_labels) {
  cat("\n=== Finding common variants across all populations ===\n")
  
  # Get variants from all LD matrices
  all_ld_variants <- list()
  for (i in seq_along(ld_files)) {
    vars_file <- paste0(ld_files[i], ".vars")
    vars_info <- read_vars_file(vars_file)
    all_ld_variants[[ancestry_labels[i]]] <- vars_info$variant_ids
    cat(sprintf("Population %s LD matrix: %d variants\n", ancestry_labels[i], length(vars_info$variant_ids)))
  }
  
  # Get variants from all GWAS files
  all_gwas_variants <- list()
  for (i in seq_along(gwas_files)) {
    # Detect columns
    detected_cols <- detect_columns(gwas_files[i], verbose = FALSE)
    
    # Read GWAS file
    gwas_df <- fread(gwas_files[i])
    snp_col <- detected_cols$snp
    if (is.null(snp_col) || !snp_col %in% colnames(gwas_df)) {
      stop(sprintf("ERROR: Could not find SNP column in GWAS file: %s", gwas_files[i]))
    }
    
    gwas_variants <- gwas_df[[snp_col]]
    all_gwas_variants[[ancestry_labels[i]]] <- gwas_variants
    cat(sprintf("Population %s GWAS: %d variants\n", ancestry_labels[i], length(gwas_variants)))
  }
  
  # Find intersection of all LD variants
  common_ld_variants <- all_ld_variants[[1]]
  for (i in 2:length(all_ld_variants)) {
    common_ld_variants <- intersect(common_ld_variants, all_ld_variants[[i]])
  }
  cat(sprintf("Common LD variants across all populations: %d\n", length(common_ld_variants)))
  
  # Find intersection of all GWAS variants  
  common_gwas_variants <- all_gwas_variants[[1]]
  for (i in 2:length(all_gwas_variants)) {
    common_gwas_variants <- intersect(common_gwas_variants, all_gwas_variants[[i]])
  }
  cat(sprintf("Common GWAS variants across all populations: %d\n", length(common_gwas_variants)))
  
  # Find final intersection of LD and GWAS variants
  final_common_variants <- intersect(common_ld_variants, common_gwas_variants)
  cat(sprintf("Final common variants (LD ∩ GWAS): %d\n", length(final_common_variants)))
  
  if (length(final_common_variants) == 0) {
    stop("ERROR: No variants are common across all LD matrices and GWAS files")
  }
  
  return(final_common_variants)
}

# Detect columns in GWAS file
detect_columns <- function(gwas_file, verbose = TRUE) {
  if (verbose) {
    cat(sprintf(" Detecting columns in %s\n", basename(gwas_file)))
  }
  
  # Try reading the header
  tryCatch({
    if (grepl("\\.gz$", gwas_file)) {
      df <- fread(gwas_file, nrows = 0)
    } else {
      df <- fread(gwas_file, nrows = 0)
    }
    headers <- colnames(df)
    
    if (verbose) {
      cat(sprintf(" Headers: %s\n", paste(headers, collapse = ", ")))
    }
    
    # Column patterns for detection
    column_patterns <- list(
      chr = c('CHR', 'chr', 'chromosome', '#CHROM', 'chrom'),
      snp = c('SNP', 'snp', 'rsid', 'RSID', 'ID', 'id', 'variant_id', 'VARIANT_ID', 'MarkerName', 'rs_number', 'MarkerID'),
      bp = c('BP', 'bp', 'pos', 'POS', 'position', 'POSITION', 'base_pair_location'),
      a1 = c('A1', 'a1', 'effect_allele', 'EFFECT_ALLELE', 'ALT', 'alt', 'ea', 'Allele1'),
      a2 = c('A2', 'a2', 'other_allele', 'OTHER_ALLELE', 'REF', 'ref', 'nea', 'Allele2'),
      beta = c('BETA', 'beta', 'Beta', 'OR', 'or', 'effect', 'EFFECT', 'b'),
      se = c('SE', 'se', 'standard_error', 'STANDARD_ERROR', 'stderr', 'STDERR'),
      pval = c('P', 'p', 'pval', 'PVAL', 'p_value', 'P_VALUE', 'pvalue', 'PVALUE', 'p.value'),
      n = c('N', 'n', 'sample_size', 'SAMPLE_SIZE', 'nobs', 'NOBS')
    )
    
    detected_cols <- list()
    
    # Detect each column type
    for (col_type in names(column_patterns)) {
      col_name <- NULL
      patterns <- column_patterns[[col_type]]
      
      for (header in headers) {
        if (header %in% patterns) {
          col_name <- header
          break
        }
      }
      
      if (is.null(col_name)) {
        if (verbose) {
          cat(sprintf(" WARNING: Could not detect %s column\n", col_type))
        }
      }
      
      detected_cols[[col_type]] <- col_name
    }
    
    return(detected_cols)
    
  }, error = function(e) {
    stop(sprintf("ERROR: Failed to read GWAS file %s: %s", gwas_file, e$message))
  })
}

# Load and process GWAS data (now filtered to common variants)
load_gwas_data <- function(gwas_file, detected_cols, common_variants, n_gwas) {
  cat(sprintf(" Loading GWAS data from %s\n", basename(gwas_file)))
  
  # Read the full GWAS file
  gwas_df <- fread(gwas_file)
  
  # Filter to common variants
  snp_col <- detected_cols$snp
  if (is.null(snp_col) || !snp_col %in% colnames(gwas_df)) {
    stop("ERROR: Could not find SNP column in GWAS file")
  }
  
  gwas_filtered <- gwas_df[get(snp_col) %in% common_variants]
  
  if (nrow(gwas_filtered) == 0) {
    stop("ERROR: No overlapping variants found between GWAS and common variant set")
  }
  
  cat(sprintf(" Found %d overlapping variants\n", nrow(gwas_filtered)))
  
  # Order by common_variants to ensure consistent ordering
  gwas_filtered <- gwas_filtered[match(common_variants, gwas_filtered[[snp_col]])]
  gwas_filtered <- gwas_filtered[!is.na(gwas_filtered[[snp_col]])]
  
  # Create standard format data frame
  result_df <- data.frame(
    SNP = gwas_filtered[[detected_cols$snp]],
    stringsAsFactors = FALSE
  )
  
  # Add other columns if available
  if (!is.null(detected_cols$beta)) {
    result_df$Beta <- gwas_filtered[[detected_cols$beta]]
  }
  if (!is.null(detected_cols$se)) {
    result_df$Se <- gwas_filtered[[detected_cols$se]]
  }
  if (!is.null(detected_cols$pval)) {
    result_df$P <- gwas_filtered[[detected_cols$pval]]
  }
  
  # Calculate Z-scores if Beta and Se are available
  if (!is.null(detected_cols$beta) && !is.null(detected_cols$se)) {
    result_df$Z <- result_df$Beta / result_df$Se
  }
  
  # Add sample size
  if (!is.null(detected_cols$n)) {
    result_df$N <- gwas_filtered[[detected_cols$n]]
  } else {
    result_df$N <- n_gwas
  }
  
  return(result_df)
}

# Load LD matrix (now pruned to common variants)
load_ld_matrix <- function(ld_file, vars_file, common_variants) {
  cat(sprintf(" Loading LD matrix from %s\n", basename(ld_file)))
  
  if (!file.exists(ld_file)) {
    stop(sprintf("ERROR: LD matrix file not found: %s", ld_file))
  }
  
  # Read the full LD matrix and variant IDs
  ld_matrix <- tryCatch({
    as.matrix(fread(ld_file, header = FALSE))
  }, error = function(e) {
    stop(sprintf("ERROR: Could not read LD matrix from %s: %s", ld_file, e$message))
  })
  
  # Read variant IDs from .vars file
  all_variant_ids <- readLines(vars_file)
  all_variant_ids <- all_variant_ids[all_variant_ids != ""]
  
  # Check original dimensions
  n_all_variants <- length(all_variant_ids)
  if (nrow(ld_matrix) != n_all_variants || ncol(ld_matrix) != n_all_variants) {
    stop(sprintf("ERROR: LD matrix dimensions (%d x %d) do not match .vars file variants (%d)", 
                 nrow(ld_matrix), ncol(ld_matrix), n_all_variants))
  }
  
  # Find indices of common variants in the full matrix
  common_indices <- which(all_variant_ids %in% common_variants)
  if (length(common_indices) == 0) {
    stop("ERROR: No common variants found in LD matrix")
  }
  
  cat(sprintf(" Pruning LD matrix from %d to %d variants\n", n_all_variants, length(common_indices)))
  
  # Subset the LD matrix to common variants
  ld_matrix_pruned <- ld_matrix[common_indices, common_indices, drop = FALSE]
  common_variant_ids <- all_variant_ids[common_indices]
  
  # Set row and column names
  rownames(ld_matrix_pruned) <- common_variant_ids
  colnames(ld_matrix_pruned) <- common_variant_ids
  
  # Reorder to match common_variants order
  variant_order <- match(common_variants, common_variant_ids)
  variant_order <- variant_order[!is.na(variant_order)]
  ld_matrix_final <- ld_matrix_pruned[variant_order, variant_order, drop = FALSE]
  
  return(ld_matrix_final)
}

# Validate inputs
validate_inputs <- function(args) {
  cat("=== Validating inputs ===\n")
  
  gwas_files <- strsplit(args$gwas_files, ",")[[1]]
  ld_files <- strsplit(args$ld_files, ",")[[1]]
  n_gwas_list <- as.numeric(strsplit(args$n_gwas, ",")[[1]])
  ancestry_labels <- strsplit(args$ancestry_labels, ",")[[1]]
  
  n_pops <- length(gwas_files)
  
  if (length(ld_files) != n_pops || length(n_gwas_list) != n_pops || length(ancestry_labels) != n_pops) {
    stop("ERROR: Number of populations must match across all inputs")
  }
  
  # Check GWAS files exist
  for (i in seq_along(gwas_files)) {
    if (!file.exists(gwas_files[i])) {
      stop(sprintf("ERROR: GWAS file not found: %s", gwas_files[i]))
    }
    cat(sprintf("Population %d (%s) GWAS file: %s ✓\n", i, ancestry_labels[i], gwas_files[i]))
  }
  
  # Check LD files exist
  for (i in seq_along(ld_files)) {
    vars_file <- paste0(ld_files[i], ".vars")
    if (!file.exists(ld_files[i]) || !file.exists(vars_file)) {
      stop(sprintf("ERROR: LD matrix (%s) or .vars file (%s) not found", ld_files[i], vars_file))
    }
    cat(sprintf("Population %d (%s) LD matrix: %s, .vars ✓\n", i, ancestry_labels[i], ld_files[i]))
  }
  
  # Create output directory if needed
  output_dir <- dirname(args$output_file)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat(sprintf("Validated input files for %d populations.\n", n_pops))
  
  return(list(
    gwas_files = gwas_files,
    ld_files = ld_files,
    n_gwas_list = n_gwas_list,
    ancestry_labels = ancestry_labels,
    n_pops = n_pops
  ))
}

# Run MESuSiE analysis
run_mesusie <- function(validated_inputs, args) {
  cat("\n=== Running MESuSiE fine-mapping ===\n")
  
  gwas_files <- validated_inputs$gwas_files
  ld_files <- validated_inputs$ld_files
  n_gwas_list <- validated_inputs$n_gwas_list
  ancestry_labels <- validated_inputs$ancestry_labels
  
  # Find common variants across all populations
  common_variants <- find_common_variants(ld_files, gwas_files, ancestry_labels)
  
  # Get chromosome and BP info from first .vars file
  vars_info <- read_vars_file(paste0(ld_files[1], ".vars"))
  cat(sprintf("Chromosome: %s, BP range: %d-%d\n", vars_info$chrom, vars_info$min_bp, vars_info$max_bp))
  
  # Process each population with common variants
  summ_stat_list <- list()
  LD_list <- list()
  
  for (i in seq_along(gwas_files)) {
    ancestry <- ancestry_labels[i]
    cat(sprintf("\nProcessing population %d (%s):\n", i, ancestry))
    
    # Detect columns
    detected_cols <- detect_columns(gwas_files[i])
    
    # Load GWAS data (filtered to common variants)
    gwas_data <- load_gwas_data(gwas_files[i], detected_cols, common_variants, n_gwas_list[i])
    
    # Load LD matrix (pruned to common variants)
    vars_file <- paste0(ld_files[i], ".vars")
    ld_matrix <- load_ld_matrix(ld_files[i], vars_file, common_variants)
    
    # Store in lists with ancestry labels
    summ_stat_list[[ancestry]] <- gwas_data
    LD_list[[ancestry]] <- ld_matrix
  }
  
  # Organize data for MESuSiE
  cat("\n=== Organizing data for MESuSiE ===\n")
  organized_gwas <- organize_gwas_list(summ_stat_list, ancestry_labels)
  organized_ld <- organize_ld_list(LD_list, organized_gwas)
  
  # Parse ancestry weights if provided
  ancestry_weight <- NULL
  if (!is.null(args$ancestry_weight)) {
    ancestry_weight <- as.numeric(strsplit(args$ancestry_weight, ",")[[1]])
    if (length(ancestry_weight) != 3) {
      stop("ERROR: ancestry_weight must have exactly 3 values (anc1, anc2, shared)")
    }
  }
  
  # Run MESuSiE
  cat("\n=== Running MESuSiE core analysis ===\n")
  mesusie_result <- meSuSie_core(
    LD_list = organized_ld,
    summ_stat_list = organized_gwas,
    L = args$L,
    coverage = args$coverage,
    estimate_residual_variance = args$estimate_residual_variance,
    cor_method = args$cor_method,
    cor_threshold = args$cor_threshold,
    prior_weights = args$prior_weights,
    ancestry_weight = ancestry_weight
  )
  
  return(mesusie_result)
}

# Organize GWAS data for MESuSiE
organize_gwas_list <- function(summ_stat_list, ancestry_labels) {
  # Use MESuSiE's organize_gwas function for two populations
  if (length(summ_stat_list) == 2) {
    organized <- organize_gwas(
      summ_stat_list[[1]], 
      summ_stat_list[[2]], 
      ancestry_labels
    )
  } else {
    # For more than 2 populations, we'd need to extend this
    stop("Currently only supports 2 populations. For more populations, extend the organize_gwas functionality.")
  }
  
  return(organized)
}

# Organize LD data for MESuSiE
organize_ld_list <- function(LD_list, organized_gwas) {
  # Use MESuSiE's organize_ld function for two populations
  if (length(LD_list) == 2) {
    organized <- organize_ld(
      LD_list[[1]], 
      LD_list[[2]], 
      organized_gwas
    )
  } else {
    # For more than 2 populations, we'd need to extend this
    stop("Currently only supports 2 populations. For more populations, extend the organize_ld functionality.")
  }
  
  return(organized)
}

# Save results
save_results <- function(mesusie_result, output_file, args) {
  cat("\n=== Saving results ===\n")
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the complete MESuSiE result object
  saveRDS(mesusie_result, file = paste0(output_file, ".rds"))
  cat(sprintf("Saved complete results to: %s.rds\n", output_file))
  
  # Extract and save summary results
  summary_results <- list(
    pip = mesusie_result$pip,
    pip_config = mesusie_result$pip_config,
    credible_sets = mesusie_result$cs,
    parameters = list(
      L = args$L,
      coverage = args$coverage,
      cor_method = args$cor_method,
      cor_threshold = args$cor_threshold
    )
  )
  
  # Save summary as CSV
  if (!is.null(mesusie_result$pip)) {
    pip_df <- data.frame(
      SNP = names(mesusie_result$pip),
      PIP = mesusie_result$pip,
      stringsAsFactors = FALSE
    )
    
    # Add ancestry-specific PIPs if available
    if (!is.null(mesusie_result$pip_config)) {
      pip_config_df <- as.data.frame(mesusie_result$pip_config)
      pip_df <- cbind(pip_df, pip_config_df)
    }
    
    write.csv(pip_df, file = paste0(output_file, "_pip.csv"), row.names = FALSE)
    cat(sprintf("Saved PIP results to: %s_pip.csv\n", output_file))
  }
  
  # Save credible sets information
  if (!is.null(mesusie_result$cs)) {
    cs_info <- data.frame(
      Effect = names(mesusie_result$cs$cs),
      Category = mesusie_result$cs$cs_category,
      Size = sapply(mesusie_result$cs$cs, length),
      Coverage = mesusie_result$cs$coverage,
      Min_Abs_Corr = mesusie_result$cs$purity[, "min.abs.corr"],
      stringsAsFactors = FALSE
    )
    
    write.csv(cs_info, file = paste0(output_file, "_credible_sets.csv"), row.names = FALSE)
    cat(sprintf("Saved credible sets info to: %s_credible_sets.csv\n", output_file))
  }
  
  cat("Results saved successfully!\n")
}

# Main function
main <- function() {
  cat("MESuSiE Multi-Population Fine-Mapping: Auto chr/bp from .vars\n")
  
  # Parse arguments
  args <- parse_arguments()
  
  tryCatch({
    # Validate inputs
    validated_inputs <- validate_inputs(args)
    
    # Run MESuSiE analysis
    mesusie_result <- run_mesusie(validated_inputs, args)
    
    # Save results
    save_results(mesusie_result, args$output_file, args)
    
    cat("\n=== Analysis completed successfully! ===\n")
    
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    quit(status = 1)
  })
}

# Run if called as script
if (!interactive()) {
  main()
}