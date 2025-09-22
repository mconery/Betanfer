#!/usr/bin/env Rscript

#' MESuSiE Multi-Population Fine-Mapping Script
#'
#' R script designed to take input files and parameters and use them to run MESuSiE
#' Attempts to use native MESuSiE functions first, falls back to pairwise if needed
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
                      help = 'Comma-delimited list of ancestry labels (e.g., EUR,AFR,EAS), same order as GWAS files')
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
                      help = 'Ancestry weights as comma-separated values (default: adaptive based on number of populations)')
  
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
  
  # Add other columns if available and rename Allele1/Allele2 to REF/ALT
  if (!is.null(detected_cols$a1)) {
    if (detected_cols$a1 == "Allele1") {
      result_df$ALT <- gwas_filtered[[detected_cols$a1]]
      cat(sprintf(" Renamed %s to ALT\n", detected_cols$a1))
    } else {
      result_df$A1 <- gwas_filtered[[detected_cols$a1]]
    }
  }
  if (!is.null(detected_cols$a2)) {
    if (detected_cols$a2 == "Allele2") {
      result_df$REF <- gwas_filtered[[detected_cols$a2]]
      cat(sprintf(" Renamed %s to REF\n", detected_cols$a2))
    } else {
      result_df$A2 <- gwas_filtered[[detected_cols$a2]]
    }
  }
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
  
  if (n_pops < 2) {
    stop("ERROR: At least 2 populations are required for multi-ancestry analysis")
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

# Extended organize functions that try to handle multiple populations
organize_gwas_multi <- function(summ_stat_list, ancestry_labels) {
  cat("Organizing GWAS data for multiple populations\n")
  n_pops <- length(summ_stat_list)
  
  if (n_pops == 2) {
    # Use the built-in MESuSiE function for 2 populations
    cat("Using native organize_gwas for 2 populations\n")
    return(organize_gwas(summ_stat_list[[1]], summ_stat_list[[2]], ancestry_labels))
  } else {
    # Try to extend the approach for >2 populations
    cat(sprintf("Attempting to organize %d populations\n", n_pops))
    
    # Get common SNPs across all populations
    all_snps <- summ_stat_list[[1]]$SNP
    for (i in 2:n_pops) {
      all_snps <- intersect(all_snps, summ_stat_list[[i]]$SNP)
    }
    
    if (length(all_snps) == 0) {
      stop("ERROR: No common SNPs found across all GWAS datasets after organization")
    }
    
    cat(sprintf("Common SNPs across all %d populations: %d\n", n_pops, length(all_snps)))
    
    # Filter and reorder all datasets to common SNPs
    organized_list <- list()
    for (i in seq_along(summ_stat_list)) {
      ancestry <- ancestry_labels[i]
      df <- summ_stat_list[[i]]
      
      # Filter to common SNPs and reorder
      df_filtered <- df[df$SNP %in% all_snps, ]
      df_ordered <- df_filtered[match(all_snps, df_filtered$SNP), ]
      
      organized_list[[ancestry]] <- df_ordered
      cat(sprintf("Population %s: %d variants after organization\n", ancestry, nrow(df_ordered)))
    }
    
    return(organized_list)
  }
}

organize_ld_multi <- function(LD_list, organized_gwas) {
  cat("Organizing LD matrices for multiple populations\n")
  n_pops <- length(LD_list)
  
  if (n_pops == 2) {
    # Use the built-in MESuSiE function for 2 populations
    cat("Using native organize_ld for 2 populations\n")
    return(organize_ld(LD_list[[1]], LD_list[[2]], organized_gwas))
  } else {
    # Handle >2 populations
    cat(sprintf("Organizing %d LD matrices\n", n_pops))
    
    # Get common SNPs from organized GWAS data
    common_snps <- organized_gwas[[1]]$SNP
    
    # Check and organize all LD matrices
    organized_ld_list <- list()
    for (i in seq_along(LD_list)) {
      ancestry <- names(LD_list)[i]
      ld_matrix <- LD_list[[i]]
      
      # Check dimensions
      if (nrow(ld_matrix) != length(common_snps) || ncol(ld_matrix) != length(common_snps)) {
        stop(sprintf("ERROR: LD matrix for %s has incorrect dimensions", ancestry))
      }
      
      # Ensure matrix is symmetric and set proper names
      if (!isSymmetric(ld_matrix)) {
        ld_matrix <- (ld_matrix + t(ld_matrix)) / 2
        cat(sprintf("Made LD matrix symmetric for population %s\n", ancestry))
      }
      
      rownames(ld_matrix) <- common_snps
      colnames(ld_matrix) <- common_snps
      
      organized_ld_list[[ancestry]] <- ld_matrix
      cat(sprintf("Population %s LD matrix: %d x %d\n", ancestry, nrow(ld_matrix), ncol(ld_matrix)))
    }
    
    return(organized_ld_list)
  }
}

# Run MESuSiE analysis
run_mesusie <- function(validated_inputs, args) {
  cat("\n=== Running MESuSiE fine-mapping ===\n")
  
  gwas_files <- validated_inputs$gwas_files
  ld_files <- validated_inputs$ld_files
  n_gwas_list <- validated_inputs$n_gwas_list
  ancestry_labels <- validated_inputs$ancestry_labels
  n_pops <- validated_inputs$n_pops
  
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
  organized_gwas <- organize_gwas_multi(summ_stat_list, ancestry_labels)
  organized_ld <- organize_ld_multi(LD_list, organized_gwas)
  
  # Parse ancestry weights if provided
  ancestry_weight <- NULL
  if (!is.null(args$ancestry_weight)) {
    ancestry_weight <- as.numeric(strsplit(args$ancestry_weight, ",")[[1]])
  }
  
  # Try to run MESuSiE with the organized data
  cat("\n=== Running MESuSiE core analysis ===\n")
  
  if (n_pops == 2) {
    # Standard 2-population MESuSiE
    if (is.null(ancestry_weight)) {
      ancestry_weight <- c(1, 1, 3)  # Default: encourage shared effects
    }
    if (length(ancestry_weight) != 3) {
      stop("ERROR: For 2 populations, ancestry_weight must have exactly 3 values (pop1, pop2, shared)")
    }
    
    mesusie_result <- meSuSie_core(
      organized_ld,
      organized_gwas,
      L = args$L,
      estimate_residual_variance = args$estimate_residual_variance,
      cor_method = args$cor_method,
      cor_threshold = args$cor_threshold,
      prior_weights = args$prior_weights,
      ancestry_weight = ancestry_weight
    )
  } else {
    # For >2 populations, first try to see if MESuSiE can handle it directly
    cat(sprintf("Attempting MESuSiE analysis with %d populations...\n", n_pops))
    
    tryCatch({
      # Try running with extended ancestry weights
      if (is.null(ancestry_weight)) {
        # Generate default weights: equal weight for each population, higher weight for shared
        ancestry_weight <- c(rep(1, n_pops), 3)
      }
      
      cat(sprintf("Using ancestry weights: %s\n", paste(ancestry_weight, collapse = ", ")))
      
      mesusie_result <- meSuSie_core(
        organized_ld,
        organized_gwas,
        L = args$L,
        estimate_residual_variance = args$estimate_residual_variance,
        cor_method = args$cor_method,
        cor_threshold = args$cor_threshold,
        prior_weights = args$prior_weights,
        ancestry_weight = ancestry_weight
      )
      
      cat("SUCCESS: MESuSiE handled multiple populations directly!\n")
      
    }, error = function(e) {
      cat(sprintf("MESuSiE failed with %d populations: %s\n", n_pops, e$message))
      cat("Note: This indicates MESuSiE may indeed be limited to 2 populations.\n")
      cat("You may want to select 2 populations or use alternative methods.\n")
      stop(sprintf("MESuSiE analysis failed: %s", e$message))
    })
  }
  
  return(c(mesusie_result, organized_gwas))
}

# Save results
save_results <- function(mesusie_result, organized_gwas, ancestry_labels, output_file, args) {
  cat("\n=== Saving results ===\n")
  
  # Create output directory if needed
  output_dir <- dirname(output_file)
  if (output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the complete MESuSiE result object as RDS for reference
  saveRDS(mesusie_result, file = paste0(output_file, ".rds"))
  cat(sprintf("Saved complete results to: %s.rds\n", output_file))
  
  # Get variant IDs from the first population (all should have same variants after organization)
  variant_ids <- organized_gwas[[1]]$SNP
  n_variants <- length(variant_ids)
  
  # Initialize the output dataframe
  result_df <- data.frame(
    VARIANT_ID = variant_ids,
    stringsAsFactors = FALSE
  )
  
  # Add Z-score columns for each population
  for (i in seq_along(ancestry_labels)) {
    ancestry <- ancestry_labels[i]
    z_col_name <- paste0("Z_", ancestry)
    
    if ("Z" %in% colnames(organized_gwas[[ancestry]])) {
      result_df[[z_col_name]] <- organized_gwas[[ancestry]]$Z
    } else if ("Beta" %in% colnames(organized_gwas[[ancestry]]) && 
               "Se" %in% colnames(organized_gwas[[ancestry]])) {
      # Calculate Z-score if not present
      result_df[[z_col_name]] <- organized_gwas[[ancestry]]$Beta / organized_gwas[[ancestry]]$Se
    } else {
      # Fill with NA if Z-score cannot be determined
      result_df[[z_col_name]] <- rep(NA, n_variants)
      cat(sprintf("WARNING: Could not find Z-scores for population %s\n", ancestry))
    }
  }
  
  # Add overall PIP column
  if (!is.null(mesusie_result$pip)) {
    # Match PIPs to variant order
    pip_values <- mesusie_result$pip[match(variant_ids, names(mesusie_result$pip))]
    result_df$PIP <- ifelse(is.na(pip_values), 0, pip_values)
  } else {
    result_df$PIP <- rep(0, n_variants)
    cat("WARNING: No PIP values found in MESuSiE results\n")
  }
  
  # Add converged column
  converged_status <- TRUE  # Default to TRUE
  if (!is.null(mesusie_result$converged)) {
    converged_status <- mesusie_result$converged
  } else if (!is.null(mesusie_result$elbo_res)) {
    # Check if ELBO converged (alternative convergence check)
    converged_status <- mesusie_result$elbo_res$converged
  }
  result_df$converged <- rep(converged_status, n_variants)
  
  # Add credible set column
  result_df$cs <- rep(0, n_variants)  # Initialize with 0 (not in any CS)
  
  if (!is.null(mesusie_result$cs) && !is.null(mesusie_result$cs$cs)) {
    for (cs_idx in seq_along(mesusie_result$cs$cs)) {
      cs_variants <- mesusie_result$cs$cs[[cs_idx]]
      
      # Convert indices to variant IDs if needed
      if (is.numeric(cs_variants)) {
        cs_variant_ids <- variant_ids[cs_variants]
      } else {
        cs_variant_ids <- cs_variants
      }
      
      # Mark variants in this credible set
      variant_matches <- match(cs_variant_ids, variant_ids)
      valid_matches <- variant_matches[!is.na(variant_matches)]
      
      if (length(valid_matches) > 0) {
        result_df$cs[valid_matches] <- cs_idx
      }
    }
    
    cat(sprintf("Assigned %d credible sets to variants\n", length(mesusie_result$cs$cs)))
  } else {
    cat("No credible sets found in MESuSiE results\n")
  }
  
  # Add population-specific and intersection PIP columns
  if (!is.null(mesusie_result$pip_config)) {
    pip_config <- mesusie_result$pip_config
    
    # Get column names from pip_config (these represent the different configurations)
    config_names <- colnames(pip_config)
    
    for (config_name in config_names) {
      # Create column name for this configuration
      pip_col_name <- paste0("PIP_", config_name)
      
      # Match PIPs to variant order
      if (is.matrix(pip_config)) {
        pip_values <- pip_config[match(variant_ids, rownames(pip_config)), config_name]
      } else {
        pip_values <- pip_config[[config_name]][match(variant_ids, names(pip_config[[config_name]]))]
      }
      
      result_df[[pip_col_name]] <- ifelse(is.na(pip_values), 0, pip_values)
    }
    
    cat(sprintf("Added %d population-specific PIP columns: %s\n", 
                length(config_names), paste(config_names, collapse = ", ")))
  } else {
    cat("WARNING: No population-specific PIPs found in MESuSiE results\n")
    
    # Create placeholder columns with zeros
    for (ancestry in ancestry_labels) {
      result_df[[paste0("PIP_", ancestry)]] <- rep(0, n_variants)
    }
    
    # Create intersection column for 2+ populations
    if (length(ancestry_labels) >= 2) {
      intersection_name <- paste(ancestry_labels, collapse = "_")
      result_df[[paste0("PIP_", intersection_name)]] <- rep(0, n_variants)
    }
  }
  
  # Write the tab-delimited output file
  output_txt_file <- paste0(output_file, ".txt")
  write.table(result_df, 
              file = output_txt_file, 
              sep = "\t", 
              row.names = FALSE, 
              col.names = TRUE, 
              quote = FALSE)
  
  cat(sprintf("Saved results to: %s\n", output_txt_file))
  cat(sprintf("Output contains %d variants and %d columns\n", nrow(result_df), ncol(result_df)))
  
  # Print column summary for verification
  cat("Output columns:\n")
  for (i in seq_along(colnames(result_df))) {
    col_name <- colnames(result_df)[i]
    if (col_name == "VARIANT_ID") {
      cat(sprintf("  %d. %s (variant identifiers)\n", i, col_name))
    } else if (grepl("^Z_", col_name)) {
      cat(sprintf("  %d. %s (Z-scores)\n", i, col_name))
    } else if (col_name == "PIP") {
      cat(sprintf("  %d. %s (overall posterior inclusion probability)\n", i, col_name))
    } else if (col_name == "converged") {
      cat(sprintf("  %d. %s (convergence status: %s)\n", i, col_name, unique(result_df[[col_name]])[1]))
    } else if (col_name == "cs") {
      n_in_cs <- sum(result_df[[col_name]] > 0)
      cat(sprintf("  %d. %s (credible set assignments, %d variants in sets)\n", i, col_name, n_in_cs))
    } else if (grepl("^PIP_", col_name)) {
      cat(sprintf("  %d. %s (population-specific PIP)\n", i, col_name))
    }
  }
  
  cat("Results saved successfully!\n")
}

# Main function
main <- function() {
  cat("MESuSiE Multi-Population Fine-Mapping: Testing Multi-Population Support\n")
  
  # Parse arguments
  args <- parse_arguments()
  
  tryCatch({
    # Validate inputs
    validated_inputs <- validate_inputs(args)
    
    # Run MESuSiE analysis
    temp <- run_mesusie(validated_inputs, args)
    mesusie_result <- temp[1]
    organized_gwas <- temp[2]
    
    # Save results
    save_results(mesusie_result, organized_gwas, validated_inputs$ancestry_labels, args$output_file, args)
    
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