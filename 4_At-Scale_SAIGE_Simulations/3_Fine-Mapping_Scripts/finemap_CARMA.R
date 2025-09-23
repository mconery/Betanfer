################################################################################

#finemap_CARMA.R

#This script is designed to execute CARMA fine-mapping given a fixed region

#The code should be run in the carma conda environment

################################################################################

# 0) Load Needed Libraries ====

library(CARMA)
library(magrittr)
library(dplyr)
library(stringr)
library(R.utils)

# 1) Read in Needed Arguments from Command ====

#Need to read in the following (matching SuSiE argument order):
# 1) summary stats file (not pruned - will be filtered using variant file)
# 2) ld matrix file in .vcor1 or similar format
# 3) output file location
# 4) column name for variant ID in summary stats [OPTIONAL; default is "markerid"]
# 5) column name for effect size (beta) in summary stats [OPTIONAL; default is "beta"]
# 6) column name for standard error in summary stats [OPTIONAL; default is "se"]
# 7) column name for sample size in summary stats [OPTIONAL; default is "n"]
# 8) Toggle for Checking Outliers [OPTIONAL; default is TRUE]
# 9) Lambda for Poisson Regularization [OPTIONAL; default is 1]
# 10) Rho for credible set confidence [OPTIONAL; default is 0.95]
# 11) random seed [OPTIONAL; default is 5]

args <- commandArgs(trailingOnly = TRUE);

if (length(args) >= 3) {
  # Required parameters (matching SuSiE)
  summary_stats_file <- args[1]
  ld_file <- args[2]
  out_file <- args[3]
  # Optional parameters matching SuSiE first
  markerid_col <- if (length(args) >= 4) args[4] else "markerid"
  beta_col <- if (length(args) >= 5) args[5] else "beta"
  se_col <- if (length(args) >= 6) args[6] else "se"
  n_col <- if (length(args) >= 7) args[7] else "n"
  # CARMA-specific optional parameters
  outlier_check_toggle <- if (length(args) >= 8) args[8] else "TRUE"
  lambda <- if (length(args) >= 9) as.numeric(args[9]) else 1
  rho <- if (length(args) >= 10) as.numeric(args[10]) else 0.95
  random_seed <- if (length(args) >= 11) as.integer(args[11]) else 5
  print("Necessary Parameters Input")
} else {
  print("Usage: Rscript finemap_CARMA.R summary_stats_file ld_file out_file [markerid_col] [beta_col] [se_col] [n_col] [outlier_check_toggle] [lambda] [rho] [random_seed]");
  print("Example: Rscript finemap_CARMA.R sumstats.txt matrix.vcor1 output.txt");
  quit(save="no");
}

# 2) Check Input Arguments ====

#Set error indicator

error_ind = FALSE

#Make file_prefix

file_prefix <- str_replace(basename(summary_stats_file),"_out_snps_sim_gwas.tsv", "")

#Check that files exist

if (!(file.exists(summary_stats_file))) {
  
  error_ind = TRUE
  
  print(paste0("ERROR: Invalid summary stats file given for ", file_prefix))
  
}

if (!(file.exists(ld_file))) {
  
  error_ind = TRUE
  
  print(paste0("ERROR: Invalid ld matrix given for ", file_prefix))
  
}

if (!(dir.exists(dirname(out_file)))) {
  
  error_ind = TRUE
  
  print(paste0("ERROR: Invalid directory given for ", file_prefix, " output file"))
  
}

#Check outlier toggle

if (typeof(outlier_check_toggle) == "character") {
  
  first_letter = tolower(substr(outlier_check_toggle, 1,1))
  
  if(first_letter == "t" || first_letter == "y"){
    
    outlier_check_toggle = TRUE
    
  }
  
  else if(first_letter == "f" || first_letter == "n") {
    
    outlier_check_toggle = FALSE
    
  }
  
  else{
    
    print(paste0("ERROR: Unrecognized outlier check option for ", file_prefix))
    
    error_ind = TRUE
    
  }
  
} else if (typeof(outlier_check_toggle) != "logical"){
  
  print(paste0("ERROR: Unrecognized outlier check option for ", file_prefix))
  
  error_ind = TRUE
  
}

#Recast lambda value

if (typeof(lambda) != "numeric") {
  
  lambda = as.numeric(lambda)
  
}

#Recast rho value

if (typeof(rho) != "numeric") {
  
  rho = as.numeric(rho)
  
}

#Quit if error indicator set to true

if (error_ind) {
  
  quit(save="no");
  
}

#Set working directory to output directory

setwd(dirname(out_file))

# 3) Prep for executing CARMA ====

### Read in input files and extract CARMA inputs ###

# Read in variant list from .vars file

vars_file <- paste0(ld_file, ".vars")

if (!file.exists(vars_file)) {
  
  print(paste0("ERROR: Variant file ", vars_file, " not found."))
  
  quit(save="no");
  
}

var_R <- read.csv(vars_file, sep = "\t", header = FALSE)[,1]

# Read in LD matrix

R <- as.matrix(read.csv(ld_file, sep = "\t", header = FALSE))

# Set row and column names to variants

rownames(R) <- var_R

colnames(R) <- var_R

# Read in summary stats

summary_stats_raw <- read.table(summary_stats_file, header = TRUE)

# Check if required columns exist

required_cols <- c(markerid_col, beta_col, se_col, n_col)

missing_cols <- required_cols[!required_cols %in% colnames(summary_stats_raw)]

if (length(missing_cols) > 0) {
  
  print(paste0("ERROR: Missing columns in summary stats file: ", paste(missing_cols, collapse=", ")))
  
  quit(save="no");
  
}

# Filter summary stats for variants in the locus

summary_stats_raw <- summary_stats_raw[summary_stats_raw[[markerid_col]] %in% var_R,]

if (nrow(summary_stats_raw) == 0) {
  
  print("ERROR: No variants from the LD matrix found in summary statistics.")
  
  quit(save="no");
  
}

# Reorder summary stats to match variant order in LD matrix

summary_stats_raw <- summary_stats_raw[match(var_R, summary_stats_raw[[markerid_col]]),]

# Check for any missing variants

missing_vars <- is.na(summary_stats_raw[[markerid_col]])

if (any(missing_vars)) {
  
  print(paste0("WARNING: ", sum(missing_vars), " variants from LD matrix not found in summary stats"))
  
  # Remove missing variants from both LD matrix and summary stats
  
  keep_idx <- !missing_vars
  
  summary_stats_raw <- summary_stats_raw[keep_idx,]
  
  R <- R[keep_idx, keep_idx]
  
  var_R <- var_R[keep_idx]
  
}

# Extract required columns

n <- mean(summary_stats_raw[[n_col]], na.rm = TRUE)

bhat <- summary_stats_raw[[beta_col]]

shat <- summary_stats_raw[[se_col]]

#Pack everything into lists for CARMA

z.list = list(bhat/shat)

ld.list = list(as.matrix(R))

lambda.list = list(lambda)

# 4) Create CARMA Execution function ====

#Check for the outlier toggle and define a function for calling CARMA accordingly

if (outlier_check_toggle == FALSE){
  
  call_carma <- function(z.list.inp, ld.list.inp, lambda.list.inp, rho.index.inp){
    
    return(suppressWarnings(CARMA(z.list.inp, ld.list.inp, lambda.list=lambda.list.inp, outlier.switch=F, rho.index = rho.index.inp, BF.index = 1000000, Max.Model.Dim = 10000000)))
    
  }
  
} else {
  
  #Run CARMA with check
  
  call_carma <- function(z.list.inp, ld.list.inp, lambda.list.inp, rho.index.inp){
    
    return(suppressWarnings(CARMA(z.list.inp, ld.list.inp, lambda.list=lambda.list.inp, outlier.switch=T, rho.index = rho.index.inp, BF.index = 1000000, Max.Model.Dim = 10000000)))
    
  }
  
}

# 5) Execute CARMA ====

### We're going to repeatedly try to fit the CARMA model using different random seeds in case of an error

#Set variables needed to track progress through loop

attempts <- 0

max_attempts <- 10

success <- FALSE

# Initialize timing variables
carma_time_seconds <- NA

# Loop until success or maximum attempts reached

while (!success && attempts < max_attempts) {
  
  tryCatch({
    
    set.seed(random_seed) #Nothing up until here has been random
    
    # Time the call_carma function specifically
    carma_start_time <- Sys.time()
    
    # Call the CARMA function here
    
    carma.results <- call_carma(z.list.inp = z.list, ld.list.inp = ld.list, lambda.list.inp = lambda.list, rho.index.inp = rho)
    
    carma_end_time <- Sys.time()
    
    # Calculate timing in seconds
    carma_time_seconds <- as.numeric(difftime(carma_end_time, carma_start_time, units = "secs"))
    
    success <- TRUE # Mark success if no error occurred
    
    break #Break the loop
    
  }, error = function(e) {
    
    attempts <<- attempts + 1
    
    random_seed <<- random_seed * 2 # Double the seed
    
    print(paste0("WARNING: Unable to fit ", file_prefix, " on attempt ", attempts, ". Retrying with seed: ", random_seed))
    
  })
  
}

if (!success) {
  
  # Handle failed runs - write output with NA values
  print(paste0("ERROR: Maximum attempts reached without success for ", file_prefix))
  
  # Write a partially empty output file and quit
  out_data <- cbind.data.frame(
    VARIANT_ID=summary_stats_raw[[markerid_col]],
    Z=bhat/shat,
    PIP=rep(NA, nrow(summary_stats_raw)),
    iterations_z=attempts,
    cs=rep(0, nrow(summary_stats_raw))
  )
  
  write.table(out_data, file = out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
  
  # Write timing file with NA for failed runs
  timing_dir <- file.path(dirname(out_file), "timing")
  dir.create(timing_dir, recursive = TRUE, showWarnings = FALSE)
  timing_file <- file.path(timing_dir, basename(out_file))
  write.table(data.frame(carma_time_seconds = NA), file = timing_file,
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
  
  quit(save="no");
  
}

cat("SUCCESS: CARMA function ran successfully after", attempts, "attempts.\n")
print(paste0("CARMA execution completed in ", round(carma_time_seconds, 2), " seconds"))

# 6) Extract Needed Info and Export Output File ====

print("Extracting results and writing output...")

#Extract Info

out_data <- cbind.data.frame(
  
  VARIANT_ID=summary_stats_raw[[markerid_col]],
  
  Z=bhat/shat,
  
  PIP=carma.results[[1]]$PIPs
  
)

#Get credible sets

cred_sets <- rep(0, nrow(out_data))

for (set in seq(1,length(carma.results[[1]]$'Credible set'[[2]]))){
  
  cred_sets[carma.results[[1]]$'Credible set'[[2]][[set]]] <- as.numeric(set)
  
}

out_data <- cbind.data.frame(out_data, iterations_z=attempts, cs=cred_sets)

#Save output dataframe

write.table(out_data, file = out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

# 7) Write Timing File ====

# Set timing subdirectory
timing_dir <- file.path(dirname(out_file), "timing")


# Create timing file with same name as output file
timing_file <- file.path(timing_dir, basename(out_file))

# Write timing data
timing_data <- data.frame(carma_time_seconds = carma_time_seconds)
write.table(timing_data, file = timing_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

print(paste0("Analysis complete. Results written to: ", out_file))
print(paste0("Timing information written to: ", timing_file))