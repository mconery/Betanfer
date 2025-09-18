################################################################################

#finemap_SuSiE.R

#This script is designed to execute SuSiE fine-mapping given a fixed region

#The code should be run in the susie conda environment

################################################################################

# 0) Load Needed Libraries ====

library(susieR)
library(stringr)

# 1) Read in Needed Arguments from Command ====

#Need to read in the following:
# 1) summary stats file (not pruned - will be filtered using variant file)
# 2) comma-delimited list of phenotype files (each must have column named "PHENO")
# 3) ld matrix file in .vcor1 or similar format
# 4) output file location
# 5) column name for variant ID in summary stats [OPTIONAL; default is "markerid"]
# 6) column name for effect size (beta) in summary stats [OPTIONAL; default is "beta"]
# 7) column name for standard error in summary stats [OPTIONAL; default is "se"]
# 8) column name for sample size in summary stats [OPTIONAL; default is "n"]
# 9) flag indicating whether in-sample fine-mapping is being used [OPTIONAL; default is TRUE]
# 10) random seed [OPTIONAL; default is 5]

args <- commandArgs(trailingOnly = TRUE);

if (length(args) >= 3) {
  # Required parameters
  summary_stats_file <- args[1]
  ld_file <- args[2]
  out_file <- args[3]
  # Optional parameters with defaults
  markerid_col <- if (length(args) >= 4) args[4] else "markerid"
  beta_col <- if (length(args) >= 5) args[5] else "beta"
  se_col <- if (length(args) >= 6) args[6] else "se"
  n_col <- if (length(args) >= 7) args[7] else "n"
  in_sample_flag <- if (length(args) >= 8) as.logical(args[8]) else TRUE
  random_seed <- if (length(args) >= 9) as.integer(args[9]) else 5
  
  print(paste0("Necessary Parameters Input. Using ", length(pheno_files), " phenotype files."))
  
} else {
  print("Usage: Rscript finemap_SuSiE.R summary_stats_file ld_file out_file [markerid_col] [beta_col] [se_col] [n_col] [in_sample_flag] [random_seed]");
  print("Example: Rscript finemap_SuSiE.R sumstats.txt matrix.vcor1 output.txt");
  quit(save="no");
}

# 2) Prep for executing SuSiE ====

#Check output location
if (file.exists(dirname(out_file)) == FALSE) {
  print("ERROR: Invalid output folder given.")
  quit(save="no");
}

#Check input locations
if (file.exists(summary_stats_file) == FALSE || file.exists(ld_file) == FALSE) {
  print("ERROR: Summary stats file or LD matrix file not found.")
  quit(save="no");
}

### Read in input files and extract susie inputs ###

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

# 3) Create Error Handling function ====

check <- function(expression){
  withCallingHandlers(expression,
                      warning = function(w){
                        # Handle warnings if needed
                      },
                      error = function(e){
                        message(paste0("ERROR: ", str_replace(basename(summary_stats_file), pattern = "_out_snps_sim_gwas.tsv", replacement = ""), " has failed prior estimation."))
                        
                        #Write a partially empty output file and quit
                        out_data <- cbind.data.frame(
                          VARIANT_ID=summary_stats_raw[[markerid_col]], 
                          Z=bhat/shat, 
                          PIP=rep(NA, nrow(summary_stats_raw)), 
                          converged=rep(NA, nrow(summary_stats_raw))
                        )
                        out_data <- cbind.data.frame(out_data, cs=rep(0, nrow(out_data)))
                        write.table(out_data, file = out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
                        quit(save="no");
                      },
                      finally = {
                        # Cleanup if needed
                      }
  )
}

# 4) Execute SuSiE ====

set.seed(random_seed) #Nothing up until here has been random

#Fit susie model
print("Running SuSiE fine-mapping...")
fitted_rss1 <- check(susie_rss(bhat = bhat, shat = shat, n = n, R = R, L = 10,
                               estimate_residual_variance = in_sample_flag))

# 5) Extract Needed Info and Export Output File ====

print("Extracting results and writing output...")

#Extract Info
out_data <- cbind.data.frame(
  VARIANT_ID=summary_stats_raw[[markerid_col]], 
  Z=bhat/shat, 
  PIP=fitted_rss1$pip, 
  converged=rep(fitted_rss1$converged, nrow(summary_stats_raw))
)

#Get credible sets
cred_sets <- rep(0, nrow(out_data))
for (set in names(fitted_rss1$sets$cs)){
  cred_sets[fitted_rss1$sets$cs[[set]]] <- as.numeric(str_replace(set, pattern = "L", replacement = ""))
}

out_data <- cbind.data.frame(out_data, cs=cred_sets)

#Save output dataframe
write.table(out_data, file = out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

print(paste0("Analysis complete. Results written to: ", out_file))
