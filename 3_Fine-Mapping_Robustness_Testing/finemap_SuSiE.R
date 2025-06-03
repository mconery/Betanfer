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
  # 1) summary stats file pruned for locus
  # 2) file of phenotype values
  # 3) ld matrix file in .vcor1 or similar format
  # 4) output file location
  # 5) flag indicating whether in-sample fine-mapping is being used [OPTIONAL; default is TRUE]
  # 6) random seed [OPTIONAL; default is FALSE]

args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 6) {
  #Print confirmation statement
  print("Necessary Parameters Input")
  #Set variables from the arguments in this case
  summary_stats_file <- args[1]
  pheno_file <- args[2]
  ld_file <- args[3]
  out_file <- args[4]
  in_sample_flag <- as.logical(args[5])
  random_seed <- as.integer(args[6])
} else if(length(args) == 5) {
  #Set variables from the arguments in this case
  summary_stats_file <- args[1]
  pheno_file <- args[2]
  ld_file <- args[3]
  out_file <- args[4]
  in_sample_flag <- as.logical(args[5])
  random_seed <- 5
} else if(length(args) == 4) {
  #Set variables from the arguments in this case
  summary_stats_file <- args[1]
  pheno_file <- args[2]
  ld_file <- args[3]
  out_file <- args[4]
  in_sample_flag <- TRUE
  random_seed <- 5
} else {
  print("Usage: %> Rscript finemap_SuSiE.R summary_stats_file pheno_file ld_file out_file [in_sample_flag] [random_seed]");
  quit(save="no");
}

# 2) Prep for executing SuSiE ====

#Check output location
if (file.exists(dirname(out_file)) == FALSE) {
  print("ERROR: Invalid output folder given.")
  quit(save="no");
}

#Check input locations
if (file.exists(summary_stats_file) == FALSE || file.exists(pheno_file) == FALSE || file.exists(ld_file) == FALSE) {
  print("ERROR: One or more invalid input files given.")
  quit(save="no");
}

### Read in input files and extract susie inputs ###
#Start with summary stats
summary_stats_raw <- read.table(summary_stats_file, header = TRUE)
n <- mean(summary_stats_raw$n)
bhat <- summary_stats_raw$beta
shat <- summary_stats_raw$se
#Next read in pheno file
pheno_raw <- read.table(pheno_file, header = FALSE)
var_y <- var(pheno_raw[,3])
remove(pheno_raw)
#Lastly read in correlation matrix
R <- read.csv(ld_file, sep = "\t", header = FALSE)

# 3) Execute SuSiE ====
set.seed(random_seed) #Nothing up until here has been random
#Fit susie model
fitted_rss1 <- susie_rss(bhat = bhat, shat = shat, n = n, R = R, var_y = var_y, L = 10,
                         estimate_residual_variance = in_sample_flag)

# 4) Extract Needed Info and Export Output File ====

#Extract Info
out_data <- cbind.data.frame(RSID=summary_stats_raw$rsid, Z=bhat/shat, PIP=fitted_rss1$pip, converged=rep(fitted_rss1$converged, nrow(summary_stats_raw)))
#Get credible sets
cred_sets <- rep(0, nrow(out_data))
for (set in names(fitted_rss1$sets$cs)){
  cred_sets[fitted_rss1$sets$cs[[set]]] <- as.numeric(str_replace(set, pattern = "L", replacement = ""))
}
out_data <- cbind.data.frame(out_data, cs=cred_sets)

#Save output dataframe
write.table(out_data, file = out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

