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

#Need to read in the following: 
  # 1) summary stats file pruned for locus
  # 2) ld matrix file in .vcor1 or similar format
  # 3) output file location
  # 4) Toggle for Checking Outliers OPTIONAL (default is TRUE)
  # 5) Lambda for Poisson Regularization OPTIONAL (default is 1)
  # 6) Rho for credible set confidence OPTIONAL (default is 0.95)
  # 7) random seed [OPTIONAL; default is FALSE]

args <- commandArgs(trailingOnly = TRUE);
if (length(args) == 7) {
  #Print confirmation statement
  print("Necessary Parameters Input")
  #Set variables from the arguments in this case
  summary_stats_file <- args[1]
  ld_file <- args[2]
  out_file <- args[3]
  outlier_check_toggle <- args[4]
  lambda <- as.numeric(args[5])
  rho <- as.numeric(args[6])
  random_seed <- as.integer(args[7])
} else if(length(args) == 6) {
  #Set variables from the arguments in this case
  summary_stats_file <- args[1]
  ld_file <- args[2]
  out_file <- args[3]
  outlier_check_toggle <- args[4]
  lambda <- as.numeric(args[5])
  rho <- as.numeric(args[6])
  random_seed <- 5
} else if(length(args) == 4) {
  #Set variables from the arguments in this case
  summary_stats_file <- args[1]
  ld_file <- args[2]
  out_file <- args[3]
  outlier_check_toggle <- TRUE
  lambda <- 1
  rho <- 0.95
  random_seed <- 5
} else {
  print("Usage: %> Rscript finemap_CARMA.R summary_stats_file  ld_file out_file [outlier_check_toggle] [lambda] [rho] [random_seed]");
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

### Read in input files and extract susie inputs ###
#Start with summary stats
summary_stats_raw <- read.table(summary_stats_file, header = TRUE)
n <- mean(summary_stats_raw$n)
bhat <- summary_stats_raw$beta
shat <- summary_stats_raw$se
#Lastly read in correlation matrix
R <- read.csv(ld_file, sep = "\t", header = FALSE)
rownames(R) <- summary_stats_raw$RSID
colnames(R) <- summary_stats_raw$RSID

#Pack everything into lists
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
set.seed(random_seed) #Nothing up until here has been random
#Fit susie model
carma.results <- call_carma(z.list.inp = z.list, ld.list.inp = ld.list, lambda.list.inp = lambda.list, rho.index.inp = rho)

# 6) Extract Needed Info and Export Output File ====

#Extract Info
out_data <- cbind.data.frame(RSID=summary_stats_raw$rsid, Z=bhat/shat, PIP=carma.results[[1]]$PIPs)
#Get credible sets and outliers
outliers <- rep(NA, nrow(out_data))
outliers[carma.results[[1]]$'Outliers'[['Index']]] <- carma.results[[1]]$'Outliers'[['Z']]
cred_sets <- rep(0, nrow(out_data))
for (set in seq(1,length(carma.results[[1]]$'Credible set'[[2]]))){
  cred_sets[carma.results[[1]]$'Credible set'[[2]][[set]]] <- as.numeric(set)
}
out_data <- cbind.data.frame(out_data, outlier_est_z=outliers, cs=cred_sets)

#Save output dataframe
write.table(out_data, file = out_file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")

