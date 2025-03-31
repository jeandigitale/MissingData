# F31_MissingData_3a_Impute_all
# Call all R scripts to impute 

# Libraries ----
library(tidyverse)
library(openxlsx)
library(data.table)
library(zoo)
library(rtemis)
library(missRanger)
library(mice)

rm(list=ls(all=TRUE))

# Toggle on check ----
check <- 1

# Set-up and load data ----
# Folder
folder <- "filepath"

# File date
file_date <- "date"
today <- format(Sys.Date(), "%Y-%m-%d")

# Filename in
infolder_bp <- file.path(folder, "MissingData_Result/2_Ampute", file_date, "BP")
infolder_ext <- file.path(folder, "MissingData_Result/2_Ampute", file_date, "Ext")

# Filename out
if(check==0) {
  out_date <- today
}
if(check==1) {
  out_date <- "date"
}

# Define the subfolders and methods
subfolders <- c("BP", "Ext")
methods <- c("mean", "locf", "mice_rf/raw", "mice_lasso/raw", "mice_pmm/raw")

# Loop through subfolders and methods
for (subfolder in subfolders) {
  for (method in methods) {
    # Remove "raw" from method for variable naming
    method_clean <- gsub("/raw", "", method)
    
    # Construct the folder path
    outfolder <- file.path(folder, "MissingData_Result/3_Impute", out_date, subfolder, method)
    
    # Create the directory
    dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
    
    # Assign the path to the global environment with a dynamic variable name
    assign(paste0("outfolder_", method_clean, "_", tolower(subfolder)), outfolder, envir = .GlobalEnv)
  }
}

# Create log directories and assign paths
for (subfolder in subfolders) {
  log_dir <- file.path(folder, "MissingData_Log/3_Impute", today, subfolder)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Assign the log directory path to the global environment
  assign(paste0("log_dir_", tolower(subfolder)), log_dir, envir = .GlobalEnv)
}

progressr::handlers(global = TRUE)
options(progressr.handlers = progressr::handler_cli)

# File names
amp_names <- c("mcar_0", "mcar_1", "mcar_2", "mar_0", "mar_1", "mar_2", "mnar_w_0", "mnar_w_1", "mnar_w_2", "mnar_m_0", "mnar_m_1", "mnar_m_2", "mnar_s_0", "mnar_s_1", "mnar_s_2")

amp_names_ext <- paste0(amp_names, "_ext")
amp_names_mnar_ext_noy <- paste0(grep("mnar", amp_names_ext, value = TRUE), "_noy")
amp_names_ext <- c(amp_names_ext, amp_names_mnar_ext_noy)

amp_names_bp <- paste0(amp_names, "_bp")
amp_names_mnar_bp_noy <- paste0(grep("mnar", amp_names_bp, value = TRUE), "_noy")
amp_names_bp <- c(amp_names_bp, amp_names_mnar_bp_noy)

dt_amp_names_ext <- paste0("dt_", amp_names_ext)
dt_amp_names_bp <- paste0("dt_", amp_names_bp)

dt_amp_names_tt_ext <- gsub("_ext", "_tt_ext", dt_amp_names_ext)
dt_amp_names_tt_bp <- gsub("_bp", "_tt_bp", dt_amp_names_bp)

list_j <- 20 # number of datasets per list

# Seeds ----
# Set a seed to fix the seeds
set.seed(382)

# Initialize lists to store the seed sets for ext and bp
seeds_mice_rf_ext <- list()
seeds_mice_lasso_ext <- list()
seeds_mice_pmm_ext <- list()

seeds_mice_rf_bp <- list()
seeds_mice_lasso_bp <- list()
seeds_mice_pmm_bp <- list()

# Number of datasets and seeds per dataset
num_sets_ext <- length(dt_amp_names_tt_ext)
num_sets_bp <- length(dt_amp_names_tt_bp)
seeds_per_set <- 20

# Generate seed sets for ext
for (i in 1:num_sets_ext) {
  seeds_mice_rf_ext[[i]] <- sample(1:10000, seeds_per_set, replace = TRUE)
  seeds_mice_lasso_ext[[i]] <- sample(1:10000, seeds_per_set, replace = TRUE)
  seeds_mice_pmm_ext[[i]] <- sample(1:10000, seeds_per_set, replace = TRUE)
}

# Generate seed sets for bp
for (i in 1:num_sets_bp) {
  seeds_mice_rf_bp[[i]] <- sample(1:10000, seeds_per_set, replace = TRUE)
  seeds_mice_lasso_bp[[i]] <- sample(1:10000, seeds_per_set, replace = TRUE)
  seeds_mice_pmm_bp[[i]] <- sample(1:10000, seeds_per_set, replace = TRUE)
}

# Run all imputation scripts
if(check==0) {
  # Run imputation scripts
  source(file.path(folder, "F31_MissingData_3b_Impute_mean.R"))
  source(file.path(folder, "F31_MissingData_3b_Impute_locf.R"))
  source(file.path(folder, "F31_MissingData_3b_Impute_mice_rf.R"))
  source(file.path(folder, "F31_MissingData_3b_Impute_mice_pmm.R"))
  source(file.path(folder, "F31_MissingData_3b_Impute_mice_lasso.R"))
}

# Check whether files imputed and store indices of ones that did not
if(check==1) {
  # Initialize a list to store the dataframes for each suffix
  missing_indices_ext <- list()
  missing_indices_bp <- list()

  # Suffixes for imputations
  suffixes <- c("mice_rf", "mice_lasso", "mice_pmm")

  # Function to check missing files and store indices
  check_missing_files <- function(suffixes, subfolder, dt_amp_names_tt, list_j) {
    missing_indices <- list()
    
    for (s in seq_along(suffixes)) {
        imp_suffix <- suffixes[s]
        
        # Dynamically get the outfolder path from the global environment
        folder_suff <- get(paste0("outfolder_", imp_suffix, "_", tolower(subfolder)))
        list_names <- paste0(dt_amp_names_tt, "_", imp_suffix)
        
        # Initialize vectors to store i and j values for non-existent files
        i_values <- c()
        j_values <- c()
        
        for (i in seq_along(list_names)) {
            list_name <- list_names[i]
            
            for (j in 1:list_j) {
                file_name <- file.path(folder_suff, paste0(list_name, "_", j, "_raw.RDS"))

                # Check if the file exists
                if (!file.exists(file_name)) {
                    # Store the i and j values
                    i_values <- c(i_values, i)
                    j_values <- c(j_values, j)
                }
            }
        }

        # Create a dataframe for the current suffix and add it to the list
        df_current_suffix <- data.frame(i = i_values, j = j_values)
        missing_indices[[imp_suffix]] <- df_current_suffix
    }
    
    return(missing_indices)
  }

  # Check for missing files for ext and bp
  missing_indices_ext <- check_missing_files(suffixes = suffixes, subfolder = "Ext", dt_amp_names_tt = dt_amp_names_tt_ext, list_j = list_j)
  missing_indices_bp <- check_missing_files(suffixes = suffixes, subfolder = "BP", dt_amp_names_tt = dt_amp_names_tt_bp, list_j = list_j)

  # Run imputation scripts
  #source(file.path(folder, "F31_MissingData_3b_Impute_mice_rf.R"))
  #source(file.path(folder, "F31_MissingData_3b_Impute_mice_pmm.R"))
  #source(file.path(folder, "F31_MissingData_3b_Impute_mice_lasso.R"))
  }
