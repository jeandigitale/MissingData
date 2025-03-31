# F31_MissingData_4a_Model
# Load 1 complete and multiple imputed datasets
# Run outcome models
# Submits one job per dataset list

# Set up ----
rm(list=ls(all=TRUE))
progressr::handlers(global = TRUE)
options(progressr.handlers = progressr::handler_cli)

# Toggle check on and off
check <- 1

# Libraries
library(rtemis)

# Folder
folder <- "filepath"

# File date
file_date <- "date"
today <- format(Sys.Date(), "%Y-%m-%d")
model_date <- ifelse(check==1, "date", today)

# Define dataset suffixes and their corresponding folder names
dataset_suffixes <- c("bp", "ext")
dataset_folders <- c("BP", "Ext")

# File path in
for (i in seq_along(dataset_suffixes)) {
  dataset_suffix <- dataset_suffixes[i]
  dataset_folder <- dataset_folders[i]
  
  assign(paste0("infolder_comp_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "comp"))
  assign(paste0("infolder_none_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "none"))
  assign(paste0("infolder_mean_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mean"))
  assign(paste0("infolder_locf_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "locf"))
  assign(paste0("infolder_mice_rf_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mice_rf", "tt"))
  assign(paste0("infolder_mice_lasso_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mice_lasso", "tt"))
  assign(paste0("infolder_mice_pmm_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mice_pmm", "tt"))
  assign(paste0("infolder_mice_rf_av_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mice_rf", "av"))
  assign(paste0("infolder_mice_lasso_av_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mice_lasso", "av"))
  assign(paste0("infolder_mice_pmm_av_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mice_pmm", "av"))
}

# File path out
# List of suffixes
suffixes <- c("comp", "none", "mean", "locf", "mice_rf", "mice_lasso", "mice_pmm", "mice_rf_av", "mice_lasso_av", "mice_pmm_av")

# List of model names
model_names <- c("gbm", "glmnet")

# Loop through BP/Ext
for (i in seq_along(dataset_suffixes)) {
  dataset_suffix <- dataset_suffixes[i]
  dataset_folder <- dataset_folders[i]

  # Loop through model names
  for (model_name in model_names) {
    # Create a folder path based on today's date and model name
    model_folder <- file.path(folder, "MissingData_Result/4_Model", model_date, dataset_folder, model_name)
    dir.create(model_folder, recursive = TRUE, showWarnings = FALSE)
    
    # Loop through suffixes
    for (suffix in suffixes) {
      # Create a folder path based on the current suffix
      outfolder <- file.path(model_folder, suffix)
      dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
      
      # Create a variable name based on the model and suffix
      variable_name <- paste("outfolder", model_name, suffix, dataset_suffix, sep = "_")
      
      # Assign the folder path to the global environment
      assign(variable_name, outfolder, envir = .GlobalEnv)
    }
  }
}

# File path log directory
for (i in seq_along(dataset_suffixes)) {
  dataset_suffix <- dataset_suffixes[i]
  dataset_folder <- dataset_folders[i]

  log_dir <- file.path(folder, "MissingData_Log/4_Model", today, dataset_folder)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Assign the log directory path to the global environment
  assign(paste0("log_dir_", tolower(dataset_folder)), log_dir, envir = .GlobalEnv)
}

# Number of datasets per missingness scenario
n_datasets <- 20

# File names
amp_names <- c("mcar_0", "mcar_1", "mcar_2", "mar_0", "mar_1", "mar_2", "mnar_w_0", "mnar_w_1", "mnar_w_2", "mnar_m_0", "mnar_m_1", "mnar_m_2", "mnar_s_0", "mnar_s_1", "mnar_s_2")

suffixes <- c("mean", "locf", "mice_rf", "mice_lasso", "mice_pmm", "mice_rf_av", "mice_lasso_av", "mice_pmm_av") # exclude comp and none here based on model functions

# Loop over dataset_suffixes
for (dataset_suffix in dataset_suffixes) {
  
  # Create extended names for the current dataset suffix
  amp_names_current <- paste0(amp_names, "_", dataset_suffix)
  
  # Handle the "_noy" addition for mnar variables
  amp_names_mnar_noy <- paste0(grep("mnar", amp_names_current, value = TRUE), "_noy")
  amp_names_current <- c(amp_names_current, amp_names_mnar_noy)
  
  # Create dt_amp_names for the current dataset suffix (without imputation suffix)
  dt_amp_names_tt <- gsub(paste0("_", dataset_suffix), 
                          paste0("_tt_", dataset_suffix), 
                          paste0("dt_", amp_names_current))
  
  # Assign dt_amp_names_tt_ext and dt_amp_names_tt_bp without suffix
  assign(paste0("dt_amp_names_tt_", dataset_suffix), dt_amp_names_tt, envir = .GlobalEnv)
  
  # Loop through each imputation suffix (locf, mean, etc.)
  for (suffix in suffixes) {
    
    # Create the final vector name
    vector_name <- paste0("dt_amp_names_tt_", dataset_suffix, "_", suffix)
    
    # Generate the vector with transformed names and make sure "_noy" comes first
    dt_amp_names_tt_with_suffix <- gsub(paste0("_", dataset_suffix), 
                                        paste0("_tt_", dataset_suffix), 
                                        paste0("dt_", amp_names_current))
    
    # Place the suffix after "_noy" if it exists
    dt_amp_names_tt_with_suffix <- gsub("_noy", paste0("_noy_", suffix), dt_amp_names_tt_with_suffix)
    
    # For non-mnar variables, append the suffix at the end
    dt_amp_names_tt_with_suffix <- ifelse(grepl("_noy", dt_amp_names_tt_with_suffix), 
                                          dt_amp_names_tt_with_suffix, 
                                          paste0(dt_amp_names_tt_with_suffix, "_", suffix))
    
    # Assign the vector with suffix to the global environment
    assign(vector_name, dt_amp_names_tt_with_suffix, envir = .GlobalEnv)
  }
}

# Check whether models ran for MI only and store indices of ones that did not ----
# Initialize a list to store the dataframes for each suffix
if(check==1) {

  # Function to find missing indices for specified model type and out_suffix
  find_missing_indices <- function(
    model_type, 
    n_datasets_f = n_datasets, 
    mice_suffixes = c("mice_rf", "mice_lasso", "mice_pmm"), 
    out_suffix) {
  
    # Initialize a list to store the missing indices for each suffix
    missing_indices <- list()

    for (imp_suffix in mice_suffixes) {
        folder_suff <- get(paste0("outfolder_", model_type, "_", imp_suffix, "_", out_suffix))

        list_names <- get(paste0("dt_amp_names_tt_", out_suffix, "_",imp_suffix))
        list_names <- paste0(gsub("_tt", "", gsub("dt", model_type, list_names)))
        
        # Initialize vectors to store i and j values for non-existent files
        i_values <- c()
        j_values <- c()
        
        for (list_name in list_names) {
            index <- which(list_names == list_name)
            print(list_name)
            print(index)

            for (j in 1:n_datasets_f) {
                file_name <- file.path(folder_suff, paste0(list_name, "_", j, ".RDS"))

                print(j)

                # Check if the file exists
                if (!file.exists(file_name)) {
                    # Store the i and j values
                    i_values <- c(i_values, index)
                    j_values <- c(j_values, j)
                }
            }
        }

        # Store the data frame for the current suffix in the missing_indices list
        df_current_suffix <- data.frame(i = i_values, j = j_values)
        missing_indices[[imp_suffix]] <- df_current_suffix
    }
    
    return(missing_indices)

  }

# Check for missing files for ext and bp using the updated function
  missing_indices_gbm_ext <- find_missing_indices(model_type = "gbm", out_suffix = "ext")
  missing_indices_gbm_bp <- find_missing_indices(model_type = "gbm", out_suffix = "bp")

  missing_indices_glmnet_ext <- find_missing_indices(model_type = "glmnet", out_suffix = "ext")
  missing_indices_glmnet_bp <- find_missing_indices(model_type = "glmnet", out_suffix = "bp")

}

# Run model scripts
#source(file.path(folder, "F31_MissingData_4b_Model_gbm.R"))
#source(file.path(folder, "F31_MissingData_4b_Model_glmnet.R"))
