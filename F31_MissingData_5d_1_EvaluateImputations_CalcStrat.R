# F31_MissingData_5d_1_EvaluateImputations
# Calculate MSE and classification error per column and overall
# Stratify based on whether cells are missing or not in original data
# This will enable us to check whether LOCF has an unfair advantage

# Set up ----
rm(list=ls(all=TRUE))

# Libraries
library(rtemis)

# Folder
folder <- "filepath"

# File date
file_date <- "date"
today <- format(Sys.Date(), "%Y-%m-%d")

# File path in
# Define dataset suffixes and their corresponding folder names
dataset_suffixes <- c("bp", "ext")
dataset_folders <- c("BP", "Ext")

# File path in
for (i in seq_along(dataset_suffixes)) {
  dataset_suffix <- dataset_suffixes[i]
  dataset_folder <- dataset_folders[i]
  assign(paste0("infolder_orig_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "orig"))
  assign(paste0("infolder_comp_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "comp"))
  assign(paste0("infolder_none_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "none"))
  assign(paste0("infolder_mean_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mean"))
  assign(paste0("infolder_locf_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "locf"))
  assign(paste0("infolder_mice_rf_av_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mice_rf", "av"))
  assign(paste0("infolder_mice_lasso_av_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mice_lasso", "av"))
  assign(paste0("infolder_mice_pmm_av_", dataset_suffix), file.path(folder, "MissingData_Result/3_Impute", file_date, dataset_folder, "mice_pmm", "av"))
}

# File path out
for (dataset_folder in dataset_folders) {
 
  outfolder <- file.path(folder, "MissingData_Result/5a_EvaluateImputations", today, dataset_folder)
  dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
  
  assign(paste0("outfolder_", tolower(dataset_folder)), outfolder, envir = .GlobalEnv)

  rm(outfolder)
}

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

# Function to calculate MSE and classification error for cells that were imputed
# Checked function in excel 12-23-23
calculate_errors_fun_strat <- function(comp, amp, imp, orig) {
  # Check if dimensions of all data frames are the same
  if (ncol(imp) != ncol(comp) || nrow(imp) != nrow(comp) || ncol(imp) != ncol(amp) || nrow(imp) != nrow(amp) || ncol(imp) != ncol(orig) || nrow(imp) != nrow(orig)) {
    stop("Dimensions of imp, comp, amp, and orig must be the same.")
  }

  # Initialize variables for numeric columns
  numeric_diffs_missing_in_orig <- c()
  numeric_diffs_not_missing_in_orig <- c()
  numeric_mse_list_missing_in_orig <- list()
  numeric_mse_list_not_missing_in_orig <- list()

  # Initialize variables for classification error
  classification_diffs_missing_in_orig <- c()
  classification_diffs_not_missing_in_orig <- c()
  classification_error_list_missing_in_orig <- list()
  classification_error_list_not_missing_in_orig <- list()

  # Initialize counters for cells
  count_missing_in_both_numeric <- 0
  count_missing_in_amp_only_numeric <- 0
  count_missing_in_both_classification <- 0
  count_missing_in_amp_only_classification <- 0

  # Loop through each column
  for (i in seq_along(imp)) {
    # Extract the i-th column from each data frame
    imp_column <- imp[[i]]
    comp_column <- comp[[i]]
    amp_column <- amp[[i]]
    orig_column <- orig[[i]]

    # Identify indices with missing values in amp and orig
    missing_indices_amp <- is.na(amp_column)
    missing_indices_orig <- is.na(orig_column)

    # Numeric Columns
    if (is.numeric(imp_column) && is.numeric(comp_column)) {
      # Calculate standard deviation of the comp column for numeric columns
      sd_comp <- sd(comp_column, na.rm = TRUE)

      if (any(missing_indices_amp)) {
        # Calculate counts for numeric columns
        count_missing_in_both_numeric <- count_missing_in_both_numeric + sum(missing_indices_amp & missing_indices_orig)
        count_missing_in_amp_only_numeric <- count_missing_in_amp_only_numeric + sum(missing_indices_amp & !missing_indices_orig)

        # Separate based on missingness in orig and standardize by SD
        imp_missing_in_orig <- imp_column[missing_indices_amp & missing_indices_orig] / sd_comp
        comp_missing_in_orig <- comp_column[missing_indices_amp & missing_indices_orig] / sd_comp
        imp_not_missing_in_orig <- imp_column[missing_indices_amp & !missing_indices_orig] / sd_comp
        comp_not_missing_in_orig <- comp_column[missing_indices_amp & !missing_indices_orig] / sd_comp

        # Calculate squared differences
        column_diffs_missing_in_orig <- comp_missing_in_orig - imp_missing_in_orig
        column_diffs_not_missing_in_orig <- comp_not_missing_in_orig - imp_not_missing_in_orig

        # Add to overall numeric differences and calculate MSE
        numeric_diffs_missing_in_orig <- c(numeric_diffs_missing_in_orig, column_diffs_missing_in_orig)
        numeric_diffs_not_missing_in_orig <- c(numeric_diffs_not_missing_in_orig, column_diffs_not_missing_in_orig)
        numeric_mse_list_missing_in_orig[[names(imp)[i]]] <- mean(column_diffs_missing_in_orig^2, na.rm = TRUE)
        numeric_mse_list_not_missing_in_orig[[names(imp)[i]]] <- mean(column_diffs_not_missing_in_orig^2, na.rm = TRUE)
      }
    } 
    
    # Classification Columns (Factor Columns)
    else if (is.factor(imp_column) && is.factor(comp_column)) {
      if (any(missing_indices_amp)) {
        # Calculate counts for classification columns
        count_missing_in_both_classification <- count_missing_in_both_classification + sum(missing_indices_amp & missing_indices_orig)
        count_missing_in_amp_only_classification <- count_missing_in_amp_only_classification + sum(missing_indices_amp & !missing_indices_orig)

        # Separate based on missingness in orig
        imp_missing_in_orig <- imp_column[missing_indices_amp & missing_indices_orig]
        comp_missing_in_orig <- comp_column[missing_indices_amp & missing_indices_orig]
        imp_not_missing_in_orig <- imp_column[missing_indices_amp & !missing_indices_orig]
        comp_not_missing_in_orig <- comp_column[missing_indices_amp & !missing_indices_orig]

        # Calculate classification errors
        column_errors_missing_in_orig <- comp_missing_in_orig != imp_missing_in_orig
        column_errors_not_missing_in_orig <- comp_not_missing_in_orig != imp_not_missing_in_orig

        # Add to overall classification error differences
        classification_diffs_missing_in_orig <- c(classification_diffs_missing_in_orig, column_errors_missing_in_orig)
        classification_diffs_not_missing_in_orig <- c(classification_diffs_not_missing_in_orig, column_errors_not_missing_in_orig)
        classification_error_list_missing_in_orig[[names(imp)[i]]] <- mean(column_errors_missing_in_orig, na.rm = TRUE)
        classification_error_list_not_missing_in_orig[[names(imp)[i]]] <- mean(column_errors_not_missing_in_orig, na.rm = TRUE)
      }
    } else {
      # Warning for non-matching column types
      warning(paste0("Column ", i, " does not have matching types. Skipping."))
    }
  }

  # Calculate overall MSE for numeric columns
  mse_missing_in_orig <- if (length(numeric_diffs_missing_in_orig) > 0) {
    mean(numeric_diffs_missing_in_orig^2, na.rm = TRUE)
  } else {
    NA
  }
  mse_not_missing_in_orig <- if (length(numeric_diffs_not_missing_in_orig) > 0) {
    mean(numeric_diffs_not_missing_in_orig^2, na.rm = TRUE)
  } else {
    NA
  }

  # Calculate overall classification error for categories
  classification_error_missing_in_orig <- if (length(classification_diffs_missing_in_orig) > 0) {
    mean(classification_diffs_missing_in_orig, na.rm = TRUE)
  } else {
    NA
  }
  classification_error_not_missing_in_orig <- if (length(classification_diffs_not_missing_in_orig) > 0) {
    mean(classification_diffs_not_missing_in_orig, na.rm = TRUE)
  } else {
    NA
  }

  # Return a list containing metrics and counts
  return(list(
    mse_missing_in_orig = mse_missing_in_orig,
    mse_not_missing_in_orig = mse_not_missing_in_orig,
    classification_error_missing_in_orig = classification_error_missing_in_orig,
    classification_error_not_missing_in_orig = classification_error_not_missing_in_orig,
    count_missing_in_both_numeric = count_missing_in_both_numeric,
    count_missing_in_amp_only_numeric = count_missing_in_amp_only_numeric,
    count_missing_in_both_classification = count_missing_in_both_classification,
    count_missing_in_amp_only_classification = count_missing_in_amp_only_classification,
    numeric_mse_per_column_missing_in_orig = numeric_mse_list_missing_in_orig,
    numeric_mse_per_column_not_missing_in_orig = numeric_mse_list_not_missing_in_orig,
    classification_error_per_column_missing_in_orig = classification_error_list_missing_in_orig,
    classification_error_per_column_not_missing_in_orig = classification_error_list_not_missing_in_orig
  ))
}

# Create function for missing and imputed datasets ----
error_imp_fun <- function(imp_suffix, out_suffix, traintest) {

  # i: the index of the 15 missingness scenarios
  # tt_list: the list of 20 datasets of one of the 15 missingness scenarios
  # j: the index of the 20 datasets within each scenario

  print(imp_suffix)
  print(out_suffix)
  print(traintest)

  file_path_comp <- file.path(get(paste0("infolder_comp_", out_suffix)), paste0("dt_tw_comp_tt_", out_suffix,".RDS"))
  tt_list_comp <- readRDS(file_path_comp)

  file_path_orig <- file.path(get(paste0("infolder_orig_", out_suffix)), paste0("dt_tw_orig_tt_", out_suffix,".RDS"))
  tt_list_orig <- readRDS(file_path_orig)

  list_names_amp <- get(paste0("dt_amp_names_tt_", out_suffix))

  list_names_imp <- paste0(list_names_amp, "_", imp_suffix)

  infolder_imp <- get(paste0("infolder_", imp_suffix, "_", out_suffix))

  outfolder <- get(paste0("outfolder_", out_suffix))

  list_out <- list()

  for (i in seq_along(list_names_imp)) {

    file_path_amp <- file.path(get(paste0("infolder_none_", out_suffix)), paste0(list_names_amp[i], ".RDS"))
    tt_list_amp <- readRDS(file_path_amp)

    list_name_imp <- list_names_imp[i]
    print(list_name_imp)

    file_path_imp <- file.path(infolder_imp, paste0(list_names_imp[i], ".RDS"))

    tt_list_imp <- readRDS(file_path_imp)
      
    list_out_i <- list()

    for (j in seq_along(tt_list_imp)) {
      
      if (is.null(tt_list_imp[[j]])) { # In case imputation fails, leave null placeholder so the remaining models can run without erroring out

        list_out_i[[j]] <- NULL

      } else {
      
      if (traintest=="train") {
       
        print(j)  
        orig <- tt_list_orig$x_train
        comp <- tt_list_comp$x_train
        amp <- tt_list_amp[[j]]$x_train
        imp <- tt_list_imp[[j]]$x_train

      }

      if (traintest=="test") {
       
        print(j)  
        orig <- tt_list_orig$x_test
        comp <- tt_list_comp$x_test
        amp <- tt_list_amp[[j]]$x_test
        imp <- tt_list_imp[[j]]$x_test

      }
                
      list_out_i[[j]] <- calculate_errors_fun_strat(comp, amp, imp, orig)

      }
    }

    list_out[[i]] <- list_out_i
    names(list_out)[i] <- get(paste0("amp_names_", out_suffix))[i]

  }

  file_path_out <<- file.path(outfolder, paste0("imp_perf_strat_", imp_suffix, "_", traintest, "_", out_suffix, ".RDS"))
  saveRDS(list_out, file_path_out)

}

# Calculate errors for imputed datasets ----
out_suffixes <- c("ext", "bp")

for (out_suffix in out_suffixes) {
  error_imp_fun("locf", out_suffix, "train")
  error_imp_fun("mean", out_suffix, "train")
  error_imp_fun("mice_rf_av", out_suffix, "train")
  error_imp_fun("mice_pmm_av", out_suffix, "train")
  error_imp_fun("mice_lasso_av", out_suffix, "train")

  error_imp_fun("locf", out_suffix, "test")
  error_imp_fun("mean", out_suffix, "test")
  error_imp_fun("mice_rf_av", out_suffix, "test")
  error_imp_fun("mice_pmm_av", out_suffix, "test")
  error_imp_fun("mice_lasso_av", out_suffix, "test")
}
