# F31_MissingData_5a_1_EvaluateImputations
# Calculate MSE and classification error per column and overall

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

# Count number of numeric vs. factor columns
file_path_comp_bp <- file.path(infolder_comp_bp, "dt_tw_comp_tt_bp.RDS")
tt_list_comp_bp <- readRDS(file_path_comp_bp)
num_numeric_cols_bp <- sum(sapply(tt_list_comp_bp$x_train, is.numeric))
num_factor_cols_bp <- sum(sapply(tt_list_comp_bp$x_train, is.factor))
print(num_numeric_cols_bp)
print(num_factor_cols_bp)

file_path_comp_ext <- file.path(infolder_comp_ext, "dt_tw_comp_tt_ext.RDS")
tt_list_comp_ext <- readRDS(file_path_comp_ext)
num_numeric_cols_ext <- sum(sapply(tt_list_comp_ext$x_train, is.numeric))
num_factor_cols_ext <- sum(sapply(tt_list_comp_ext$x_train, is.factor))
print(num_numeric_cols_ext)
print(num_factor_cols_ext)

# Function to calculate MSE and classification error for cells that were imputed
calculate_errors_fun <- function(comp, amp, imp) {
  # Check if dimensions of all data frames are the same
  if (ncol(imp) != ncol(comp) || nrow(imp) != nrow(comp) || ncol(imp) != ncol(amp) || nrow(imp) != nrow(amp)) {
    stop("Dimensions of imp, comp, and amp must be the same.")
  }

  # Initialize variables to hold overall differences for numeric and factor columns
  numeric_diffs <- c()
  factor_diffs <- c()

  # Initialize lists to store individual column-wise MSE and classification errors
  numeric_mse_list <- list()
  classification_error_list <- list()

  # Loop through each column
  for (i in seq_along(imp)) {
    # Extract the i-th column from each data frame
    imp_column <- imp[[i]]  # imputed data
    comp_column <- comp[[i]]  # ground truth
    amp_column <- amp[[i]]  # column indicating missing values

    # Calculate standard deviation of the comp column for numeric columns
    if (is.numeric(comp_column)) {
      sd_comp_bp <- sd(comp_column, na.rm = TRUE)
    }

    # Identify which indices have missing values in amp
    missing_indices <- is.na(amp_column)

    # If the column has missing values in amp, proceed
    if (any(missing_indices)) {
      # Only keep entries where the amp column has missing values
      imp_column_missing <- imp_column[missing_indices]
      comp_column_missing <- comp_column[missing_indices]

      # Case for Numeric Columns
      if (is.numeric(imp_column) && is.numeric(comp_column)) {
        # Standardize columns by the standard deviation
        imp_column_standardized <- imp_column_missing / sd_comp_bp
        comp_column_standardized <- comp_column_missing / sd_comp_bp

        # Calculate squared differences (truth - imputed)^2 for this column
        column_diffs <- comp_column_standardized - imp_column_standardized
        # Add to overall numeric differences
        numeric_diffs <- c(numeric_diffs, column_diffs)
        # Calculate and store the MSE for this specific column
        numeric_mse_list[[names(imp)[i]]] <- mean(column_diffs^2, na.rm = TRUE)

      # Case for Factor Columns
      } else if (is.factor(imp_column) && is.factor(comp_column)) {
        # Calculate classification error (truth != imputed) for this column
        column_errors <- comp_column_missing != imp_column_missing
        # Add to overall factor differences
        factor_diffs <- c(factor_diffs, column_errors)
        # Calculate and store the classification error for this specific column
        classification_error_list[[names(imp)[i]]] <- mean(column_errors, na.rm = TRUE)

      } else {
        # Issue a warning if the column types do not match between imp and comp
        warning(paste0("Column ", i, " does not have matching types. Skipping."))
      }
    }
  }

  # Calculate overall MSE for numeric columns if applicable
  if (length(numeric_diffs) > 0) {
    mse <- mean(numeric_diffs^2, na.rm = TRUE)
  } else {
    mse <- NA
  }

  # Calculate overall classification error for factor columns if applicable
  if (length(factor_diffs) > 0) {
    classification_error <- mean(factor_diffs, na.rm = TRUE)
  } else {
    classification_error <- NA
  }

  # Return a list containing overall and column-wise error metrics
  return(list(
    mse = mse,
    classification_error = classification_error,
    numeric_mse_per_column = numeric_mse_list,
    classification_error_per_column = classification_error_list
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

  file_path_comp <- file.path(get(paste0("infolder_comp_", out_suffix)), paste0("dt_tw_comp_tt_", out_suffix, ".RDS"))
  tt_list_comp <- readRDS(file_path_comp)

  names_amp <- get(paste0("amp_names_", out_suffix))

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
        comp <- tt_list_comp$x_train
        amp <- tt_list_amp[[j]]$x_train
        imp <- tt_list_imp[[j]]$x_train

      }

      if (traintest=="test") {
       
        print(j)  
        comp <- tt_list_comp$x_test
        amp <- tt_list_amp[[j]]$x_test
        imp <- tt_list_imp[[j]]$x_test

      }
                
      list_out_i[[j]] <- calculate_errors_fun(comp, amp, imp)

      }
    }

    list_out[[i]] <- list_out_i
    names(list_out)[i] <- names_amp[i]

  }

  file_path_out <<- file.path(outfolder, paste0("imp_perf_", imp_suffix, "_", traintest, "_", out_suffix, ".RDS"))
  saveRDS(list_out, file_path_out)

}

# Calculate errors for imputed datasets ----
# Define suffixes to loop over
out_suffixes <- c("ext", "bp")

# Loop through out_suffix and call error_imp_fun for each method
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
