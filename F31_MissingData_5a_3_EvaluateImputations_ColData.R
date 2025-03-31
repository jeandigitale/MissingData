# F31_MissingData_5a_3_EvaluateImputations
# Create dataset of column errors to feed into F31_MissingData_5c_ModelsofModels

# Set up ----
rm(list=ls(all=TRUE))

# Libraries
library(rtemis)
library(data.table)
library(ggplot2)
library(gridExtra)

# Folder
folder <- "filepath"
today <- format(Sys.Date(), "%Y-%m-%d")

# File path in: imputation performance
file_date <- "date"

# Define dataset suffixes and their corresponding folder names
dataset_suffixes <- c("bp", "ext")
dataset_folders <- c("BP", "Ext")

for (dataset_folder in dataset_folders) {
 
  outfolder <- file.path(folder, "MissingData_Result/5a_EvaluateImputations", file_date, dataset_folder)
  dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
  
  assign(paste0("infolder_", tolower(dataset_folder)), outfolder, envir = .GlobalEnv)

  rm(outfolder)
}

# Out suffixes
out_suffixes <- c("ext", "bp")

# List of suffixes
suffixes <- c("comp", "none", "mean", "locf", "mice_rf_av",  "mice_rf",  "mice_pmm_av",  "mice_pmm", "mice_lasso_av",  "mice_lasso")

## File names ----
amp_names <- c("mcar_0", "mcar_1", "mcar_2", "mar_0", "mar_1", "mar_2", "mnar_w_0", "mnar_w_1", "mnar_w_2", "mnar_m_0", "mnar_m_1", "mnar_m_2", "mnar_s_0", "mnar_s_1", "mnar_s_2")

amp_names_ext <- paste0(amp_names, "_ext")
amp_names_mnar_ext_noy <- paste0(grep("mnar", amp_names_ext, value = TRUE), "_noy")
amp_names_ext <- c(amp_names_ext, amp_names_mnar_ext_noy)

amp_names_bp <- paste0(amp_names, "_bp")
amp_names_mnar_bp_noy <- paste0(grep("mnar", amp_names_bp, value = TRUE), "_noy")
amp_names_bp <- c(amp_names_bp, amp_names_mnar_bp_noy)

col_perf_dt_fun <- function(imp_suffix, out_suffix) {

  infolder_imp <- get(paste0("infolder_", out_suffix))

  # Read in the data
  imp_train_list <- readRDS(file.path(infolder_imp, paste0("imp_perf_", imp_suffix, "_train_", out_suffix, ".RDS")))
  imp_test_list <- readRDS(file.path(infolder_imp, paste0("imp_perf_", imp_suffix, "_test_", out_suffix, ".RDS")))

  # Initialize empty vectors to store the results
  imp_vec <- c()
  amp_vec <- c()
  variable_vec <- c()
  metric_type_vec <- c()
  error_value_vec <- c()
  set_type_vec <- c()
  j_index_vec <- c()  # New vector to store the value of j

  # Function to populate vectors with metric data
  populate_vectors <- function(imp_list, set_type, amp_name, imp_suffix) {
    for(j in seq_along(imp_list)) {
      # Check if it's a list and has the necessary components
      if(is.list(imp_list[[j]])) {
        if (!is.null(imp_list[[j]]$numeric_mse_per_column)) {
          for(var_name in names(imp_list[[j]]$numeric_mse_per_column)) {
            imp_vec <<- c(imp_vec, imp_suffix)
            amp_vec <<- c(amp_vec, amp_name)
            variable_vec <<- c(variable_vec, var_name)
            metric_type_vec <<- c(metric_type_vec, "numeric")
            error_value_vec <<- c(error_value_vec, imp_list[[j]]$numeric_mse_per_column[[var_name]])
            set_type_vec <<- c(set_type_vec, set_type)
            j_index_vec <<- c(j_index_vec, j)  # Add j to the vector
          }
        }
        if (!is.null(imp_list[[j]]$classification_error_per_column)) {
          for(var_name in names(imp_list[[j]]$classification_error_per_column)) {
            imp_vec <<- c(imp_vec, imp_suffix)
            amp_vec <<- c(amp_vec, amp_name)
            variable_vec <<- c(variable_vec, var_name)
            metric_type_vec <<- c(metric_type_vec, "factor")
            error_value_vec <<- c(error_value_vec, imp_list[[j]]$classification_error_per_column[[var_name]])
            set_type_vec <<- c(set_type_vec, set_type)
            j_index_vec <<- c(j_index_vec, j)  # Add j to the vector
          }
        }
      }
    }
  }

  # Loop through train and test lists to populate vectors
  for(i in seq_along(imp_train_list)) {
    amp_name <- names(imp_train_list)[i]

    # Remove any occurrence of "bp" or "ext"
    amp_name <- gsub("bp_|ext_|_bp|_ext", "", amp_name)

    # Remove the last suffix (from the imputation types)
    amp_name <- sub(paste0("(_", paste(suffixes, collapse = "|_"), ")$"), "", amp_name)

    populate_vectors(imp_train_list[[i]], "train", amp_name, imp_suffix)
  }

  for(i in seq_along(imp_test_list)) {
    amp_name <- names(imp_test_list)[i]

    # Remove any occurrence of "bp" or "ext"
    amp_name <- gsub("bp_|ext_|_bp|_ext", "", amp_name)

    # Remove the last suffix (from the imputation types)
    amp_name <- sub(paste0("(_", paste(suffixes, collapse = "|_"), ")$"), "", amp_name)

    populate_vectors(imp_test_list[[i]], "test", amp_name, imp_suffix)
  }

  # Combine all the vectors into a data frame
  final_df <- data.frame(
    imp = imp_vec,
    miss = amp_vec,
    variable = variable_vec,
    metric_type = metric_type_vec,
    error_value = error_value_vec,
    set_type = set_type_vec,
    index = j_index_vec  # Include j_index in the final dataframe
  )
  
  return(final_df)
}

# Combine imputation performance into one dataset ----
imp_suffix_all <- c("locf", "mean", "mice_rf_av", "mice_pmm_av", "mice_lasso_av")

# Use lapply to apply the function to each imp_suffix and store results in a list
col_perf_list_ext <- lapply(imp_suffix_all, col_perf_dt_fun, out_suffix = "ext")
col_perf_list_bp <- lapply(imp_suffix_all, col_perf_dt_fun, out_suffix = "bp")

# Concatenate all data frames in the list into one large data frame
col_perf_ext <- as.data.table(do.call(rbind, col_perf_list_ext))
col_perf_bp <- as.data.table(do.call(rbind, col_perf_list_bp))

# Order the imputation variable
col_perf_ext$imp <- factor(col_perf_ext$imp, levels = imp_suffix_all)
col_perf_bp$imp <- factor(col_perf_bp$imp, levels = imp_suffix_all)

# Create new columns for type and times missing (0,1,2) ----
col_perf_ext[, c("miss_type", "p_miss") := .(sub("_\\d$", "", miss), as.numeric(sub("^.*_", "", miss)))]
col_perf_ext[miss=="comp", miss_type := "comp"]
col_perf_ext$miss_type <- factor(col_perf_ext$miss_type, 
                             levels = c("mcar", "mar", "mnar_w", "mnar_m", "mnar_s", "comp"))
col_perf_ext$p_miss <- factor(col_perf_ext$p_miss, levels = c(0, 1, 2), labels = c("0.5", "1", "2"))
col_perf_ext[miss=="comp", p_miss := "0"]
col_perf_ext$p_miss <- factor(col_perf_ext$p_miss, levels = c(0.5, 1, 2, 0))

col_perf_bp[, c("miss_type", "p_miss") := .(sub("_\\d$", "", miss), as.numeric(sub("^.*_", "", miss)))]
col_perf_bp[miss=="comp", miss_type := "comp"]
col_perf_bp$miss_type <- factor(col_perf_bp$miss_type, 
                             levels = c("mcar", "mar", "mnar_w", "mnar_m", "mnar_s", "comp"))
col_perf_bp$p_miss <- factor(col_perf_bp$p_miss, levels = c(0, 1, 2), labels = c("0.5", "1", "2"))
col_perf_bp[miss=="comp", p_miss := "0"]
col_perf_bp$p_miss <- factor(col_perf_bp$p_miss, levels = c(0.5, 1, 2, 0))

# Save
saveRDS(col_perf_bp, paste0(infolder_bp, "/col_perf_bp.RDS"))
saveRDS(col_perf_ext, paste0(infolder_ext, "/col_perf_ext.RDS"))
