# F31_MissingData_5b_1_AverageMI
# For each of the 20 datasets per missingness scenario, average the predicted train and test probabilities of the multiple imputation results
# Calculate performance metrics
# Output list of length 20 per missingness scenario

# Set up ----
rm(list=ls(all=TRUE))

# Libraries
library(rtemis)
library(data.table)
library(ggplot2)
library(gridExtra)

# Folders ----
folder <- "filepath"

today <- format(Sys.Date(), "%Y-%m-%d")

# File path in
file_date <- "date"

# List of suffixes
suffixes <- c("mice_rf", "mice_pmm", "mice_lasso")

# List of model names
model_names <- c("gbm", "glmnet")

# Dataset suffixes
dataset_suffixes <- c("bp", "ext")
dataset_folders  <- c("BP", "Ext")

# Number of datasets
n_datasets <- 20

# Initialize an empty list to hold the infolders
infolder_list <- list()

for (model_name in model_names) {
  for (i in seq_along(dataset_suffixes)) {
    dataset_suffix <- dataset_suffixes[i]  # e.g., "bp" or "ext"
    dataset_folder <- dataset_folders[i]   # e.g., "BP" or "Ext"

    # Build the path to the model folder
    model_folder <- file.path(folder, "MissingData_Result/4_Model", 
                              file_date, dataset_folder, model_name)

    for (suffix in suffixes) {
      # The directory containing each imputation method
      infolder <- file.path(model_folder, suffix)

      # We'll store the variable name using the dataset_suffix in lowercase
      # e.g., "infolder_gbm_mice_rf_bp"
      variable_name <- paste("infolder", model_name, suffix, dataset_suffix, sep = "_")

      infolder_list[[variable_name]] <- infolder
    }
  }
}

## File names ----
amp_names <- c("mcar_0", "mcar_1", "mcar_2", "mar_0", "mar_1", "mar_2", "mnar_w_0", "mnar_w_1", "mnar_w_2", "mnar_m_0", "mnar_m_1", "mnar_m_2", "mnar_s_0", "mnar_s_1", "mnar_s_2")

amp_names_ext <- paste0(amp_names, "_ext")
amp_names_mnar_ext_noy <- paste0(grep("mnar", amp_names_ext, value = TRUE), "_noy")
amp_names_ext <- c(amp_names_ext, amp_names_mnar_ext_noy)

amp_names_bp <- paste0(amp_names, "_bp")
amp_names_mnar_bp_noy <- paste0(grep("mnar", amp_names_bp, value = TRUE), "_noy")
amp_names_bp <- c(amp_names_bp, amp_names_mnar_bp_noy)

# Initialize an empty list to hold the results
filename_list_bp <- list()
filename_list_ext <- list()

# List of suffixes for filenames
amp_suffixes <- c("_mice_rf", "_mice_pmm", "_mice_lasso")

# Nested loops to fill filename_list_bp and filename_list_ext
for (model_name in model_names) {
  for (amp_suffix in amp_suffixes) {
    # ----- For BP -----
    # Construct a name to store in filename_list_bp
    list_name_bp <- paste0(model_name, "_amp_names", amp_suffix, "_bp")
    # Construct the values using amp_names_bp
    list_value_bp <- paste0(model_name, "_", amp_names_bp, amp_suffix)
    # Store in filename_list_bp
    filename_list_bp[[list_name_bp]] <- list_value_bp

    # ----- For Ext -----
    # Construct a name to store in filename_list_ext
    list_name_ext <- paste0(model_name, "_amp_names", amp_suffix, "_ext")
    # Construct the values using amp_names_ext
    list_value_ext <- paste0(model_name, "_", amp_names_ext, amp_suffix)
    # Store in filename_list_ext
    filename_list_ext[[list_name_ext]] <- list_value_ext
  }
}

# Calculate train and test error for multiple imputation methods: classification ----
calculate_class_errors <- function(mod) {
    # Helper function to calculate error
    calculateError <- function(probField, outcomeField) {
        # Excluding null models
        mod_valid <- Filter(function(x) !is.null(x), mod)
        
        # Extracting probability vectors
        probs <- lapply(mod_valid, function(x) x[[probField]])

        # Check if all vectors are of the same length
        if (length(unique(sapply(probs, length))) != 1) {
            stop("Not all probability vectors are of the same length.")
        }

        # Calculating the mean for each corresponding element
        average_vector <- sapply(1:length(probs[[1]]), function(i) {
            mean(sapply(probs, function(x) x[i]))
        })

        # Converting to binary vector
        binary_vector <- ifelse(average_vector > 0.5, 1, 0)

        # Converting the numeric vector to a factor with labels
        factor_vector <- factor(binary_vector, levels = c(1, 0), labels = c("1", "0"))

        # Calculating error using the outcome from the first model
        error <- class_error(true = mod_valid[[1]][[outcomeField]], estimated = factor_vector, estimated.prob = average_vector)
        return(error)
    }

    # Calculate training error
    training_error <- calculateError("fitted.prob", "y.train")

    # Calculate testing error
    testing_error <- calculateError("predicted.prob", "y.test")

    # Return a list of errors
    return(list(error.train = training_error, error.test = testing_error))
}

# Calculate train and test error for multiple imputation methods: continuous ----
calculate_cont_errors <- function(mod) {
    # Helper function to calculate error
    calculateError <- function(predField, outcomeField) {
        # Excluding null models
        mod_valid <- Filter(function(x) !is.null(x), mod)
        
        # Extracting predicted vectors
        pred <- lapply(mod_valid, function(x) x[[predField]])

        # Check if all vectors are of the same length
        if (length(unique(sapply(pred, length))) != 1) {
            stop("Not all prediction vectors are of the same length.")
        }

        # Calculating the mean for each corresponding element
        average_vector <- sapply(1:length(pred[[1]]), function(i) {
            mean(sapply(pred, function(x) x[i]))
        })

        # Calculating error using the outcome from the first model
        error <- mod_error(true = mod_valid[[1]][[outcomeField]], estimated = average_vector)
        return(error)
    }

    # Calculate training error
    training_error <- calculateError("fitted", "y.train")

    # Calculate testing error
    testing_error <- calculateError("predicted", "y.test")

    # Return a list of errors
    return(list(error.train = training_error, error.test = testing_error))
}

# Function to apply the calculateError function to all mice models
mice_error_list_fun <- function(imp_suffix, model_name, dataset_suffix, n_datasets = n_datasets) {
  
  # Choose which filename_list to use based on dataset_suffix
  if (dataset_suffix == "bp") {
    list_name <- paste0(model_name, "_amp_names_", imp_suffix, "_bp")
    list_names <- filename_list_bp[[list_name]]
  } else {
    list_name <- paste0(model_name, "_amp_names_", imp_suffix, "_ext")
    list_names <- filename_list_ext[[list_name]]
  }

  # Get infolder filepath
  infolder <- infolder_list[[paste0("infolder_", model_name, "_", imp_suffix, "_", dataset_suffix)]]

  for (i in seq_along(list_names)) {
    error_out_list <- list()
    for (j in 1:n_datasets) {
      file_path <- paste0(infolder, "/", list_names[i], "_", j, ".RDS")
      
      # For tracking
      print(file_path)

      # Check if the file exists
      if (file.exists(file_path)) {

        mod <- readRDS(file_path)

        if (dataset_suffix == "bp") {
          error_out_list[[j]] <- calculate_cont_errors(mod)
        } else {
          error_out_list[[j]] <- calculate_class_errors(mod)
        }

      } else {
        # Handle the case where the file does not exist
        error_out_list[[j]] <- NULL
      }
    }
    filename_out <- paste0(infolder, "/", list_names[i], ".RDS")
    saveRDS(error_out_list, file = filename_out)
  }
}

# Call the function
for (model_name in model_names) {
  for (dataset_suffix in dataset_suffixes) {
    for (imp_suffix in suffixes) {
      mice_error_list_fun(imp_suffix, model_name, dataset_suffix, n_datasets)
    }
  }
}
