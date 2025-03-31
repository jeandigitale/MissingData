# F31_MissingData_5b_EvaluateModel
# Compare performance of models on imputed datasets with truth

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
suffixes <- c("comp", "none", "mean", "locf", "mice_rf_av",  "mice_rf",  "mice_pmm_av",  "mice_pmm", "mice_lasso_av",  "mice_lasso")

# List of suffixes for filenames
amp_suffixes <- c("", "_mean", "_locf", "_mice_rf_av",  "_mice_rf",  "_mice_pmm_av",  "_mice_pmm", "_mice_lasso_av",  "_mice_lasso")

# List of model names
model_names <- c("gbm", "glmnet")

# Number of datasets
n_datasets <- 20

# Initialize an empty list to hold the infolders
infolder_list <- list()

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

# File path out
outfolder <- file.path(folder, "MissingData_Result/5b_EvaluateModel", today)
dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)

# Create subfolders for BP and Ext
for (i in seq_along(dataset_suffixes)) {
  dataset_suffix <- dataset_suffixes[i]  # e.g., "bp" or "ext"
  dataset_folder <- dataset_folders[i]   # e.g., "BP" or "Ext"
  
  # Create a variable name based on the model and suffix
  variable_name <- paste("outfolder", dataset_suffix, sep = "_")

  # Create a folder path based on the current suffix
  outfolder_path <- file.path(outfolder, dataset_folder)
  dir.create(outfolder_path, recursive = TRUE, showWarnings = FALSE)

  # Assign the folder path to the global environment
  assign(variable_name, outfolder_path, envir = .GlobalEnv)
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

# Function to extract model performance into a list ----
mod_perf_list_fun <- function(imp_suffix, out_suffix, model_name) {
  filename_list <- get(paste0("filename_list_", out_suffix))

  # Get filename for the models
  if (imp_suffix != "none") {
    list_names <- filename_list[[paste0(model_name, "_amp_names_", imp_suffix, "_", out_suffix)]]
  }
  if (imp_suffix == "none") {
    list_names <- filename_list[[paste0(model_name, "_amp_names_", out_suffix)]]
  }

  # Get infolder filepath
  infolder <- infolder_list[[paste0("infolder_", model_name, "_", imp_suffix, "_", out_suffix)]]

  # Create an empty list to hold each combined_df
  all_model_perf <- list()

  if(out_suffix=="ext") {
    col_names <- c("Sensitivity", "Specificity", "Balanced Accuracy", "PPV", "NPV", "F1", "Accuracy", "AUC", "Log loss")
  }
  if(out_suffix=="bp"){
    col_names <- c("Mean Absolute Error","Mean Squared Error", "Root Mean Squared Error", "Normalized Root Mean Squared Error", "Pearson’s Correlation Coefficient", "p-value", "SSE", "SSR", "SST", "R-squared", "stderr")
  }

  # Loop through each of the list names (1-15 missingness scenarios)
  for (i in seq_along(list_names)) {
    list_name <- list_names[i]

    # 1) Remove the first prefix (up to the first underscore)
    amp_name <- sub("^[^_]*_", "", list_name)

    # 2) Remove any occurrence of "bp" or "ext"
    amp_name <- gsub("bp_|ext_|_bp|_ext", "", amp_name)

    # 3) Remove the last suffix (from the imputation types)
    amp_name <- sub(paste0("(_", paste(suffixes, collapse = "|_"), ")$"), "", amp_name)

    file_path <- file.path(infolder, paste0(list_names[i], ".RDS"))
    mod_list <- readRDS(file_path)

    # Loop over the models in mod_list to extract and aggregate metrics
    # Extract and rbind the dataframes together
    combined_df <- do.call(rbind, lapply(seq_along(mod_list), function(j) {
      x <- mod_list[[j]] # j is 1-20 for each missingness type
      if (is.null(x)) { # If null, insert row of NAs so the function doesn't error out
        print(paste("Model is NA for", list_name, ": ", j))
        df_train <- data.frame(matrix(ncol = length(col_names), nrow = 1))
        df_test <- data.frame(matrix(ncol = length(col_names), nrow = 1))
        names(df_train) <- paste0("train_", col_names)
        names(df_test) <- paste0("test_", col_names)
      } else { # Extract train and test errors
        if(out_suffix=="ext") {
          df_train <- x$error.train[[2]]
          df_test <- x$error.test[[2]]
        }
        if(out_suffix=="bp") {
          df_train <- x$error.train
          df_test <- x$error.test
        }
        names(df_train) <- paste0("train_", names(df_train))
        names(df_test) <- paste0("test_", names(df_test))
      }
      df <- cbind(df_train, df_test)
      df$imp <- imp_suffix # Fill in imputation type for all rows
      df$miss <- amp_name # Fill in missingness type for all rows
      setcolorder(df, c("imp", "miss"))
      return(df)
    }))

    # Append to list of all model_perf
    all_model_perf[[i]] <- combined_df
  }

  # rbind all the model_perf together after the loop: returns one dataframe per imp_suffix
  final_combined_df <- do.call(rbind, all_model_perf)
  
  return(final_combined_df)
}

# Combine model performance into one dataset ----
# Exclude complete here because it is loaded separately
imp_suffix_all <- setdiff(suffixes, "comp")
# Exclude none for models that require complete data
imp_suffix_comp <- setdiff(suffixes, c("comp", "none")) # LASSO requires complete data

# Use lapply to apply the function to each imp_suffix and store results in a list
# EXT
mod_perf_list_gbm_ext <- lapply(imp_suffix_all, mod_perf_list_fun, 
                                out_suffix="ext", model_name="gbm")
names(mod_perf_list_gbm_ext) <- imp_suffix_all

mod_perf_list_glmnet_ext <- lapply(imp_suffix_comp, mod_perf_list_fun, 
                                   out_suffix="ext", model_name="glmnet")
names(mod_perf_list_glmnet_ext) <- imp_suffix_comp

# BP
mod_perf_list_gbm_bp <- lapply(imp_suffix_all, mod_perf_list_fun, 
                               out_suffix="bp", model_name="gbm")
names(mod_perf_list_gbm_bp) <- imp_suffix_all

mod_perf_list_glmnet_bp <- lapply(imp_suffix_comp, mod_perf_list_fun, 
                                  out_suffix="bp", model_name="glmnet")
names(mod_perf_list_glmnet_bp) <- imp_suffix_comp

# Extract model performance for complete data
mod_perf_dt_fun <- function(model_name, out_suffix, imp_suffix) {
  mod_perf_list <- get(paste0("mod_perf_list_", model_name, "_", out_suffix))
  infolder_comp <- infolder_list[[paste0("infolder_", model_name, "_comp_", out_suffix)]]
  mod_comp <- readRDS(file.path(infolder_comp, paste0(model_name, "_comp_", out_suffix, ".RDS")))
  if(out_suffix=="ext") {
    df_comp_train <- mod_comp$error.train[[2]]
    df_comp_test <- mod_comp$error.test[[2]]
  }
  if(out_suffix=="bp") {
    df_comp_train <- mod_comp$error.train
    df_comp_test <- mod_comp$error.test
  }
  names(df_comp_train) <- paste0("train_", names(df_comp_train))
  names(df_comp_test) <- paste0("test_", names(df_comp_test))
  df_comp <- cbind(df_comp_train, df_comp_test)
  df_comp$imp <- "comp"
  df_comp$miss <- "comp"
  mod_perf_list[[length(mod_perf_list)+1]] <- df_comp

  # Concatenate all data frames in the list into one large data frame
  mod_perf <- as.data.table(do.call(rbind, mod_perf_list))
  names(mod_perf) <- tolower(gsub(" ", "_", names(mod_perf)))

  # Order the imputation variable
  mod_perf$imp <- factor(mod_perf$imp, levels = c(imp_suffix, "comp"))

  # Take the difference between the test and train performance ----
  # Identify the unique metrics, e.g., 'sensitivity', 'specificity', etc.
  unique_metrics <- unique(gsub("^train_|^test_", "", grep("^(train|test)_", names(mod_perf), value = TRUE)))

  # Create new variables for the differences
  for (metric in unique_metrics) {
    train_col <- paste("train_", metric, sep = "")
    test_col <- paste("test_", metric, sep = "")
    diff_col <- paste("diff_", metric, sep = "")
    
    # Calculate the difference and assign it to a new column
    mod_perf[, (diff_col) := get(train_col) - get(test_col)]
  }

  return(mod_perf)

}

# Apply the function to create dataset
mod_perf_gbm_bp <- mod_perf_dt_fun(imp_suffix=imp_suffix_all, model_name = "gbm", out_suffix = "bp")
mod_perf_glmnet_bp <- mod_perf_dt_fun(imp_suffix=imp_suffix_comp, model_name = "glmnet", out_suffix = "bp")

mod_perf_gbm_ext <- mod_perf_dt_fun(imp_suffix=imp_suffix_all, model_name = "gbm", out_suffix = "ext")
mod_perf_glmnet_ext <- mod_perf_dt_fun(imp_suffix=imp_suffix_comp, model_name = "glmnet", out_suffix = "ext")

# Create one large dataset
# Add a new column to each data.table with the model name
mod_perf_gbm_bp[, model := "gbm"]
mod_perf_glmnet_bp[, model := "glmnet"]

mod_perf_gbm_ext[, model := "gbm"]
mod_perf_glmnet_ext[, model := "glmnet"]

# Combine into one dataset
mod_perf_bp <- rbind(mod_perf_glmnet_bp, mod_perf_gbm_bp)
mod_perf_ext <- rbind(mod_perf_glmnet_ext, mod_perf_gbm_ext)

# Calculate ranges of differences from complete dataset for abstract ----
library(data.table)

calculate_metric_difference <- function(model_type, outcome, out_suffix) {
  # Select the appropriate dataset
  dataset_name <- paste0("mod_perf_", model_type, "_", out_suffix)
  mod_perf_data <- get(dataset_name, envir = .GlobalEnv) # Retrieve dataset dynamically

  # Define the column name dynamically
  outcome_col <- paste0("test_", outcome)

  # Get complete dataset value
  complete_value <- mod_perf_data[imp == "comp", get(outcome_col)]

  # Initialize a list to store results
  metric_min_max <- list()

  # Loop through each unique 'imp' value except 'comp'
  for (unique_imp in setdiff(unique(mod_perf_data$imp), "comp")) {
    # Filter and calculate the mean outcome for each imp, grouped by 'miss'
    current_agg <- mod_perf_data[imp == unique_imp & !grepl("noy", miss), 
                                 .(mean_outcome = mean(get(outcome_col))), by = miss]

    # Perform the absolute difference
    diff_vector <- complete_value - current_agg$mean_outcome

    # Compute the mean, min, and max absolute differences
    mean_val <- mean(diff_vector)
    min_val <- min(diff_vector)
    max_val <- max(diff_vector)

    # Compute percentage differences
    percent_diff_vector <- (diff_vector / complete_value) * 100
    mean_percent <- mean(percent_diff_vector)
    min_percent <- min(percent_diff_vector)
    max_percent <- max(percent_diff_vector)

    # Store the results in the list
    metric_min_max[[unique_imp]] <- list(
      mean = mean_val, min = min_val, max = max_val,
      mean_percent = mean_percent, min_percent = min_percent, max_percent = max_percent
    )
  }

  # Convert the list into a data.table
  metric_min_max_dt <- rbindlist(lapply(names(metric_min_max), function(imp) {
    data.table(
      imp = imp, 
      mean_diff = metric_min_max[[imp]]$mean, 
      min_diff = metric_min_max[[imp]]$min, 
      max_diff = metric_min_max[[imp]]$max,
      mean_percent_diff = metric_min_max[[imp]]$mean_percent, 
      min_percent_diff = metric_min_max[[imp]]$min_percent, 
      max_percent_diff = metric_min_max[[imp]]$max_percent
    )
  }))

  # Define output file name
  output_file <- paste0(outfolder, "/", model_type, "_", outcome, "_", out_suffix, "_min_max.csv")

  # Export to CSV
  fwrite(metric_min_max_dt, output_file)

  # Print the data.table to check
  print(metric_min_max_dt)

  # Return the filename for reference
  return(output_file)
}

# Example usage:
calculate_metric_difference("gbm", "balanced_accuracy", "ext")
calculate_metric_difference("gbm", "auc", "ext")
calculate_metric_difference("gbm", "mse", "bp")
calculate_metric_difference("glmnet", "balanced_accuracy", "ext")
calculate_metric_difference("glmnet", "auc", "ext")
calculate_metric_difference("glmnet", "mse", "bp")

# Coefficient of variation
cv <- function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

cv(mod_perf_gbm_ext[!grepl("noy", miss), test_balanced_accuracy])
# [1] 0.04230932
cv(mod_perf_gbm_ext[!grepl("noy", miss), test_auc])
# 0.0123901
cv(mod_perf_gbm_bp[!grepl("noy", miss), test_mse])
# [1] 0.001160136

cv(mod_perf_glmnet_ext[!grepl("noy", miss), test_balanced_accuracy])
# [1] 0.02912393
cv(mod_perf_glmnet_ext[!grepl("noy", miss), test_auc])
# [1] 0.01105491
cv(mod_perf_glmnet_bp[!grepl("noy", miss), test_mse])
# [1] 0.0007804555

# Summarize model performance by imputation ----
cols_perf_all_ext <- setdiff(names(mod_perf_ext), c("miss", "imp", "model"))
cols_perf_all_bp <- grep("_mae|_mse|_rmse|_rsq", names(mod_perf_bp), value = TRUE)
#mod_perf_imp <- mod_perf[, lapply(.SD, mean, na.rm=TRUE), by = .(imp, model), .SDcols = cols_perf_all_ext]
#mod_perf_imp_miss <- mod_perf[, lapply(.SD, mean), by = .(miss, imp, model)]

# Create new columns for type and times missing (0,1,2) ----
mod_perf_ext[, c("miss_type", "p_miss") := .(sub("_\\d$", "", miss), as.numeric(sub("^.*_", "", miss)))]
mod_perf_ext[miss=="comp", miss_type := "comp"]
mod_perf_ext$miss_type <- factor(mod_perf_ext$miss_type, 
                             levels = c("mcar", "mar", "mnar_w", "mnar_m", "mnar_s", "comp"))
mod_perf_ext$p_miss <- factor(mod_perf_ext$p_miss, levels = c(0, 1, 2), labels = c("0.5", "1", "2"))
mod_perf_ext[miss=="comp", p_miss := "0"]
mod_perf_ext$p_miss <- factor(mod_perf_ext$p_miss, levels = c(0.5, 1, 2, 0))

mod_perf_bp[, c("miss_type", "p_miss") := .(sub("_\\d$", "", miss), as.numeric(sub("^.*_", "", miss)))]
mod_perf_bp[miss=="comp", miss_type := "comp"]
mod_perf_bp$miss_type <- factor(mod_perf_bp$miss_type, 
                             levels = c("mcar", "mar", "mnar_w", "mnar_m", "mnar_s", "comp"))
mod_perf_bp$p_miss <- factor(mod_perf_bp$p_miss, levels = c(0, 1, 2), labels = c("0.5", "1", "2"))
mod_perf_bp[miss=="comp", p_miss := "0"]
mod_perf_bp$p_miss <- factor(mod_perf_bp$p_miss, levels = c(0.5, 1, 2, 0))

mod_perf_bp_imp_p_miss <- mod_perf_bp[, lapply(.SD, mean), by = .(p_miss, imp, model), .SD=cols_perf_all_bp]
write.csv(mod_perf_bp_imp_p_miss, paste0(outfolder, "/mod_perf_bp_imp_p_miss.csv"))

mod_perf_ext_imp_p_miss <- mod_perf_ext[, lapply(.SD, mean), by = .(p_miss, imp, model), .SD=cols_perf_all_ext]
write.csv(mod_perf_ext_imp_p_miss, paste0(outfolder, "/mod_perf_ext_imp_p_miss.csv"))

# Save
saveRDS(mod_perf_ext, paste0(outfolder_ext, "/mod_perf_ext.RDS"))
saveRDS(mod_perf_bp, paste0(outfolder_bp, "/mod_perf_bp.RDS"))

# Create box plots ----
cbPalette <- c("#E69F00", "#6b62af", "#1c8042", "#0072B2", "#D55E00", "#c6659a", "#3a86b1", "#d17d3c", "#c38aa9", "#6e98b0", "#d49565", "#caaebe")

# Define your custom colors
custom_colors <- c("None" = cbPalette[2], "Mean" = cbPalette[1], "LOCF" = cbPalette[3], 
                  "Random Forest MI" = cbPalette[4], "Bayesian/PMM MI" = cbPalette[5], "LASSO MI" = cbPalette[6],
                  "Random Forest AV" = cbPalette[7], "Bayesian/PMM AV" = cbPalette[8], "LASSO AV" = cbPalette[9])

# Loop through each column and create a ggplot
plot_performance_metrics <- function(data, cols_perf, title, model_name, y_limits, n_col=3) {
  
  # Prep data for graphs
  mod_perf_model <- data[model==model_name]

  # Separate the 'comp' rows and then remove them from the original data
  comp_perf <- mod_perf_model[miss == "comp"]
  mod_perf_graph <- mod_perf_model[miss != "comp"]

  # Update the labels of miss to label graphs
  miss_labels_old <- c("mcar_0", "mcar_1", "mcar_2",
                       "mar_0", "mar_1", "mar_2",
                       "mnar_w_0", "mnar_w_0_noy",
                       "mnar_w_1", "mnar_w_1_noy",
                       "mnar_w_2", "mnar_w_2_noy",
                       "mnar_m_0", "mnar_m_0_noy",
                       "mnar_m_1", "mnar_m_1_noy", 
                       "mnar_m_2", "mnar_m_2_noy",
                       "mnar_s_0", "mnar_s_0_noy", 
                       "mnar_s_1", "mnar_s_1_noy", 
                       "mnar_s_2", "mnar_s_2_noy")

  miss_labels_new <- c(
    "MCAR: 0.5x", "MCAR: 1x", "MCAR: 2x",         # for mcar_0, mcar_1, mcar_2
    "MAR: 0.5x", "MAR: 1x", "MAR: 2x",              # for mar_0, mar_1, mar_2
    "MNAR (weak): 0.5x", "MNAR no Y (weak): 0.5x",   # for mnar_w_0, mnar_w_0_noy
    "MNAR (weak): 1x", "MNAR no Y (weak): 1x",        # for mnar_w_1, mnar_w_1_noy
    "MNAR (weak): 2x", "MNAR no Y (weak): 2x",        # for mnar_w_2, mnar_w_2_noy
    "MNAR (moderate): 0.5x", "MNAR no Y (moderate): 0.5x",  # for mnar_m_0, mnar_m_0_noy
    "MNAR (moderate): 1x", "MNAR no Y (moderate): 1x",       # for mnar_m_1, mnar_m_1_noy
    "MNAR (moderate): 2x", "MNAR no Y (moderate): 2x",       # for mnar_m_2, mnar_m_2_noy
    "MNAR (strong): 0.5x", "MNAR no Y (strong): 0.5x",       # for mnar_s_0, mnar_s_0_noy
    "MNAR (strong): 1x", "MNAR no Y (strong): 1x",           # for mnar_s_1, mnar_s_1_noy
    "MNAR (strong): 2x", "MNAR no Y (strong): 2x"            # for mnar_s_2, mnar_s_2_noy
  )

  mod_perf_graph[, miss := factor(miss, levels = miss_labels_old, labels = miss_labels_new)]

  # Update the labels of imp to label graphs
  mod_perf_graph$imp <- droplevels(mod_perf_graph$imp)
  imp_labels_old <- levels(mod_perf_graph$imp)
  imp_labels_new <- c("Mean", "LOCF", "Random Forest AV", "Random Forest MI", "Bayesian/PMM AV", "Bayesian/PMM MI", "LASSO AV", "LASSO MI")
  if (model_name=="gbm") {
    imp_labels_new <- c(imp_labels_new, "None")
  }  
  mod_perf_graph[, imp := factor(imp, levels = imp_labels_old, labels = imp_labels_new)]

  # Loop using 'i' as a counter, ranging from 1 to the length of 'cols_perf'
  for (i in 1:length(cols_perf)) {
    # Get the column name using 'i'
    col <- cols_perf[i]

    # Prepare column title
    # Remove prefixes like 'test_', 'diff_', 'train_'
    col_title <- gsub("^(test_|diff_|train_)", "", col)
    # Replace underscores with spaces
    col_title <- gsub("_", " ", col_title)
    # Capitalize first letter
    if(col_title == "ppv" | col_title == "npv" | col_title == "auc" | col_title == "mae" | col_title == "mse" | col_title == "rmse") {
      col_title <- toupper(col_title)
    } else if(col_title == "log loss") {
      col_title <- "Log Loss"
    } else if(col_title == "rsq") {
      col_title <- "R-squared"
    } else {
      col_title <- tools::toTitleCase(col_title)
    }

    # Get the value of the 'comp' row for the current column
    comp_value <- comp_perf[[col]][1]

    p <-  ggplot(mod_perf_graph, aes(x = imp, y = !!sym(col))) +
          geom_boxplot(aes(color = imp), outlier.size = 0.75) +
          geom_jitter(aes(color = imp), width = 0.2, height = 0, size = 0.5, alpha = 0.5) +
          geom_hline(aes(yintercept = comp_value, linetype = "complete data"), color = "black", linewidth = 1) +
          scale_linetype_manual(values = c("complete data" = "dashed"),
                                name = "",
                                guide = guide_legend(override.aes = list(color = "black", linewidth = 1))) +
          scale_fill_manual(values = custom_colors) +  # Manually set fill colors
          scale_color_manual(values = custom_colors) + # Manually set color (for geom_jitter)
          facet_wrap(~ miss, ncol = n_col) +
          labs(title = paste0(title, col_title),
              x = "Imputation",
              y = "Value", 
              fill = "Imputation",
              color = "Imputation" 
              ) +
          theme_minimal() +
          theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.key.width = unit(1, "cm"),
            strip.background = element_rect(fill = "lightgrey", color = "white"),
            axis.text.x = element_blank()  # Updated line to remove x-axis text
          ) +
          theme(axis.text.y.left = element_text(),
                axis.ticks.y.left = element_line()) +
          ylim(y_limits[[col]][1], y_limits[[col]][2])


    # Print the plot to the PDF device (adding a new page to the PDF)
    print(p)
  }
}

# Create plots of model perfomance
cols_perf_train_ext <- grep("train", cols_perf_all_ext, value=TRUE)
cols_perf_test_ext <- grep("test", cols_perf_all_ext, value=TRUE)
cols_perf_diff_ext <- grep("diff", cols_perf_all_ext, value=TRUE)

cols_perf_train_bp <- grep("train", cols_perf_all_bp, value=TRUE)
cols_perf_test_bp <- grep("test", cols_perf_all_bp, value=TRUE)
cols_perf_diff_bp <- grep("diff", cols_perf_all_bp, value=TRUE)

# Function to calculate y-axis limits for consistent scaling
calculate_y_limits <- function(data, cols_perf) {
  y_limits <- list()
  for (col in cols_perf) {
    y_limits[[col]] <- range(data[[col]], na.rm = TRUE)
  }
  return(y_limits)
}

# Calculate y-axis limits for performance metrics
ext_y_limits <- calculate_y_limits(mod_perf_ext[!grepl("noy", miss)], cols_perf_all_ext)
bp_y_limits <- calculate_y_limits(mod_perf_bp[!grepl("noy", miss)], cols_perf_all_bp)

# GBM
## Extubation
pdf(file.path(outfolder_ext, "Mod_Perf_GBM_Ext.pdf"))

plot_performance_metrics(data = mod_perf_ext[!grepl("noy", miss)],
                         cols_perf = cols_perf_train_ext, 
                         title = "Extubation Gradient Boosted Model Train Performance:\n",
                         model_name = "gbm",
                         y_limits = ext_y_limits)

plot_performance_metrics(data = mod_perf_ext[!grepl("noy", miss)],
                         cols_perf = cols_perf_test_ext,
                         title = "Extubation Gradient Boosted Model Test Performance:\n",
                         model_name = "gbm",
                         y_limits = ext_y_limits)

plot_performance_metrics(data = mod_perf_ext[!grepl("noy", miss)],
                         cols_perf = cols_perf_diff_ext,
                         title = "Extubation Gradient Boosted Model Difference Between Train and Test:\n",
                         model_name = "gbm",
                         y_limits = ext_y_limits)                       

dev.off()

## Blood Pressure
pdf(file.path(outfolder_bp, "Mod_Perf_GBM_BP.pdf"))

plot_performance_metrics(data = mod_perf_bp[!grepl("noy", miss)],
                         cols_perf = cols_perf_train_bp, 
                         title = "Blood Pressure Gradient Boosted Model Train Performance:\n",
                         model_name = "gbm",
                         y_limits = bp_y_limits)                       

plot_performance_metrics(data = mod_perf_bp[!grepl("noy", miss)],
                         cols_perf = cols_perf_test_bp,
                         title = "Blood Pressure Gradient Boosted Model Test Performance:\n",
                         model_name = "gbm",
                         y_limits = bp_y_limits)   

plot_performance_metrics(data = mod_perf_bp[!grepl("noy", miss)],
                         cols_perf = cols_perf_diff_bp,
                         title = "Blood Pressure Gradient Boosted Model Difference Between Train and Test:\n",
                         model_name = "gbm",
                         y_limits = bp_y_limits)                            

dev.off()

# LASSO
## Extubation
pdf(file.path(outfolder_ext, "Mod_Perf_LASSO_Ext.pdf"))

plot_performance_metrics(data = mod_perf_ext[!grepl("noy", miss)],
                         cols_perf = cols_perf_train_ext, 
                          title = "Extubation LASSO Train Performance:\n",
                          model_name = "glmnet",
                         y_limits = ext_y_limits)

plot_performance_metrics(data = mod_perf_ext[!grepl("noy", miss)],
                         cols_perf = cols_perf_test_ext,
                         title = "Extubation LASSO Test Performance:\n",
                         model_name = "glmnet",
                         y_limits = ext_y_limits)

plot_performance_metrics(data = mod_perf_ext[!grepl("noy", miss)],
                         cols_perf = cols_perf_diff_ext,
                         title = "Extubation LASSO Difference Between Train and Test:\n",
                         model_name = "glmnet",
                         y_limits = ext_y_limits)                             

dev.off()

## Blood Pressure
pdf(file.path(outfolder_bp, "Mod_Perf_LASSO_BP.pdf"))

plot_performance_metrics(data = mod_perf_bp[!grepl("noy", miss)],
                         cols_perf = cols_perf_train_bp, 
                          title = "Blood Pressure LASSO Train Performance:\n",
                          model_name = "glmnet",
                         y_limits = bp_y_limits)   

plot_performance_metrics(data = mod_perf_bp[!grepl("noy", miss)],
                         cols_perf = cols_perf_test_bp,
                         title = "Blood Pressure LASSO Test Performance:\n",
                         model_name = "glmnet",
                         y_limits = bp_y_limits)   

plot_performance_metrics(data = mod_perf_bp[!grepl("noy", miss)],
                         cols_perf = cols_perf_diff_bp,
                         title = "Blood Pressure LASSO Difference Between Train and Test:\n",
                         model_name = "glmnet",
                         y_limits = bp_y_limits)                                

dev.off()

# Comparing MNAR with Y to no Y
mod_perf_ext_mnar_y <-  mod_perf_ext[grepl("mnar|comp", miss)]
mod_perf_bp_mnar_y <-  mod_perf_bp[grepl("mnar|comp", miss)]

# GBM
## Extubation
pdf(file.path(outfolder_ext, "Mod_Perf_GBM_Ext_MNAR.pdf"))

plot_performance_metrics(data = mod_perf_ext_mnar_y,
                         cols_perf = cols_perf_test_ext,
                         title = "Extubation Gradient Boosted Model Test Performance:\n",
                         model_name = "gbm",
                         n_col = 2,
                         y_limits = ext_y_limits)                

dev.off()

## Blood Pressure
pdf(file.path(outfolder_bp, "Mod_Perf_GBM_BP_MNAR.pdf"))

plot_performance_metrics(data = mod_perf_bp_mnar_y,
                         cols_perf = cols_perf_test_bp,
                         title = "Blood Pressure Gradient Boosted Model Test Performance:\n",
                         model_name = "gbm",
                         n_col = 2,
                         y_limits = bp_y_limits)   
                  
dev.off()

# LASSO
## Extubation
pdf(file.path(outfolder_ext, "Mod_Perf_LASSO_Ext_MNAR.pdf"))

plot_performance_metrics(data = mod_perf_ext_mnar_y,
                         cols_perf = cols_perf_test_ext,
                         title = "Extubation LASSO Test Performance:\n",
                         model_name = "glmnet",
                         n_col = 2,
                         y_limits = ext_y_limits)                            

dev.off()

## Blood Pressure
pdf(file.path(outfolder_bp, "Mod_Perf_LASSO_BP_MNAR.pdf"))

plot_performance_metrics(data = mod_perf_bp_mnar_y,
                         cols_perf = cols_perf_test_bp,
                         title = "Blood Pressure LASSO Test Performance:\n",
                         model_name = "glmnet",
                         n_col = 2,
                         y_limits = bp_y_limits)                                

dev.off()

## Standard deviations
sd(mod_perf_ext[!grepl("noy", miss) & model=="gbm", test_balanced_accuracy])
sd(mod_perf_ext[!grepl("noy", miss) & model=="glmnet", test_balanced_accuracy])

sd(mod_perf_ext[!grepl("noy", miss) & model=="gbm", test_auc])
sd(mod_perf_ext[!grepl("noy", miss) & model=="glmnet", test_auc])

sd(mod_perf_bp[!grepl("noy", miss) & model=="gbm", test_mse])
sd(mod_perf_bp[!grepl("noy", miss) & model=="glmnet", test_mse])

# Build level models for results
# All interactions between imputation method, missingness type, and amount missing
#mod_perf$imp <- factor(mod_perf$imp, levels = c("comp", imp_suffix_all))
#mod_perf$miss_type <- factor(mod_perf$miss_type, 
#                             levels = c("comp", "mcar", "mar", "mnar_w", "mnar_m", "mnar_s"))
#mod_perf$p_miss <- factor(mod_perf$p_miss, levels = c(0, 0.5, 1, 2))
#model <- lm(test_balanced_accuracy ~ imp * miss_type + imp * p_miss, data = mod_perf)
# Summarize the model
#summary(model)

# # All interactions between imputation method, missingness type, and amount missing
# model <- lm(test_balanced_accuracy ~ imp * miss_type + imp * p_miss, data = mod_perf[imp != "comp"])
# # Summarize the model
# summary(model)
