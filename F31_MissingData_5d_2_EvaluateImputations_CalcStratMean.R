# F31_MissingData_5d_2_EvaluateImputations
# This file turns results calculated in 5d_1 into a data frame and creates graphs

# Set up ----
rm(list=ls(all=TRUE))

# Libraries
library(rtemis)
library(tidyverse)
library(data.table)
library(ggplot2)

# Folder
folder <- "filepath"

# File date
file_date <- "date"
today <- format(Sys.Date(), "%Y-%m-%d")

# File path in
# Define dataset suffixes and their corresponding folder names
dataset_suffixes <- c("bp", "ext")
dataset_folders <- c("BP", "Ext")

for (dataset_folder in dataset_folders) {
 
  infolder <- file.path(folder, "MissingData_Result/5a_EvaluateImputations", file_date, dataset_folder)
  
  assign(paste0("infolder_", tolower(dataset_folder)), infolder, envir = .GlobalEnv)

  rm(infolder)
}

outfolder <- file.path(folder, "MissingData_Result/5a_EvaluateImputations", file_date)

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

# Imputation suffixes
imp_suffixes <- c("mean", "locf", "mice_rf_av", "mice_pmm_av", "mice_lasso_av")

# Function to make a data frame
imp_perf_strat_fun <- function(imp_suffix, out_suffix) {

  infolder <- get(paste0("infolder_", out_suffix))

  # Read in the training and testing data
  train <- readRDS(paste0(infolder, "/imp_perf_strat_", imp_suffix, "_train_", out_suffix, ".RDS"))
  test <- readRDS(paste0(infolder, "/imp_perf_strat_", imp_suffix, "_test_", out_suffix, ".RDS"))

  # Function to process each dataset (train or test)
  process_data <- function(data, prefix) {
    # Initialize lists to store the variables
    imp_list <- list()
    miss_list <- list()
    j_list <- list()
    mse_missing_in_orig_list <- list()
    mse_not_missing_in_orig_list <- list()
    classification_error_missing_in_orig_list <- list()
    classification_error_not_missing_in_orig_list <- list()
    count_missing_in_both_numeric_list <- list()
    count_missing_in_amp_only_numeric_list <- list()
    count_missing_in_both_classification_list <- list()
    count_missing_in_amp_only_classification_list <- list()

    # Loop through the data and extract the relevant variables
    for (i in seq_along(data)) {
        miss_name <- names(data)[i]

        # Remove any occurrence of "bp" or "ext"
        miss_name <- gsub("bp_|ext_|_bp|_ext", "", miss_name)

        # Remove the last suffix (from the imputation types)
        miss_name <- sub(paste0("(_", paste(imp_suffixes, collapse = "|_"), ")$"), "", miss_name)

        for (j in seq_along(data[[i]])) {
            item <- data[[i]][[j]]

            imp_list[[length(imp_list) + 1]] <- imp_suffix
            miss_list[[length(miss_list) + 1]] <- miss_name
            j_list[[length(j_list) + 1]] <- j
            mse_missing_in_orig_list[[length(mse_missing_in_orig_list) + 1]] <- item$mse_missing_in_orig %||% NA
            mse_not_missing_in_orig_list[[length(mse_not_missing_in_orig_list) + 1]] <- item$mse_not_missing_in_orig %||% NA
            classification_error_missing_in_orig_list[[length(classification_error_missing_in_orig_list) + 1]] <- item$classification_error_missing_in_orig %||% NA
            classification_error_not_missing_in_orig_list[[length(classification_error_not_missing_in_orig_list) + 1]] <- item$classification_error_not_missing_in_orig %||% NA
            count_missing_in_both_numeric_list[[length(count_missing_in_both_numeric_list) + 1]] <- item$count_missing_in_both_numeric %||% NA
            count_missing_in_amp_only_numeric_list[[length(count_missing_in_amp_only_numeric_list) + 1]] <- item$count_missing_in_amp_only_numeric %||% NA
            count_missing_in_both_classification_list[[length(count_missing_in_both_classification_list) + 1]] <- item$count_missing_in_both_classification %||% NA
            count_missing_in_amp_only_classification_list[[length(count_missing_in_amp_only_classification_list) + 1]] <- item$count_missing_in_amp_only_classification %||% NA
        }
    }

    # Combine the lists into a data frame
    df <- data.frame(
        imp = unlist(imp_list),
        miss = unlist(miss_list),
        j = unlist(j_list),
        mse_missing_in_orig = unlist(mse_missing_in_orig_list),
        mse_not_missing_in_orig = unlist(mse_not_missing_in_orig_list),
        classification_error_missing_in_orig = unlist(classification_error_missing_in_orig_list),
        classification_error_not_missing_in_orig = unlist(classification_error_not_missing_in_orig_list),
        count_missing_in_both_numeric = unlist(count_missing_in_both_numeric_list),
        count_missing_in_amp_only_numeric = unlist(count_missing_in_amp_only_numeric_list),
        count_missing_in_both_classification = unlist(count_missing_in_both_classification_list),
        count_missing_in_amp_only_classification = unlist(count_missing_in_amp_only_classification_list),
        stringsAsFactors = FALSE
    )

    # Add prefix to column names, except for 'imp', 'miss', and 'j'
    colnames(df)[-c(1, 2, 3)] <- paste0(prefix, colnames(df)[-c(1, 2, 3)])

    return(df)
}

  # Process train and test data
  train_df <- process_data(train, "train_")
  test_df <- process_data(test, "test_")

  # Combine train and test dataframes using imp, miss, and j for merging
  final_df <- merge(train_df, test_df, by = c("imp", "miss", "j"))

  return(final_df)
}

# Apply the function to each element in imp_suffixes
list_of_dfs_ext <- lapply(imp_suffixes, imp_perf_strat_fun, out_suffix = "ext")
list_of_dfs_bp <- lapply(imp_suffixes, imp_perf_strat_fun, out_suffix = "bp")

# Combine all the data frames into one
imp_perf_strat_ext <- setDT(do.call(rbind, list_of_dfs_ext))
imp_perf_strat_bp <- setDT(do.call(rbind, list_of_dfs_bp))

# Create new columns for type and times missing (0,1,2) ----
imp_perf_strat_ext[, c("miss_type", "p_miss") := .(sub("_\\d$", "", miss), as.numeric(sub("^.*_", "", miss)))]
imp_perf_strat_ext$miss_type <- factor(imp_perf_strat_ext$miss_type, 
                             levels = c("mcar", "mar", "mnar_w", "mnar_m", "mnar_s"))
imp_perf_strat_ext$p_miss <- factor(imp_perf_strat_ext$p_miss, levels = c(0, 1, 2), labels = c("0.5", "1", "2"))
imp_perf_strat_ext$p_miss <- factor(imp_perf_strat_ext$p_miss, levels = c(0.5, 1, 2))

imp_perf_strat_bp[, c("miss_type", "p_miss") := .(sub("_\\d$", "", miss), as.numeric(sub("^.*_", "", miss)))]
imp_perf_strat_bp$miss_type <- factor(imp_perf_strat_bp$miss_type, 
                             levels = c("mcar", "mar", "mnar_w", "mnar_m", "mnar_s"))
imp_perf_strat_bp$p_miss <- factor(imp_perf_strat_bp$p_miss, levels = c(0, 1, 2), labels = c("0.5", "1", "2"))
imp_perf_strat_bp$p_miss <- factor(imp_perf_strat_bp$p_miss, levels = c(0.5, 1, 2))

# Function to reshape columns to match imp_perf data
# Need one row per missing in orig/not and one column per performance metric
# Should have double the rows of imp_perf
# Function to reshape columns
reshape_columns <- function(data, col1, col2, new_col, missing_status_values, id_cols) {
  # Subset for the first column and add the missing_status
  part1 <- data[, c(.SD, list(missing_status = missing_status_values[1])), .SDcols = c(id_cols, col1)]
  setnames(part1, col1, "value")

  # Subset for the second column and add the missing_status
  part2 <- data[, c(.SD, list(missing_status = missing_status_values[2])), .SDcols = c(id_cols, col2)]
  setnames(part2, col2, "value")

  # Combine both parts, rename the column, and remove the temporary 'value' column
  rbind(part1, part2)[, (new_col) := value][, value := NULL]
}

# ID columns
id_cols <- c("imp", "miss", "j", "miss_type", "p_miss")

# Apply function to train and test columns
## Ext
train_mse_ext <- reshape_columns(imp_perf_strat_ext,  "train_mse_missing_in_orig", 
                             "train_mse_not_missing_in_orig", "train_mse", 
                             c("missing_in_orig", "not_missing_in_orig"), id_cols)

train_classification_error_ext <- reshape_columns(imp_perf_strat_ext, 
                                              "train_classification_error_missing_in_orig", 
                                              "train_classification_error_not_missing_in_orig", 
                                              "train_classification_error", 
                                              c("missing_in_orig", "not_missing_in_orig"), id_cols)

test_mse_ext <- reshape_columns(imp_perf_strat_ext, "test_mse_missing_in_orig", 
                            "test_mse_not_missing_in_orig", "test_mse", 
                            c("missing_in_orig", "not_missing_in_orig"), id_cols)

test_classification_error_ext <- reshape_columns(imp_perf_strat_ext, 
                                             "test_classification_error_missing_in_orig", 
                                             "test_classification_error_not_missing_in_orig", 
                                             "test_classification_error", 
                                             c("missing_in_orig", "not_missing_in_orig"), id_cols)

# Apply function to train and test columns
## BP
train_mse_bp <- reshape_columns(imp_perf_strat_bp,  "train_mse_missing_in_orig", 
                             "train_mse_not_missing_in_orig", "train_mse", 
                             c("missing_in_orig", "not_missing_in_orig"), id_cols)

train_classification_error_bp <- reshape_columns(imp_perf_strat_bp, 
                                              "train_classification_error_missing_in_orig", 
                                              "train_classification_error_not_missing_in_orig", 
                                              "train_classification_error", 
                                              c("missing_in_orig", "not_missing_in_orig"), id_cols)

test_mse_bp <- reshape_columns(imp_perf_strat_bp, "test_mse_missing_in_orig", 
                            "test_mse_not_missing_in_orig", "test_mse", 
                            c("missing_in_orig", "not_missing_in_orig"), id_cols)

test_classification_error_bp <- reshape_columns(imp_perf_strat_bp, 
                                             "test_classification_error_missing_in_orig", 
                                             "test_classification_error_not_missing_in_orig", 
                                             "test_classification_error", 
                                             c("missing_in_orig", "not_missing_in_orig"), id_cols)                                             

# Update the ID columns for merge to include the 'missing_status' column
id_cols_for_merge <- c(id_cols, "missing_status")

# Merge the reshaped subsets
imp_perf_strat_long_ext <- merge(train_mse_ext, train_classification_error_ext, by = id_cols_for_merge, all = TRUE)
imp_perf_strat_long_ext <- merge(imp_perf_strat_long_ext, test_mse_ext, by = id_cols_for_merge, all = TRUE)
imp_perf_strat_long_ext <- merge(imp_perf_strat_long_ext, test_classification_error_ext, by = id_cols_for_merge, all = TRUE)

imp_perf_strat_long_bp <- merge(train_mse_bp, train_classification_error_bp, by = id_cols_for_merge, all = TRUE)
imp_perf_strat_long_bp <- merge(imp_perf_strat_long_bp, test_mse_bp, by = id_cols_for_merge, all = TRUE)
imp_perf_strat_long_bp <- merge(imp_perf_strat_long_bp, test_classification_error_bp, by = id_cols_for_merge, all = TRUE)

# Create new factor variable: imputation type and missing in orig
imp_labels_new <- c("Mean", "LOCF", "Random Forest AV", "Bayesian/PMM AV", "LASSO AV")
imp_perf_strat_long_ext[, imp := factor(imp, levels = c("mean", "locf", "mice_rf_av", "mice_pmm_av", "mice_lasso_av"), labels = imp_labels_new)]
imp_perf_strat_long_bp[, imp := factor(imp, levels = c("mean", "locf", "mice_rf_av", "mice_pmm_av", "mice_lasso_av"), labels = imp_labels_new)]

imp_perf_strat_long_ext[, imp_missing_status := paste0(imp, "_", missing_status)]
imp_perf_strat_long_bp[, imp_missing_status := paste0(imp, "_", missing_status)]

# Graph
imp_perf_strat_graph_ext <- copy(imp_perf_strat_long_ext)
imp_perf_strat_graph_bp <- copy(imp_perf_strat_long_bp)

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

imp_perf_strat_graph_ext[, miss := factor(miss, levels = miss_labels_old, labels = miss_labels_new)]
imp_perf_strat_graph_bp[, miss := factor(miss, levels = miss_labels_old, labels = miss_labels_new)]

# Use gsub to replace parts of imp_missing_status
imp_perf_strat_graph_ext[, imp_missing_status := gsub("_not_missing_in_orig$", ": no", imp_missing_status)]
imp_perf_strat_graph_ext[, imp_missing_status := gsub("_missing_in_orig$", ": yes", imp_missing_status)]

imp_perf_strat_graph_bp[, imp_missing_status := gsub("_not_missing_in_orig$", ": no", imp_missing_status)]
imp_perf_strat_graph_bp[, imp_missing_status := gsub("_missing_in_orig$", ": yes", imp_missing_status)]

# Convert imp_missing_status to a factor with new labels
imp_perf_strat_graph_ext[, imp_missing_status := factor(imp_missing_status, levels = unique(imp_missing_status))]
imp_perf_strat_graph_bp[, imp_missing_status := factor(imp_missing_status, levels = unique(imp_missing_status))]

# Define the desired order of levels
order_levels <- c("Mean: yes", "Mean: no", "LOCF: yes", "LOCF: no", 
                   "Random Forest AV: yes", "Random Forest AV: no", 
                   "Bayesian/PMM AV: yes", "Bayesian/PMM AV: no", 
                   "LASSO AV: yes", "LASSO AV: no")

# Reorder the levels of imp_missing_status
imp_perf_strat_graph_ext$imp_missing_status <- factor(imp_perf_strat_graph_ext$imp_missing_status, levels = order_levels)
imp_perf_strat_graph_bp$imp_missing_status <- factor(imp_perf_strat_graph_bp$imp_missing_status, levels = order_levels)

# Check the table to see if the factor levels are set correctly
table(imp_perf_strat_graph_ext$imp_missing_status)
table(imp_perf_strat_graph_bp$imp_missing_status)

# Color palette ----
cbPalette <- c("#E69F00", "#6b62af", "#1c8042", "#0072B2", "#D55E00", "#c6659a", "#3a86b1", "#d17d3c", "#c38aa9", "#6e98b0", "#d49565", "#caaebe")

cbPaletteStrat <- c("#f0bf55", "#807ab1", "#54a271", "#0072B2", "#D55E00", "#c6659a", "#5baad7", "#ed9653", "#dbabc5", "#96c7e3", "#f8bc8e", "#fae3f0")

custom_colors_imp_strat <- c("Mean: yes" = cbPalette[1],           
                  "Mean: no" = cbPaletteStrat[1], 
                  "LOCF: yes" = cbPalette[3], 
                  "LOCF: no" = cbPaletteStrat[3], 
                  "Random Forest AV: yes" = cbPalette[7], 
                  "Random Forest AV: no" = cbPaletteStrat[7], 
                  "Bayesian/PMM AV: yes" = cbPalette[8], 
                  "Bayesian/PMM AV: no" = cbPaletteStrat[8],
                  "LASSO AV: yes" = cbPalette[9],
                  "LASSO AV: no" = cbPaletteStrat[9]
                  )

# Create box plot function ----
# Loop through each column and create a ggplot
plot_performance_metrics <- function(data, cols_perf, title, perf_names) {
  
  # Loop using 'i' as a counter, ranging from 1 to the length of 'cols_perf'
  for (i in 1:length(cols_perf)) {
    # Get the column name using 'i'
    col <- cols_perf[i]

    p <- ggplot(data, aes(x = imp_missing_status, y = !!sym(col))) +
  geom_boxplot(aes(color = imp_missing_status), outlier.size = 0.75) +
  geom_jitter(aes(color = imp_missing_status), width = 0.2, height = 0, size = 0.5, alpha = 0.5) +
  scale_fill_manual(name = "Imputation:\nMissing in Original Data", values = custom_colors_imp_strat) +  # Manually set fill colors
  scale_color_manual(name = "Imputation:\nMissing in Original Data", values = custom_colors_imp_strat) + # Manually set color (for geom_jitter)
  facet_wrap(~ miss, ncol = 3) +
  labs(title = paste(title, "\n", perf_names[i]),
       x = "",
       y = "Value",
       fill = "Imputation:\nMissing in Original Data") +  # Set legend title with a newline
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Remove x-axis text
    axis.title.x = element_blank(), # Remove x-axis title
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.key.width = unit(1, "cm"),
    strip.background = element_rect(fill = "lightgrey", color = "white"),
    axis.text.y.left = element_text(),
    axis.ticks.y.left = element_line(),
    legend.text = element_text(size = 8),  # Adjust the size of the legend text
    legend.title = element_text(size = 9), # Adjust the size of the legend title
    legend.key.size = unit(0.5, "cm"),      # Adjust the size of the legend keys
    plot.title = element_text(hjust = 0.5, size = 12) # Center the plot title
  )

# Print the plot to the PDF device (adding a new page to the PDF)
print(p)

  }
}
# Create plots of overall imputation performance ----
perf_names <- c("Mean Squared Error", "Classification Error")

# Extubation
pdf(file.path(outfolder, "Imp_Perf_Strat_Ext.pdf"))

cols_perf_test_ext <- grep("test", names(imp_perf_strat_graph_ext), value=TRUE)
plot_performance_metrics(data = imp_perf_strat_graph_ext, 
                         cols_perf = cols_perf_test_ext, 
                         title = "Extubation Imputation Test Performance Metrics\nby Missingness in Original Data:\n",
                         perf_names = perf_names)

cols_perf_train_ext <- grep("train", names(imp_perf_strat_graph_ext), value=TRUE)
plot_performance_metrics(data = imp_perf_strat_graph_ext, 
                         cols_perf = cols_perf_train_ext, 
                         title = "Extubation Imputation Train Performance Metrics\nby Missingness in Original Data:\n",
                         perf_names = perf_names)                         

dev.off()

# BP
pdf(file.path(outfolder, "Imp_Perf_Strat_BP.pdf"))

cols_perf_test_bp <- grep("test", names(imp_perf_strat_graph_bp), value=TRUE)
plot_performance_metrics(data = imp_perf_strat_graph_bp, 
                         cols_perf = cols_perf_test_bp, 
                         title = "Blood Pressure Imputation Test Performance Metrics\nby Missingness in Original Data:\n",
                         perf_names = perf_names)

cols_perf_train_bp <- grep("train", names(imp_perf_strat_graph_bp), value=TRUE)
plot_performance_metrics(data = imp_perf_strat_graph_bp, 
                         cols_perf = cols_perf_train_bp, 
                         title = "Blood Pressure Imputation Train Performance Metrics\nby Missingness in Original Data:\n",
                         perf_names = perf_names)                         

dev.off()
