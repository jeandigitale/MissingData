# F31_MissingData_5a_2_EvaluateImputations
# Turn MSE/classification error into datasets (model level)
# Graph performance of MSE overall and by column

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
 
  infolder <- file.path(folder, "MissingData_Result/5a_EvaluateImputations", file_date, dataset_folder)
  dir.create(infolder, recursive = TRUE, showWarnings = FALSE)
  
  assign(paste0("infolder_", tolower(dataset_folder)), infolder, envir = .GlobalEnv)

  rm(infolder)
}

# File path out
outfolder <- file.path(folder, "MissingData_Result/5a_EvaluateImputations", file_date)
dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)

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

# Function to extract imputation performance into a dataset ----
imp_perf_dt_fun <- function(imp_suffix, out_suffix) {
  infolder_imp <- get(paste0("infolder_", out_suffix))

  # Read in the data
  imp_train_list <- readRDS(file.path(infolder_imp, paste0("imp_perf_", imp_suffix, "_train_", out_suffix, ".RDS")))
  imp_test_list <- readRDS(file.path(infolder_imp, paste0("imp_perf_", imp_suffix, "_test_", out_suffix, ".RDS")))

  # Initialize empty vectors to store the results
  imp_vec <- c()
  amp_vec <- c()
  train_mse_vec <- c()
  train_class_error_vec <- c()
  test_mse_vec <- c()
  test_class_error_vec <- c()

  # Loop through the lists and extract the relevant variables
  for(i in seq_along(imp_train_list)) {
    amp_name <- names(imp_train_list)[i]

    # Remove any occurrence of "bp" or "ext"
    amp_name <- gsub("bp_|ext_|_bp|_ext", "", amp_name)

    # Remove the last suffix (from the imputation types)
    amp_name <- sub(paste0("(_", paste(suffixes, collapse = "|_"), ")$"), "", amp_name)

    for(j in seq_along(imp_train_list[[i]])) {
      imp_vec <- c(imp_vec, imp_suffix)
      amp_vec <- c(amp_vec, amp_name)

      # Check if it's a list and has the necessary components
      if(is.list(imp_train_list[[i]][[j]]) && !is.null(imp_train_list[[i]][[j]]$mse)) {
        train_mse_vec <- c(train_mse_vec, imp_train_list[[i]][[j]]$mse)
        train_class_error_vec <- c(train_class_error_vec, imp_train_list[[i]][[j]]$classification_error)
      } else {
        train_mse_vec <- c(train_mse_vec, NA)
        train_class_error_vec <- c(train_class_error_vec, NA)
      }

      if(is.list(imp_test_list[[i]][[j]]) && !is.null(imp_test_list[[i]][[j]]$mse)) {
        test_mse_vec <- c(test_mse_vec, imp_test_list[[i]][[j]]$mse)
        test_class_error_vec <- c(test_class_error_vec, imp_test_list[[i]][[j]]$classification_error)
      } else {
        test_mse_vec <- c(test_mse_vec, NA)
        test_class_error_vec <- c(test_class_error_vec, NA)
      }
    }
  }

  # Combine all the vectors into a data frame
  final_df <- data.frame(
    imp = imp_vec,
    miss = amp_vec,
    train_mse = train_mse_vec,
    train_class_error = train_class_error_vec,
    test_mse = test_mse_vec,
    test_class_error = test_class_error_vec
  )
  
  return(final_df)
}

# Combine imputation performance into one dataset ----
imp_suffix_all <- c("mean", "locf", "mice_rf_av", "mice_pmm_av", "mice_lasso_av")

# Use lapply to apply the function to each imp_suffix and store results in a list
imp_perf_list_bp <- lapply(imp_suffix_all, imp_perf_dt_fun, out_suffix = "bp")
imp_perf_list_ext <- lapply(imp_suffix_all, imp_perf_dt_fun, out_suffix = "ext")

names(imp_perf_list_bp) <- imp_suffix_all
names(imp_perf_list_ext) <- imp_suffix_all

# Concatenate all data frames in the list into one large data frame
imp_perf_bp <- as.data.table(do.call(rbind, imp_perf_list_bp))
imp_perf_ext <- as.data.table(do.call(rbind, imp_perf_list_ext))

# Order the imputation variable
imp_perf_bp$imp <- factor(imp_perf_bp$imp, levels = imp_suffix_all)
imp_perf_ext$imp <- factor(imp_perf_ext$imp, levels = imp_suffix_all)

# Take the difference between the test and train performance ----
unique_metrics <- c("mse", "class_error")

# Create new variables for the differences
for (metric in unique_metrics) {
  train_col <- paste("train_", metric, sep = "")
  test_col <- paste("test_", metric, sep = "")
  diff_col <- paste("diff_", metric, sep = "")
  
  # Calculate the difference and assign it to a new column
  imp_perf_bp[, (diff_col) := get(train_col) - get(test_col)]
  imp_perf_ext[, (diff_col) := get(train_col) - get(test_col)]
}

# Summarize model performance by imputation ----
cols_perf_all <- setdiff(names(imp_perf_bp), c("miss", "imp"))

# Take mean by missingness and imputation type
#imp_perf_imp_miss <- imp_perf[, lapply(.SD, mean, na.rm=TRUE), by = .(miss, imp)]

# Create new columns for type and times missing (0,1,2)
imp_perf_bp[, c("miss_type", "p_miss") := .(sub("_\\d$", "", miss), as.numeric(sub("^.*_", "", miss)))]
imp_perf_bp$miss_type <- factor(imp_perf_bp$miss_type, 
                             levels = c("mcar", "mar", "mnar_w", "mnar_m", "mnar_s"))
imp_perf_bp$p_miss <- factor(imp_perf_bp$p_miss, levels = c(0, 1, 2), labels = c("0.5", "1", "2"))

imp_perf_ext[, c("miss_type", "p_miss") := .(sub("_\\d$", "", miss), as.numeric(sub("^.*_", "", miss)))]
imp_perf_ext$miss_type <- factor(imp_perf_ext$miss_type, 
                             levels = c("mcar", "mar", "mnar_w", "mnar_m", "mnar_s"))
imp_perf_ext$p_miss <- factor(imp_perf_ext$p_miss, levels = c(0, 1, 2), labels = c("0.5", "1", "2"))

# Take mean by imputation type
#imp_perf_imp <- imp_perf[, lapply(.SD, mean, na.rm=TRUE), by = .(imp), .SDcols = cols_perf_all]

# Take mean by missingness type (15 types)
#imp_perf_miss <- imp_perf[, lapply(.SD, mean, na.rm=TRUE), by = .(miss), .SDcols = cols_perf_all]

# Take mean by missingness type (5 types) 
#imp_perf_miss_type <- imp_perf[, lapply(.SD, mean, na.rm=TRUE), by = .(miss_type), .SDcols = cols_perf_all]

# Take mean by missingness type (5 types) only for mice imputations
#imp_perf_miss_type_mice <- imp_perf[grepl("mice", imp), lapply(.SD, mean, na.rm=TRUE), by = .(miss_type), .SDcols = cols_perf_all]

# Take mean by missingness type (5 types) and imputation type
#imp_perf_imp_miss_type <- imp_perf[, lapply(.SD, mean, na.rm=TRUE), by = .(imp, miss_type), .SDcols = cols_perf_all]

#imp_perf_imp_p_miss <- imp_perf[, lapply(.SD, mean), by = .(p_miss, imp), .SD=cols_perf_all]
#write.csv(imp_perf_imp_p_miss, paste0(outfolder, "/imp_perf_imp_p_miss.csv"))

# Prep data for graphs ----
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

# Update the labels of imp to label graphs
imp_labels_old <- levels(imp_perf_bp$imp)
imp_labels_new <- c("Mean", "LOCF", "Random Forest AV", "Bayesian/PMM AV", "LASSO AV")

imp_perf_bp_graph <- copy(imp_perf_bp)
imp_perf_bp_graph[, miss := factor(miss, levels = miss_labels_old, labels = miss_labels_new)]
imp_perf_bp_graph[, imp := factor(imp, levels = imp_labels_old, labels = imp_labels_new)]

imp_labels_old <- levels(imp_perf_ext$imp)
imp_labels_new <- c("Mean", "LOCF", "Random Forest AV", "Bayesian/PMM AV", "LASSO AV")

imp_perf_ext_graph <- copy(imp_perf_ext)
imp_perf_ext_graph[, miss := factor(miss, levels = miss_labels_old, labels = miss_labels_new)]
imp_perf_ext_graph[, imp := factor(imp, levels = imp_labels_old, labels = imp_labels_new)]

# Color palette ----
cbPalette <- c("#E69F00", "#6b62af", "#1c8042", "#0072B2", "#D55E00", "#c6659a", "#3a86b1", "#d17d3c", "#c38aa9", "#6e98b0", "#d49565", "#caaebe")

# Define your custom colors
custom_colors_imp <- c("Mean" = cbPalette[1], "LOCF" = cbPalette[3], 
                  "Random Forest AV" = cbPalette[7], "Bayesian/PMM AV" = cbPalette[8], "LASSO AV" = cbPalette[9])


# Create box plot function ----
# Loop through each column and create a ggplot
plot_performance_metrics <- function(data, cols_perf, title, perf_names, n_col = 3) {
  
  # Loop using 'i' as a counter, ranging from 1 to the length of 'cols_perf'
  for (i in 1:length(cols_perf)) {
    # Get the column name using 'i'
    col <- cols_perf[i]

    p <-  ggplot(data, aes(x = imp, y = !!sym(col))) +
          geom_boxplot(aes(color = imp), outlier.size = 0.75) +
          geom_jitter(aes(color = imp), width = 0.2, height = 0, size = 0.5, alpha = 0.5) +
          scale_fill_manual(values = custom_colors_imp) +  # Manually set fill colors
          scale_color_manual(values = custom_colors_imp) + # Manually set color (for geom_jitter)
          facet_wrap(~ miss, ncol = n_col) +
          labs(title = paste0(title, perf_names[i]),
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
                axis.ticks.y.left = element_line())

    # Print the plot to the PDF device (adding a new page to the PDF)
    print(p)
  }
}
# Create plots of overall imputation performance ----
# Initialize the PDF device
pdf(file.path(outfolder, "Imp_Perf_BP.pdf"))

perf_names <- c("Mean Squared Error", "Classification Error")

cols_perf_test <- grep("test", cols_perf_all, value=TRUE)
plot_performance_metrics(data = imp_perf_bp_graph[!grepl("no Y", miss)], 
                         cols_perf = cols_perf_test, 
                         title = "Blood Pressure Imputation Test Performance Metrics:\n",
                         perf_names = perf_names)

cols_perf_train <- grep("train", cols_perf_all, value=TRUE)
plot_performance_metrics(data = imp_perf_bp_graph[!grepl("no Y", miss)], 
                         cols_perf = cols_perf_train, 
                         title = "Blood Pressure Imputation Train Performance Metrics:\n",
                         perf_names = perf_names)                         

cols_perf_diff <- grep("diff", cols_perf_all, value=TRUE)
plot_performance_metrics(data = imp_perf_bp_graph[!grepl("no Y", miss)], 
                         cols_perf = cols_perf_diff, 
                         title = "Blood Pressure Imputation Performance Difference Between Train and Test:\n",
                         perf_names = perf_names)                         

# Close the PDF device
dev.off()

# Initialize the PDF device
pdf(file.path(outfolder, "Imp_Perf_Ext.pdf"))

perf_names <- c("Mean Squared Error", "Classification Error")

cols_perf_test <- grep("test", cols_perf_all, value=TRUE)
plot_performance_metrics(data = imp_perf_ext_graph[!grepl("no Y", miss)], 
                         cols_perf = cols_perf_test, 
                         title = "Extubation Imputation Test Performance Metrics:\n",
                         perf_names = perf_names)

cols_perf_train <- grep("train", cols_perf_all, value=TRUE)
plot_performance_metrics(data = imp_perf_ext_graph[!grepl("no Y", miss)], 
                         cols_perf = cols_perf_train, 
                         title = "Extubation Imputation Train Performance Metrics:\n",
                         perf_names = perf_names)                         

cols_perf_diff <- grep("diff", cols_perf_all, value=TRUE)
plot_performance_metrics(data = imp_perf_ext_graph[!grepl("no Y", miss)], 
                         cols_perf = cols_perf_diff, 
                         title = "Extubation Imputation Performance Difference Between Train and Test:\n",
                         perf_names = perf_names)                         

# Close the PDF device
dev.off()

# Initialize the PDF device
pdf(file.path(outfolder, "Imp_Perf_BP_MNAR.pdf"))

perf_names <- c("Mean Squared Error", "Classification Error")

cols_perf_test <- grep("test", cols_perf_all, value=TRUE)
plot_performance_metrics(data = imp_perf_bp_graph[grepl("MNAR", miss)], 
                         cols_perf = cols_perf_test, 
                         title = "Blood Pressure Imputation Test Performance Metrics:\n",
                         perf_names = perf_names,
                         n_col = 2)

# Close the PDF device
dev.off()

# Initialize the PDF device
pdf(file.path(outfolder, "Imp_Perf_Ext_MNAR.pdf"))

perf_names <- c("Mean Squared Error", "Classification Error")

cols_perf_test <- grep("test", cols_perf_all, value=TRUE)
plot_performance_metrics(data = imp_perf_ext_graph[grepl("MNAR", miss)], 
                         cols_perf = cols_perf_test, 
                         title = "Extubation Imputation Test Performance Metrics:\n",
                         perf_names = perf_names,
                         n_col = 2)

# Close the PDF device
dev.off()

# Create function to calculate difference between mean MSE and other method MSE
calc_diff <- function(out_suffix) {
  # Select the correct dataset based on out_suffix
  if (out_suffix == "ext") {
    imp_perf <- imp_perf_ext
  } else if (out_suffix == "bp") {
    imp_perf <- imp_perf_bp
  } else {
    stop("Invalid out_suffix. Choose 'bp' or 'ext'.")
  }

  min_max_results <- list()

  # Subset where imp == "Mean" and aggregate by miss
  mean_agg <- imp_perf[imp == "mean" & !grepl("noy", miss), .(mean_test_mse = mean(test_mse)), by = miss]

  # Loop through each unique imp value except "Mean"
  for (unique_imp in unique(imp_perf$imp)) {
    if (unique_imp != "mean") {
      # Subset by current imp value and aggregate by miss
      current_agg <- imp_perf[imp == unique_imp & !grepl("noy", miss), .(mean_test_mse = mean(test_mse)), by = miss]
      
      # Merge the mean and current aggregated data on miss
      merged_data <- merge(mean_agg, current_agg, by = "miss", suffixes = c("_mean", paste0("_", unique_imp)))
      
      # Rename columns for clarity
      setnames(merged_data, old = paste0("mean_test_mse_", unique_imp), new = "mean_test_mse_current")
      
      # Create a new column for the difference
      merged_data[, diff := mean_test_mse_mean - mean_test_mse_current]
      
      # Compute the mean, min, and max of the difference
      mean_val <- mean(merged_data$diff)
      min_val <- min(merged_data$diff)
      max_val <- max(merged_data$diff)
      
      # Store the results in the list
      min_max_results[[unique_imp]] <- list(mean = mean_val, min = min_val, max = max_val)
    }
  }

  # Convert the list into a data.table
  min_max_results_dt <- rbindlist(lapply(names(min_max_results), function(imp) {
    data.table(imp = imp, 
               mean_diff = min_max_results[[imp]]$mean, 
               min_diff = min_max_results[[imp]]$min, 
               max_diff = min_max_results[[imp]]$max)
  }))

  # Define output file name
  output_file <- paste0(outfolder, "/min_max_diff_", out_suffix, ".csv")

  # Export to CSV
  fwrite(min_max_results_dt, output_file)

  # Print the data.table to check
  print(min_max_results_dt)

}

# Calculate the results
calc_diff("ext")
calc_diff("bp")
