# LOCF ----
# Submits entire LOCF for each outcome as one job 
locf_fun <- function(data_list) {

  data_list <- lapply(seq_along(data_list), function(index) {
    
    tt <- data_list[[index]]
    
    # Print the current position and the corresponding list name
    print(paste("Processing position:", index))
    
    x_train <- tt$x_train
    x_test <- tt$x_test

    # Sort rows to ensure they are in correct order by pat_mrn_id
    setorder(x_train, pat_mrn_id, placement_instant, date_time_start) 
    setorder(x_test, pat_mrn_id, placement_instant, date_time_start) 

    # Calculate means for numeric variables in training data
    numeric_means <- x_train[, lapply(.SD, mean, na.rm = TRUE), .SDcols = sapply(x_train, is.numeric)]

    # Calculate modes for factor variables in training data
    factor_modes <- x_train[, lapply(.SD, get_mode), .SDcols = sapply(x_train, is.factor)]

    # Combine the means and modes into a single list
    train_mean_mode <- c(numeric_means, factor_modes)

    locf_data_fun <- function(data) {

      # Carry forward lagged variables on same row
      column_names <- names(data)
      lag2_names <- grep("_lag2$", column_names, value = TRUE)
      lag1_names <- grep("_lag1$", column_names, value = TRUE)
      lag0_names <- sub("_lag1$", "", lag1_names)

      # Loop through lagged variables
      # Carry forward lag2 variables to fill in lag1
      for (lag1_col in lag1_names) {
        lag2_col <- sub("_lag1$", "_lag2", lag1_col)
        
        # Check for NULL columns and print relevant info
        if (is.null(data[[lag2_col]]) || is.null(data[[lag1_col]])) {
          print(paste("Column is NULL:", lag2_col, "or", lag1_col, "at position", index))
        } else {
          data[, (lag1_col) := case_when(
            is.na(data[[lag1_col]]) ~ data[[lag2_col]],
            TRUE ~ data[[lag1_col]])]
        }
      }

      # Carry forward lag1 variables to fill in unlagged variables
      for (lag0_col in lag0_names) {
        lag1_col <- sub("_lag1$", "", lag0_col)
        
        # Check for NULL columns and print relevant info
        if (is.null(data[[lag1_col]]) || is.null(data[[lag0_col]])) {
          print(paste("Column is NULL:", lag1_col, "or", lag0_col, "at position", index))
        } else {
          data[, (lag0_col) := case_when(
            is.na(data[[lag0_col]]) ~ data[[lag1_col]],
            TRUE ~ data[[lag0_col]])]
        }
      }
      
      # Selecting columns that are not ID/time variables
      cols <- setdiff(column_names, c("pat_mrn_id", "placement_instant", "date_time_start", "date_time_end"))

      data[, (cols) := lapply(.SD, function(x) zoo::na.locf(x, na.rm = FALSE)), .SDcols = cols, by=.(pat_mrn_id, placement_instant)]

      # Number of missing cells after LOCF excluding ID/time columns
      n_cell <- nrow(data) * length(cols)
      n_cell_miss <- sum(is.na(data[, ..cols]))
      p_miss_afterLOCF <- n_cell_miss / n_cell

      ## Fill in missing data with mean (numeric) and mode (factor) of training set ----
      for (col_name in names(train_mean_mode)) {
            data[is.na(data[[col_name]]), (col_name) := train_mean_mode[[col_name]]]
      }

      output_list_data <- list(x = data, p = p_miss_afterLOCF)

      return(output_list_data)

    } 

    x_train_locf <- locf_data_fun(x_train)
    x_test_locf <- locf_data_fun(x_test)

    # Sort rows to ensure they are in descending order by pat_mrn_id to match y vectors/ID
    setorder(x_train_locf$x, pat_mrn_id, placement_instant, -date_time_start)
    setorder(x_test_locf$x, pat_mrn_id, placement_instant, -date_time_start) 

    ## Remove pat_mrn_id, placement_instant, date_time_start, date_time_end before outcome model
    x_train_locf$x[, c("pat_mrn_id", "placement_instant", "date_time_start", "date_time_end") := NULL]
    x_test_locf$x[, c("pat_mrn_id", "placement_instant", "date_time_start", "date_time_end") := NULL]

    output_list <- list(x_train = x_train_locf$x, 
                        x_test = x_test_locf$x, 
                        y_train = tt$y_train, 
                        y_test = tt$y_test, 
                        id_train = tt$id_train,
                        p_miss_afterLOCF_train = x_train_locf$p, 
                        p_miss_afterLOCF_test = x_test_locf$p)

    return(output_list)

  })

  return(data_list)

} 

# Loop through each list and apply locf_fun, saving each result as an .RDS file
# Extubation
sge_submit(
  {
    for (list_name in dt_amp_names_tt_ext) {
      print(list_name)
      assign(list_name, readRDS(file.path(infolder_ext, paste0(list_name, ".RDS"))))
      list_obj <- get(list_name)  # Get the list object by its name
      new_name <- paste0(list_name, "_locf")
      assign(new_name, locf_fun(list_obj))
      saveRDS(get(new_name), file = file.path(outfolder_locf_ext, paste0(new_name, ".RDS")))
    }
  },
  obj_names = c("dt_amp_names_tt_ext", "infolder_ext", "locf_fun", "outfolder_locf_ext"),
  packages = c("rtemis", "data.table", "tidyverse"),
  n_threads = 1,
  sge_out = log_dir_ext,
  sge_opts = "#$ -cwd -N 'locf'",
  h_rt = "10:00:00",
  system_command = "module load CBI r"
)

# Blood pressure
sge_submit(
  {
    for (list_name in dt_amp_names_tt_bp) {
      print(list_name)
      assign(list_name, readRDS(file.path(infolder_bp, paste0(list_name, ".RDS"))))
      list_obj <- get(list_name)  # Get the list object by its name
      new_name <- paste0(list_name, "_locf")
      assign(new_name, locf_fun(list_obj))
      saveRDS(get(new_name), file = file.path(outfolder_locf_bp, paste0(new_name, ".RDS")))
    }
  },
  obj_names = c("dt_amp_names_tt_bp", "infolder_bp", "locf_fun", "outfolder_locf_bp"),
  packages = c("rtemis", "data.table", "tidyverse"),
  n_threads = 1,
  sge_out = log_dir_bp,
  sge_opts = "#$ -cwd -N 'locf'",
  h_rt = "10:00:00",
  system_command = "module load CBI r"
)

