# Mean/mode ----
# Submits entire mean/mode for each outcome as one job 

mm_fun <- function(data_list) {

  data_list <- lapply(data_list, function(tt) {
  
    x_train <- tt$x_train
    x_test <- tt$x_test

    # Calculate means for numeric variables in training data
    numeric_means <- x_train[, lapply(.SD, mean, na.rm = TRUE), .SDcols = sapply(x_train, is.numeric)]

    # Calculate modes for factor variables in training data
    factor_modes <- x_train[, lapply(.SD, get_mode), .SDcols = sapply(x_train, is.factor)]

    # Combine the means and modes into a single list
    train_mean_mode <- c(numeric_means, factor_modes)

    ## Fill in missing data with mean (numeric) and mode (factor) of training set ----
    for (col_name in names(train_mean_mode)) {
        x_train[is.na(x_train[[col_name]]), (col_name) := train_mean_mode[[col_name]]]
        x_test[is.na(x_test[[col_name]]), (col_name) := train_mean_mode[[col_name]]]
    }

    ## Remove pat_mrn_id, placement_instant, date_time_start, date_time_end before outcome model
    x_train[, c("pat_mrn_id", "placement_instant", "date_time_start", "date_time_end") := NULL]
    x_test[, c("pat_mrn_id", "placement_instant", "date_time_start", "date_time_end") := NULL]

    output_list <- list(x_train = x_train, 
                        x_test = x_test, 
                        y_train = tt$y_train, 
                        y_test = tt$y_test,
                        id_train = tt$id_train)

    return(output_list)

  })

  return(data_list)

} 

# Loop through each list and apply mm_fun, saving each result as an .RDS file
# Extubation
sge_submit(
  {
    for (list_name in dt_amp_names_tt_ext) {
      assign(list_name, readRDS(file.path(infolder_ext, paste0(list_name, ".RDS"))))
      list_obj <- get(list_name)  # Get the list object by its name
      new_name <- paste0(list_name, "_mean")
      assign(new_name, mm_fun(list_obj))
      saveRDS(get(new_name), file = file.path(outfolder_mean_ext, paste0(new_name, ".RDS")))
    }
  },
  obj_names = c("dt_amp_names_tt_ext", "infolder_ext", "mm_fun", "outfolder_mean_ext"),
  packages = c("rtemis", "data.table", "tidyverse"),
  n_threads = 1,
  sge_out = log_dir_ext,
  sge_opts = "#$ -cwd -N 'mean'",
  h_rt = "10:00:00",
  system_command = "module load CBI r"
)

# Blood pressure
sge_submit(
  {
    for (list_name in dt_amp_names_tt_bp) {
      assign(list_name, readRDS(file.path(infolder_bp, paste0(list_name, ".RDS"))))
      list_obj <- get(list_name)  # Get the list object by its name
      new_name <- paste0(list_name, "_mean")
      assign(new_name, mm_fun(list_obj))
      saveRDS(get(new_name), file = file.path(outfolder_mean_bp, paste0(new_name, ".RDS")))
    }
  },
  obj_names = c("dt_amp_names_tt_bp", "infolder_bp", "mm_fun", "outfolder_mean_bp"),
  packages = c("rtemis", "data.table", "tidyverse"),
  n_threads = 1,
  sge_out = log_dir_bp,
  sge_opts = "#$ -cwd -N 'mean'",
  h_rt = "10:00:00",
  system_command = "module load CBI r"
)

