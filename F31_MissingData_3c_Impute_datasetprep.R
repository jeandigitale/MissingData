# Final data preparation for rtemis models
# Complete dataset: remove ID variables
# No imputation datasets: remove ID variables
# LOCF: calculate mean percent missing after LOCF, nothing needed on data
# Mean: no action needed
# mice: take one raw file per dataset post-imputation -
  # 1 - take averages, save into one file per missingness scenario (length 20) 
  # 2 - combine into one list of 30 imputations, output 20 files (one per dataset with 30 imputations) per missingness scenario; remove ID variables

# Libraries
library(mice)
library(data.table)
library(rtemis)

rm(list=ls(all=TRUE))

# Folder
folder <- "filepath"

# File date
file_date_amp <- "date"
file_date <- "date"
today <- format(Sys.Date(), "%Y-%m-%d")

# Folders
infolder_amp_bp <- file.path(folder, "MissingData_Result/2_Ampute", file_date_amp, "BP")
infolder_amp_ext <- file.path(folder, "MissingData_Result/2_Ampute", file_date_amp, "Ext")

# Define the subfolders and methods
subfolders <- c("BP", "Ext")
methods <- c("mean", "locf", "mice_rf", "mice_lasso", "mice_pmm", "comp", "none", "orig")

# Loop through subfolders and methods
for (subfolder in subfolders) {
  for (method in methods) {
    # Remove "raw" from method for variable naming
    method_clean <- gsub("/raw", "", method)
    
    # Construct the folder path
    outfolder <- file.path(folder, "MissingData_Result/3_Impute", file_date, subfolder, method)
    
    # Create the directory
    dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)
    
    # Assign the path to the global environment with a dynamic variable name
    assign(paste0("folder_", method_clean, "_", tolower(subfolder)), outfolder, envir = .GlobalEnv)
  }
}

# File path log directory
# Create log directories and assign paths
for (subfolder in subfolders) {
  log_dir <- file.path(folder, "MissingData_Log/3_Impute", today, subfolder, "Prep")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Assign the log directory path to the global environment
  assign(paste0("log_dir_", tolower(subfolder)), log_dir, envir = .GlobalEnv)
}

# Number of datasets
n_datasets <- 20

# ID vars to delete 
id_vars <- c("pat_mrn_id", "placement_instant", "date_time_start", "date_time_end")

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

# LOCF: calculate % missing after LOCF ----
dt_amp_names_tt_ext_locf <- paste0(dt_amp_names_tt_ext, "_locf")
dt_amp_names_tt_bp_locf <- paste0(dt_amp_names_tt_bp, "_locf")

process_locf_list <- function(dt_amp_names_tt_locf, folder_locf, amp_names) {
  
  locf_list <- list()
  
  for (i in seq_along(dt_amp_names_tt_locf)) {
      file_path <- file.path(folder_locf, paste0(dt_amp_names_tt_locf[i], ".RDS"))
      locf_list[[i]] <- readRDS(file_path)
  }
  names(locf_list) <- amp_names

  # Initialize vectors to store the averages for train and test separately
  averages_train <- numeric(length(locf_list))
  averages_test <- numeric(length(locf_list))

  # Loop over each missing data scenario in locf_list
  for (miss_scen_index in seq_along(locf_list)) {
    # Extract the datasets for the current missing data scenario
    datasets <- locf_list[[miss_scen_index]]
    
    # Calculate the train averages for each dataset within the current missing data scenario
    train_averages_dataset <- sapply(datasets, function(dataset) {
      if (!is.null(dataset$p_miss_afterLOCF_train)) {
        dataset$p_miss_afterLOCF_train
      } else {
        NA # Return NA if the value is missing
      }
    })
    
    # Calculate the test averages for each dataset within the current missing data scenario
    test_averages_dataset <- sapply(datasets, function(dataset) {
      if (!is.null(dataset$p_miss_afterLOCF_test)) {
        dataset$p_miss_afterLOCF_test
      } else {
        NA # Return NA if the value is missing
      }
    })
    
    # Store the mean of all dataset averages for train and test separately, excluding NAs
    averages_train[miss_scen_index] <- mean(train_averages_dataset, na.rm = TRUE)
    averages_test[miss_scen_index] <- mean(test_averages_dataset, na.rm = TRUE)
  }

  # Calculate the overall average across all missing data scenarios for train and test separately
  overall_average_train <- mean(averages_train, na.rm = TRUE)
  overall_average_test <- mean(averages_test, na.rm = TRUE)

  # Return the overall averages
  return(list(overall_average_train = overall_average_train, overall_average_test = overall_average_test))
}

locf_av_ext <- process_locf_list(dt_amp_names_tt_locf = dt_amp_names_tt_ext_locf, folder_locf = folder_locf_ext, amp_names = amp_names_ext)

locf_av_bp <- process_locf_list(dt_amp_names_tt_locf = dt_amp_names_tt_bp_locf, folder_locf = folder_locf_bp, amp_names = amp_names_bp)

# Printing the overall averages
print(paste("Ext Overall Average Train:", locf_av_ext$overall_average_train))
print(paste("Ext Overall Average Test:", locf_av_ext$overall_average_test))

print(paste("BP Overall Average Train:", locf_av_bp$overall_average_train))
print(paste("BP Overall Average Test:", locf_av_bp$overall_average_test))

# Prep non-mice data ----
prep_data <- function(suffix) {
  # Choose the correct folder based on the suffix
  infolder_amp <- if (suffix == "bp") {
    infolder_amp_bp
  } else if (suffix == "ext") {
    infolder_amp_ext
  }

  folder_comp <- if (suffix == "bp") {
    folder_comp_bp
  } else if (suffix == "ext") {
    folder_comp_ext
  }

  folder_orig <- if (suffix == "bp") {
    folder_orig_bp
  } else if (suffix == "ext") {
    folder_orig_ext
  }

  folder_none <- if (suffix == "bp") {
    folder_none_bp
  } else if (suffix == "ext") {
    folder_none_ext
  }

  dt_amp_names_tt <- if (suffix == "bp") {
    dt_amp_names_tt_bp
  } else if (suffix == "ext") {
    dt_amp_names_tt_ext
  }
  
  # Load the complete dataset
  comp <- readRDS(file.path(infolder_amp, paste0("dt_tw_comp_tt_", suffix, ".RDS")))

  # Check number of time windows per intubation in analytic dataset
  tw <- rbind(comp$x_train[, ..id_vars], comp$x_test[, ..id_vars])
  nrow <- tw[, .(Count = .N), by = .(pat_mrn_id, placement_instant)]
  summary(nrow)
  sd(nrow$Count)

  print("number of time windows per intubation in analytic dataset")
  print(suffix)
  print(nrow(tw))
  print(summary(nrow))
  print(sd(nrow$Count))

  # Set ID vars in complete data to NULL
  comp$x_train[, (id_vars) := NULL]
  comp$x_test[, (id_vars) := NULL]
  saveRDS(comp, file.path(folder_comp, paste0("dt_tw_comp_tt_", suffix, ".RDS")))

  # Set ID vars in original data to NULL so dimensions match
  orig <- readRDS(file.path(infolder_amp, paste0("dt_tw_orig_tt_", suffix, ".RDS")))
  orig$x_train[, (id_vars) := NULL]
  orig$x_test[, (id_vars) := NULL]
  saveRDS(orig, file.path(folder_orig, paste0("dt_tw_orig_tt_", suffix, ".RDS")))

  # No imputation datasets: prep data
  for (i in seq_along(dt_amp_names_tt)) {
    file_path <- file.path(infolder_amp, paste0(dt_amp_names_tt[i], ".RDS"))
    tt_list <- readRDS(file_path)
    
    for (j in seq_along(tt_list)) {
        tt_list[[j]]$x_train[, (id_vars) := NULL]
        tt_list[[j]]$x_test[, (id_vars) := NULL]
    }

    saveRDS(tt_list, file.path(folder_none, paste0(dt_amp_names_tt[i], ".RDS")))
  }
}

# Call prep_data: ext ----
sge_opts_note <<- paste0("#$ -cwd -N '", "impute_prep_ext", "'")
sge_submit(
  { 
  starttime <- proc.time() 

  prep_data("ext")

  endtime <- proc.time() 
  totaltime <- endtime - starttime
  print(totaltime)
  print(totaltime[3]/60)
  },

  obj_names = c("infolder_amp_bp", "infolder_amp_ext", "folder_comp_bp", "folder_comp_ext", "folder_orig_bp", "folder_orig_ext", "folder_none_bp", "folder_none_ext", "dt_amp_names_tt_bp", "dt_amp_names_tt_ext", "id_vars", "prep_data"),
  packages = c("rtemis", "data.table"),
  n_threads = 1,
  sge_out = log_dir_ext,
  sge_opts = sge_opts_note,
  h_rt = "60:00:00",
  system_command = "module load CBI r"
  )

# Call prep_data: bp ----
sge_opts_note <<- paste0("#$ -cwd -N '", "impute_prep_bp", "'")
sge_submit(
  { 
  starttime <- proc.time() 

  prep_data("bp")

  endtime <- proc.time() 
  totaltime <- endtime - starttime
  print(totaltime)
  print(totaltime[3]/60)
  },

  obj_names = c("infolder_amp_bp", "infolder_amp_ext", "folder_comp_bp", "folder_comp_ext", "folder_orig_bp", "folder_orig_ext", "folder_none_bp", "folder_none_ext", "dt_amp_names_tt_bp", "dt_amp_names_tt_ext", "id_vars", "prep_data"),
  packages = c("rtemis", "data.table"),
  n_threads = 1,
  sge_out = log_dir_bp,
  sge_opts = sge_opts_note,
  h_rt = "60:00:00",
  system_command = "module load CBI r"
  )

# Mice: function to prep data
mice_suffix <- c("mice_rf", "mice_pmm", "mice_lasso")

for (suffix in c("bp", "ext")) { 
  for (i in seq_along(mice_suffix)) {
    imp_suffix <- mice_suffix[i]
    folder <<- get(paste0("folder_", imp_suffix, "_", suffix))
    list_names <<- if(suffix=="bp"){
        paste0(dt_amp_names_tt_bp, "_", imp_suffix)
      } else if(suffix=="ext"){
        paste0(dt_amp_names_tt_ext, "_", imp_suffix)
      }
    
    folder_raw <<- file.path(folder, "raw")
    folder_av <<- file.path(folder, "av")
    dir.create(folder_av, recursive = TRUE, showWarnings = FALSE)
    folder_tt <<- file.path(folder, "tt")
    dir.create(folder_tt, recursive = TRUE, showWarnings = FALSE)

    for(l in seq_along(list_names)) {
      list_name <<- list_names[[l]]
      
      sge_opts_note <<- if(suffix=="bp"){
          paste0("#$ -cwd -N '", "imp_prep_bp_", imp_suffix, "_", l, "'")
        } else if(suffix=="ext"){
          paste0("#$ -cwd -N '", "imp_prep_ext_", imp_suffix, "_", l, "'")
        }

      sge_submit(
      { 
      
      # Initialize list of averages
      list_av <- list()

      file_name_av <- file.path(folder_av, paste0(list_name, "_av.RDS"))

      # Loop through each of the 20 files and read them into the list
      for(j in 1:n_datasets) {
        
          file_name <- file.path(folder_raw, paste0(list_name, "_", j, "_raw.RDS"))
          print(file_name)

          file_name_tt <- file.path(folder_tt, paste0(list_name, "_", j, ".RDS"))

          if(file.exists(file_name)) {
            mice_imp <- readRDS(file_name)
            
            # Averaged imputation data ----
            # Create complete dataset by stacking all n_imputations into one long dataset
            x_train_imp <- as.data.table(complete(mice_imp$x_train, c(1:mice_imp$x_train$m)))
            x_test_imp <- as.data.table(complete(mice_imp$x_test, c(1:mice_imp$x_test$m)))
            
            mean_mode <- function(x) {
              if (is.numeric(x)) {
                return(mean(x))
              } else if (is.factor(x) || is.character(x)) {
                return(as.factor(get_mode(x)))
              } else {
                stop("Unsupported column type")
              }
            }

            # Take the mean of numeric columns and the mode of factor columns by pat_mrn_id/placement_instant/date time start/date time end
            x_train_imp <- x_train_imp[, lapply(.SD, mean_mode), by = id_vars, .SDcols = setdiff(names(x_train_imp), id_vars)]
            x_test_imp <- x_test_imp[, lapply(.SD, mean_mode), by = id_vars, .SDcols = setdiff(names(x_test_imp), id_vars)]

            # Remove pat_mrn_id, placement_instant, date_time_start, date_time_end 
            x_train_imp[, (id_vars) := NULL]
            # Remove pat_mrn_id, placement_instant, date_time_start, date_time_end 
            x_test_imp[, (id_vars) := NULL]

            av_imp <- list(x_train_imp, x_test_imp, mice_imp$y_train, mice_imp$y_test, mice_imp$id_train)
            names(av_imp) <- c("x_train", "x_test", "y_train", "y_test", "id_train")

            # Save into list
            list_av[[j]] <- av_imp

            # Save list of imputations ----
            # Initialize list to save imputed data (length=m)
            num_imp <- mice_imp$x_train$m
            imp_list <- list()

            for (m in 1:num_imp) { # Save 1-final imputations into MI dataset
              x_train <- setDT(complete(mice_imp$x_train, m))
              x_test <- setDT(complete(mice_imp$x_test, m))

              # Remove pat_mrn_id, placement_instant, date_time_start, date_time_end 
              x_train[, (id_vars) := NULL]
              # Remove pat_mrn_id, placement_instant, date_time_start, date_time_end 
              x_test[, (id_vars) := NULL]

              imp_list[[m]] <- list(x_train, x_test, mice_imp$y_train, mice_imp$y_test, mice_imp$id_train)
              names(imp_list[[m]]) <- c("x_train", "x_test", "y_train", "y_test", "id_train")
            }

          # Save the list as an RDS file with the prefix name
          saveRDS(imp_list, file_name_tt)

          } else {
            message(paste0("File ", file_name, " does not exist"))
          }
        }
      
      saveRDS(list_av, file_name_av)

      }
        ,

          obj_names = c("folder", "folder_raw", "folder_av", "folder_tt", "amp_names", "id_vars", "n_datasets", "list_name"),
          packages = c("rtemis", "data.table", "mice"),
          n_threads = 1,
          sge_out = log_dir,
          sge_opts = sge_opts_note,
          h_rt = "60:00:00",
          system_command = "module load CBI r"
          )

      } 
    }
}
