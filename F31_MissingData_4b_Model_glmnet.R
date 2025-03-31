# GLMNET models
run_glmnet <- function(tt) {

  options(future.globals.maxSize = 1 * 1024^3)

  s_GLMNET(
            x = tt$x_train, 
            y = tt$y_train,
            x.test = tt$x_test,
            y.test = tt$y_test,
            n.cores = parallelly::availableCores(),
            alpha = 1,
            grid.resample.params = setup.resample(
              resampler = "kfold", n.resamples = 5, 
              id.strat = tt$id_train
            )) # check the grid search results to see tuned lambda, automatically tuning
}

# Model for no imputation and non-MI imputed data ----
glmnet_fun <- function(imp_suffix, out_suffix) {

  # i: the index of the 15 missingness scenarios
  # tt_list: the list of 20 datasets of one of the 15 missingness scenarios
  # j: the index of the 20 datasets within each scenario

  if(imp_suffix != "none") {
    list_names <<- get(paste0("dt_amp_names_tt_", out_suffix, "_",imp_suffix))
  }
  if(imp_suffix == "none") {
    list_names <<- get(paste0("dt_amp_names_tt_", out_suffix))
  }

  infolder <<- get(paste0("infolder_", imp_suffix, "_", out_suffix))
  outfolder <<- get(paste0("outfolder_glmnet_", imp_suffix, "_", out_suffix))
  log_dir <<- get(paste0("log_dir_", out_suffix))  # Dynamically set the log directory

  for (i in seq_along(list_names)) {
    i <<- i
    sge_opts_note <<- paste0("#$ -cwd -N '", "glmnet_", imp_suffix, "_", out_suffix, "_", i, "'")

    sge_submit(
      { 

      starttime <- proc.time()

      list_name <- list_names[i]
      print(list_name)
      file_path <- file.path(infolder, paste0(list_names[i], ".RDS"))
      tt_list <- readRDS(file_path)
        
      list_out <- list()

      for (j in seq_along(tt_list)) {
        print(j)  
        tt <<- tt_list[[j]] # Test/train dataset

        if (is.null(tt)) { # In case imputation fails, leave null placeholder so the remaining models can run without erroring out

          list_out[[j]] <- NULL

        } else {

          # Special case for i=12 and j=18: Fill NAs in wbc_mean_lag2 with 0
          if (i == 12 && j == 18) {
            message("Replacing NA values in wbc_mean_lag2 with 0 for i=12, j=18")
            tt$x_train[is.na(tt$x_train$wbc_mean_lag2), wbc_mean_lag2 := 0]
          }

          mod_glmnet <- run_glmnet(tt)

          list_out[[j]] <- mod_glmnet

        }
      }

      mod_path <- file.path(outfolder, paste0(gsub("_tt", "", gsub("dt", "glmnet", list_name)), ".RDS"))
      saveRDS(list_out, mod_path)

      endtime <- proc.time() 
      totaltime <- endtime - starttime
      print(totaltime)
      print(totaltime[3]/60)
      },
        obj_names = c("list_names", "infolder", "outfolder", "i", "run_glmnet"),
        packages = c("rtemis", "data.table"),
        n_threads = 12,
        sge_out = log_dir,
        sge_opts = sge_opts_note,
        h_rt = "60:00:00",
        system_command = "module load CBI r"
        )
  }
}

# Model for MI imputed data ----
glmnet_fun_mi <- function(imp_suffix, out_suffix) {
  # i: the index of the 15 missingness scenarios
  # tt_list_imp: the list of m imputations of one of the 20 datasets of one of the 15 missingness scenarios
  # j: the index of the 20 datasets within each scenario
  # m: the index of the m imputed datasets

  # Check if missing_indices_glmnet exists and is relevant for the current imp_suffix
  # This will allow only missing models to re-run
  if (exists(paste0("missing_indices_glmnet_", out_suffix))) {
    
      missing_indices_glmnet <- get(paste0("missing_indices_glmnet_", out_suffix))

      if (is.null(missing_indices_glmnet[[imp_suffix]]) || nrow(missing_indices_glmnet[[imp_suffix]]) == 0) {
        # If the imp_suffix is not in missing_indices_glmnet or is empty, skip this imp_suffix
        return()
      } else {
        # If imp_suffix is in missing_indices_glmnet and is not empty, set run_missing_indices to TRUE
        run_missing_indices <- TRUE
        missing_indices <- missing_indices_glmnet[[imp_suffix]]
      }
    } else {
      # If missing_indices_glmnet does not exist, run all models
      run_missing_indices <- FALSE
    }

  out_suffix <<- out_suffix
  infolder <<- get(paste0("infolder_", imp_suffix, "_", out_suffix))
  outfolder <<- get(paste0("outfolder_glmnet_", imp_suffix, "_", out_suffix))
  log_dir <<- get(paste0("log_dir_", out_suffix))  # Dynamically set the log directory
  list_names <<- get(paste0("dt_amp_names_tt_", out_suffix, "_",imp_suffix))

  for (i in seq_along(list_names)) {
    i <<- i

    for (j in 1:n_datasets) {
        j <<- j

      # Skip this iteration if missing_indices exist and this (i, j) is not in it
      if (run_missing_indices && !any(missing_indices$i == i & missing_indices$j == j)) {
        next
      }

      sge_opts_note <<- paste0("#$ -cwd -N '", "glmnet_", imp_suffix, "_", "_", out_suffix, "_", i, "_", j, "'")

      sge_submit(
        { 

        starttime <- proc.time()
        list_name <- list_names[i]
        file_name <- paste0(list_name, "_", j, ".RDS")
        print(file_name)
        file_path <- file.path(infolder, file_name)
        tt_list_imp <- readRDS(file_path)
          
        list_out <- list()

        for (m in seq_along(tt_list_imp)) {  

          print(m)
          tt <<- tt_list_imp[[m]] # Test/train dataset

          if (is.null(tt)) { # In case imputation fails, leave null placeholder so the remaining models can run without erroring out

            list_out[[m]] <- NULL

          } else {

            # Special case for i=12 and j=18: Fill NAs in wbc_mean_lag2 with 0
            if (i == 12 && j == 18) {
              message("Replacing NA values in wbc_mean_lag2 with 0 for i=12, j=18, m=", m)
              tt$x_train[is.na(tt$x_train$wbc_mean_lag2), wbc_mean_lag2 := 0]
            }

            mod_glmnet <- run_glmnet(tt)

            list_out[[m]] <- mod_glmnet

          }
        }

        mod_path <- file.path(outfolder, paste0(gsub("_tt", "", gsub("dt", "glmnet", file_name))))
        saveRDS(list_out, mod_path)

        endtime <- proc.time() 
        totaltime <- endtime - starttime
        print(totaltime)
        print(totaltime[3]/60)
        },

        obj_names = c("list_names", "infolder", "outfolder", "i", "j", "run_glmnet", "out_suffix"),
        packages = c("rtemis", "data.table"),
        n_threads = 12,
        sge_out = log_dir,
        sge_opts = sge_opts_note,
        h_rt = "60:00:00",
        system_command = "module load CBI r"
        )
    }
  }
}

if (!exists("missing_indices_glmnet_bp") && !exists("missing_indices_glmnet_ext")){ # Run if NOT running just models that did not run prior (missing indices)

# Model for completed ("truth") data ----
glmnet_comp_fun <- function(out_suffix) {
  
  out_suffix <<- out_suffix
  infolder_comp <<- get(paste0("infolder_comp_", out_suffix))
  outfolder_glmnet_comp <<- get(paste0("outfolder_glmnet_comp_", out_suffix))
  log_dir <<- get(paste0("log_dir_", out_suffix))  # Dynamically set the log directory

  sge_opts_note <<- paste0("#$ -cwd -N '", "glmnet_comp_", out_suffix, "'")

  sge_submit(
        { 
        starttime <- proc.time()

        tt <- readRDS(file.path(infolder_comp, paste0("dt_tw_comp_tt_", out_suffix, ".RDS"))) 
        mod_glmnet <- run_glmnet(tt)
        tt_path <<- file.path(outfolder_glmnet_comp, paste0("glmnet_comp_", out_suffix, ".RDS"))
        saveRDS(mod_glmnet, tt_path)

        endtime <- proc.time() 
        totaltime <- endtime - starttime
        print(totaltime)
        print(totaltime[3]/60)
        },
        obj_names = c("infolder_comp", "outfolder_glmnet_comp", "run_glmnet", "out_suffix"),
        packages = c("rtemis", "data.table"),
        n_threads = 12,
        sge_out = log_dir,
        sge_opts = sge_opts_note,
        h_rt = "60:00:00",
        system_command = "module load CBI r"
        )
}

# Run models for completed dataset ----
set.seed(2649)
glmnet_comp_fun("ext")

set.seed(993)
glmnet_comp_fun("bp")

# Run models for missing and imputed datasets ----
set.seed(3498)
glmnet_fun("none", "ext")

set.seed(682)
glmnet_fun("none", "bp")

set.seed(8172)
glmnet_fun("locf", "ext")

set.seed(475)
glmnet_fun("locf", "bp")

set.seed(3502)
glmnet_fun("mean", "ext")

set.seed(9844)
glmnet_fun("mean", "bp")

set.seed(357)
glmnet_fun("mice_rf_av", "ext")

set.seed(6245)
glmnet_fun("mice_rf_av", "bp")

set.seed(1673)
glmnet_fun("mice_pmm_av", "ext")

set.seed(5046)
glmnet_fun("mice_pmm_av", "bp")

set.seed(2649)
glmnet_fun("mice_lasso_av", "ext")

set.seed(6940)
glmnet_fun("mice_lasso_av", "bp")

}

# Models for multiply imputed datasets ----
set.seed(5931)
glmnet_fun_mi("mice_rf", "ext")

set.seed(9031)
glmnet_fun_mi("mice_rf", "bp")

set.seed(2457)
glmnet_fun_mi("mice_pmm", "ext")

set.seed(6474)
glmnet_fun_mi("mice_pmm", "bp")

set.seed(9096)
glmnet_fun_mi("mice_lasso", "ext")

set.seed(5553)
glmnet_fun_mi("mice_lasso", "bp")

###########################################################################

# Added 2-14-2025 because of multicollinearity 
# For i=12, j=18 in the BP outcome datasets, wbc_mean_lag2 did not impute because of multicollinearity

# Demonstrate collinearity in amputed dataset
i=12
j=18

list_name <<- dt_amp_names_tt_bp[i]
print(list_name)
file_path <- "filepath/MissingData_Result/2_Ampute/2024-09-22/BP/dt_mnar_m_2_tt_bp.RDS"
tt_list <- readRDS(file_path)

tt <<- tt_list[[j]] # Test/train dataset

x_train <- tt$x_train
mice:::find.collinear(x_train)

# Identify numeric columns
numeric_mask <- sapply(x_train, is.numeric)

# Select only numeric columns using data.table syntax
numeric_vars <- x_train[, ..numeric_mask]

# Compute correlation matrix
cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")

# Extract correlations of 'wbc_mean_lag2'
cor_wbc <- cor_matrix["wbc_mean_lag2", ]

# Variables highly correlated are the lagged vars: wbc_mean, wbc_mean_lag1
# Find variables highly correlated with 'wbc_mean_lag2' (|correlation| > 0.8)
names(cor_wbc[abs(cor_wbc) > 0.8 & names(cor_wbc) != "wbc_mean_lag2"])

# Demonstrate that in every other GLMNET model for glmnet_mnar_m_2_bp_mice_rf, wbc_mean_lag2 has a coefficient of 0
# Define the base path and file pattern
base_path <- "filepathMissingData_Result/4_Model/2024-12-17/BP/glmnet/mice_rf/"

# File indices to loop through (excluding 18)
file_indices <- setdiff(1:20, 18)

# Variable to check
target_var <- "wbc_mean_lag2"

# Initialize a flag to track if wbc_mean_lag2 is ever non-zero
is_nonzero <- FALSE

# Loop through the files
for (i in file_indices) {
  # Construct file name
  file_name <- paste0(base_path, "glmnet_mnar_m_2_bp_mice_rf_", i, ".RDS")
  
  # Read the RDS file
  test <- readRDS(file_name)
  
  # Loop through the 30 objects inside
  for (j in 1:30) {
    # Extract varimp
    varimp <- test[[j]]$varimp
    
    # Check if wbc_mean_lag2 exists and is non-zero
    if (target_var %in% names(varimp) && varimp[target_var] != 0) {
      is_nonzero <- TRUE
      message("wbc_mean_lag2 is non-zero in file: ", file_name, ", object: ", j)
    }
  }
}

# Final result
if (!is_nonzero) {
  message("wbc_mean_lag2 is zero in all files and objects.")
}

# Solution: inserted code above to make wbc_mean_lag2 zero only for i=12 and j=18 (which mice did not impute). It will drop out and not be included in the final models.
