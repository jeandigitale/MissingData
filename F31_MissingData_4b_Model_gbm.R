# GBM models
run_lightgbm <- function(tt) {

  options(future.globals.maxSize = 1 * 1024^3)

  s_LightGBM(
    x = tt$x_train, 
    y = tt$y_train,
    x.test = tt$x_test,
    y.test = tt$y_test,
    n.cores = parallelly::availableCores(),
    num_leaves = 2^c(1, 2, 3), 
    learning.rate = c(0.001, 0.01, 0.1),
    lambda_l1 = 0, # specific for LightGBM
    lambda_l2 = 0,
    grid.resample.params = setup.resample(
      resampler = "kfold", n.resamples = 5, 
      id.strat = tt$id_train
    )
  )
}

# Model for no imputation and non-MI imputed data ----
gbm_fun <- function(imp_suffix, out_suffix) {

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
  outfolder <<- get(paste0("outfolder_gbm_", imp_suffix, "_", out_suffix))
  log_dir <<- get(paste0("log_dir_", out_suffix))  # Dynamically set the log directory

  for (i in seq_along(list_names)) {
    i <<- i
    sge_opts_note <<- paste0("#$ -cwd -N '", "gbm_", imp_suffix, "_", out_suffix, "_", i, "'")

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

          mod_gbm <- run_lightgbm(tt)

          list_out[[j]] <- mod_gbm

        }
      }

      mod_path <- file.path(outfolder, paste0(gsub("_tt", "", gsub("dt", "gbm", list_name)), ".RDS"))
      saveRDS(list_out, mod_path)

      endtime <- proc.time() 
      totaltime <- endtime - starttime
      print(totaltime)
      print(totaltime[3]/60)
      },
        obj_names = c("list_names", "infolder", "outfolder", "i", "run_lightgbm"),
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
gbm_fun_mi <- function(imp_suffix, out_suffix) {
  # i: the index of the 15 missingness scenarios
  # tt_list_imp: the list of m imputations of one of the 20 datasets of one of the 15 missingness scenarios
  # j: the index of the 20 datasets within each scenario
  # m: the index of the m imputed datasets

  # Check if missing_indices_gbm exists and is relevant for the current imp_suffix
  # This will allow only missing models to re-run
  if (exists(paste0("missing_indices_gbm_", out_suffix))) {
    
      missing_indices_gbm <- get(paste0("missing_indices_gbm_", out_suffix))

      if (is.null(missing_indices_gbm[[imp_suffix]]) || nrow(missing_indices_gbm[[imp_suffix]]) == 0) {
        # If the imp_suffix is not in missing_indices_gbm or is empty, skip this imp_suffix
        return()
      } else {
        # If imp_suffix is in missing_indices_gbm and is not empty, set run_missing_indices to TRUE
        run_missing_indices <- TRUE
        missing_indices <- missing_indices_gbm[[imp_suffix]]
      }
    } else {
      # If missing_indices_gbm does not exist, run all models
      run_missing_indices <- FALSE
    }

  out_suffix <<- out_suffix
  infolder <<- get(paste0("infolder_", imp_suffix, "_", out_suffix))
  outfolder <<- get(paste0("outfolder_gbm_", imp_suffix, "_", out_suffix))
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

      sge_opts_note <<- paste0("#$ -cwd -N '", "gbm_", imp_suffix, "_", "_", out_suffix, "_", i, "_", j, "'")

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
          tt <<- tt_list_imp[[m]] # Test/train dataset

          if (is.null(tt)) { # In case imputation fails, leave null placeholder so the remaining models can run without erroring out

            list_out[[m]] <- NULL

          } else {

            mod_gbm <- run_lightgbm(tt)

            list_out[[m]] <- mod_gbm

          }
        }

        mod_path <- file.path(outfolder, paste0(gsub("_tt", "", gsub("dt", "gbm", file_name))))
        saveRDS(list_out, mod_path)

        endtime <- proc.time() 
        totaltime <- endtime - starttime
        print(totaltime)
        print(totaltime[3]/60)
        },

        obj_names = c("list_names", "infolder", "outfolder", "i", "j", "run_lightgbm", "out_suffix"),
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

if (!exists("missing_indices_gbm_bp") && !exists("missing_indices_gbm_ext")) { # Run if NOT running just models that did not run prior (missing indices)

# Model for completed ("truth") data ----
gbm_comp_fun <- function(out_suffix) {
  
  out_suffix <<- out_suffix
  infolder_comp <<- get(paste0("infolder_comp_", out_suffix))
  outfolder_gbm_comp <<- get(paste0("outfolder_gbm_comp_", out_suffix))
  log_dir <<- get(paste0("log_dir_", out_suffix))  # Dynamically set the log directory

  sge_opts_note <<- paste0("#$ -cwd -N '", "gbm_comp_", out_suffix, "'")

  sge_submit(
        { 
        starttime <- proc.time()

        tt <- readRDS(file.path(infolder_comp, paste0("dt_tw_comp_tt_", out_suffix, ".RDS"))) 
        mod_gbm <- run_lightgbm(tt)
        tt_path <<- file.path(outfolder_gbm_comp, paste0("gbm_comp_", out_suffix, ".RDS"))
        saveRDS(mod_gbm, tt_path)

        endtime <- proc.time() 
        totaltime <- endtime - starttime
        print(totaltime)
        print(totaltime[3]/60)
        },
        obj_names = c("infolder_comp", "outfolder_gbm_comp", "run_lightgbm", "out_suffix"),
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
gbm_comp_fun("ext")

set.seed(993)
gbm_comp_fun("bp")

# Run models for missing and imputed datasets ----
set.seed(3498)
gbm_fun("none", "ext")

set.seed(682)
gbm_fun("none", "bp")

set.seed(8172)
gbm_fun("locf", "ext")

set.seed(475)
gbm_fun("locf", "bp")

set.seed(3502)
gbm_fun("mean", "ext")

set.seed(9844)
gbm_fun("mean", "bp")

set.seed(357)
gbm_fun("mice_rf_av", "ext")

set.seed(6245)
gbm_fun("mice_rf_av", "bp")

set.seed(1673)
gbm_fun("mice_pmm_av", "ext")

set.seed(5046)
gbm_fun("mice_pmm_av", "bp")

set.seed(2649)
gbm_fun("mice_lasso_av", "ext")

set.seed(6940)
gbm_fun("mice_lasso_av", "bp")

}

# Models for multiply imputed datasets ----
set.seed(5931)
gbm_fun_mi("mice_rf", "ext")

set.seed(9031)
gbm_fun_mi("mice_rf", "bp")

set.seed(2457)
gbm_fun_mi("mice_pmm", "ext")

set.seed(6474)
gbm_fun_mi("mice_pmm", "bp")

set.seed(9096)
gbm_fun_mi("mice_lasso", "ext")

set.seed(5553)
gbm_fun_mi("mice_lasso", "bp")
