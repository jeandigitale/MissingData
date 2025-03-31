# mice RF ----
# Submits one job per dataset (e.g. 20 jobs per missingness list)

run_mice_rf <- function(dt_amp_names_tt, 
                        infolder, 
                        outfolder_mice_rf, 
                        seeds_mice_rf, 
                        log_dir, 
                        missing_indices_mice_rf = NULL) {

dt_amp_names_tt <<- dt_amp_names_tt   
infolder <<- infolder
outfolder_mice_rf <<- outfolder_mice_rf
seeds_mice_rf <<- seeds_mice_rf                       

for (i in seq_along(dt_amp_names_tt)) {
    
    i <<- i

    for (j in 1:list_j) {

      # Check if we're running only missing indices and skip this iteration if not needed
        if (!is.null(missing_indices_mice_rf) && 
            !any(missing_indices_mice_rf$i == i & missing_indices_mice_rf$j == j)) {
          next
        }

      j <<- j
      sge_opts_note <<- paste0("#$ -cwd -N '", "mice_rf_", i, "_", j, "'")

      print(paste0("mice_rf_", i, "_", j, "'"))

      sge_submit(
      { 

        starttime <- proc.time()

        list_name <<- dt_amp_names_tt[i]
        print(list_name)
        file_path <- file.path(infolder, paste0(dt_amp_names_tt[i], ".RDS"))
        tt_list <- readRDS(file_path)

        print(j)  
        tt <<- tt_list[[j]] # Test/train dataset
        print(names(tt))

        x_train <- tt$x_train
        x_test <- tt$x_test
        print(dim(x_train))

        # Create a matrix specifying the predictor relationships
        # Start with a default predictor matrix
        predictor_matrix <- make.predictorMatrix(x_train)
        print(dim(predictor_matrix))

        # Set the columns corresponding to the variables to exclude to 0
        predictor_matrix[,"pat_mrn_id"] <- 0
        predictor_matrix[,"placement_instant"] <- 0
        predictor_matrix[,"date_time_start"] <- 0
        predictor_matrix[,"date_time_end"] <- 0

        # Get seed
        seed_mice <- seeds_mice_rf[[i]][[j]]

        # Run the mice function with the custom predictor matrix on the training data
        starttime_train <- proc.time()

        mice_train <- futuremice(x_train, 
                            m = 30,
                            predictorMatrix = predictor_matrix, 
                            method="rf", 
                            parallelseed=seed_mice,
                            visitSequence = "monotone")
        endtime_train <- proc.time() 
        totaltime_train <- endtime_train - starttime_train
        print("Time to impute on train data")
        print(totaltime_train)
        print("Minutes")
        print(totaltime_train[3]/60)

        # Impute on the test data
        starttime_test <- proc.time()
        mice_test <- mice.mids(mice_train, newdata = x_test)
        endtime_test <- proc.time() 
        totaltime_test <- endtime_test - starttime_test
        print("Time to impute on test data: all")
        print(totaltime_test)
        print("Minutes")
        print(totaltime_test[3]/60)

        # Impute on one row of test data
        # starttime_test <- proc.time()
        # mice_test_single1 <- mice.mids(mice_train, newdata = x_test[3967])
        # endtime_test <- proc.time() 
        # totaltime_test <- endtime_test - starttime_test
        # print("Time to impute on test data: single row")
        # print(totaltime_test)
        # print("Minutes")
        # print(totaltime_test[3]/60)

        # starttime_test <- proc.time()
        # mice_test_single2 <- mice.mids(mice_train, newdata = x_test[29])
        # endtime_test <- proc.time() 
        # totaltime_test <- endtime_test - starttime_test
        # print("Time to impute on test data: single row")
        # print(totaltime_test)
        # print("Minutes")
        # print(totaltime_test[3]/60)

        # starttime_test <- proc.time()
        # mice_test_single3 <- mice.mids(mice_train, newdata = x_test[470])
        # endtime_test <- proc.time() 
        # totaltime_test <- endtime_test - starttime_test
        # print("Time to impute on test data: single row")
        # print(totaltime_test)
        # print("Minutes")
        # print(totaltime_test[3]/60)

        # Create list of outputs needed for model
        impute_mice_rf <- list(x_train = mice_train, 
                            x_test = mice_test, 
                            y_train = tt$y_train, 
                            y_test = tt$y_test, 
                            id_train = tt$id_train, 
                            seed = seed_mice)

      tt_path <<- file.path(outfolder_mice_rf, paste0(list_name, "_mice_rf_", j, "_raw.RDS"))
      saveRDS(impute_mice_rf, tt_path)

    endtime <- proc.time() 
    totaltime <- endtime - starttime
    print("Total time")
    print(totaltime)
    print("Minutes")
    print(totaltime[3]/60)
    },
      obj_names = c("dt_amp_names_tt", "infolder", "outfolder_mice_rf", "seeds_mice_rf", "i", "j"),
      packages = c("rtemis", "data.table", "mice"),
      n_threads = 12,
      sge_out = log_dir,
      sge_opts = sge_opts_note,
      h_rt = "200:00:00",
      system_command = "module load CBI r"
      )
  }
}
}

# Run all or based on missing indices
run_if_needed <- function(missing_indices, dt_amp_names_tt, infolder, outfolder_mice_rf, seeds_mice_rf, log_dir) {

if (exists(missing_indices, envir = .GlobalEnv)) {
    indices <- get(missing_indices, envir = .GlobalEnv)
    if (nrow(indices$mice_rf) > 0) {
            # Run the function with missing_indices_mice_rf
            run_mice_rf(dt_amp_names_tt = dt_amp_names_tt, 
                        infolder = infolder,
                        outfolder_mice_rf = outfolder_mice_rf, 
                        seeds_mice_rf = seeds_mice_rf,
                        log_dir = log_dir,
                        missing_indices_mice_rf = indices$mice_rf)
        } else {
            # Don't run the function and log a message
            print(paste("No missing files, skipping run_mice_rf for", missing_indices))
        }
    } else {
        # Run the function with missing_indices_mice_rf = NULL
        run_mice_rf(dt_amp_names_tt = dt_amp_names_tt, 
                    infolder = infolder,
                    outfolder_mice_rf = outfolder_mice_rf, 
                    seeds_mice_rf = seeds_mice_rf,
                    log_dir = log_dir,
                    missing_indices_mice_rf = NULL)
    }
}

# Call the function for Ext
run_if_needed("missing_indices_ext", dt_amp_names_tt_ext, infolder_ext, outfolder_mice_rf_ext, seeds_mice_rf_ext, log_dir_ext)

# Call the function for BP
run_if_needed("missing_indices_bp", dt_amp_names_tt_bp, infolder_bp, outfolder_mice_rf_bp, seeds_mice_rf_bp, log_dir_bp)
