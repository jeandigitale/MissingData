# mice lasso ----
# Submits one job per dataset per missingness list

run_mice_lasso <- function(dt_amp_names_tt, 
                        infolder, 
                        outfolder_mice_lasso, 
                        seeds_mice_lasso, 
                        log_dir, 
                        missing_indices_mice_lasso = NULL) {

dt_amp_names_tt <<- dt_amp_names_tt   
infolder <<- infolder
outfolder_mice_lasso <<- outfolder_mice_lasso
seeds_mice_lasso <<- seeds_mice_lasso                       

for (i in seq_along(dt_amp_names_tt)) {
    
    i <<- i

    for (j in 1:list_j) {

      # Check if we're running only missing indices and skip this iteration if not needed
        if (!is.null(missing_indices_mice_lasso) && 
            !any(missing_indices_mice_lasso$i == i & missing_indices_mice_lasso$j == j)) {
          next
        }

      j <<- j
      sge_opts_note <<- paste0("#$ -cwd -N '", "mice_lasso_", i, "_", j, "'")

      print(paste0("mice_lasso_", i, "_", j, "'"))

      sge_submit(
      { 
        
        starttime <- proc.time()  

        list_name <<- dt_amp_names_tt[i]
        print(list_name)
        file_path <- file.path(infolder, paste0(dt_amp_names_tt[i], ".RDS"))
        tt_list <- readRDS(file_path)

        print(j)  
        tt <<- tt_list[[j]] # Test/train dataset
    
        x_train <- tt$x_train
        x_test <- tt$x_test

        # Create a matrix specifying the predictor relationships
        # Start with a default predictor matrix
        predictor_matrix <- make.predictorMatrix(x_train)

        # Set the columns corresponding to the variables to exclude to 0
        predictor_matrix[,"pat_mrn_id"] <- 0
        predictor_matrix[,"placement_instant"] <- 0
        predictor_matrix[,"date_time_start"] <- 0
        predictor_matrix[,"date_time_end"] <- 0

        # Get seed
        seed_mice <- seeds_mice_lasso[[i]][[j]]

        # Run the mice function with the custom predictor matrix on the training data
        mice_train <- futuremice(x_train, 
                            m = 30, 
                            predictorMatrix = predictor_matrix, 
                            defaultMethod = c("lasso.norm", "lasso.logreg", "cart", "cart"),
                            parallelseed = seed_mice,
                            visitSequence = "monotone")
        # defaultMethod: numeric, binary factor (1 var), unordered factor (9 vars), ordered factor (3 vars)
        # NOTE: NO binary factors missing, 3 unordered factors missing (oi_osi), 3 ordered factors missing (secretion_amount_mode)

        # Impute on the test data
        mice_test <- mice.mids(mice_train, newdata = x_test)

        # Create list of outputs needed for model
        impute_mice_lasso <- list(x_train = mice_train, 
                            x_test = mice_test, 
                            y_train = tt$y_train, 
                            y_test = tt$y_test, 
                            id_train = tt$id_train, 
                            seed = seed_mice)

      tt_path <<- file.path(outfolder_mice_lasso, paste0(list_name, "_mice_lasso_", j, "_raw.RDS"))
      saveRDS(impute_mice_lasso, tt_path)

      endtime <- proc.time() 
      totaltime <- endtime - starttime
      print(totaltime)
      print(totaltime[3]/60)
    },
      obj_names = c("dt_amp_names_tt", "infolder", "outfolder_mice_lasso", "seeds_mice_lasso", "i", "j"),
      packages = c("rtemis", "data.table", "mice"),
      n_threads = 12,
      sge_out = log_dir,
      sge_opts = sge_opts_note,
      h_rt = "100:00:00",
      system_command = "module load CBI r"
      )
  }
}
}

# Run all or based on missing indices
run_if_needed <- function(missing_indices, dt_amp_names_tt, infolder, outfolder_mice_lasso, seeds_mice_lasso, log_dir) {

if (exists(missing_indices, envir = .GlobalEnv)) {
    indices <- get(missing_indices, envir = .GlobalEnv)
    if (nrow(indices$mice_lasso) > 0) {
            # Run the function with missing_indices_mice_lasso
            run_mice_lasso(dt_amp_names_tt = dt_amp_names_tt, 
                        infolder = infolder,
                        outfolder_mice_lasso = outfolder_mice_lasso, 
                        seeds_mice_lasso = seeds_mice_lasso,
                        log_dir = log_dir,
                        missing_indices_mice_lasso = indices$mice_lasso)
        } else {
            # Don't run the function and log a message
            print(paste("No missing files, skipping run_mice_lasso for", missing_indices))
        }
    } else {
        # Run the function with missing_indices_mice_lasso = NULL
        run_mice_lasso(dt_amp_names_tt = dt_amp_names_tt, 
                    infolder = infolder,
                    outfolder_mice_lasso = outfolder_mice_lasso, 
                    seeds_mice_lasso = seeds_mice_lasso,
                    log_dir = log_dir,
                    missing_indices_mice_lasso = NULL)
    }
}

# Call the function for Ext
run_if_needed("missing_indices_ext", dt_amp_names_tt_ext, infolder_ext, outfolder_mice_lasso_ext, seeds_mice_lasso_ext, log_dir_ext)

# Call the function for BP
run_if_needed("missing_indices_bp", dt_amp_names_tt_bp, infolder_bp, outfolder_mice_lasso_bp, seeds_mice_lasso_bp, log_dir_bp)

