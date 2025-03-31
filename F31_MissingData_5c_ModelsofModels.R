# F31_MissingData_5c_ModelofModels
# Models of Models to evaluate imputation and outcome model performance ----
library(data.table)
library(lme4)
library(lmerTest)
library(broom)
library(broom.mixed)
library(emmeans)
library(ggplot2)
library(openxlsx)
rm(list=ls(all=TRUE))

emm_options(pbkrtest.limit = 300000)
emm_options(lmerTest.limit = 300000)

# Folder
folder <- "filepath"

# File date
today <- format(Sys.Date(), "%Y-%m-%d")
file_date_col <- "date"
file_date_mod <- "date"

# File path out
outfolder <- file.path(folder, "MissingData_Result/5c_ModelofModels", today)
dir.create(outfolder, recursive = TRUE, showWarnings = FALSE)

# Create a new workbook
wb <- createWorkbook()

# File paths
col_perf_filepath_ext <- file.path(folder, "MissingData_Result", "5a_EvaluateImputations", file_date_col, "Ext", "col_perf_ext.RDS")
col_perf_filepath_bp <- file.path(folder, "MissingData_Result", "5a_EvaluateImputations", file_date_col, "BP", "col_perf_bp.RDS")

mod_perf_filepath_ext <- file.path(folder, "MissingData_Result", "5b_EvaluateModel",file_date_mod, "Ext", "mod_perf_ext.RDS")
mod_perf_filepath_bp <- file.path(folder, "MissingData_Result", "5b_EvaluateModel",file_date_mod, "BP", "mod_perf_bp.RDS")

# Stepwise interaction function ----
stepwise_interaction_fun <- function(data, full_formula, interactions, model_type, threshold = 0.05) {
  # Fit the full model based on the specified model type
  if (model_type == "lmer") {
    full_model <- lmer(as.formula(full_formula), data = data, REML = FALSE)
  } else if (model_type == "lm") {
    full_model <- lm(as.formula(full_formula), data = data)
  } else {
    stop("Invalid model type. Choose 'lmer' or 'lm'.")
  }
  
  print(summary(full_model))
  current_model <- full_model
  final_formula <- full_formula
  results <- list()
  
  # Loop through each interaction
  for (interaction in interactions) {
    # Create the formula for the reduced model (without the current interaction)
    reduced_formula <- gsub(paste("\\+", interaction), "", final_formula)
    reduced_formula <- gsub(paste(interaction, "\\+"), "", reduced_formula)
    reduced_formula <- gsub("\\++\\+", "+", reduced_formula)
    
    # Fit the reduced model
    if (model_type == "lmer") {
      reduced_model <- lmer(as.formula(reduced_formula), data = data, REML = FALSE)
    } else {
      reduced_model <- lm(as.formula(reduced_formula), data = data)
    }
    
    print(summary(reduced_model))
    
    # Conduct ANOVA between the current model and the reduced model
    anova_result <- anova(current_model, reduced_model)
    print(final_formula)
    print(reduced_formula)
    print(anova_result)
    
    # Extract the p-value based on model type
    p_value <- if (model_type == "lmer") {
      anova_result$"Pr(>Chisq)"[2]
    } else {
      tail(anova_result$"Pr(>F)", 1)
    }
    
    # Handle NA or missing p-values
    if (is.na(p_value) || length(p_value) == 0) {
      results[[interaction]] <- "Kept; p-value = NA"
      # Note if p-value is null, the reduced model is equal to the current model because the reduced model discarded a two-way interaction when a three-way interaction that subsumes it is still included; thus, the two-way interaction should remain in the model
    } else {
      # Decide to keep or drop based on the p-value
      if (p_value > threshold) {
        final_formula <- reduced_formula
        current_model <- reduced_model
        results[[interaction]] <- paste("Dropped; p-value =", p_value)
      } else {
        results[[interaction]] <- paste("Kept; p-value =", p_value)
      }
    }
  }
  
  return(list("Final Formula" = final_formula, "Results" = results))
}

# Function to run stepwise interaction and store results ----
run_stepwise_fun <- function(datasets, full_formula, all_interactions, model_type, outcome=NULL) {
  # Initialize lists to store results
  stepwise_results <- list()
  final_models <- list()
  final_summaries <- list()
  final_tidy_models <- list()
  
  for (dataset_name in names(datasets)) {
    # Modify sheet names based on the presence of 'outcome'
    sheet_name_suffix <- if (!is.null(outcome)) {
      paste(outcome, dataset_name, sep = "_")
    } else {
      dataset_name
    }

    # Run the stepwise interaction test
    results <- stepwise_interaction_fun(data = datasets[[dataset_name]], full_formula, all_interactions, model_type = model_type)
    stepwise_results[[dataset_name]] <- results$Results
    
    # Fit the final model based on specified model type
    if (model_type == "lmer") {
      final_model <- lmer(as.formula(results$`Final Formula`), data = datasets[[dataset_name]], REML = FALSE)
    } else if (model_type == "lm") {
      final_model <- lm(as.formula(results$`Final Formula`), data = datasets[[dataset_name]])
    } else {
      stop("Invalid model type. Use 'lm' or 'lmer'.")
    }
    
    # Store the final model, its summary, and tidy results
    final_models[[dataset_name]] <- final_model
    final_summaries[[dataset_name]] <- summary(final_model)
    if (model_type == "lmer") {
      final_tidy_models[[dataset_name]] <- broom.mixed::tidy(final_model)
    } else if (model_type == "lm") {
      final_tidy_models[[dataset_name]] <- broom::tidy(final_model)
    }
    
    # Add a sheet with tidy model results to the workbook
    addWorksheet(wb, paste0(sheet_name_suffix, "_m"))
    writeData(wb, paste0(sheet_name_suffix, "_m"), final_tidy_models[[dataset_name]])
    
    # Convert stepwise_results to a matrix and include row names
    stepwise_matrix <- as.matrix(stepwise_results[[dataset_name]])
    addWorksheet(wb, paste0(sheet_name_suffix, "_p"))
    writeData(wb, paste0(sheet_name_suffix, "_p"), stepwise_matrix, withFilter = TRUE, rowNames = TRUE)
  }
  
  # Combine and return the final results as a list
  return(list(stepwise_results = stepwise_results,
              final_models = final_models, 
              final_summaries = final_summaries, 
              final_tidy_models = final_tidy_models))
}

# Imputation performance results ----
# Multilevel model
col_perf_ext <- as.data.table(readRDS(col_perf_filepath_ext))
col_perf_bp <- as.data.table(readRDS(col_perf_filepath_bp))

# Exclude noy models
col_perf_ext <- col_perf_ext[!grepl("noy", miss)]
col_perf_bp <- col_perf_bp[!grepl("noy", miss)]

# Function for col ext and bp ----
stepwise_col_fun <- function(out_suffix) {

  col_perf <- get(paste0("col_perf_", out_suffix))

  # Create group variable
  # Define groups
  # 1: Always observed
  #group_1 <- grep("ext_success|age_int_days|alternative_gases|bp|cuffed|dx|etco2|ett_cuff_leak|ett_placed|ett_size|fio2|hrs_int|intake_output|map|med_|present_on_hosp|pulse|resp_|sbt|sex|spo2|temp|wfaz_bl|weight_change_p|race_c", unique(col_perf$variable), value=TRUE) 
  # 2: Ventilator parameters/measures
  group_2 <- grep("exhaled_vt|mean_airway_pressure|peep|pip|pressure_support|rr_set|total_rr|vent_aprv", unique(col_perf$variable), value=TRUE) 
  # 3: Blood gases/OI/OSI
  group_3 <- grep("base_excess|hgb|lactate|oi_osi|pco2|ph|sample_type", unique(col_perf$variable), value=TRUE) 
  # 4: Resp assessment
  group_4 <- grep("breath_sounds|respiratory_pattern|secretion_|cough", unique(col_perf$variable), value=TRUE)
  # 5: Other clinical parameters (Most frequently missing)
  group_5 <- grep("cvp|glasgow|hfaz_bl|motor|state_behavioral|wbc|wflz_bmiz_bl", unique(col_perf$variable), value=TRUE)

  # Define a function to find the group number for a given column
  find_group <- function(column_name) {
    if (column_name %in% group_2) return(2)
    if (column_name %in% group_3) return(3)
    if (column_name %in% group_4) return(4)
    if (column_name %in% group_5) return(5)
    return(NA)  # Return NA if the column doesn't fit any group
  }

  # Apply this function to create groups of variables
  col_perf[, group := factor(sapply(variable, find_group))]

  # Create a unique ID by concatenating imp, miss_type, p_miss, and index
  col_perf[, id := paste(imp, miss_type, p_miss, index, sep = "_")]

  # Limit to test
  col_perf_test <- col_perf[set_type=="test"]

  # Subset by type
  col_perf_test_mse <- col_perf_test[metric_type=="numeric"]
  col_perf_test_class <- col_perf_test[metric_type=="factor"]

  # Count how many variables in each model
  # MSE
  print(paste(out_suffix, ":number of vars in MSE model:", length(unique(col_perf_test_mse$variable))))
  table(col_perf_test_mse$group)

  # Classification error
  print(paste(out_suffix, ":number of vars in class error model:", length(unique(col_perf_test_class$variable))))
  table(col_perf_test_class$group)

  # Apply droplevels() to all factor variables in the data.table
  droplevels_fun <- function(data) {
    data[, names(data) := lapply(.SD, function(x) if(is.factor(x)) droplevels(x) else x)]
    return(data)
  }
  col_perf_test_mse <- droplevels_fun(col_perf_test_mse)
  col_perf_test_class <- droplevels_fun(col_perf_test_class)

  # Interactions
  col_three_way_interactions <- c("imp:miss_type:p_miss", "imp:miss_type:group", 
                              "imp:p_miss:group", "miss_type:p_miss:group")

  col_two_way_interactions <- c("imp:miss_type", "imp:p_miss", "imp:group", 
                            "miss_type:p_miss", "miss_type:group", "p_miss:group")

  col_all_interactions <- c(col_three_way_interactions, col_two_way_interactions)

  # Base part of the formula
  col_base_formula <- "error_value ~ imp + miss_type + p_miss + group" # + (1 | id)"

  # Add interaction terms to the base formula
  col_full_formula <- paste(col_base_formula, paste(col_all_interactions, collapse = " + "), sep = " + ")

  # List of datasets
  col_datasets <- list(
    col_perf_test_mse = col_perf_test_mse, 
    col_perf_test_class = col_perf_test_class
  )

  col_output <- run_stepwise_fun(datasets = col_datasets, 
                      full_formula = col_full_formula, 
                      all_interactions = col_all_interactions, 
                      model_type = "lm", #"lmer", 
                      outcome = out_suffix)

  return(col_output)
}

col_output_ext <- stepwise_col_fun("ext")
col_output_bp <- stepwise_col_fun("bp")

# Imputation performance models: Margins plots ----
# Define the function to plot EMMs
plot_emm_interaction <- function(final_model, model_name, outfolder, error_type) {
  # Compute the EMMs for the interaction of interest
  emm <- emmeans(final_model, ~ imp * p_miss * group, weights = "proportional")
is  
  # You can transform the EMMs into a data frame for plotting
  emm_df <- as.data.frame(emm)

  # Update the labels
  levels(emm_df$imp) <- c("LOCF", "Mean", "Random Forest AV", "Bayesian/PMM AV", "LASSO AV")

  # Color palette
  cbPalette <- c("#E69F00", "#6b62af", "#1c8042", "#0072B2", "#D55E00", "#c6659a", "#3a86b1", "#d17d3c", "#c38aa9", "#6e98b0", "#d49565", "#caaebe")

  # Define your custom colors
  custom_colors_imp <- c("Mean" = cbPalette[1], "LOCF" = cbPalette[3], 
                    "Random Forest AV" = cbPalette[7], "Bayesian/PMM AV" = cbPalette[8], "LASSO AV" = cbPalette[9])
  
  # Plot title
  if(grepl("ext", model_name)) {
    plot_title <- paste("Imputation performance (extubation):", error_type, "\nMarginal means for interaction between imputation method, proportion missing, and variable group")
  }
  if(grepl("bp", model_name)) {
    plot_title <- paste("Imputation performance (blood pressure):", error_type, "\nMarginal means for interaction between imputation method, proportion missing, and variable group")
  }

  # Plot the EMMs
  plot_emm <- ggplot(emm_df, aes(x = p_miss, y = emmean, color = imp, group = imp)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = custom_colors_imp) +
    facet_wrap(~ group) +
    labs(x = 'Proportion Missing', y = 'Estimated Marginal Mean', color = 'Imputation Method') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(plot_title)

  # Save plot to the specified output folder
  # Create the plot title based on the model_name
  error_title <- ifelse(grepl("mse", model_name), "Mean squared error", "Classification error")


  file_name <- paste0(outfolder, "/", model_name, "_", tolower(gsub(" ", "_", error_type)), ".pdf")
  ggsave(file_name, plot = plot_emm, width = 11, height = 8.5, units = "in", device = 'pdf')
}

# Assuming 'col_output$final_models' contains your models and 'outfolder' is defined
## Extubation
# Call the function for mean squared error
plot_emm_interaction(col_output_ext$final_models$col_perf_test_mse, "col_perf_test_mse_ext", outfolder, "Mean squared error")

# Call the function for classification error
plot_emm_interaction(col_output_ext$final_models$col_perf_test_class, "col_perf_test_class_ext", outfolder, "Classification error")

## BP
# Call the function for mean squared error
plot_emm_interaction(col_output_bp$final_models$col_perf_test_mse, "col_perf_test_mse_bp", outfolder, "Mean squared error")

# Call the function for classification error
plot_emm_interaction(col_output_bp$final_models$col_perf_test_class, "col_perf_test_class_bp", outfolder, "Classification error")

# Outcome model performance results ----
# NOT a multilevel model
mod_perf_ext <- as.data.table(readRDS(mod_perf_filepath_ext))
mod_perf_bp <- as.data.table(readRDS(mod_perf_filepath_bp))

# Extubation
# Exclude complete data
mod_perf_ext <- mod_perf_ext[imp != "comp"]
mod_perf_ext$imp <- relevel(mod_perf_ext$imp, ref = "locf") # Set 'locf' as reference level
mod_perf_ext <- mod_perf_ext[!grepl("noy", miss)]

# Interactions
mod_three_way_interactions <- c("imp:miss_type:p_miss")

mod_two_way_interactions <- c("imp:miss_type", "imp:p_miss", "miss_type:p_miss")

mod_all_interactions <- c(mod_three_way_interactions, mod_two_way_interactions)

# Subset data by outcome model: GBM and glmnet
mod_perf_ext_gbm <- mod_perf_ext[model == "gbm"]
mod_perf_ext_glmnet <- mod_perf_ext[model == "glmnet"]
# Apply droplevels() to all factor variables in the data.table
mod_perf_ext_gbm[, (names(mod_perf_ext_gbm)) := lapply(.SD, function(x) if(is.factor(x)) droplevels(x) else x)]
mod_perf_ext_glmnet[, (names(mod_perf_ext_glmnet)) := lapply(.SD, function(x) if(is.factor(x)) droplevels(x) else x)]

# Multiply balanced accuracy and AUC by 100
mod_perf_ext_gbm[, test_balanced_accuracy_100 := test_balanced_accuracy*100]
mod_perf_ext_glmnet[, test_balanced_accuracy_100 := test_balanced_accuracy*100]

mod_perf_ext_gbm[, test_auc_100 := test_auc*100]
mod_perf_ext_glmnet[, test_auc_100 := test_auc*100]

# List of datasets
mod_datasets <- list(
  mod_perf_ext_gbm = mod_perf_ext_gbm,
  mod_perf_ext_glmnet = mod_perf_ext_glmnet
)

# Balanced accuracy
# Base part of the formula
mod_base_formula <- "test_balanced_accuracy_100 ~ imp + miss_type + p_miss"

# Add interaction terms to the base formula
mod_full_formula <- paste(mod_base_formula, paste(mod_all_interactions, collapse = " + "), sep = " + ")

mod_output_balacc_100 <- run_stepwise_fun(datasets = mod_datasets, 
                     full_formula = mod_full_formula, 
                     all_interactions = mod_all_interactions, 
                     model_type = "lm",
                     outcome = "ba_100")

# AUC
# Base part of the formula
mod_base_formula <- "test_auc_100 ~ imp + miss_type + p_miss"

# Add interaction terms to the base formula
mod_full_formula <- paste(mod_base_formula, paste(mod_all_interactions, collapse = " + "), sep = " + ")

mod_output_auc_100 <- run_stepwise_fun(datasets = mod_datasets, 
                     full_formula = mod_full_formula, 
                     all_interactions = mod_all_interactions, 
                     model_type = "lm",
                     outcome = "auc_100")

# Blood Pressure
# Exclude complete data
mod_perf_bp <- mod_perf_bp[imp != "comp"]
mod_perf_bp$imp <- relevel(mod_perf_bp$imp, ref = "locf") # Set 'locf' as reference level
mod_perf_bp <- mod_perf_bp[!grepl("noy", miss)]

# Interactions
mod_three_way_interactions <- c("imp:miss_type:p_miss")

mod_two_way_interactions <- c("imp:miss_type", "imp:p_miss", "miss_type:p_miss")

mod_all_interactions <- c(mod_three_way_interactions, mod_two_way_interactions)

# Subset data by outcome model: GBM and glmnet
mod_perf_bp_gbm <- mod_perf_bp[model == "gbm"]
mod_perf_bp_glmnet <- mod_perf_bp[model == "glmnet"]
# Apply droplevels() to all factor variables in the data.table
mod_perf_bp_gbm[, (names(mod_perf_bp_gbm)) := lapply(.SD, function(x) if(is.factor(x)) droplevels(x) else x)]
mod_perf_bp_glmnet[, (names(mod_perf_bp_glmnet)) := lapply(.SD, function(x) if(is.factor(x)) droplevels(x) else x)]

# Multiply MSE by 100
mod_perf_bp_gbm[, test_mse_100 := test_mse*100]
mod_perf_bp_glmnet[, test_mse_100 := test_mse*100]

# List of datasets
mod_datasets <- list(
  mod_perf_bp_gbm = mod_perf_bp_gbm,
  mod_perf_bp_glmnet = mod_perf_bp_glmnet
)

# MSE
# Base part of the formula
mod_base_formula <- "test_mse_100 ~ imp + miss_type + p_miss"

# Add interaction terms to the base formula
mod_full_formula <- paste(mod_base_formula, paste(mod_all_interactions, collapse = " + "), sep = " + ")

mod_output_mse_100 <- run_stepwise_fun(datasets = mod_datasets, 
                     full_formula = mod_full_formula, 
                     all_interactions = mod_all_interactions, 
                     model_type = "lm",
                     outcome = "mse_100")          

# Save the excel output ----
saveWorkbook(wb, paste0(outfolder, "/Col_Mod_Results.xlsx"), overwrite = TRUE)
