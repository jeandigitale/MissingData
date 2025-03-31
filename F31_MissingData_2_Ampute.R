# F31_MissingData_2_Ampute
## Make variables numeric and remove ID variables to prepare for ampute
## Group variables for ampute
## Create pattern and weight matrices for ampute
## Ampute
## Convert variables back to factors
## Characterize missingness in amputed data compared with original and export results to excel
## Create lagged variables (not including baseline variables due to computational burden for imputation)
## Limit dataset to time windows ≥12 hours after intubation through time window ending 12 hours prior to intubation (because outcome is lagged by 12 hours)
## Test/train split
## Preprocess for rtemis
## Save data in Amp_Data

# Libraries ----
library(tidyverse)
library(openxlsx)
library(data.table)
library(knitr)
library(mice)
library(zoo)
library(rtemis)

# Set-up ----
rm(list=ls(all=TRUE))

today <- format(Sys.Date(), "%Y-%m-%d")

starttime <- proc.time()

# Set number of hours per time window
timewindow_hrs <- 4

# Set characteristics of lag
lag_hrs <- 12

# Folder
folder <- "filepath"

# File date
file_date <- "date"

# Filename in
filename_in <- paste0(folder, "MissingData_Result/1_CompleteData/dt_tw_1_", file_date, ".RDS")
filename_in_comp <- paste0(folder, "MissingData_Result/1_CompleteData/dt_tw_comp_1_", file_date, ".RDS")

# Filename out
subfolder <- file.path(folder, "MissingData_Result/2_Ampute", today)
dir.create(subfolder, showWarnings = FALSE)
subfolder_ext <- file.path(subfolder, "Ext")
subfolder_bp <- file.path(subfolder, "BP")
dir.create(subfolder_ext, showWarnings = FALSE)
dir.create(subfolder_bp, showWarnings = FALSE)

# Load data ----
dt_tw <- readRDS(filename_in)
dt_tw_comp <- readRDS(filename_in_comp)

# Sort rows to ensure they are in correct order by pat_mrn_id
setorder(dt_tw, pat_mrn_id, placement_instant, date_time_start)
setorder(dt_tw_comp, pat_mrn_id, placement_instant, date_time_start)

# Create systolic BP outcome as mean value one time window ahead
dt_tw[, y_sys_bpcentile_mean_lag_4 := shift(sys_bpcentile_mean, type = "lead"), 
      by = .(pat_mrn_id, placement_instant)]
dt_tw_comp[, y_sys_bpcentile_mean_lag_4 := shift(sys_bpcentile_mean, type = "lead"), 
      by = .(pat_mrn_id, placement_instant)]
#View(dt_tw_comp[, .(pat_mrn_id, date_time_start, sys_bp_mean, sys_bpcentile_mean, y_sys_bpcentile_mean_lag_4, drop_lt12h_toext, sex, hfaz_bl)])

# Filter out incomplete cases ----
# NOTE: ALL CASES COMPLETE AS THEY SHOULD BE 12-22-23
# Count number of missing per row
comp_ck <- copy(dt_tw_comp[drop_lt12h_toext==0])
comp_ck[drop_lt12h_toext==0, num_miss := rowSums(is.na(comp_ck))]
table(comp_ck$num_miss)
rm(comp_ck)

# Create vector to filter out incomplete cases
complete_cases <- complete.cases(dt_tw_comp)
dt_tw_comp <- dt_tw_comp[complete_cases]

# Original dataset: Keep only cases and variables in completed dataset ----
# Keep only complete cases
dt_tw <- dt_tw[complete_cases]

# Remove non-numeric/ID variables so dataset will match when characterizing missingness ----
dt_tw_id <- dt_tw[, c("year", "pat_mrn_id", "placement_instant", "ucsf_bcho", "date_time_start", "date_time_end", "tw_id", "drop_lt12h_toext")]
dt_tw[, c("year", "pat_mrn_id", "placement_instant", "ucsf_bcho", "date_time_start", "date_time_end", "tw_id", "drop_lt12h_toext") := NULL]

# Remove non-numeric/ID variables for ampute
dt_tw_comp_id <- dt_tw_comp[, c("year", "pat_mrn_id", "placement_instant", "ucsf_bcho", "date_time_start", "date_time_end", "tw_id", "drop_lt12h_toext")]
dt_tw_comp[, c("year", "pat_mrn_id", "placement_instant", "ucsf_bcho", "date_time_start", "date_time_end", "tw_id", "drop_lt12h_toext") := NULL]

# Check dt_tw_id == dt_tw_comp_id: TRUE 12/22/23
for (col in names(dt_tw_id)) {
   if (!all.equal(dt_tw_id[[col]], dt_tw_comp_id[[col]])) {
      cat("Differences found in column:", col, "\n")
   }
}

# Original dataset: Make variables numeric to match complete dataset ----
dt_tw <- preprocess(dt_tw, character2factor = TRUE)
dt_tw <- preprocess(dt_tw, oneHot = TRUE)

check_data(dt_tw)

# Completed dataset: Make variables numeric for ampute ----
dt_tw_comp <- preprocess(dt_tw_comp, character2factor = TRUE)
dt_tw_comp <- preprocess(dt_tw_comp, oneHot = TRUE)

check_data(dt_tw_comp)

# Drop constant variables in completed dataset
# NOTE: NO CONSTANT VARIABLES 12-21-23
constant <- which(sapply(dt_tw_comp, is_constant))
if (length(constant) > 0) {
  print(paste0("The following variables are constant:", constant))
  dt_tw_comp[, (constant) := NULL]
}

# In original data, keep only variables in completed dataset
cols_drop <- names(dt_tw)[!(names(dt_tw) %in% names(dt_tw_comp))]
if (length(cols_drop) > 0) {
  dt_tw[, (cols_drop) := NULL]
}

# Completed dataset: Group variables for ampute ----
# Groups defined based on theoretical groupings and empirical cluster analysis: /Users/jeandigitale/Dropbox/Dropbox_JD/UCSF/Extubation/Analysis/MissingData_Result/1_CompleteData/Ext_MissingData_2023-12-21_final.xlsx

# 1: Always observed
group_1 <- grep("ext_success|age_int_days|alternative_gases|bp|cuffed|dx|etco2|ett_cuff_leak|ett_placed|ett_size|fio2|hrs_int|intake_output|map|med_|present_on_hosp|pulse|resp_|sbt|sex|spo2|temp|wfaz_bl|weight_change_p", names(dt_tw_comp), value=TRUE) 
# 2: Ventilator parameters/measures
group_2 <- grep("exhaled_vt|mean_airway_pressure|peep|pip|pressure_support|rr_set|total_rr|vent_aprv", names(dt_tw_comp), value=TRUE) 
# 3: Blood gases/OI/OSI
group_3 <- grep("base_excess|hgb|lactate|oi_osi|pco2|ph|sample_type", names(dt_tw_comp), value=TRUE) 
# 4: Resp assessment
group_4 <- grep("breath_sounds|respiratory_pattern|secretion_|cough", names(dt_tw_comp), value=TRUE)
# 5: Other clinical parameters (Most frequently missing)
group_5 <- grep("cvp|glasgow|hfaz_bl|motor|state_behavioral|wbc|wflz_bmiz_bl", names(dt_tw_comp), value=TRUE)

# Check exclusions
length(names(dt_tw)[names(dt_tw) %in% c(group_1, group_2, group_3, group_4, group_5)])
names(dt_tw)[!(names(dt_tw) %in% c(group_1, group_2, group_3, group_4, group_5))]
# Number variables in each group
n1 <- length(group_1)
n2 <- length(group_2)
n3 <- length(group_3)
n4 <- length(group_4)
n5 <- length(group_5)

# Ampute pattern matrix ----
# Pattern matrix: 0 = missing, 1 = observed (15 combinations)
miss_pattern <- rbind(as.numeric(!(names(dt_tw_comp) %in% group_2)), # missing group 2
                      as.numeric(!(names(dt_tw_comp) %in% group_3)), # missing group 3
                      as.numeric(!(names(dt_tw_comp) %in% group_4)), # missing group 4
                      as.numeric(!(names(dt_tw_comp) %in% group_5)), # missing group 5
                      as.numeric(!(names(dt_tw_comp) %in% c(group_2, group_3))), # missing group 2 and 3
                      as.numeric(!(names(dt_tw_comp) %in% c(group_2, group_4))), # missing group 2 and 4
                      as.numeric(!(names(dt_tw_comp) %in% c(group_2, group_5))), # missing group 2 and 5
                      as.numeric(!(names(dt_tw_comp) %in% c(group_3, group_4))), # missing group 3 and 4
                      as.numeric(!(names(dt_tw_comp) %in% c(group_3, group_5))), # missing group 3 and 5
                      as.numeric(!(names(dt_tw_comp) %in% c(group_4, group_4))), # missing group 4 and 5
                      as.numeric(!(names(dt_tw_comp) %in% c(group_2, group_3, group_4))), # missing group 2, 3, 4
                      as.numeric(!(names(dt_tw_comp) %in% c(group_2, group_3, group_5))), # missing group 2, 3, 5
                      as.numeric(!(names(dt_tw_comp) %in% c(group_2, group_4, group_5))), # missing group 2, 4, 5
                      as.numeric(!(names(dt_tw_comp) %in% c(group_3, group_4, group_5))), # missing group 3, 4, 5
                      as.numeric(!(names(dt_tw_comp) %in% c(group_2, group_3, group_4, group_5)))) # missing group 2, 3, 4, 5
colnames(miss_pattern) <- names(dt_tw_comp)

# Ampute weights matrices ----
# Capture prefixes of unique concepts for weights
unique_pre <- unique(sub("_mean.*|_median.*|_max.*|_min.*|_first.*|_last.*|_mode.*|\\.x.*", "", names(dt_tw_comp)))
unique_pre[unique_pre=="resp"] <- "resp_" # make resp unique from "respiratory"
unique_pre[unique_pre=="dia_bp"] <- "dia_bp_" # make bp unique from "bpcentile"
unique_pre[unique_pre=="sys_bp"] <- "sys_bp_" # make bp unique from "bpcentile"

# Make multiple select prefixes unique
unique_pre <- unique(sub("(cough).*", "\\1", unique_pre))
unique_pre <- unique(sub("(_sounds).*", "\\1", unique_pre))
unique_pre <- unique(sub("(_color).*", "\\1", unique_pre))
unique_pre <- unique(sub("(_consistency).*", "\\1", unique_pre))
unique_pre <- unique(sub("(_pattern).*", "\\1", unique_pre))

# Weights matrices

# MAR: observed variables = 1/n_obs/y (including outcome); missing variables = 0
# weak MNAR (y-related): observed variables = 1/n_obs/y; missing variables = 0.5/n_miss/y (including outcome)
# moderate MNAR (y-related): observed variables = 1/n_obs/y; missing variables = 1/n_miss/y (including outcome)
# strong MNAR (y-related): observed variables = 0; missing variables = 1/n_miss/y (including outcome)
# weak MNAR (NOT y-related): observed variables = 1/n_obs/y; missing variables = 0.5/n_miss/y
# moderate MNAR (NOT y-related): observed variables = 1/n_obs/y; missing variables = 1/n_miss/y
# strong MNAR (NOT y-related): observed variables = 0; missing variables = 1/n_miss/y

# Observed weights:
  # NOTE: y_sys_bpcentile_mean_lag_4 and ext_success_lag_12 are automatically included in observed weights because they are in group 1 (always observed variables)
  # For MAR extubation ampute - need to exclude y_sys_bpcentile_mean_lag_4
  # For MAR bp ampute - need to exclude ext_success_lag_12
  # For all MNAR - need exclude both outcomes because the outcome will be included in the missing weights
# Missing weights:
  # For MNAR with outcome - need to add in relevant outcome
  # For MNAR without outcome - do not add in anything

weight_mat_fun_ex <- function(obs, miss, outcome, include_y_miss = 0) {
  
  # Create empty matrix
  weight_mat <- array(dim = dim(miss_pattern))
  colnames(weight_mat) <- names(dt_tw_comp)
  
  # For MAR
  if(outcome == "ext_success_lag_12" & miss == 0 ){ # miss = 0 is MAR
    exclude_var_obs <- "y_sys_bpcentile_mean_lag_4"
  }
  if(outcome == "y_sys_bpcentile_mean_lag_4" & miss == 0) { # miss = 0 is MAR
    exclude_var_obs <- "ext_success_lag_12"
  }

  # For MNAR
  if(miss > 0) { # miss > 0 is MNAR
    exclude_var_obs <- c("ext_success_lag_12", "y_sys_bpcentile_mean_lag_4")
  }

  # For MNAR with outcome
  include_var_miss <- ""

  if(include_y_miss==1) {
    include_var_miss <- outcome
  }

  for (i in 1:nrow(weight_mat)) {
    # Capture names of non-missing variables, excluding those in 'exclude_var_obs'
    names_obs <- names(dt_tw_comp)[miss_pattern[i, ] != 0 & !names(dt_tw_comp) %in% exclude_var_obs] 
    # Check how many non-missing variables there are per unique concept (e.g. pulse, resp, etc are each 1 unique variable)
    unique_obs_n <- sapply(unique_pre, function(x) sum(grepl(paste0("^", x), names_obs)))
    # Count how many non-missing unique concepts there are
    unique_obs_n_tot <- sum(unique_obs_n > 0)
    # Create weights for non-missing variables
    unique_obs_weights <- obs / unique_obs_n_tot / unique_obs_n[unique_obs_n > 0]
    
    # Capture names of missing variables, including those in 'include_var_miss'
    names_miss <- names(dt_tw_comp)[miss_pattern[i, ] == 0] 
    if(include_var_miss!=""){
      names_miss <- c(names_miss, include_var_miss)
    }   
    # Check how many missing variables there are per unique concept (e.g. pulse, resp, etc are each 1 unique variable)
    unique_miss_n <- sapply(unique_pre, function(x) sum(grepl(paste0("^", x),  names_miss)))
    # Count how many missing unique concepts there are
    unique_miss_n_tot <- sum(unique_miss_n > 0)
    # Create weights for missing variables
    unique_miss_weights <- miss / unique_miss_n_tot / unique_miss_n[unique_miss_n > 0]
    
    # Create weight vector combining observed and missing weights
    unique_weights <- c(unique_obs_weights, unique_miss_weights)
  
    for (j in 1:length(unique_weights)) { # loop through each unique weight
      # Capture which variables match the unique prefix and should get the corresponding weight
      matching_names <- colnames(weight_mat)[grep(paste0("^", names(unique_weights)[j]), colnames(weight_mat))]
      # Set the corresponding cells in the result matrix to the matching weights
      weight_mat[i, which(colnames(weight_mat) %in% matching_names)] <- unique_weights[j]
    }
    
    # Set all NA values in weight_mat to 0
    weight_mat[is.na(weight_mat)] <- 0

  }
  return(weight_mat)
}

# Weight matrices - Extubation ----
# Missing at random weight matrix - extubation (weight assigned to ext_success_lag_12)
  ## Observed = 1 / unique OBSERVED concepts / number variables per concept (e.g. pulse is a unique concept, 6 pulse variables, weight is 1/n concepts/6)
  ## Missing = 0
weight_mar_ext <- weight_mat_fun_ex(obs=1, miss=0, outcome="ext_success_lag_12")

# Missing not at random weight matrix - including Y in weights
# Weak MNAR 
  ## Observed = 1 / unique OBSERVED concepts / number variables per concept
  ## Missing = 0.5 / unique MISSING concepts / number variables per concept
weight_mnar_w_ext <- weight_mat_fun_ex(obs=1, miss=0.5, outcome="ext_success_lag_12", include_y_miss = 1)
# Moderate MNAR
  ## Observed = 1 / unique OBSERVED concepts / number variables per concept
  ## Missing = 1 / unique MISSING concepts / number variables per concept
weight_mnar_m_ext <- weight_mat_fun_ex(obs=1, miss=1, outcome="ext_success_lag_12", include_y_miss = 1)
# Strong MNAR
  ## Observed = 0
  ## Missing = 1 / unique MISSING concepts / number variables per concept
weight_mnar_s_ext <- weight_mat_fun_ex(obs=0, miss=1, outcome="ext_success_lag_12", include_y_miss = 1)

# Missing not at random weight matrix - excluding Y from weights
# Weak MNAR 
weight_mnar_w_ext_noy <- weight_mat_fun_ex(obs=1, miss=0.5, outcome="ext_success_lag_12")
# Moderate MNAR
weight_mnar_m_ext_noy <- weight_mat_fun_ex(obs=1, miss=1, outcome="ext_success_lag_12")
# Strong MNAR
weight_mnar_s_ext_noy <- weight_mat_fun_ex(obs=0, miss=1, outcome="ext_success_lag_12")

# Weight matrices - BP ----
# Missing at random weight matrix
weight_mar_bp <- weight_mat_fun_ex(obs=1, miss=0,  outcome="y_sys_bpcentile_mean_lag_4")
# Missing not at random weight matrix - including Y in weights
# Weak MNAR 
weight_mnar_w_bp <- weight_mat_fun_ex(obs=1, miss=0.5,  outcome="y_sys_bpcentile_mean_lag_4", include_y_miss = 1)
# Moderate MNAR
weight_mnar_m_bp <- weight_mat_fun_ex(obs=1, miss=1, outcome="y_sys_bpcentile_mean_lag_4", include_y_miss = 1)
# Strong MNAR
weight_mnar_s_bp <- weight_mat_fun_ex(obs=0, miss=1, outcome="y_sys_bpcentile_mean_lag_4", include_y_miss = 1)

# Missing not at random weight matrix - excluding Y from weights
# Weak MNAR 
weight_mnar_w_bp_noy <- weight_mat_fun_ex(obs=1, miss=0.5,  outcome="y_sys_bpcentile_mean_lag_4")
# Moderate MNAR
weight_mnar_m_bp_noy <- weight_mat_fun_ex(obs=1, miss=1, outcome="y_sys_bpcentile_mean_lag_4")
# Strong MNAR
weight_mnar_s_bp_noy <- weight_mat_fun_ex(obs=0, miss=1, outcome="y_sys_bpcentile_mean_lag_4")

# Check weights
# weight_ck <- weight_mnar_w_ext
# exclude_miss_var <- c("y_sys_bpcentile_mean_lag_4")
# exclude_obs_var <- c("ext_success_lag_12", "y_sys_bpcentile_mean_lag_4")
# outcome_miss <- "ext_success_lag_12"
# weight_ck[, "ext_success_lag_12"]
# weight_ck[, "y_sys_bpcentile_mean_lag_4"]

# # Check: Sum all the cells that are 0 per row (missing), excluding the columns in exclude_var
# rowSums(weight_ck * ((miss_pattern == 0 & !colnames(miss_pattern)[col(miss_pattern)] %in% exclude_miss_var)| 
#                     colnames(miss_pattern)[col(miss_pattern)] == outcome_miss))

# # Check: Sum all the cells that are 1 per row (observed)
# rowSums(weight_ck * (miss_pattern == 1 & !colnames(miss_pattern)[col(miss_pattern)] %in% exclude_obs_var))

# Define frequencies for missingness patterns for three levels of missingness ----
# Frequency patterns for 1/2x percent missing in original data
pat_freq_0 <- c(0.05,  # missing group 2
              0.18,   # missing group 3
              0.09,   # missing group 4
              0.31,  # missing group 5
              0.08,  # missing group 2 and 3
              0.02,  # missing group 2 and 4
              0.01,  # missing group 2 and 5
              0.1,   # missing group 3 and 4
              0.08,   # missing group 3 and 5
              0.03,   # missing group 4 and 5
              0.01,   # missing group 2, 3, 4
              0.01,  # missing group 2, 3, 5
              0.01,  # missing group 2, 4, 5
              0.01,  # missing group 3, 4, 5
              0.01)  # missing group 2, 3, 4, 5

# Frequency patterns for 1x percent missing in original data
pat_freq_1 <- c(0.01,  # missing group 2
              0.2,   # missing group 3
              0.02,   # missing group 4
              0.26,  # missing group 5
              0.06,  # missing group 2 and 3
              0.01,  # missing group 2 and 4
              0.03,  # missing group 2 and 5
              0.12,   # missing group 3 and 4
              0.08,   # missing group 3 and 5
              0.05,   # missing group 4 and 5
              0.03,   # missing group 2, 3, 4
              0.05,  # missing group 2, 3, 5
              0.02,  # missing group 2, 4, 5
              0.05,  # missing group 3, 4, 5
              0.01)  # missing group 2, 3, 4, 5

# Frequency patterns for 2x percent missing in original data
pat_freq_2 <- c(0.01,  # missing group 2
              0.01,   # missing group 3
              0.01,   # missing group 4
              0.01,  # missing group 5
              0.01,  # missing group 2 and 3
              0.01,  # missing group 2 and 4
              0.01,  # missing group 2 and 5
              0.07,   # missing group 3 and 4
              0.26,   # missing group 3 and 5
              0.03,   # missing group 4 and 5
              0.01,   # missing group 2, 3, 4
              0.14,  # missing group 2, 3, 5
              0.08,  # missing group 2, 4, 5
              0.2,  # missing group 3, 4, 5
              0.14)  # missing group 2, 3, 4, 5

# Set proportion for levels of missingness
prop_0 <- 0.5
prop_1 <- 0.83
prop_2 <- 0.97

# Ampute function ----
# Function to generate ampute outputs with random types
ampute_reps_fun <- function(data=dt_tw_comp, mech, prop, freq, patterns=miss_pattern, weights=NULL, num_replications=20) {
  # Generate random types for each replication
  random_types <- replicate(num_replications, {
    random_num <- runif(1)
    if (random_num < 0.25) {
      type <- "LEFT"
    } else if (random_num < 0.5) {
      type <- "MID"
    } else if (random_num < 0.75) {
      type <- "RIGHT"
    } else {
      type <- "TAIL"
    }
    type
  })
  
  # Call ampute() for each random type using lapply
  output_list <- lapply(random_types, function(type) {
    ampute(
      data = data,
      mech = mech,
      prop = prop,
      type = type,
      freq = freq,
      patterns = patterns, 
      weights = weights
    )
  })
  
  # Return the output list
  return(output_list)
}

# Create amputed datasets ----
# Set the seed for reproducibility
set.seed(406)
# Extubation
# MCAR: 1/2x missingness
mcar_0_ext <- ampute_reps_fun(mech="MCAR", 
                      prop=prop_0, 
                      freq=pat_freq_0)

# MCAR: 1x missingness
mcar_1_ext <- ampute_reps_fun(mech="MCAR", 
                      prop=prop_1, 
                      freq=pat_freq_1)

# MCAR: 2x missingness
mcar_2_ext <- ampute_reps_fun(mech="MCAR", 
                      prop=prop_2, 
                      freq=pat_freq_2)

# MAR: 1/2x missingness
mar_0_ext <- ampute_reps_fun(mech="MAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mar_ext)

# MAR: 1x missingness
mar_1_ext <- ampute_reps_fun(mech="MAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mar_ext)

# MAR: 2x missingness
mar_2_ext <- ampute_reps_fun(mech="MAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mar_ext)

# Weak MNAR: 1/2x missingness
mnar_w_0_ext <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_w_ext)

# Weak MNAR: 1x missingness
mnar_w_1_ext <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_w_ext)

# Weak MNAR: 2x missingness
mnar_w_2_ext <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_w_ext)

# Moderate MNAR: 1/2x missingness
mnar_m_0_ext <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_m_ext)

# Moderate MNAR: 1x missingness
mnar_m_1_ext <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_m_ext)

# Moderate MNAR: 2x missingness
mnar_m_2_ext <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_m_ext)

# Strong MNAR: 1/2x missingness
mnar_s_0_ext <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_s_ext)

# Strong MNAR: 1x missingness
mnar_s_1_ext <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_s_ext)

# Strong MNAR: 2x missingness
mnar_s_2_ext <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_s_ext)

# MNAR: No Y ----
# MNAR with no ext_success_lag_12
# Weak MNAR: 1/2x missingness
mnar_w_0_ext_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_w_ext_noy)

# Weak MNAR: 1x missingness
mnar_w_1_ext_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_w_ext_noy)

# Weak MNAR: 2x missingness
mnar_w_2_ext_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_w_ext_noy)

# Moderate MNAR: 1/2x missingness
mnar_m_0_ext_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_m_ext_noy)

# Moderate MNAR: 1x missingness
mnar_m_1_ext_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_m_ext_noy)

# Moderate MNAR: 2x missingness
mnar_m_2_ext_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_m_ext_noy)

# Strong MNAR: 1/2x missingness
mnar_s_0_ext_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_s_ext_noy)

# Strong MNAR: 1x missingness
mnar_s_1_ext_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_s_ext_noy)

# Strong MNAR: 2x missingness
mnar_s_2_ext_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_s_ext_noy)

# Set the seed for reproducibility
set.seed(394)
# BP
# MCAR: 1/2x missingness
mcar_0_bp <- ampute_reps_fun(mech="MCAR", 
                      prop=prop_0, 
                      freq=pat_freq_0)

# MCAR: 1x missingness
mcar_1_bp <- ampute_reps_fun(mech="MCAR", 
                      prop=prop_1, 
                      freq=pat_freq_1)

# MCAR: 2x missingness
mcar_2_bp <- ampute_reps_fun(mech="MCAR", 
                      prop=prop_2, 
                      freq=pat_freq_2)

# MAR: 1/2x missingness
mar_0_bp <- ampute_reps_fun(mech="MAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mar_bp)

# MAR: 1x missingness
mar_1_bp <- ampute_reps_fun(mech="MAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mar_bp)

# MAR: 2x missingness
mar_2_bp <- ampute_reps_fun(mech="MAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mar_bp)

# Weak MNAR: 1/2x missingness
mnar_w_0_bp <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_w_bp)

# Weak MNAR: 1x missingness
mnar_w_1_bp <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_w_bp)

# Weak MNAR: 2x missingness
mnar_w_2_bp <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_w_bp)

# Moderate MNAR: 1/2x missingness
mnar_m_0_bp <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_m_bp)

# Moderate MNAR: 1x missingness
mnar_m_1_bp <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_m_bp)

# Moderate MNAR: 2x missingness
mnar_m_2_bp <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_m_bp)

# Strong MNAR: 1/2x missingness
mnar_s_0_bp <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_s_bp)

# Strong MNAR: 1x missingness
mnar_s_1_bp <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_s_bp)

# Strong MNAR: 2x missingness
mnar_s_2_bp <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_s_bp)

# MNAR: No Y ----
# MNAR with no bp_success_lag_12
# Weak MNAR: 1/2x missingness
mnar_w_0_bp_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_w_bp_noy)

# Weak MNAR: 1x missingness
mnar_w_1_bp_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_w_bp_noy)

# Weak MNAR: 2x missingness
mnar_w_2_bp_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_w_bp_noy)

# Moderate MNAR: 1/2x missingness
mnar_m_0_bp_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_m_bp_noy)

# Moderate MNAR: 1x missingness
mnar_m_1_bp_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_m_bp_noy)

# Moderate MNAR: 2x missingness
mnar_m_2_bp_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_m_bp_noy)

# Strong MNAR: 1/2x missingness
mnar_s_0_bp_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_0, 
                      freq=pat_freq_0,
                      weights=weight_mnar_s_bp_noy)

# Strong MNAR: 1x missingness
mnar_s_1_bp_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_1, 
                      freq=pat_freq_1,
                      weights=weight_mnar_s_bp_noy)

# Strong MNAR: 2x missingness
mnar_s_2_bp_noy <- ampute_reps_fun(mech="MNAR", 
                      prop=prop_2, 
                      freq=pat_freq_2,
                      weights=weight_mnar_s_bp_noy)

# Extract and save just the amputed datasets ----
extract_amp_dataframes <- function(amp_name) {
  list_object <- get(amp_name)
  
  amp_dataframes <- vector("list", length = length(list_object))
  
  for (i in seq_along(list_object)) {
    amp_dataframes[[i]] <- as.data.table(list_object[[i]]$amp)
  }
  
  result <- setNames(amp_dataframes, paste0("dt_", amp_name))
  
  return(result)
}

amp_names <- c("mcar_0", "mcar_1", "mcar_2", "mar_0", "mar_1", "mar_2", "mnar_w_0", "mnar_w_1", "mnar_w_2", "mnar_m_0", "mnar_m_1", "mnar_m_2", "mnar_s_0", "mnar_s_1", "mnar_s_2")

amp_names_ext <- paste0(amp_names, "_ext")
amp_names_mnar_ext_noy <- paste0(grep("mnar", amp_names_ext, value = TRUE), "_noy")
amp_names_ext <- c(amp_names_ext, amp_names_mnar_ext_noy)

amp_names_bp <- paste0(amp_names, "_bp")
amp_names_mnar_bp_noy <- paste0(grep("mnar", amp_names_bp, value = TRUE), "_noy")
amp_names_bp <- c(amp_names_bp, amp_names_mnar_bp_noy)

amp_names_all <- c(amp_names_ext, amp_names_bp)

for (amp_name in amp_names_all) {
  result <- extract_amp_dataframes(amp_name)
  assign(paste0("dt_", amp_name), result)
}

dt_amp_names_ext <- paste0("dt_", amp_names_ext)
dt_amp_names_bp <- paste0("dt_", amp_names_bp)
dt_amp_names_all <- c(dt_amp_names_ext, dt_amp_names_bp)

# Convert variables back to factors ----
# Variables that need to be converted back to factors:
  # secretion_amount_mode, secretion_amount_first, secretion_amount_last
  # ett_cuff_leak_pressure_cat_last
  # sex
  # oi_osi                     
  # sbt_cat    
factor_cols <- c("secretion_amount_mode", "ett_cuff_leak_pressure_cat_last", "sex", "oi_osi", "sbt_cat")

create_factor_fun <- function(dt) {
  for (factor_col in factor_cols) {
    # Convert indicator_cols to character vector
    indicator_cols <- grep(factor_col, colnames(dt), value = TRUE)
  
    # Modify the indicator column names to create factor levels
    factor_levels <- gsub(paste0(factor_col, ".x_"), "", indicator_cols)
  
    # Create the factor variable using interaction of indicator variables
    dt[, (factor_col) := interaction(.SD, drop = TRUE), .SDcols = indicator_cols]
  
    # Check if all indicator variables are 0 and set the factor variable to NA
    dt[rowSums(dt[, .SD, .SDcols = indicator_cols]) == 0, (factor_col) := NA]
  
    # Drop unused levels using set()
    set(dt, j = factor_col, value = droplevels(dt[[factor_col]]))
    # Set the factor levels using modified indicator column names
    setattr(dt[[factor_col]], "levels", factor_levels)
    
    # Delete the indicator columns
    dt[, (indicator_cols) := NULL]
  }
  
  # Define the desired order of levels
  secretion_amount_levels <- c("None", "Scant", "Small", "Moderate", "Large", "Copious")

  # Make secretion amount an ordered factor
  dt[, c("secretion_amount_mode") := lapply(.SD, function(x) factor(x, levels = secretion_amount_levels, ordered = TRUE)), .SDcols = c("secretion_amount_mode")]

  # Return the modified data.table
  return(dt)
}

# Create the factor variables in original data
dt_tw <- create_factor_fun(dt_tw)

# Create the factor variables in completed data
dt_tw_comp <- create_factor_fun(dt_tw_comp)

# Apply the create_factor_fun to each amputed dataset
for (i in 1:length(dt_amp_names_all)) {
  amp_name <- dt_amp_names_all[i]  # Get the list name
  current_list <- get(amp_name)  # Get the current list
  
  # Apply the create_factor_fun to each dataset within the current list
  for (j in 1:length(current_list)) {
    current_list[[j]] <- create_factor_fun(current_list[[j]])
  }
  
  # Store the updated current list back in the global environment
  assign(amp_name, current_list)
}

# Re-define groups of variables with names of factors ----
# 1: Always observed
group_1 <- grep("ext_success|age_int_days|alternative_gases|bp|cuffed|dx|etco2|ett_cuff_leak|ett_placed|ett_size|fio2|hrs_int|intake_output|map|med_|present_on_hosp|pulse|resp_|sbt|sex|spo2|temp|wfaz_bl|weight_change_p", names(dt_tw), value=TRUE) 
# 2: Ventilator parameters/measures
group_2 <- grep("exhaled_vt|mean_airway_pressure|peep|pip|pressure_support|rr_set|total_rr|vent_aprv", names(dt_tw), value=TRUE) 
# 3: Blood gases/OI/OSI
group_3 <- grep("base_excess|hgb|lactate|oi_osi|pco2|ph|sample_type", names(dt_tw), value=TRUE) 
# 4: Resp assessment
group_4 <- grep("breath_sounds|respiratory_pattern|secretion_|cough", names(dt_tw), value=TRUE)
# 5: Other clinical parameters (Most frequently missing)
group_5 <- grep("cvp|glasgow|hfaz_bl|motor|state_behavioral|wbc|wflz_bmiz_bl", names(dt_tw), value=TRUE)

# Check exclusions
length(names(dt_tw)[names(dt_tw) %in% c(group_1, group_2, group_3, group_4, group_5)])
names(dt_tw)[!(names(dt_tw) %in% c(group_1, group_2, group_3, group_4, group_5))]

# Function to characterize missingness ----
missing_ck_fun <- function(data) {
  # Number of missing cells
  n_cell <- nrow(data)* ncol(data)
  n_cell_miss <- sum(is.na(data))
  p_miss_o <- n_cell_miss/n_cell

  # Count number of missing per row
  num_miss_o <- rowSums(is.na(data))
  cumsum(prop.table(table(num_miss_o)))
  sum_num_miss <- summary(num_miss_o)
  quant_num_miss <- quantile(num_miss_o, probs = seq(0,1,0.1))

  # Number of missing cells by group in original data
  # Create vector tagging which group each variable is in
  group_o <- rep(NA, length(names(data)))
  group_o[(names(data) %in% group_1)] <- 1 
  group_o[(names(data) %in% group_2)] <- 2 
  group_o[(names(data) %in% group_3)] <- 3 
  group_o[(names(data) %in% group_4)] <- 4 
  group_o[(names(data) %in% group_5)] <- 5 
  #group_o[is.na(group_o)] <- 0

  # Sum missing per column
  na_var_n_o <- colSums(is.na(data))

  # Sum missing per group
  na_group_n_o <- tapply(na_var_n_o, group_o, sum)

  # Total cells per group = cols per group * total rows
  tot_group_n_o <- tapply(group_o, group_o, length)*nrow(data)

  # Percent missing per group
  na_group_p_o <- na_group_n_o/tot_group_n_o

  # Create a named vector of percent missing
  p_miss <- c(p_miss_o, na_group_p_o)
  names(p_miss) <- c("total_p_missing", "group_1_p_missing", "group_2_p_missing", "group_3_p_missing", "group_4_p_missing", "group_5_p_missing")

  return(list(
    p_miss = p_miss,
    sum_num_miss = sum_num_miss,
    quant_num_miss = quant_num_miss
  ))
}

# Apply function to check missingness to amputed datasets ----
missing_ck_list <- lapply(dt_amp_names_all, function(amp_name) {
  list_object <- get(amp_name)
  
  result <- lapply(list_object, function(df) {
    result <- missing_ck_fun(df)
    return(result)
  })
  
  return(result)
})
names(missing_ck_list) <- amp_names_all

# Apply function to check missingness to original dataset for comparison
missing_ck_dt_tw <- missing_ck_fun(dt_tw)

# Function to create matrix of missingness outcomes to compare ----
missing_ck_matrix_fun <- function(missing_ck_list, element_index) {
  # Get the number of lists within missing_ck_list
  num_lists <- length(missing_ck_list)

  # Get the number of datasets within each list
  num_datasets <- length(missing_ck_list[[1]])

  # Calculate the number of columns in the result matrix
  num_cols <- length(missing_ck_list[[1]][[1]][[element_index]]) + 3

  # Create an empty matrix to store the results
  result_matrix <- matrix(nrow = num_lists * num_datasets, ncol = num_cols)

  # Loop through each list
  for (i in 1:num_lists) {
    amp_name <- dt_amp_names_all[i]  # Get the list name
    current_list <- missing_ck_list[[i]]  # Get the current list

    # Extract the suffix from the dataset name
    suffix <- gsub("[^0-9]", "", amp_name)

    # Loop through each dataset within the current list
    for (j in 1:num_datasets) {
      element <- current_list[[j]][[element_index]]  # Get the specified element of the current list
      row_index <- (i - 1) * num_datasets + j  # Calculate the row index

      # Store the list name, suffix (0, 1, 2 for 1/2x, 1x, or 2x missingness), element number, and element value in the result matrix
      result_matrix[row_index, 1] <- amp_name
      result_matrix[row_index, 2] <- suffix
      result_matrix[row_index, 3] <- j
      result_matrix[row_index, 4:num_cols] <- as.numeric(element)
    }
  }

  return(result_matrix)
}

# Create matrices to explore missingness and compare to original dataset ----
# Percent missing by group
missing_ck_p <- as.data.table(rbind(c("dt_tw", "NA", "NA", missing_ck_dt_tw[[1]]), missing_ck_matrix_fun(missing_ck_list, 1)))

# Number variables missing by group
missing_ck_num <- as.data.table(rbind(c("dt_tw", "NA", "NA", missing_ck_dt_tw[[2]]), missing_ck_matrix_fun(missing_ck_list, 2)))

# Quantiles of number variables missing by group
missing_ck_quant <- as.data.table(rbind(c("dt_tw", "NA", "NA", missing_ck_dt_tw[[3]]), missing_ck_matrix_fun(missing_ck_list, 3)))

# Prepare the missingness check data to summarize it ----
missing_ck_prep_sum_fun <- function(data) {
  # Set column names
  setnames(data, c("V1", "V2", "V3"), c("dt_amp", "missing_x", "dt_amp_num"))

  # Get the column names excluding "dt_amp" and "dt_amp_num"
  cols_to_convert <- setdiff(names(data), c("dt_amp", "dt_amp_num"))
  
  # Convert the selected columns to numeric using lapply and set
  lapply(cols_to_convert, function(col) set(data, j = col, value = as.numeric(data[[col]])))
  
  # Return the modified data.table
  return(data)
}

# Apply function to prep data
missing_ck_p <- missing_ck_prep_sum_fun(missing_ck_p)
missing_ck_num <- missing_ck_prep_sum_fun(missing_ck_num)
missing_ck_quant <- missing_ck_prep_sum_fun(missing_ck_quant)
# THIS CODE GENERATES NA WARNINGS FOR COMPLETE DATA, OK 07/19/24

# Average all missingness checks by missingness scenario ----
missing_ck_p_av <- missing_ck_p[, lapply(.SD, mean), by = dt_amp, .SDcols = !c("dt_amp_num")]
missing_ck_num_av <- missing_ck_num[, lapply(.SD, mean), by = dt_amp, .SDcols = !c("dt_amp_num")]
missing_ck_quant_av <- missing_ck_quant[, lapply(.SD, mean), by = dt_amp, .SDcols = !c("dt_amp_num")]

# Put data to output into a list and write excel
missingck_tables <- list()
missingck_tables[[1]] <- missing_ck_p_av
missingck_tables[[2]] <- missing_ck_num_av
missingck_tables[[3]] <- missing_ck_quant_av
missingck_tables[[4]] <- missing_ck_p
missingck_tables[[5]] <- missing_ck_num
missingck_tables[[6]] <- missing_ck_quant
names(missingck_tables) <- c("p_av", "num_av", "quant_av", "p", "num", "quant")
write.xlsx(missingck_tables, file = paste0(subfolder, "/Ext_AmpMissingCk_", today, ".xlsx"))

# Create variables to identify baseline values and start of observations in dataset ----
# Tag all time windows where date_time_end is ≥12 hours after intubation
dt_tw_comp_id[, hour_12 := fifelse(as.numeric(difftime(date_time_end, placement_instant, units = "hour"))>=12, 1, 0)]
#table(dt_tw_comp_id[hour_12==1, tw_id])
# Mark first time window where date_time_end is ≥12 hours after intubation 
# Get row numbers of the first row by patient/placement_instant
indx <- dt_tw_comp_id[hour_12==1, .I[1L], by=.(pat_mrn_id, placement_instant)]$V1
dt_tw_comp_id[indx, hour_12_1st := 1]
#table(dt_tw_comp_id[hour_12_1st==1, tw_id])
#View(dt_tw_comp_id[, .(pat_mrn_id, placement_instant, date_time_start, date_time_end, tw_id, hour_12, hour_12_1st)])

# Test/train split ----
dt_tw_comp_traintest <- copy(unique(dt_tw_comp_id[, .(pat_mrn_id, placement_instant)]))
# Determine the first placement_instant for each pat_mrn_id
first_instant <- dt_tw_comp_traintest[, .(first_placement_instant = min(placement_instant)), by = pat_mrn_id]
# Sort the pat_mrn_id by the first placement_instant
setorder(first_instant, first_placement_instant)
# Determine the split sizes
n_patients <- nrow(first_instant)
n_train <- round(n_patients * 0.75)
n_test <- n_patients - n_train
# Split pat_mrn_id into train and test sets
train_ids <- first_instant[1:n_train, pat_mrn_id]
test_ids <- first_instant[(n_train + 1):n_patients, pat_mrn_id]
# Create a data.table to map pat_mrn_id to train/test
train_test_mapping <- data.table(pat_mrn_id = first_instant$pat_mrn_id)
train_test_mapping[, traintest := ifelse(pat_mrn_id %in% train_ids, "train", "test")]
# Merge the mapping back with the original dataset
dt_tw_comp_traintest <- merge(dt_tw_comp_traintest, train_test_mapping, by = "pat_mrn_id", all.x = TRUE)
# Merge train/test variable into ID data
dt_tw_comp_id[dt_tw_comp_traintest, traintest := traintest, on = .(pat_mrn_id, placement_instant)]

print("dt_tw_id$pat_mrn_id==dt_tw_comp_id$pat_mrn_id")
all.equal(dt_tw_id$pat_mrn_id, dt_tw_comp_id$pat_mrn_id)
print("dt_tw_id$placement_instant==dt_tw_comp_id$placement_instant")
all.equal(dt_tw_id$placement_instant, dt_tw_comp_id$placement_instant)
print("dt_tw_id$date_time_start==dt_tw_comp_id$date_time_start")
all.equal(dt_tw_id$date_time_start, dt_tw_comp_id$date_time_start)

# ADD ARGUMENT TO EXCLUDE OTHER VARIABLE
# Create function to add baseline and lagged variables and preprocess for rtemis ----
preprocess_fun <- function(id_vars = dt_tw_comp_id, data, y, exclude_var) {
  # Merge in ID variables
  data <- cbind(id_vars, data)
  
  # Make outcome factor variable where first level is positive case (for rtemis)
  data[, ext_success_lag_12 := factor(ext_success_lag_12, levels=c(1, 0))]

  # Exclude variables as needed
  data[, (exclude_var) := NULL]

# Lagged variables
  # Variables to lag: summary variables, medications, intake/output, oi/osi, spontaneous breathing test
  cols_tolag <- grep("first|last|mean|mode|min|max|median|med|intake_output|oi_osi|sbt_cat", colnames(data), value=TRUE) 
    #setdiff(colnames(dt_tw), cols_tolag)
    # Note: not lagging diagnoses, arterial line presence, weight change % as these are not expected to change frequently between time windows; their most recent measure should be sufficient
  cols_tolag <- cols_tolag[!cols_tolag %in% "y_sys_bpcentile_mean_lag_4"]
  # Lag 1
  cols_lagged <- paste(cols_tolag, "lag1", sep="_")
  data[, (cols_lagged) := shift(.SD), by=c("pat_mrn_id", "placement_instant"), .SDcols=cols_tolag]
  # Lag 2
  cols_lagged <- paste(cols_tolag, "lag2", sep="_")
  data[, (cols_lagged) := shift(.SD, n=2), by=c("pat_mrn_id", "placement_instant"), .SDcols=cols_tolag]

  # Keep only time windows where date_time_end is ≥ 12 hours after intubation
  data <- data[hour_12 == 1]
  
  # Remove rows < 12 hours before extubation (keep drop_lt12h_toext==0, drop drop_lt12h_toext==1)
  data <- data[drop_lt12h_toext==0]
  
  # Remove time variables (keep ID vars: pat_mrn_id and placement_instant to ensure all end up in the same fold, keep date_time_start/end for imputation code)
  data[, c("tw_id", "hour_12", "hour_12_1st", "year") := NULL]

  # Make train and test datasets
  x_train <- data[traintest=="train", -("traintest")]
  x_test <-  data[traintest=="test", -("traintest")]

  # Create ID vector
  # Sort by descending date_time_start so the last row for each ID is first 
  # In rtemis resampling, the outcome is stratified on and the first value per ID is used for stratification
  setorder(x_train, pat_mrn_id, placement_instant, -date_time_start)
  x_train[ , id:=.GRP, by = c("pat_mrn_id", "placement_instant")]
  id_train <- x_train$id

  # Sort x_test the same way so it's the same as x_train to minimize confusion
  setorder(x_test, pat_mrn_id, placement_instant, -date_time_start)
  
  y_train <- x_train[[y]]
  x_train[, (y) := NULL]

  y_test <- x_test[[y]]
  x_test[, (y) := NULL]

  # Remove ucsf_bcho, id
  # Keep pat_mrn_id and placement_instant for imputation
  x_train[, c("ucsf_bcho", "id", "drop_lt12h_toext") := NULL]
  x_test[, c("ucsf_bcho", "drop_lt12h_toext") := NULL]
  
  data_list <- list(x_train, y_train, x_test, y_test, id_train, y)
  names(data_list) <- c("x_train", "y_train", "x_test", "y_test", "id_train", "y")

  return(data_list)
}

# Apply the processing function to all the datasets ----
# Apply to completed dataset which will serve as reference
dt_tw_comp_p <- preprocess_fun(data = dt_tw_comp, y = "ext_success_lag_12", exclude_var = "y_sys_bpcentile_mean_lag_4")
saveRDS(dt_tw_comp_p, file = file.path(subfolder_ext, "dt_tw_comp_tt_ext.RDS"))

dt_tw_comp_p <- preprocess_fun(data = dt_tw_comp, y = "y_sys_bpcentile_mean_lag_4", exclude_var = "ext_success_lag_12")
saveRDS(dt_tw_comp_p, file = file.path(subfolder_bp, "dt_tw_comp_tt_bp.RDS"))

# Apply to original dataset (needed to check LOCF % later)
# Checked dt_tw_id and dt_tw_comp_id are equal above
dt_tw_o_p <- preprocess_fun(data = dt_tw, y = "ext_success_lag_12", exclude_var = "y_sys_bpcentile_mean_lag_4")
saveRDS(dt_tw_o_p, file = file.path(subfolder_ext, "dt_tw_orig_tt_ext.RDS"))

dt_tw_o_p <- preprocess_fun(data = dt_tw, y = "y_sys_bpcentile_mean_lag_4", exclude_var = "ext_success_lag_12")
saveRDS(dt_tw_o_p, file = file.path(subfolder_bp, "dt_tw_orig_tt_bp.RDS"))

# Loop through each list in dt_names_amp and store the processed datasets as separate objects in the global environment
# Loop for extubation
for (i in 1:length(dt_amp_names_ext)) {
  # Get the list name
  list_name <- dt_amp_names_ext[i]
  
  # Access the list using the name
  current_list <- get(list_name)
  
  # Apply the process_dataset function to each dataset in the current list
  processed_list <- lapply(current_list, function(x) preprocess_fun(data = x, y = "ext_success_lag_12", exclude_var = "y_sys_bpcentile_mean_lag_4"))
  
  # Create a unique object name for the processed list
  object_name <- gsub("_ext", "_tt_ext", dt_amp_names_ext[i])
  
  # Assign the processed list as a separate R object in the global environment
  assign(object_name, processed_list, envir = .GlobalEnv)

  # Create the filename for the RDS file in the subfolder
  filename <- file.path(subfolder_ext, paste(object_name, ".RDS", sep = ""))
    
  # Save the object as an RDS file
  saveRDS(get(object_name), file = filename)
}

# Loop for BP
for (i in 1:length(dt_amp_names_bp)) {
  # Get the list name
  list_name <- dt_amp_names_bp[i]
  
  # Access the list using the name
  current_list <- get(list_name)
  
  # Apply the process_dataset function to each dataset in the current list
  processed_list <- lapply(current_list, function(x) preprocess_fun(data = x, y = "y_sys_bpcentile_mean_lag_4", exclude_var = "ext_success_lag_12"))
  
  # Create a unique object name for the processed list
  object_name <- gsub("_bp", "_tt_bp", dt_amp_names_bp[i])
  
  # Assign the processed list as a separate R object in the global environment
  assign(object_name, processed_list, envir = .GlobalEnv)

  # Create the filename for the RDS file in the subfolder
  filename <- file.path(subfolder_bp, paste(object_name, ".RDS", sep = ""))
    
  # Save the object as an RDS file
  saveRDS(get(object_name), file = filename)
}
