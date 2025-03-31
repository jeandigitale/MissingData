# F31_MissingData_1_CompleteData
## Explore missingness in original data
## Create completed dataset
## Cluster analysis to group variables
## Make variables numeric and remove ID variables to prepare for ampute
## Returns:
    # dt_tw_1: original dataset with rows and variables to match completed dataset
    # dt_tw_comp_1: completed dataset ready for ampute function

# Libraries ----
library(tidyverse)
library(openxlsx)
library(data.table)
library(knitr)
library(mice)
library(zoo)
library(rtemis)

# Set-up -----
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
filename_all <- paste0(folder, "Data_Clean/", "Ext_tw", timewindow_hrs, "_lag", lag_hrs, "_", file_date)
filename_dt_tw <- paste0(folder, "Data_Clean/", "Ext_tw", timewindow_hrs, "_lag", lag_hrs, "_dt_tw_", file_date)
filename_cohort <- paste0(folder, "Data_Clean/", "Ext_tw", timewindow_hrs, "_lag", lag_hrs, "_cohort_", file_date)

# Filename out
filename_out <- paste0(folder, "MissingData_Result/1_CompleteData/dt_tw_1_", today, ".RDS")
filename_out_comp <- paste0(folder, "MissingData_Result/1_CompleteData/dt_tw_comp_1_", today, ".RDS")

# File path log directory
log_dir <- file.path(folder, "MissingData_Log/1_CompleteData", today)
dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)

# Load data ----
load(paste0(filename_dt_tw,"_missing.RData"))
load(paste0(filename_cohort,"_missing.RData"))

# Characterize cohort
# N intubations
nrow(cohort)
# N patients
length(unique(cohort$pat_mrn_id))
# Sex
# Keeping the first row for each 'pat_mrn_id'
first_rows <- cohort[, .SD[1], by = "pat_mrn_id"]
# Create a frequency table for 'sex' from the first rows
table(first_rows$sex)
prop.table(table(first_rows$sex))
# Outcome
table(cohort$ext_outcome)
prop.table(table(cohort$ext_outcome))
# Age
summary(cohort$age_int_days) # days
summary(cohort$age_int_days/365.25) # years
sd(cohort$age_int_days/365.25)
# Length of intubation
summary(cohort$int_length_days)
sd(cohort$int_length_days)

# Sort rows to ensure they are in correct order by pat_mrn_id
setorder(dt_tw, pat_mrn_id, placement_instant, date_time_start)

# Remove unlagged outcome variable (keep only lagged outcome)
dt_tw <- dt_tw[, ext_success:=NULL]
# Create a complete variable that tags when ext_success_lag_12 is missing - these observations should be dropped
dt_tw[, drop_lt12h_toext := ifelse(is.na(ext_success_lag_12), 1, 0)]

# Remove variables to decrease the number so mice will run ----
# Ultimately, this code proceeds in steps and only keeps MEAN/MODE variables in the end
# Remove race - low variable importance, poor justification for this model, high amount of missingness
dt_tw[, race_c := NULL]

# Remove variables with minimal variation within time windows
# >90% of min/max/first/last = mean: /Users/jeandigitale/Dropbox/Dropbox_JD/UCSF/Extubation/Analysis/mice_subtests.xlsx
# Kept ventilator settings/measures even when little variation due to weaning being potentially important and to keep categories consistent
# Some respiratory pattern and breath sounds indicator variables also higher variation but deleted variables other than mean to keep consistent for all such indicator variables 
prefixes <- c("respiratory_pattern", "vent_aprv_hfov", "breath_sounds", "secretion_color", "secretion_consistency", "cough", "wbc", "le_motor", "ue_motor")

# Initialize an empty vector to store the columns to remove
cols_to_remove <- c()

# Loop through each prefix to identify columns that are not _mean
for (prefix in prefixes) {
  not_mean_cols <- grep(paste0("^", prefix, ".*_(min|max|first|last|median)$"), names(dt_tw), value = TRUE)
  cols_to_remove <- c(cols_to_remove, not_mean_cols)
}

# Remove the unwanted columns from the data.table
dt_tw[, (cols_to_remove) := NULL]

# Identify columns containing the word median, first, last
cols_to_remove <- grep("median|first|last", names(dt_tw), value = TRUE)
# Columns to keep
cols_to_keep <- c("ett_cuff_leak_pressure_cat_last", "alternative_gases_last_locf")
# Update the list by removing the columns to keep
cols_to_remove <- setdiff(cols_to_remove, cols_to_keep)
# Remove the unwanted columns from the data.table
dt_tw[, (cols_to_remove) := NULL]

# Identify columns containing _min or _max
cols_to_remove <- grep("_min|_max", names(dt_tw), value = TRUE)
# Remove the unwanted columns from the data.table
dt_tw[, (cols_to_remove) := NULL]

# Original data: Explore missing data by site} ----
dt_tw_sum <- dt_tw[, lapply(.SD, function(x) sum(is.na(x))), by=.(ucsf_bcho)]
dt_tw_sum <- data.table::transpose(dt_tw_sum, keep.names = "var", make.names = "ucsf_bcho")
ucsf_tot_n <- sum(dt_tw$ucsf_bcho=="UCSF")
bcho_tot_n <- sum(dt_tw$ucsf_bcho=="BCHO")
names(dt_tw_sum)[2] <- "bcho_miss_n"
names(dt_tw_sum)[3] <- "ucsf_miss_n"
dt_tw_sum <- dt_tw_sum[-1,]
dt_tw_sum[, bcho_miss_n := as.numeric(bcho_miss_n)]
dt_tw_sum[, ucsf_miss_n := as.numeric(ucsf_miss_n)]
dt_tw_sum[, bcho_miss_p := bcho_miss_n/bcho_tot_n]
dt_tw_sum[, ucsf_miss_p := ucsf_miss_n/ucsf_tot_n]

# Original data: Explore missing data by year ----
# Create year variable for when ETT was placed
dt_tw[, year := year(placement_instant)]
# Sum number missing by site and year
dt_tw_sum_yr_n <- dt_tw[, lapply(.SD, function(x) sum(is.na(x))), by=.(ucsf_bcho, year)][order(year, ucsf_bcho)]
# Sum number rows by site and year for denominator, merge into missing data summary
dt_tw_sum_yr_nrow <- dt_tw[, .N, by=.(ucsf_bcho, year)][order(year, ucsf_bcho)]
dt_tw_sum_yr_p <- copy(dt_tw_sum_yr_n)
dt_tw_sum_yr_p[dt_tw_sum_yr_nrow, N := N, on=c("year", "ucsf_bcho")]
# Divide to get percent missing
cols <- setdiff(colnames(dt_tw_sum_yr_p), c("year", "ucsf_bcho", "N"))
dt_tw_sum_yr_p[, (cols) := lapply(.SD, function(d) d/get('N')), .SDcols = cols]
dt_tw_sum_yr_p[, N := NULL]
# Create variables that hold what will be column names after transposed
dt_tw_sum_yr_n[, ucsf_bcho_yr_miss_n := paste(tolower(ucsf_bcho), year, "miss_n", sep = "_")][ ,`:=`(ucsf_bcho = NULL, year = NULL)]
dt_tw_sum_yr_p[, ucsf_bcho_yr_miss_p := paste(tolower(ucsf_bcho), year, "miss_p", sep = "_")][ ,`:=`(ucsf_bcho = NULL, year = NULL)]
# Transpose, merge, make columns numeric
dt_tw_sum_yr_n <- data.table::transpose(dt_tw_sum_yr_n, keep.names = "var", make.names = "ucsf_bcho_yr_miss_n")
dt_tw_sum_yr_p <- data.table::transpose(dt_tw_sum_yr_p, keep.names = "var", make.names = "ucsf_bcho_yr_miss_p")
dt_tw_sum_yr <- merge(dt_tw_sum_yr_p, dt_tw_sum_yr_n, on="var")
cols <- setdiff(colnames(dt_tw_sum_yr), "var")
dt_tw_sum_yr[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]

# Create completed dataset with interpolation for numeric variables and nearest non-missing observation ----
# Imputation to complete initial dataset to act as "TRUTH" for missing data analysis
# Linear interpolation before and after for numerics
# Then identify nearest non-missing value and fill in for remaining missing numerics, numerics that should only be integers, and non-numerics

# Copy original dataset
dt_tw_comp <- copy(dt_tw)

# Get numeric column names, then exclude variables that should be integers so that a decimal is not interpolated
colnames_num <- dt_tw_comp[, names(.SD), .SDcols = is.numeric]
prefixes <- c('glasgow_coma_scale_score', 'state_behavioral_scale', 'cough', 'breath_sounds', 'secretion_color', 'secretion_consistency', 'respiratory_pattern', 'vent_aprv_hfov', 'ue_motor_response', 'le_motor_response', 'secretion')
suffixes <- c('first', 'last', 'min', 'max')
colnames_num <- grep(paste0("^(", paste(prefixes, collapse = "|"), ").*(", paste(suffixes, collapse = "|"), ")$"), colnames_num, value = TRUE, invert = TRUE)

# Fill in with linear interpolation for numeric variables excluding integers
dt_tw_comp[, (colnames_num) := lapply(.SD, function(x) zoo::na.approx(x, na.rm=FALSE)), .SDcols = colnames_num, by = .(pat_mrn_id, placement_instant)]

# Get colnames to apply LOCF/NOCB
colnames <- setdiff(names(dt_tw_comp), c("pat_mrn_id", "date_time_start", "date_time_end", "placement_instant", "ext_success_lag_12", "tw_id", "year", "ucsf_bcho", "drop_lt12h_toext"))

# Identify nearest non-missing value and fill in
for (x in colnames) {
  dt_tw_comp[is.na(get(x)), (x) := 
    dt_tw_comp[!is.na(get(x))][.SD, on = .(pat_mrn_id, placement_instant, date_time_start), roll = "nearest", get(x)]]
}

# Completed data: Explore missing data by site ----
dt_tw_comp_sum <- dt_tw_comp[, lapply(.SD, function(x) sum(is.na(x))), by=.(ucsf_bcho)]
dt_tw_comp_sum <- data.table::transpose(dt_tw_comp_sum, keep.names = "var", make.names = "ucsf_bcho")
ucsf_tot_n <- sum(dt_tw_comp$ucsf_bcho=="UCSF")
bcho_tot_n <- sum(dt_tw_comp$ucsf_bcho=="BCHO")
names(dt_tw_comp_sum)[2] <- "bcho_miss_n"
names(dt_tw_comp_sum)[3] <- "ucsf_miss_n"
dt_tw_comp_sum <- dt_tw_comp_sum[-1,]
dt_tw_comp_sum[, bcho_miss_n := as.numeric(bcho_miss_n)]
dt_tw_comp_sum[, ucsf_miss_n := as.numeric(ucsf_miss_n)]
dt_tw_comp_sum[, bcho_miss_p := bcho_miss_n/bcho_tot_n]
dt_tw_comp_sum[, ucsf_miss_p := ucsf_miss_n/ucsf_tot_n]

# Completed data: Explore missing data by year ----
# Create year variable for when ETT was placed
dt_tw_comp[, year := year(placement_instant)]
# Sum number missing by site and year
dt_tw_comp_sum_yr_n <- dt_tw_comp[, lapply(.SD, function(x) sum(is.na(x))), by=.(ucsf_bcho, year)][order(year, ucsf_bcho)]
# Sum number rows by site and year for denominator, merge into missing data summary
dt_tw_comp_sum_yr_nrow <- dt_tw_comp[, .N, by=.(ucsf_bcho, year)][order(year, ucsf_bcho)]
dt_tw_comp_sum_yr_p <- copy(dt_tw_comp_sum_yr_n)
dt_tw_comp_sum_yr_p[dt_tw_comp_sum_yr_nrow, N := N, on=c("year", "ucsf_bcho")]
# Divide to get percent missing
cols <- setdiff(colnames(dt_tw_comp_sum_yr_p), c("year", "ucsf_bcho", "N"))
dt_tw_comp_sum_yr_p[, (cols) := lapply(.SD, function(d) d/get('N')), .SDcols = cols]
dt_tw_comp_sum_yr_p[, N := NULL]
# Create variables that hold what will be column names after transposed
dt_tw_comp_sum_yr_n[, ucsf_bcho_yr_miss_n := paste(tolower(ucsf_bcho), year, "miss_n", sep = "_")][ ,`:=`(ucsf_bcho = NULL, year = NULL)]
dt_tw_comp_sum_yr_p[, ucsf_bcho_yr_miss_p := paste(tolower(ucsf_bcho), year, "miss_p", sep = "_")][ ,`:=`(ucsf_bcho = NULL, year = NULL)]
# Transpose, merge, make columns numeric
dt_tw_comp_sum_yr_n <- data.table::transpose(dt_tw_comp_sum_yr_n, keep.names = "var", make.names = "ucsf_bcho_yr_miss_n")
dt_tw_comp_sum_yr_p <- data.table::transpose(dt_tw_comp_sum_yr_p, keep.names = "var", make.names = "ucsf_bcho_yr_miss_p")
dt_tw_comp_sum_yr <- merge(dt_tw_comp_sum_yr_p, dt_tw_comp_sum_yr_n, on="var")
cols <- setdiff(colnames(dt_tw_comp_sum_yr), "var")
dt_tw_comp_sum_yr[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]

# REMOVE BCHO PATIENTS ----
dt_tw <- dt_tw[ucsf_bcho=="UCSF"]
dt_tw_comp <- dt_tw_comp[ucsf_bcho=="UCSF"]

# Cluster analysis of missing variables ----
# Copy dt_tw to exclude ID variables from cluster analysis
dt_tw_cluster <- copy(dt_tw)
dt_tw_cluster[, c("year", "pat_mrn_id", "placement_instant", "date_time_start", "date_time_end", "tw_id", "drop_lt12h_toext", "ucsf_bcho") := NULL]
# Create matrix copy of dt_tw where 1 if missing
missing_mat <- matrix(as.numeric(is.na(dt_tw_cluster)), nrow=nrow(dt_tw_cluster), ncol=ncol(dt_tw_cluster), dimnames=list(NULL, names(dt_tw_cluster)))
# Transpose
missing_mat <- t(missing_mat)
set.seed(730)
mod_kmeans <- c_KMeans(missing_mat, k=5)
missing_mat <- as.data.table(missing_mat, keep.rownames=TRUE)
setnames(missing_mat, "rn", "var")
missing_mat[, cluster := mod_kmeans$clusters.train]
# 1: Always observed
missing_mat[grep("ext_success|age_int_days|wfaz_bl|ett_cuff_leak|fio2|hrs_int|med_|pulse|resp_|sex|spo2|temp|intake_output|bp|map|dx", var), group := 1]
# 2: Respiratory
missing_mat[grep("alternative_gases|etco2|mean_airway_pressure|peep|pip|cuffed|ett_size|present_on_hosp|pressure_support|rr_set|total_rr|oi_osi|vent_aprv|exhaled_vt|sbt|ett_placed", var), group := 2]
# 3: Blood gas
missing_mat[grep("base_excess|hgb|lactate|pco2|ph|sample_type|wbc", var), group := 3]
# 4: Other assessment
missing_mat[grep("cough|delirium|glasgow|motor|state_behavioral|weight_change_p|hfaz_bl|wflz_bmiz_bl", var), group := 4]
# 5: Resp assessment
missing_mat[grep("breath_sounds|respiratory_pattern|secretion_", var), group := 5]
setcolorder(missing_mat, c("var", "cluster", "group"))

# Write spreadsheet ----
# Create matrix of missingness patterns in completed dataset
#md_pattern_comp <- md.pattern(dt_tw_comp)
# Put data to output into a list
missingdata_tables <- list()
missingdata_tables[[1]] <- dt_tw_sum
missingdata_tables[[2]] <- dt_tw_sum_yr
missingdata_tables[[3]] <- dt_tw_comp_sum
missingdata_tables[[4]] <- dt_tw_comp_sum_yr
#missingdata_tables[[5]] <- md_pattern_comp
missingdata_tables[[6]] <- missing_mat[, c("var", "cluster", "group")]
names(missingdata_tables) <- c("site_miss", "site_year_miss", "site_comp", "site_year_comp", "md_pattern_comp", "missing_mat_sum")
write.xlsx(missingdata_tables, file = paste0(folder, "/MissingData_Result/1_CompleteData/Ext_MissingData_", today, ".xlsx"))

# Completed dataset: Fill in variables with too much missingness and keep only complete cases ----
# Count number of missing per row
# dt_tw_comp[, num_miss := rowSums(is.na(dt_tw_comp))]
# table(dt_tw_comp$num_miss)

# Quantify missingness after filling in with interpolation and nearest non-missing value
n_cell <- nrow(dt_tw_comp)* ncol(dt_tw_comp)
n_cell_miss <- sum(is.na(dt_tw_comp))
p_miss <- n_cell_miss/n_cell

# Fill in present on hospital admission using assumptions 
# In main dataset, only filled in as NOT present on admission if placement instant was within 4 hours of hospital admission
# Here, will fill in as NOT present on admission if placement instant was >=1 hour from hospital admission time
# Will fill in as present on admission if placement instant was <1 hour from hospital admission time
# Calculate time from hospital admission to placement_instant
dt_tw_comp[cohort, placement_hosp_adsm_diff := placement_hosp_adsm_diff, on = .(pat_mrn_id, placement_instant)]
# If time from hospital admission to placement_instant >= 1 hours, then make var = 0
dt_tw_comp[placement_hosp_adsm_diff>=(1/24) & is.na(present_on_hosp_admission), present_on_hosp_admission := 0]
# If time from hospital admission to placement_instant < 1 hours, then make var = 1
dt_tw_comp[placement_hosp_adsm_diff<(1/24) & is.na(present_on_hosp_admission), present_on_hosp_admission := 1]
# Remove placement_hosp_adsm_diff
dt_tw_comp[, placement_hosp_adsm_diff := NULL]
# Fill in ETT size with formula
# https://www.uptodate.com/contents/image?imageKey=EM%2F72932&topicKey=EM%2F6316&source=see_link
# https://www.ncbi.nlm.nih.gov/books/NBK539747/#:~:text=%5B5%5D%20The%20typical%20depth%20of,somewhat%20an%20institution%20dependent%20practice
# For now assumed cuffed (95% are cuffed) if missing/fill in formula for size
# 3.5 + (age in years/4) for cuffed, 4 + (age in years/4) for uncuffed
# Max size to be filled in: 7.5 for male, 7.0 for female
# To round to nearest 0.5: round(formula/0.5)*0.5 
dt_tw_comp[is.na(cuffed), cuffed := 1]
dt_tw_comp[cuffed==1 & is.na(ett_size), ett_size_fill := round((3.5 + (age_int_days/365.25/4))/0.5)*0.5]
dt_tw_comp[cuffed==0 & is.na(ett_size), ett_size_fill := round((4 + (age_int_days/365.25/4))/0.5)*0.5]
dt_tw_comp[ett_size_fill > 7 & sex=="Female", ett_size_fill := 7]
dt_tw_comp[ett_size_fill > 7.5 & sex=="Male", ett_size_fill := 7.5]
dt_tw_comp[is.na(ett_size), ett_size := ett_size_fill]
dt_tw_comp[, ett_size_fill := NULL]
#View(cohort[is.na(cuffed)|is.na(ett_size), .(cuffed, ett_size,description)])
#View(dt_tw_comp[is.na(ett_size), .(age_int_days, ett_size_ck)])

# Fill in with random forest/PMM remaining missing data ----
# Note: this will fill in ext_success_lag_12 where it should be missing (time windows <12 hr before extubation), these are filtered out using drop_lt12h_toext in the next R script

# Save original dataset ----
saveRDS(dt_tw, filename_out)

# SGE_submit to pre-process and complete dataset ----
progressr::handlers(global = TRUE)
options(progressr.handlers = progressr::handler_cli)

sge_opts_note <<- paste0("#$ -cwd -N '", "comp_data")

mod <- sge_submit(
  {
    set.seed(2876)
    exclude_vars <- which(names(dt_tw_comp) %in% c("year", "pat_mrn_id", "placement_instant", "date_time_start", "date_time_end", "tw_id", "drop_lt12h_toext", "ucsf_bcho"))
    dt_tw_comp_1 <- preprocess(dt_tw_comp, impute=TRUE)
    saveRDS(dt_tw_comp_1, filename_out_comp)
  },
  obj_names = c("dt_tw_comp", "filename_out_comp"),
  packages = c("rtemis", "data.table"),
  n_threads = 12,
  sge_out = log_dir,
  sge_opts = sge_opts_note,
  h_rt = "20:00:00",
  system_command = "module load CBI r"
)
