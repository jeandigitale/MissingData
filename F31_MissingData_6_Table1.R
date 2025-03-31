# F31_MissingData_6_Table1
# Create Table 1 to describe sample for missing data paper

# Libraries ----
library(data.table)
library(tableone)
library(labelled)

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

# Filename in
filename_in <- "filepath"

filename_in_dt <- "filepath"

# Filename out
subfolder <- file.path(folder, "MissingData_Result/6_Table1", today)
dir.create(subfolder, recursive = TRUE, showWarnings = FALSE)

# Load cohort data ----
load(filename_in)
dt_tw_tt <- readRDS(filename_in_dt)
dt_tw <- rbind(dt_tw_tt$x_train, dt_tw_tt$x_test)

# Check which patients remain in dataset after dropping time windows ----
dt_tw_cohort <- unique(copy(dt_tw[, .(pat_mrn_id, placement_instant)]))

# Ensure both data.tables have the same key columns set for joining
setkey(cohort, pat_mrn_id, placement_instant)
setkey(dt_tw_cohort, pat_mrn_id, placement_instant)

# Perform an inner join to keep only matching rows
final_cohort <- cohort[dt_tw_cohort, nomatch=0]

# Get maximum value of diagnosis data
cohort_dx <- dt_tw[, .(max_dx_airway_comp = max(dx_airway_comp_locf, na.rm = TRUE),
                      max_dx_cardiomyopathy = max(dx_cardiomyopathy_locf, na.rm = TRUE),
                      max_dx_neuromusc = max(dx_neuromusc_locf, na.rm = TRUE),
                      max_dx_critical_airway = max(dx_critical_airway_locf, na.rm = TRUE)),
                  by = .(pat_mrn_id, placement_instant)]

# Merge cohort_dx into final_cohort
final_cohort <- final_cohort[cohort_dx, on = .(pat_mrn_id, placement_instant), nomatch = NA]

# Count number of patients ----
length(unique(final_cohort$pat_mrn_id))

# Count number of time windows in analytic dataset ----
nrow(dt_tw_tt$x_train)+nrow(dt_tw_tt$x_test)

# Convert age from days to years
final_cohort[, age_int_years := age_int_days / 365.25]

# Create new variable 'race_eth'
final_cohort[, race_eth := ifelse(ethnicity == "Hispanic or Latino", "Latinx", race_c)]

# Convert 'sex' from character to factor with "Female" as the first level
final_cohort[, sex := factor(sex, levels = c("Male", "Female"))]

# Map existing ext_success values to new labels
outcome_levels <- setNames(c("Success", "Failure", "Death", "Tracheostomy", "Transfer to another unit", "ETT change"),
                           c("success", "failure", "death", "trach", "transfer", "tubechange"))

# Convert 'ext_outcome' to a factor with the specified levels and labels
final_cohort[, ext_outcome := factor(ext_outcome, levels = names(outcome_levels), labels = outcome_levels)]

# Create patient-level dataset
# Keep only one row per 'pat_mrn_id' (keeping the first occurrence)
cohort_pt <- copy(final_cohort)
cohort_pt <- unique(cohort_pt, by = "pat_mrn_id")

# Specify variables to include in Table 1
vars_pt <- c("sex", "race_eth")
vars_int <- c("age_int_years", "ext_outcome", "int_length_days")

# Adding labels to the variables
var_label(cohort_pt$sex) <- "Sex"
var_label(cohort_pt$race_eth) <- "Race/Ethnicity"
var_label(final_cohort$age_int_years) <- "Age at intubation (years)"
var_label(final_cohort$ext_outcome) <- "Extubation outcome"
var_label(final_cohort$int_length_days) <- "Length of intubation (days)"

# Create Table 1 using the CreateTableOne function: patient level
table1_pt <- CreateTableOne(vars = vars_pt, data = cohort_pt)

# Convert Table 1 to a data frame
table1_pt_df <- print(table1_pt, printToggle = FALSE, varLabels = TRUE, nonnormal = TRUE)

# Create Table 1 using the CreateTableOne function: intubation level
table1_int <- CreateTableOne(vars = vars_int, data = final_cohort)

# Convert Table 1 to a data frame
table1_int_df <- print(table1_int, printToggle = FALSE, noSpaces = TRUE, varLabels = TRUE, contDigits = 1, nonnormal = TRUE)

# Combine tables
# Convert matrices to data.table, including row names as a new column
table1_pt_dt <- data.table(table1_pt_df, keep.rownames = "Variable")
table1_int_dt <- data.table(table1_int_df, keep.rownames = "Variable")
# Combine the data tables
table1_df <- rbindlist(list(table1_pt_dt, table1_int_dt), use.names = TRUE)

# Export to CSV
write.csv(table1_df, paste0(subfolder, "/Table1_MissingData.csv"), row.names = FALSE)
