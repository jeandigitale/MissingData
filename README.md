This repository contains code used in the analysis for the manuscript:

**"Comparing Methods to Handle Missing Longitudinal Electronic Health Record Data for Prediction Models"**  
Authors: Jean C. Digitale, PhD, MPH, RN, Deborah Franzon, MD, Mark J. Pletcher, MD, MPH, Charles E. McCulloch, PhD, Efstathios D. Gennatas, MBBS, PhD

---

## Abstract

**Objective:** We evaluated methods for handling missing data when using longitudinal electronic health record (EHR) data to build clinical prediction models in pediatric intensive care unit (PICU) patients.

**Materials and Methods:** Using EHR data containing missing values from an academic medical center PICU, we generated a synthetic complete dataset. From this, we created 300 datasets with missing data under varying mechanisms and proportions of missingness for the outcomes of 1) successful extubation (binary) and 2) blood pressure (continuous). We assessed strategies to address missing data including simple methods (e.g., last observation carried forward [LOCF]), complex methods (e.g., random forest multiple imputation), and native support for missing values in outcome prediction models.

**Results:** Across 886 patients and 1,220 intubation events, 18.2% of original data were missing. LOCF had the lowest imputation error, followed by random forest imputation (average mean squared error [MSE] improvement over mean imputation: 0.41 [range: 0.30, 0.50] and 0.33 [0.21, 0.43], respectively). LOCF generally outperformed other imputation methods across outcome metrics and models (mean improvement: 1.28% [range: -0.07%-7.2%]). Imputation methods showed more performance variability for the binary outcome (balanced accuracy coefficient of variation [CV]: 0.042) than the continuous outcome (MSE CV: 0.001).

**Conclusion:** In datasets with frequent measurements, LOCF and native support for missing values in machine learning models offer reasonable performance for handling missingness at minimal computational cost in predictive analyses.

---

## Repository Structure

- F31_MissingData_1_CompleteData.R – Create complete (“truth”) dataset using interpolation and LOCF

- F31_MissingData_2_Ampute.R – Simulate missingness mechanisms (MCAR, MAR, MNAR) by amputing complete data

- F31_MissingData_3a_Impute_all.R – Master script to run all imputation types

- F31_MissingData_3b_Impute_locf.R – LOCF imputation

- F31_MissingData_3b_Impute_mean.R – Mean/mode imputation

- F31_MissingData_3b_Impute_mice_rf.R – MICE with random forest

- F31_MissingData_3b_Impute_mice_pmm.R – MICE with predictive mean matching

- F31_MissingData_3b_Impute_mice_lasso.R – MICE with LASSO

- F31_MissingData_3c_Impute_datasetprep.R – Final prep of imputed datasets

- F31_MissingData_4a_Model_all.R – Run all outcome models across all imputation strategies

- F31_MissingData_4b_Model_gbm.R – Gradient boosted models

- F31_MissingData_4b_Model_glmnet.R – LASSO models

- F31_MissingData_5a_1_EvaluateImputations_Calc.R – Compute imputation error metrics

- F31_MissingData_5a_2_EvaluateImputations_GraphAll.R – Plot imputation performance

- F31_MissingData_5a_3_EvaluateImputations_ColData.R – Evaluate imputation error by variable

- F31_MissingData_5b_1_AverageMI.R – Average model results across imputations

- F31_MissingData_5b_2_EvaluateModel.R – Compare model performance

- F31_MissingData_5c_ModelsofModels.R – Model imputation and prediction performance metrics

- F31_MissingData_5d_1_EvaluateImputations_CalcStrat.R – Stratify imputation performance by original missingness

- F31_MissingData_5d_2_EvaluateImputations_CalcStratMean.R – Plot stratified imputation performance

- F31_MissingData_6_Table1.R – Generate Table 1 summary statistics

- README.md – This file

---

## Data Availability
The data used in this analysis include protected health information (PHI) and cannot be shared publicly.

---

## Software
This project was developed in R version 4.3.2.
