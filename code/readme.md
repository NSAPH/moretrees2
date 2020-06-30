A Bayesian multi-outcome analysis of fine particulate matter and cardiorespiratory hospitalizations
================
Emma Thomas
30/06/2020

This directory contains code for processing data and performing analyses for the paper "A Bayesian multi-outcome analysis of fine particulate matter and cardiorespiratory hospitalizations".
Scripts should be run in order of the digit at the start of the file name, i.e., scripts named 0_* should be run first, followed by 1_*, 2_*, and so on. Scripts with the same number can be run in any order.

# Files

  - `0_unique_qids_mapping.R`
      - Generates unique numeric ID for each Medicare beneficiary in the dataset.
  - `1_get_admissions.R`
      - Extracts cardiovascular and respiratory admissions for the study period from Medicare data.
  - `2_get_enviro_vars`
      - Extracts PM2.5, humidity, temperate, and ozone data for each day and ZIP code in the study period.
  - `3_merge_cvd_enviro.R`
      - Merges environmental variables with cardiovascular hospital admissions.
  - `4_merge_resp_enviro.R`
      - Merges environmental variables with respiratory hospital admissions.
  - `5_main_analysis.R`
      - Fits MOReTreeS models for main analysis (and sensitivity analysis 1).
  - `6_CLR_models.R`
      - Fits conditional logistic regression (CLR) models using groups specified at different levels of the outcome tree for comparison with MOReTreeS.
  - `6_cross_validation.R` and
      - Runs cross validation analysis comparing MOReTreeS to CLR with pre-specified groups.
  - `7_figures_and_tables_main.R`
      - Produces figures and tables presented in the main manuscript.
  - `8_sensitivityX.R`
      - Runs sensitivity analysis X were X = 2, 3, 4, 5.
  - `9_figures_supplementary.R`
      - Produces all figures shown in appendices.
  - `9_tables_sensitivityX.R`
      - Produces tables with results of sensitivity analysis X shown in appendices.
