# merge denominator and admission data
library(NSAPHutils)
set_threads()

# Make sure working directory is set to moretrees2/code
# setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(data.table)
library(fst)
library(icd)

# Get list of relevant ICD9 codes at billable level
chaps <- icd9_chapters[c("Mental Disorders", "Diseases Of The Nervous System And Sense Organs")]
codes <- c(expand_range(chaps[[1]][1], chaps[[1]][2]), 
           expand_range(chaps[[2]][1], chaps[[2]][2]))
codes <- get_leaf(codes)

enviro_columns <- c("zip", "date", "pm25_ensemble", "tmmx",
                    "humidity", "ozone")
                         #"race", "sex", "age",
                         #"dual", "dead", "pm25_ensemble", "tmmx",
                         #"whatever else you think you need")
admission_columns <- c("zipcode_R", "qid", "adate",
                       "icd9", "ccs_l1", "ccs_l2", 
                       "ccs_l3", "ccs_l4")

enviro_path <- "../data/enviro_vars/"
admissions_cvd <- "../data/admissions_cvd/"

for (year_ in 2000:2014) {
  enviro <- read_data(enviro_path, years = year_, 
                               columns = enviro_columns)
  admissions <- read_data(admissions_cvd, years = year_, 
                          columns = admission_columns)
  merged <- merge(patient_summary, admissions, all.y = T, by = "qid")
  
  # memory management, don't need the other data in memory any more
  rm(patient_summary)
  rm(admissions)
  
  merged[is.na(condition1), condition1 := 0]
  merged[is.na(condition2), condtion2 := 0]
  # etc etc
  
  write_fst(merged, paste0("../data/merged_admission_counts/denom_admission_counts_",year_,".fst")
  
  
}