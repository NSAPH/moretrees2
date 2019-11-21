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

denominator_columns <- c("zip", "year", "qid")
                         #"race", "sex", "age",
                         #"dual", "dead", "pm25_ensemble", "tmmx",
                         #"whatever else you think you need")
admission_columns <- c("qid", codes)

denom_path <- "../data/denominator/"
admission_counts <- "../data/admission_counts/"

for (year_ in 2000:2014) {
  patient_summary <- read_data(denom_path, years = year_, 
                               columns = denominator_columns)
  admissions <- read_data(admission_counts, years = year_, 
                          columns = admission_columns)
  merged <- merge(patient_summary, admissions, all.y = T, by = "qid")
  
  # aggregate by zip code
  denom <- dcast(patient_summary, zip ~ year_, fun.aggregate = length, drop = F)
  
  # memory management, don't need the other data in memory any more
  rm(patient_summary)
  rm(admissions)
  
  merged[is.na(condition1), condition1 := 0]
  merged[is.na(condition2), condtion2 := 0]
  # etc etc
  
  write_fst(merged, paste0("../data/merged_admission_counts/denom_admission_counts_",year_,".fst")
  
  
}