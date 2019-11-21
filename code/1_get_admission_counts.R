# Code to get admission counts by person, by year for medicare data

# Make sure working directory is set to moretrees2/code
# setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(NSAPHutils)
set_threads()

library(data.table)
library(fst)
library(icd)

# Get list of relevant ICD9 codes at billable level
chaps <- icd9_chapters[c("Mental Disorders", "Diseases Of The Nervous System And Sense Organs")]
codes <- c(expand_range(chaps[[1]][1], chaps[[1]][2]), 
           expand_range(chaps[[2]][1], chaps[[2]][2]))
codes <- get_leaf(codes)

admissions <- "../data/admissions"

for (year_ in 2000:2014) {
  admission_data <- read_data(admissions, years = year_, columns = c("QID", "DIAG1"))
  names(admission_data) <- tolower(names(admission_data))
  
  # Keep only relevant diagnoses
  admission_data <- subset(admission_data, diag1 %in% codes
  admission_data$diag1 <- factor(admission_data$diag1, levels = codes)
  
  # This creates new columns for each ICD9 code with counts of every time each Medicare
  # beneficiary was hospitalized with that code as principal diagnosis
  admission_data <- dcast(admission_data, qid ~ diag1, fun.aggregate = length, drop = F)
  
  # This creates annual counts of how many hospitalizations for each individual occured
  # The .SD notation is a bit confusing so please ask me about it
  
  # admission_data <- admission_data[,lapply(.SD, sum, na.rm = T), by = "qid",
  #                                  .SDcols = codes]
  # admission_data[, year := year_]
  
  write_fst(admission_data, paste0("../data/admission_counts/admission_counts_", year_, ".fst"))
}
