# Code to get admission counts by person, by year for medicare data

library(NSAPHutils)
set_threads()

library(data.table)
library(fst)
library(icd)

# Get list of relevant ICD9 codes at billable level
chaps <- icd9_chapters[c("Mental Disorders", "Diseases Of The Nervous System And Sense Organs")]
codes <- c(expand_range(chaps[[1]][1], chaps[[2]][1]))
codes <- children(codes, billable = T)

admissions <- "../data/admissions"

for (year_ in 2000:2014) {
  admission_data <- read_data(admissions, years = year_, columns = c("QID", "DIAG1"))
  names(admission_data) <- tolower(names(admission_data))
  
  # This creates a new column that is true or false for whether the admission is
  # part of the ICD code set
  sapply(codes, function(code) admission_data[, paste(code) := diag1 == code])
  
  # This creates annual counts of how many hospitalizations for each individual occured
  # The .SD notation is a bit confusing so please ask me about it
  
  admission_data <- admission_data[,lapply(.SD, sum, na.rm = T), by = "qid",
                                   .SDcols = codes]
  admission_data[, year := year_]
  
  write_fst(admission_data, paste0("../data/admission_counts/admission_counts_", year_, ".fst"))
}