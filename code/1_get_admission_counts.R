# Code to get admission counts by person, by year for medicare data

library(NSAPHutils)
set_threads()

library(data.table)
library(fst)

admissions <- "../data/admissions"

for (year_ in 2000:2014) {
  admission_data <- read_data(admissions, years = year_, columns = c("QID", "DIAG1"))
  names(admission_data) <- tolower(names(admission_data))
  
  # This creates a new column that is true or false for whether the admission is
  # part of the ICD code set
  admission_data[, made_up_condition := diag1 %in% c("code1", "code2", "you can use the icd range thing here")]
  
  # This creates annual counts of how many hospitalizations for each individual occured
  # The .SD notation is a bit confusing so please ask me about it
  
  condition_list <- c("made_up_condition1", "made_up_condition2")
  admission_data <- admission_data[,lapply(.SD, sum, na.rm = T), by = "qid",
                                   .SDcols = condition_list]
  admission_data[, year := year_]
  
  write_fst(admission_data, paste0("../data/admission_counts/admission_counts_", year_, ".fst"))
}