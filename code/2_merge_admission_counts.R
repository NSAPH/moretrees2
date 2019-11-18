# merge denominator and admission data
library(NSAPHutils)
set_threads()

library(data.table)
library(fst)


denominator_columns <- c("zip", "year", "qid", "race", "sex", "age",
                         "dual", "dead", "pm25_ensemble", "tmmx",
                         "whatever else you think you need")
admission_columns <- c("qid", "condition1", "condition2", "whatever")
conditions

denom_path <- "../data/denominator/"
admission_counts <- "../data/admission_counts/"

for (year_ in 2000:2014) {
  patient_summary <- read_data(denom_path, years = year_, 
                               columns = denominator_columns)
  admissions <- read_data(admission_counts, years = year_, 
                          columns = admission_columns)
  merged <- merge(patient_summary, admissions, all.x = T, by = "qid")
  
  # memory management, don't need the other data in memory any more
  rm(patient_summary)
  rm(admissions)
  
  merged[is.na(condition1), condition1 := 0]
  merged[is.na(condition2), condtion2 := 0]
  # etc etc
  
  write_fst(merged, paste0("../data/merged_admission_counts/denom_admission_counts_",year_,".fst")
  
  
}