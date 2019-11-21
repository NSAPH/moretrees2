# Code to get admission counts by person, by year for medicare data

# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

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
denom_path <- "../data/denominator/"

denom_columns <- c("qid", "zip", "pm25_ensemble")

for (year_ in 2000:2014) {
  admission_data <- read_data(admissions, years = year_, columns = c("QID", "DIAG1"))
  names(admission_data) <- tolower(names(admission_data))
  
  # Keep only relevant diagnoses
  admission_data <- subset(admission_data, diag1 %in% codes)
  admission_data$diag1 <- factor(admission_data$diag1, levels = codes)
  
  # Get zip codes
  patient_summary <- read_data(denom_path, years = year_, columns = denom_columns)
  
  # Merge in zip codes
  admission_data <- merge(admission_data, patient_summary, by = "qid", all.x = T)
  
  # This creates new columns for each ICD9 codes with counts of principal diagnosis by zip
  admission_data$zip <- factor(admission_data$zip, levels = unique(patient_summary$zip))
  admission_data <- dcast(admission_data, zip ~ diag1, fun.aggregate = length, drop = F) 
  
  # Merge in denominator/pm25 data
  counts <- patient_summary[ , .(denom = .N, pm25 = unique(pm25_ensemble)), by = zip]
  counts$zip <- factor(counts$zip, levels = unique(patient_summary$zip))
  counts <- merge(counts, admission_data, by = "zip", all.x = T)
  # note: NA zip codes in admissions data get dropped here
  
  # print(summary(colSums(counts[ , ..codes])))
  
  # This creates annual counts of how many hospitalizations for each individual occured
  # The .SD notation is a bit confusing so please ask me about it
  
  # admission_data <- admission_data[,lapply(.SD, sum, na.rm = T), by = "qid",
  #                                  .SDcols = codes]
  # admission_data[, year := year_]
  write_fst(counts, paste0("../data/admission_counts/admission_counts_", year_, ".fst"))
}

for (year_ in 2000:2014) {
  counts <- read_fst(paste0("../data/admission_counts/admission_counts_", year_, ".fst"))
  print(year_)
  print(summary(colSums(counts[ , codes])))
}
