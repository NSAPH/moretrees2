# Code to get admission counts by person, by year for medicare data

# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(NSAPHutils)
set_threads()

library(data.table)
library(fst)
library(icd)
library(magrittr)

# Get list of relevant ICD9 codes at billable level
codes <- icd9_chapters$"Diseases Of The Circulatory System"
codes <- expand_range(codes[1], codes[2])
codes <- get_leaf(codes)

# Get data.frame showing mapping from ICD9 to CCS
ccs_list <- icd9_map_multi_ccs[[4]]
ccs_icd9 <- data.frame(icd9 = unlist(ccs_list), 
                       stringsAsFactors = F)
row.names(ccs_icd9) <- NULL
ccs_icd9$ccs <- ccs_list %>% 
  names %>%
  sapply(FUN = function(nm) rep(nm, length(ccs_list[[nm]]))) %>%
  unlist
ccs_icd9 <- subset(ccs_icd9, icd9 %in% codes)

admissions <- "../data/admissions"
# denom_path <- "../data/denominator/"

admissions_columns <- c("QID", "DIAG1", "ADATE", "zipcode_R")
# denom_columns <- c("qid", "zip")

for (year_ in 2000:2014) {
  admission_data <- read_data(admissions, years = year_, columns = admissions_columns)
  names(admission_data) <- tolower(names(admission_data))
  
  # Keep only relevant diagnoses
  admission_data <- subset(admission_data, diag1 %in% codes)
  
  # Merge in ccs codes
  admission_data <- merge(admission_data, ccs_icd9, by.x = "diag1",
                          by.y = "icd9", all.x = T)
  
  write_fst(counts, paste0("../data/admission_counts/admission_counts_", year_, ".fst"))
}

