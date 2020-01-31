# Code to get admission counts by person, by year for medicare data

# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(NSAPHutils)
set_threads()

library(data.table)
library(fst)
library(icd)
library(magrittr)

# Get data.frame showing mapping from ICD9 to multilevel CCS
ccs_icd9 <- data.frame(icd9 = unlist(icd9_map_multi_ccs[[1]]), stringsAsFactors = F)
for (i in 1:4) {
  ccs_list <- icd9_map_multi_ccs[[i]]
  ccs_df <- data.frame(icd9 = unlist(ccs_list), 
                         stringsAsFactors = F)
  ccs_df$ccs <- ccs_list %>% 
    names %>% # names of the list entries are the CCS codes
    sapply(FUN = function(nm) rep(nm, length(ccs_list[[nm]]))) %>%
    unlist
  names(ccs_df)[2] <- paste0("ccs_l", i)
  ccs_icd9 <- merge(ccs_icd9, ccs_df, by = "icd9", all.x = T, all.y = F)
}

# Keep only diseases of the circulatory system
ccs_icd9 <- subset(ccs_icd9, ccs_l1 == "7")

admissions <- "../data/admissions"

admissions_columns <- c("QID", "DIAG1", "ADATE", "ADM_TYPE", "zipcode_R")

for (year_ in 2000:2014) {
  
  # Read in data
  admission_data <- read_data(admissions, years = year_, columns = admissions_columns)
  names(admission_data) <- tolower(names(admission_data))
  
  # Keep only relevant diagnoses
  admission_data <- subset(admission_data, diag1 %in% ccs_icd9$icd9)
  
  # Keep only urgent/emergency hospital admissions
  admission_data <- subset(admission_data, adm_type %in% c(1, 2))
  
  # Merge in ccs codes
  admission_data <- merge(admission_data, ccs_icd9, by.x = "diag1",
                          by.y = "icd9", all.x = T)
  
  # Write to file
  write_fst(admission_data, paste0("../data/admissions_cvd/admissions_cvd_", year_, ".fst"))
}

