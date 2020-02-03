# Code to get admission counts by person, by year for medicare data

# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(NSAPHutils)
set_threads()

library(data.table)
library(fst)
library(icd)
library(magrittr)

# function to reverse zip codes
zipReverse <- function(x){
  as.character(x) %>%
    strsplit(split = NULL) %>%
    lapply(FUN = rev) %>%
    sapply(FUN = paste, collapse = "") %>%
    as.integer
} 
zipReverseVec <- function(x) {
  y <- x
  y[!is.na(y)] <- zipReverse(x[!is.na(x)])
  return(y)
}

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

leftover_admissions <- NULL

for (year_ in 2015:2000) { 
  # need to go backwards (2015 to 2000) as some admission dates for previous years are
  # in future years' datasets
  # excluding 2016 as no ICD9 codes were in use by 2016
  # will only save data from 2000 to 2014 (years ICD9 was used 100% of time)
  
  # Read in data
  admission_data <- read_data(admissions, years = year_, columns = admissions_columns)
  names(admission_data) <- tolower(names(admission_data))
  
  # Keep only relevant diagnoses
  admission_data <- admission_data[diag1 %in% ccs_icd9$icd9]
  
  # Keep only urgent/emergency hospital admissions
  admission_data <- subset(admission_data, adm_type %in% c(1, 2))
  
  # Convert admission date variable to date format
  admission_data[ , adate := as.Date(adate, "%d%b%Y")]
  
  # Make sure date range is correct
  admission_data <- admission_data[adate <= as.Date("2014-12-31") &
                             adate >= as.Date("2000-01-01")]
  
  # Reverse zip codes
  admission_data[ , zip := zipReverseVec(zipcode_r)]
  admission_data[ , zipcode_r := NULL]
  
  # Append previous leftover admissions from future years
  admission_data <- rbind(admission_data, leftover_admissions)
  
  # Collect any admission dates in previous years
  leftover_admissions <- admission_data[year(adate) < year_]
  
  if (year_ %in% 2000:2014) {
    # Save only current year's admissions
    admission_data <- admission_data[year(adate) == year_]
    
    # Merge in ccs codes
    admission_data <- merge(admission_data, ccs_icd9, by.x = "diag1",
                            by.y = "icd9", all.x = T)
    
    # Set data.table keys
    setkey(admission_data, zip, adate)
    
    # Write to file
    write_fst(admission_data, paste0("../data/admissions_cvd/admissions_cvd_", year_, ".fst"))
  }
  
  # Remove dataset
  rm(admission_data)
}

