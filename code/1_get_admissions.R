# Code to get admission counts by person, by year for medicare data

# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(NSAPHutils)
set_threads()

library(data.table)
library(fst)
library(icd)
library(magrittr)
require(lubridate)

# Get data.frame showing mapping from ICD9 to multilevel CCS
ccs_icd9 <- data.table(icd9 = unlist(icd9_map_multi_ccs[[1]]))
for (i in 1:4) {
  ccs_list <- icd9_map_multi_ccs[[i]]
  ccs_dt <- data.table(icd9 = unlist(ccs_list))
  ccs_dt[ , ccs := ccs_list %>% 
           names %>% # names of the list entries are the CCS codes
           sapply(FUN = function(nm) rep(nm, length(ccs_list[[nm]]))) %>%
           unlist]
  nm <- paste0("ccs_l", i)
  setnames(ccs_dt, "ccs", nm)
  ccs_icd9 <- merge(ccs_icd9, ccs_dt, by = "icd9", all.x = T, all.y = F)
  # Fill in blanks for lower levels with previous level
  if (i > 1) {
    nm_p <- paste0("ccs_l", i - 1)
    ccs_icd9[get(nm) == " ", (nm) := .SD, .SDcols = nm_p] 
  }
}

# Keep only diseases of the circulatory system
ccs_icd9 <- subset(ccs_icd9, ccs_l1 == "7")

# Keep only 4 level code, which contains all the relevant information
ccs_icd9[ , c("ccs_l1", "ccs_l2", "ccs_l3") := NULL]
setnames(ccs_icd9, "ccs_l4", "ccs")

# Mapping from unique QIDs to integer IDs
qids <- read_fst("../data/unique_qids/qids.fst", as.data.table = T)

# Extract relevant admissions
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
  admission_data <- admission_data[adm_type %in% c(1, 2)]
  
  # Convert admission date variable to date format
  admission_data[ , adate := dmy(adate)]
  
  # Make sure date range is correct
  admission_data <- admission_data[adate <= make_date(year = 2014, month = 12, day = 31) &
                             adate >= make_date(year = 2000, month = 1, day = 1)]
  
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

    # Reverse zip codes
    zips <- data.table(zipcode_r = unique(admission_data$zipcode_r), key = "zipcode_r")
    zips <- zips[!is.na(zipcode_r)] # NA zipcodes can be left as NA
    zips[ , zip := sapply(zipcode_r, function(z) int_to_zip_str(z) %>%
                            reverse_string %>%
                            as.integer)]
    # Merge in reversed zip codes
    admission_data <- merge(admission_data, zips, by = "zipcode_r", all.x = T, all.y = F)
    admission_data[ , zipcode_r := NULL]
    
    # Merge in integer IDs to replace QIDs
    admission_data <- merge(admission_data, qids, by = "qid",
                            all.x = T, all.y = F)
    admission_data[ , qid := NULL]
    
    # Set data.table keys
    setkey(admission_data, zip, adate)
    
    # Write to file
    write_fst(admission_data, paste0("../data/admissions_cvd/admissions_cvd_", year_, ".fst"))
  }
  
  # Remove dataset
  rm(admission_data)
  
  print(year_)
}

