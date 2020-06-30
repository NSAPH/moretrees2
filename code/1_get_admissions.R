# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(NSAPHutils)
set_threads()
# devtools::install_github("emgthomas/moretrees_pkg", ref = "devel")
require(moretrees)
library(data.table)
library(fst)
library(icd)
library(magrittr)
require(lubridate)

# Get data.frame showing mapping from ICD9 to multilevel CCS for cardiovascular disease
ccs_icd_cvd <- moretrees::ccs_tree("7")$ccs_icd_mapping
names(ccs_icd_cvd) <- c("icd9", "ccs", "ccs_added_zeros")
ccs_icd_cvd$dataset <- "cvd"

# Get data.frame showing mapping from ICD9 to multilevel CCS for respiratory disease
ccs_icd_resp <- moretrees::ccs_tree("8")$ccs_icd_mapping
names(ccs_icd_resp) <- c("icd9", "ccs", "ccs_added_zeros")
ccs_icd_resp$dataset <- "resp"

# Both
ccs_icd <- rbind(ccs_icd_cvd, ccs_icd_resp)

# Mapping from unique QIDs to integer IDs
qids <- read_fst("../data/unique_qids/qids.fst", as.data.table = T)

# Extract relevant admissions
admissions <- "../data/admissions"
admissions_columns <- c("QID", "DIAG1", "ADATE", "ADM_TYPE", "zipcode_R", "SSA_STATE_CD",
                        "Race_gp", "Sex_gp", "age_gp", "Dual")
leftover_admissions <- NULL

for (year_ in 2015:2000) { 
  # need to go backwards (2015 to 2000) as some admission dates for previous years are
  # in future years' datasets
  # excluding 2016 as no ICD9 codes were in use by 2016
  # will only save data from 2000 to 2014 (years ICD9 was used 100% of time)
  
  # Read in data
  admission_data <- read_data(admissions, years = year_, columns = admissions_columns)
  names(admission_data) <- tolower(names(admission_data))
  
  # Keep only urgent/emergency hospital admissions
  admission_data <- admission_data[adm_type %in% c(1, 2)]
  
  # Keep only relevant diagnoses
  admission_data <- admission_data[diag1 %in% ccs_icd$icd9]
  
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
    admission_data <- merge(admission_data, ccs_icd, by.x = "diag1",
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
    
    # Write to file separately for CVD and respiratory disease
    write_fst(admission_data[dataset == "cvd"][ , dataset := NULL], paste0("../data/admissions_cvd/admissions_cvd_", year_, ".fst"))
    write_fst(admission_data[dataset == "resp"][ , dataset := NULL], paste0("../data/admissions_resp/admissions_resp_", year_, ".fst"))
  }
  
  # Remove dataset
  rm(admission_data)
  
  print(year_)
}

