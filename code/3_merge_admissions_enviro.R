
library(NSAPHutils)
set_threads()

# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(data.table)
library(fst)
library(icd)
library(stringr)
library(readr)

admission_columns <- c("zipcode_r", "qid", "adate",
                       "icd9", "ccs_l1", "ccs_l2", 
                       "ccs_l3", "ccs_l4")

pm_path <- "../data/daily_pm/"
temp_path <- "../data/temperature/"
ozone_path <- "../data/ozone/"

# function to reverse zip codes
zipReverse <- function(x){
  x <- as.character(x)
  l <- str_length(x)
  if (l != 5) {
    x <- paste0(c(rep("0", 5 - l), x), collapse = "")
  }
  x <- paste0(rev(strsplit(x, NULL)[[1]]), collapse = "")
  return(x)
}

set.seed(6357312)

for (year_ in 2000:2014) {
  # read in admissions
  admissions <- read_fst(paste0("../data/admissions_cvd/admissions_cvd_", year_, ".fst"))
  
  # select control days (randomly exactly one week before or after)
  admissions$cdate <- admissions$adate + sample(c(-7, 7), nrow(admissions), replace = T)
  
  # merge in lag01 PM2.5 data for case and control days
  admissions$pm25_case <- numeric(nrow(admissions))
  admissions$pm25_control <- numeric(nrow(admissions))
  for (i in 1:(nrow(admissions))) {
    case0 <- str_remove_all(as.character(admissions$adate[i]), pattern = "-")
    case1 <- str_remove_all(as.character(admissions$adate[i] - 1), pattern = "-")
    control0 <- str_remove_all(as.character(admissions$cdate[i]), pattern = "-")
    control1 <- str_remove_all(as.character(admissions$cdate[i] - 1), pattern = "-")
    zip <- zipReverse(admissions$zipcode_r[i])
    # case day PM2.5 lag01
    case_pm0 <- subset(read_rds(paste0(pm_path, case0, ".rds")), ZIP == zip)$pm25
    case_pm1 <- subset(read_rds(paste0(pm_path, case1, ".rds")), ZIP == zip)$pm25
    admissions$pm25_case[i] <- (case_pm0 + case_pm1) / 2
    # control day PM2.5 lag01
    control_pm0 <- subset(read_rds(paste0(pm_path, control0, ".rds")), ZIP == zip)$pm25
    control_pm1 <- subset(read_rds(paste0(pm_path, control1, ".rds")), ZIP == zip)$pm25
    admissions$pm25_control[i] <- (control_pm0 + control_pm1) / 2
    # print
    if (i %% 1000 == 0) print(i / nrow(admissions) * 100)
  }
  
  # memory management, don't need the other data in memory any more
  rm(admissions)
  
  write_fst(admissions, paste0("../data/merged_admissions_enviro/admissions_enviro_",year_,".fst")
  
  
}
