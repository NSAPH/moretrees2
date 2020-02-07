
library(NSAPHutils)
set_threads()

# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(data.table)
library(fst)
library(icd)
library(stringr)
library(readr)

admissions_path <- "../data/admissions_cvd/"
enviro_path <- "../data/temperature/"
merged_path <- "../data/merged_admissions_enviro/"

# set seed because we will randomly select control day
set.seed(6357312)

for (year_ in 2000:2014) {
  # read in admissions
  admissions <- read_fst(paste0("../data/admissions_cvd/admissions_cvd_", year_, ".fst"),
                        as.data.table = T, columns = c("zip", "adate", "qid", "ccs_l4"))
  
  # select control days (randomly exactly one week before or after)
  admissions$cdate <- admissions$adate + sample(c(-7, 7), nrow(admissions), replace = T)
  
  # read in PM2.5 data
  pm25 <- read_fst(paste0("../data/enviro/pm_", year_, ".fst"), 
                   as.data.table = T, columns = c("ZIP", "date", "pm25", "pm25_lag1"))
  
  # compute lag01 PM2.5
  pm25[ , pm25_lag01_case := (pm25 + pm25_lag1) / 2]
  pm25[ , c("pm25","pm25_lag1") := NULL]
  
  # merge lag01 PM2.5 for case and control days
  admissions <- merge(admissions, pm25, by.x = c("zip", "adate"), by.y = c("ZIP", "date"),
                      all.x = T, all.y = F)
  setnames(pm25, "pm25_lag01_case", "pm25_lag01_control")
  admissions <- merge(admissions, pm25, by.x = c("zip", "cdate"), by.y = c("ZIP", "date"),
                      all.x = T, all.y = F)

  # write to file
  write_fst(admissions, paste0("../data/merged_admissions_enviro/admissions_enviro_",year_,".fst")
  
}
