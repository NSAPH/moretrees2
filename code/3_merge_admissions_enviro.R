
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

admissions_all <- NULL

for (year_ in 2000:2014) {
  # read in admissions
  # NOTE: for now I am ignoring the fact that multiple hospitalizations may occur for same individual
  admissions <- read_fst(paste0("../data/admissions_cvd/admissions_cvd_", year_, ".fst"),
                        as.data.table = T, columns = c("zip", "adate", "ccs"))
  
  # select control days (randomly exactly one week before or after)
  # NOTE: could later swap this to "time-stratified" design
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
  
  # Remove PM2.5 data
  rm(pm25)
  
  # read in temperature & humidity data
  temp <- read_fst(paste0("../data/enviro/temp_", year_, ".fst"), as.data.table = T,
                          columns = c("ZIP", "date", "tmmx", "tmmx_lag1", "rmax", "rmax_lag1"))
  
  # compute lag01 temperature & humidity
  temp[ , tmmx_lag01_case := (tmmx + tmmx_lag1) / 2]
  temp[ , rmax_lag01_case := (rmax + rmax_lag1) / 2]
  temp[ , c("tmmx", "tmmx_lag1", "rmax", "rmax_lag1") := NULL]
  
  # merge temp & humidity data for case and control days
  admissions <- merge(admissions, temp, by.x = c("zip", "adate"), by.y = c("ZIP", "date"),
                      all.x = T, all.y = F)
  setnames(temp, c("tmmx_lag01_case", "rmax_lag01_case"),
                 c("tmmx_lag01_control", "rmax_lag01_control"))
  admissions <- merge(admissions, temp, by.x = c("zip", "cdate"), by.y = c("ZIP", "date"),
                      all.x = T, all.y = F)
  
  # Remove temp data
  rm(temp)
  
  # Drop unnecessary variables
  admissions[ , c("adate", "cdate", "zip") := NULL]
  
  # concatenate with other year's data
  admissions_all <- rbind(admissions_all, admissions)
  
  # Remove yearly dataset
  rm(admissions)
  
  # Keep track of progress
  print(year_)
}

# write to file
write_fst(admissions_all, "../data/merged_admissions_enviro/admissions_enviro.fst")
