
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
enviro_path <- "../data/enviro/"
merged_path <- "../data/merged_admissions_enviro/"

# function that gets potential control dates for all case dates in a range
get_control_dates <- function(begin, end) {
  dates <- seq(begin, end, by = "day")
  controls_dt <- data.table(dates = dates, control_days = rep(list(), length(dates)),
                            key = "dates")
  for (i in 1:nrow(controls_dt)) {
    date <- controls_dt$dates[i]
    start <- make_date(day = 1, month = month(date), year = year(date))
    control_days <- seq(start, start + days_in_month(start) - 1, by = "day") # match on month/year
    control_days <- control_days[wday(control_days) == wday(date)] # match on day of week
    control_days <- control_days[control_days != date] # remove date itself
    controls_dt$control_days[[i]] <- control_days
  }
  return(controls_dt)
}

admissions_all <- NULL

# set seed because we will randomly select control day
set.seed(6357312)

for (year_ in 2000:2014) {
  # read in admissions
  # NOTE: for now I am ignoring the fact that multiple hospitalizations may occur for same individual
  admissions <- read_fst(paste0("../data/admissions_cvd/admissions_cvd_", year_, ".fst"),
                        as.data.table = T, columns = c("id", "zip", "adate", "ccs"))
  
  # select control days (stratified on year, month, day of week)
  # get list of potential controls
  potential_controls <- get_control_dates(begin = make_date(day = 1, month = 1, year = year_),
                                          end = make_date(day = 31, month = 12, year = year_))
  admissions <- merge(admissions, potential_controls, by.x = "adate", by.y = "dates",
                      all.x = T, all.y = F)
  # randomly select one control per case
  admissions[ , cdate := lapply(control_days, sample, size = 1)]
  admissions[ , cdate := do.call("c", cdate)]
  admissions[ , control_days := NULL]
  
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
  
  # remove temp data
  rm(temp)

  # read in ozone data
  ozone <- read_fst(paste0("../data/enviro/ozone_", year_, ".fst"), as.data.table = T,
			columns = c("ZIP", "date", "ozone", "ozone_lag1"))

  # Compute lag01 ozone
  ozone[ , ozone_lag01_case := (ozone + ozone_lag1) / 2]
  ozone[ , c("ozone", "ozone_lag1") := NULL]

  # merge ozone data for case and control days
  admissions <- merge(admissions, ozone, by.x = c("zip", "adate"), by.y = c("ZIP", "date"),
		all.x = T, all.y = F)
  setnames(ozone, "ozone_lag01_case", "ozone_lag01_control")
  admissions <- merge(admissions, ozone, by.x = c("zip", "cdate"), by.y = c("ZIP", "date"),
		all.x = T, all.y = F)

  # Remove ozone data
  rm(ozone)
  
  # Drop unnecessary variables
  # admissions[ , c("adate", "cdate", "zip") := NULL]
  
  # concatenate with other year's data
  admissions_all <- rbind(admissions_all, admissions)
  
  # Remove yearly dataset
  rm(admissions)
  
  # Keep track of progress
  print(year_)
}

# write to file (all admissions)
write_fst(admissions_all, "../data/merged_admissions_enviro/admissions_enviro.fst")
