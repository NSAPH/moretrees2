require(readr)
require(fst)
require(magrittr)
require(reshape2)
require(data.table)

setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

# function to reverse strings
zipReverse <- function(x) as.character(x) %>%
  strsplit(split = NULL) %>%
  lapply(FUN = rev) %>%
  sapply(FUN = paste, collapse = "") %>%
  as.numeric

# ------------------------------------------------------ PM2.5 ------------------------------------------------
pm <- fread("../data/daily_pm/all_days_PM.csv", data.table = T,
            key = c("ZIP", "date"))

# Temperature
# temp <- fread("../data/temperature/temperature_daily_zipcode_combined.csv",
#              nrow = 10)

for (year_ in 2000:2014) {
  # subset to particular year
  pm_year <- pm[date <= year_ * 10000 + 9999 & 
                      date >= year_ * 10000 - 1]
  # reverse zip codes (???)
  pm_year[ , ZIP := zipReverse(ZIP)]
  # re-format date variable
  pm_year$date <- as.character(pm_year$date) %>% 
    as.Date(format = "%Y%m%d")
  # sort by zip/date
  pm_year <- pm_year[order(ZIP, date)]
  # compute lag01 PM2.5
  pm_year[ , c("pm25_lag1", "pm25_lag2") := shift(pm25, n = 1:2, type = "lag"),
                           by = ZIP]
  # write to file
  write_fst(pm_year, path = paste0("../data/enviro/pm_", year_, ".fst"))
  # remove subsetted data.table for memory purposes
  rm(pm_year)
}

# ----------------------------------------------- Temperature ------------------------------------------------

# Temperature
# temp <- fread("../data/temperature/temperature_daily_zipcode_combined.csv",
#              nrow = 10)

for (year_ in 2000:2014) {
  # subset to particular year
  temp_year <- temp[date <= year_ * 10000 + 9999 & 
                  date >= year_ * 10000 - 1]
  # reverse zip codes (???)
  temp_year[ , ZIP := zipReverse(ZIP)]
  # re-format date variable
  temp_year$date <- as.character(pm_year$date) %>% 
    as.Date(format = "%Y%m%d")
  # sort by zip/date
  pm_year <- pm_year[order(ZIP, date)]
  # compute lag01 PM2.5
  pm_year[ , c("pm25_lag1", "pm25_lag2") := shift(pm25, n = 1:2, type = "lag"),
          by = ZIP]
  # write to file
  write_fst(pm_year, file = paste0("../data/enviro/pm_", year_, ".fst"))
  # remove subsetted data.table for memory purposes
  rm(pm_year)
}

# ozone <- read_rds("../data/ozone_data/PREDICTIONGeneral2_OZONE_zipcode_SUBSET_2000_2012.rds")
