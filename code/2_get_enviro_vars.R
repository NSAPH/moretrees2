require(readr)
require(fst)
require(magrittr)
require(data.table)

setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

# ------------------------------------------------------ PM2.5 ------------------------------------------------

pm <- fread("../data/daily_pm/all_days_PM.csv", data.table = T,
            key = c("ZIP", "date"))
# re-format date variable
pm[ , date := as.Date(date, format = "%Y-%m-%d")]
# Keep only relevant years
pm <- pm[date >= as.Date("2000-01-01") & 
               date <= as.Date("2014-12-31")]

for (year_ in 2000:2014) {
  # subset to particular year
  pm_year <- pm[date <= as.Date(paste0(year_, "-12-31")) &
                  date >= as.Date(paste0(year_, "-01-01")) - 2] # need to include two days due to lag
  # reverse zip codes (???)
  pm_year[ , ZIP := zipReverse(ZIP)]
  # sort by zip/date
  pm_year <- pm_year[order(ZIP, date)]
  # compute lags 1, 2 for PM2.5
  pm_year[ , c("pm25_lag1", "pm25_lag2") := shift(pm25, n = 1:2, type = "lag"),
                           by = ZIP]
  # exclude extra lag days
  pm_year <- pm_year[date >= as.Date(paste0(year_, "-01-01"))]
  # write to file
  write_fst(pm_year, path = paste0("../data/enviro/pm_", year_, ".fst"))
  # remove subsetted data.table for memory purposes
  rm(pm_year)
}
# remove large PM dataset
rm(pm)

# ----------------------------------------------- Temperature ------------------------------------------------

# Temperature
temp <- fread("../data/temperature/temperature_daily_zipcode_combined.csv", data.table = T,
                 key = c("ZIP", "date"))
temp[ , c("pr", "year") := NULL]
temp[ , date := as.Date(date, format = "%Y-%m-%d")]
temp <- temp[date >= as.Date("2000-01-01") & 
               date <= as.Date("2014-12-31")]

for (year_ in 2000:2014) {
  # subset to particular year
  temp_year <- temp[date <= as.Date(paste0(year_, "-12-31")) &
                      date >= as.Date(paste0(year_, "-01-01")) - 3] # need to include three days due to lag
  # sort by zip/date
  temp_year <- temp_year[order(ZIP, date)]
  # compute lags 1, 2, 3 for temp
  temp_year[ , c("tmmx_lag1", "tmmx_lag2", "tmmx_lag3") := shift(tmmx, n = 1:3, type = "lag"),
          by = ZIP]
  # compute lags 1, 2 for humidity
  temp_year[ , c("rmax_lag1", "rmax_lag2") := shift(rmax, n = 1:2, type = "lag"),
            by = ZIP]
  # write to file
  write_fst(temp_year, path = paste0("../data/enviro/temp_", year_, ".fst"))
  # remove subsetted data.table for memory purposes
  rm(temp_year)
}
# remove large temp dataset
rm(temp)

# ----------------------------------------------- Ozone ---------------------------------------------------

# ozone <- read_rds("../data/ozone_data/PREDICTIONGeneral2_OZONE_zipcode_SUBSET_2000_2012.rds")
