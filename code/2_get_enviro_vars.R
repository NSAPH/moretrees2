require(readr)
require(fst)
require(magrittr)
require(data.table)
require(lubridate)

setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

# ------------------------------------------------------ PM2.5 ------------------------------------------------

pm <- fread("../data/daily_pm/all_days_PM.csv", data.table = T,
            key = c("ZIP", "date"))
# re-format date variable
pm[ , date := ymd(date)]
# Keep only relevant years
pm <- pm[date >= make_date(year = 2000, month = 1, day = 1) & 
               date <= make_date(year = 2014, month = 12, day = 31)]

for (year_ in 2000:2014) {
  # subset to particular year
  pm_year <- pm[date <= make_date(year = year_, month = 12, day = 31) &
                  date >= make_date(year = year_, month = 1, day = 1) - 2] # need to include two days due to lag
  # sort by zip/date
  pm_year <- pm_year[order(ZIP, date)]
  # compute lags 1, 2 for PM2.5
  pm_year[ , c("pm25_lag1", "pm25_lag2") := shift(pm25, n = 1:2, type = "lag"),
                           by = ZIP]
  # exclude extra lag days
  pm_year <- pm_year[date >= make_date(year = year_, month = 1, day = 1)]
  # write to file
  write_fst(pm_year, path = paste0("../data/enviro/pm_", year_, ".fst"))
  # remove subsetted data.table for memory purposes
  rm(pm_year)
  print(year_)
}
# remove large PM dataset
rm(pm)

# ----------------------------------------------- Temperature ------------------------------------------------

# Temperature
temp <- fread("../data/temperature/temperature_daily_zipcode_combined.csv", data.table = T,
                 key = c("ZIP", "date"))
temp[ , c("pr", "year") := NULL]
temp[ , date := ymd(date)]
temp <- temp[date >= make_date(year = 2000, month = 1, day = 1) & 
               date <= make_date(year = 2014, month = 12, day = 31)]

for (year_ in 2000:2014) {
  # subset to particular year
  temp_year <- temp[date <= make_date(year = year_, month = 12, day = 31) &
                      date >= make_date(year = year_, month = 1, day = 1) - 3] # need to include three days due to lag
  # sort by zip/date
  temp_year <- temp_year[order(ZIP, date)]
  # compute lags 1, 2, 3 for temp
  temp_year[ , c("tmmx_lag1", "tmmx_lag2", "tmmx_lag3") := shift(tmmx, n = 1:3, type = "lag"),
          by = ZIP]
  # compute lags 1, 2 for humidity
  temp_year[ , c("rmax_lag1", "rmax_lag2") := shift(rmax, n = 1:2, type = "lag"),
            by = ZIP]
  # exclude extra lag days
  temp_year <- temp_year[date >= make_date(year = year_, month = 1, day = 1)]
  # write to file
  write_fst(temp_year, path = paste0("../data/enviro/temp_", year_, ".fst"))
  # remove subsetted data.table for memory purposes
  rm(temp_year)
  print(year_)
}
# remove large temp dataset
rm(temp)

# ----------------------------------------------- Ozone ---------------------------------------------------

# ozone <- read_rds("../data/ozone_data/PREDICTIONGeneral2_OZONE_zipcode_SUBSET_2000_2012.rds")


