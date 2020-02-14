require(readr)
require(fst)
require(magrittr)
require(data.table)
require(lubridate)

setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

# Relevant dates
control_dist <- 0 # max number of days between case day and control day
maxlag <- 3 # maximum lag to consider
buffer <- control_dist + maxlag # total buffer
year_range <- 2000:2014 # years to include

# ------------------------------------------------------ PM2.5 ------------------------------------------------

pm <- fread("../data/daily_pm/all_days_PM.csv", data.table = T,
            key = c("ZIP", "date"))
# re-format date variable
pm[ , date := ymd(date)]
pm <- pm[date >= make_date(year = 2000, month = 1, day = 1) - buffer & 
               date <= make_date(year = 2014, month = 12, day = 31) + buffer]

for (year_ in year_range) {
  # subset to particular year
  # (need to keep buffer of one week + 3 days on either end due to lag and control days)
  pm_year <- pm[date >= make_date(year = year_, month = 1, day = 1) - buffer &
                  date <= make_date(year = year_, month = 12, day = 31) + buffer]
  
  # sort by zip/date
  pm_year <- pm_year[order(ZIP, date)]
  # compute lags 1, 2 for PM2.5
  pm_year[ , c("pm25_lag1", "pm25_lag2") := shift(pm25, n = 1:2, type = "lag"),
                           by = ZIP]
  # exclude extra lag days
  pm_year <- pm_year[date >= make_date(year = year_, month = 1, day = 1) - control_dist &
                       date <= make_date(year = year_, month = 12, day = 31) + control_dist]
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
temp <- temp[date >= make_date(year = 2000, month = 1, day = 1) - buffer & 
               date <= make_date(year = 2014, month = 12, day = 31) + buffer]

for (year_ in year_range) {
  # subset to particular year
  # (need to keep buffer of one week + 3 days on either end due to lag and control days)
  temp_year <- temp[date >= make_date(year = year_, month = 1, day = 1) - buffer &
                      date <= make_date(year = year_, month = 12, day = 31) + buffer] 
  # sort by zip/date
  temp_year <- temp_year[order(ZIP, date)]
  # compute lags 1, 2, 3 for temp
  temp_year[ , c("tmmx_lag1", "tmmx_lag2", "tmmx_lag3") := shift(tmmx, n = 1:3, type = "lag"),
          by = ZIP]
  # compute lags 1, 2 for humidity
  temp_year[ , c("rmax_lag1", "rmax_lag2") := shift(rmax, n = 1:2, type = "lag"),
            by = ZIP]
  # exclude extra lag days
  temp_year <- temp_year[date >= make_date(year = year_, month = 1, day = 1) - control_dist &
                           date <= make_date(year = year_, month = 12, day = 31) + control_dist]
  # write to file
  write_fst(temp_year, path = paste0("../data/enviro/temp_", year_, ".fst"))
  # remove subsetted data.table for memory purposes
  rm(temp_year)
  print(year_)
}
# remove large temp dataset
rm(temp)

# ----------------------------------------------- Ozone ---------------------------------------------------

for (year_ in year_range) {
  # get relevant dates
  dates <- seq(make_date(year = year_, month = 1, day = 1) - buffer,
                 make_date(year = year_, month = 12, day = 31) + buffer,
               by = "days")
  dates <- format(dates, "%Y%m%d")
  ozone_year <- data.table(ZIP = integer(0),
                           ozone = numeric(0),
                           date = integer(0))
  # read in ozone data by days
  for (date_ in dates) {
    date_file <- paste0("../data/ozone/", date_, ".csv")
    if (file.exists(date_file)) {
      ozone_date <- fread(file = date_file, data.table = T)
      ozone_year <- rbind(ozone_year, ozone_date)
    }
  }
  
  # format date variable
  ozone_year[ , date := ymd(date)]
  
  # sort by zip/date
  setkey(ozone_year, ZIP, date)
  ozone_year <- ozone_year[order(ZIP, date)]
  # compute lags 1, 2 for PM2.5
  ozone_year[ , c("ozone_lag1", "ozone_lag2") := shift(ozone, n = 1:2, type = "lag"),
          by = ZIP]
  # exclude extra lag days
  ozone_year <- ozone_year[date >= make_date(year = year_, month = 1, day = 1) - control_dist &
                       date <= make_date(year = year_, month = 12, day = 31) + control_dist]
  # write to file
  write_fst(ozone_year, path = paste0("../data/enviro/ozone_", year_, ".fst"))
  # remove subsetted data.table for memory purposes
  rm(ozone_year)
  print(year_)
}


