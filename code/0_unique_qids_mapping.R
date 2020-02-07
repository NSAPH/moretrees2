# Creating mapping from unique QIDs to integer ID for memory purposes

# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

library(NSAPHutils)
set_threads()

library(data.table)
library(fst)

# admissions path
admissions <- "../data/admissions"

qids <- character(0)

for (year_ in 2000:2015) { 
  # need to go up to 2015 as some admissions
  
  # Read in data
  qids_new <- read_data(admissions, years = year_, columns = "QID")$QID
  
  # Add unique new QIDs to list
  qids <- union(qids, qids_new)
  
  print(year_)
}

# Save as a data.table
qids <- data.table(qid = qids, id = 1:length(qids), key = "qid")
write_fst(qids, path = "../data/unique_qids/qids.fst")
