# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

# Check for updates on moretrees master branch
# devtools::install_github("emgthomas/moretrees_pkg", ref = "devel")
require(moretrees)
# note: for some updates, may have to restart R session
require(fst)

# Load zipcode/fips info
# Got this data from https://www.unitedstateszipcodes.org/zip-code-database/
zips <- read.csv("./data/zip_code_database.csv")
zips <- data.table::data.table(zips, key = "zip")
zips <- zips[, c("zip", "state")]
states_list <- c("CT", "ME", "MA", "NH", "RI", 
                 "VT", "NJ", "NY", "PA")

# Load data
dt <- read_fst("../data/merged_admissions_enviro/admissions_enviro.fst",
               as.data.table = T, 
               columns = c("id", "adate", "zip",
                           "ccs_added_zeros", "pm25_lag01_case", "pm25_lag01_control",
                           "tmmx_lag01_case", "tmmx_lag01_control",
                           "rmax_lag01_case", "rmax_lag01_control"))

# Keep only north east region
dt <- merge(dt, zips, by = "zip", all.x = T, all.y = F)
dt <- dt[state %in% states_list]

# First admission only
dt <- dt[order(id, adate)]
dt <- dt[ , .SD[1], by = id]

# Get difference between case and control
dt[ , pm25 := pm25_lag01_case - pm25_lag01_control]
dt[ , tmmx := tmmx_lag01_case - tmmx_lag01_control]
dt[ , rmax := rmax_lag01_case - rmax_lag01_control]

# Divide vars by their standard deviation
dt[ , pm25 := pm25 / sd(pm25, na.rm = T)]
dt[ , tmmx := tmmx / sd(tmmx, na.rm = T)]
dt[ , rmax := rmax / sd(rmax, na.rm = T)]

# Remove unnecessary columns
dt[ , c("id", "adate", "pm25_lag01_case", "pm25_lag01_control",
        "zip", "state",
        "tmmx_lag01_case", "tmmx_lag01_control",
        "rmax_lag01_case", "rmax_lag01_control") := NULL]

# Remove NA rows (moretrees doesn't do this automatically)
dt <- na.omit(dt)

# Get tree
require(magrittr)
require(igraph)
tr <- moretrees::ccs_tree("7")$tr # note: we have an error here.

# check all outcome codes are leaves of tree
sum(!(dt$ccs_added_zeros %in% names(V(tr))[V(tr)$leaf])) == 0

# Take a subsample (stratified on outcomes)
set.seed(24568)
# require(splitstackshape)
# dt <- stratified(indt = dt, group = "ccs_added_zeros", size = 50)

# Run MOReTreeS
# Started around 4pm Thursday
mod1 <- moretrees::moretrees(X = as.matrix(dt$pm25, ncol = 1), 
                             W = as.matrix(dt[ , c("tmmx", "rmax")]),
                             y = rep(1, nrow(dt)),
                             outcomes = dt$ccs_added_zeros,
                             max_iter = 1E4,
                             update_hyper_freq = 20,
                             tr = tr, 
                             method = "tree",
                             nrestarts = 1,
                             W_method = "shared",
                             print_freq = 1,  
                             get_ml = FALSE)

# Delete g
moretrees_results <- mod1
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL
save(moretrees_results, file = "../results/attempt1.RData")

set.seed(9384765)

# keep running
mod2 <- moretrees::moretrees(X = as.matrix(dt$pm25, ncol = 1), 
                             W = as.matrix(dt[ , c("tmmx", "rmax")]),
                             y = rep(1, nrow(dt)),
                             outcomes = dt$ccs_added_zeros,
                             initial_values = mod1$mod,
                             max_iter = 1E4,
                             update_hyper_freq = 20,
                             tr = tr, 
                             method = "tree",
                             nrestarts = 1,
                             W_method = "shared",
                             print_freq = 1,  
                             get_ml = FALSE)

# Delete g
moretrees_results <- mod2
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL
save(moretrees_results, file = "../results/attempt1_run2.RData")

# keep running
mod3 <- moretrees::moretrees(X = as.matrix(dt$pm25, ncol = 1), 
                             W = as.matrix(dt[ , c("tmmx", "rmax")]),
                             y = rep(1, nrow(dt)),
                             outcomes = dt$ccs_added_zeros,
                             initial_values = mod2$mod,
                             max_iter = 1E4,
                             update_hyper_freq = 20,
                             tr = tr, 
                             method = "tree",
                             nrestarts = 1,
                             W_method = "shared",
                             print_freq = 1,  
                             get_ml = FALSE)




