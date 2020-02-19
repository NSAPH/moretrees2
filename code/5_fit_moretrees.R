# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

# Check for updates on moretrees master branch
# devtools::install_github("emgthomas/moretrees_pkg", ref = "devel")
require(moretrees)
# note: for some updates, may have to restart R session
require(fst)

# Load data
dt <- read_fst("../data/merged_admissions_enviro/admissions_enviro.fst",
               as.data.table = T, 
               columns = c("id", "adate", "ccs", "pm25_lag01_case", "pm25_lag01_control",
                           "tmmx_lag01_case", "tmmx_lag01_control",
                           "rmax_lag01_case", "rmax_lag01_control"))
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
        "tmmx_lag01_case", "tmmx_lag01_control",
        "rmax_lag01_case", "rmax_lag01_control") := NULL]

# Remove NA rows (moretrees doesn't do this automatically)
dt <- na.omit(dt)

# Get tree
require(magrittr)
require(igraph)
tr <- moretrees::ccs_tree("7")$tr # note: we have an error here.

dt <- dt[ccs %in% names(V(tr))[V(tr)$leaf]]

# Take a subsample (stratified on outcomes)
# require(splitstackshape)
set.seed(24568)
# dt <- stratified(indt = dt, group = "ccs", size = 50)

# Run MOReTreeS
mod1 <- moretrees::moretrees(X = as.matrix(dt$pm25, ncol = 1), 
                             W = as.matrix(dt[ , c("tmmx", "rmax")]),
                             y = rep(1, nrow(dt)),
                             outcomes = dt$ccs,
                             max_iter = 50,
                             tr = tr, 
                             method = "tree",
                             nrestarts = 1,
                             W_method = "shared",
                             print_freq = 1,  
                             get_ml = FALSE)





