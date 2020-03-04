# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/")

# Check for updates on moretrees master branch
# devtools::install_github("emgthomas/moretrees_pkg", ref = "devel")
require(moretrees)
# note: for some updates, may have to restart R session
require(fst)

# Select states
states_list <- c(7, 20, 22, 30, 41,
                 47, 31, 33, 39)
# states_list <- c("CT", "ME", "MA", "NH", "RI", 
#                 "VT", "NJ", "NY", "PA")

# Load data
dt <- read_fst("./data/merged_admissions_enviro/admissions_enviro_resp.fst",
               as.data.table = T, 
               columns = c("id", "adate", "ssa_state_cd",
                           "ccs_added_zeros", "pm25_lag01_case", "pm25_lag01_control",
                           "tmmx_lag01_case", "tmmx_lag01_control",
                           "rmax_lag01_case", "rmax_lag01_control"))

# Keep only north east region
dt <- dt[ssa_state_cd %in% states_list]

# First admission only
dt <- dt[order(id, adate)]
dt <- dt[ , .SD[1], by = id]

# Split pm2.5 into above and below EPA threshold
dt[ , pm25_blw35_case := ifelse(pm25_lag01_case <= 35, pm25_lag01_case, 0)]
dt[ , pm25_abv35_case := ifelse(pm25_lag01_case > 35, pm25_lag01_case, 0)]
dt[ , pm25_blw35_control := ifelse(pm25_lag01_control <= 35, pm25_lag01_control, 0)]
dt[ , pm25_abv35_control := ifelse(pm25_lag01_control > 35, pm25_lag01_control, 0)]

# Get difference between case and control
dt[ , pm25_blw35 := pm25_blw35_case - pm25_blw35_control]
dt[ , pm25_abv35 := pm25_abv35_case - pm25_abv35_control]
dt[ , tmmx := tmmx_lag01_case - tmmx_lag01_control]
dt[ , rmax := rmax_lag01_case - rmax_lag01_control]

# Divide covariates by their standard deviation
sd_tmmx <- sd(dt$tmmx, na.rm = T)
sd_rmax <- sd(dt$rmax, na.rm = T)
dt[ , tmmx := tmmx / sd_tmmx]
dt[ , rmax := rmax / sd_rmax]

# Remove unnecessary variables
dt[ , c("id", "adate", "ssa_state_cd", "pm25_lag01_case", "pm25_lag01_control",
        "pm25_blw35_case", "pm25_abv35_case", "pm25_blw35_control", "pm25_abv35_control") := NULL]

# Remove NA rows (moretrees doesn't do this automatically)
dt <- na.omit(dt)

# Get tree
require(magrittr)
require(igraph)
tr <- moretrees::ccs_tree("8")$tr 

# check all outcome codes are leaves of tree
setequal(unique(dt$ccs_added_zeros), names(V(tr))[V(tr)$leaf])

# Some leaves are not outcomes, so take subtree
vids <- unique(dt$ccs_added_zeros)
vids <- ego(tr, order = 100, nodes = vids, mode = "in")
vids <- Reduce(union, vids)
tr <- induced_subgraph(graph = tr, vids = vids)
# check again
setequal(unique(dt$ccs_added_zeros), names(V(tr))[V(tr)$leaf])

# Model 1: no covariate control ------------------------------------------------------------------------------------------
set.seed(6239502)
mod1 <- moretrees::moretrees(X = as.matrix(dt[, c("pm25_blw35", "pm25_abv35")]), 
                             W = NULL,
                             y = rep(1, nrow(dt)),
                             outcomes = dt$ccs_added_zeros,
                             max_iter = 1E5,
                             update_hyper_freq = 20,
                             tr = tr, 
                             method = "tree",
                             nrestarts = 1,
                             W_method = "shared",
                             print_freq = 1,  
                             get_ml = TRUE)

# Delete g
moretrees_results <- mod1
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL
save(moretrees_results, file = "./results/mod1_split35_northEast_resp.RData")

mod1_run2 <- moretrees::moretrees(X = as.matrix(dt[, c("pm25_blw35", "pm25_abv35")]), 
                             W = NULL,
                             y = rep(1, nrow(dt)),
                             outcomes = dt$ccs_added_zeros,
                             max_iter = 1E5,
                             update_hyper_freq = 20,
                             tr = tr, 
                             initial_values = mod1$mod,
                             method = "tree",
                             nrestarts = 1,
                             W_method = "shared",
                             print_freq = 1,  
                             get_ml = TRUE)

# Delete g
moretrees_results <- mod1_run2
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL

# Get obs counts by group
obs_counts <- sapply(moretrees_results$beta_moretrees$outcomes,
                     function(out, dat) sum(dat %in% out),
                     dat = dt$ccs_added_zeros)
moretrees_results$beta_moretrees$n_obs <- obs_counts

# save
save(moretrees_results, file = "./results/mod1_split35_northEast_resp_run2.RData")
# rm(mod1, mod1_run2)

# Model 2: linear covariate control ------------------------------------------------------------------------------------
mod2 <- moretrees::moretrees(X = as.matrix(dt[, c("pm25_blw35", "pm25_abv35")]), 
                             W = as.matrix(dt[ , c("tmmx", "rmax")]),
                             y = rep(1, nrow(dt)),
                             initial_values = mod1_run2$mod,
                             outcomes = dt$ccs_added_zeros,
                             max_iter = 1E5,
                             update_hyper_freq = 20,
                             tr = tr, 
                             method = "tree",
                             nrestarts = 1,
                             W_method = "shared",
                             print_freq = 1,  
                             get_ml = TRUE)

# Delete g
moretrees_results <- mod2
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL

# save
save(moretrees_results, sd_tmmx, sd_rmax, file = "./results/mod2_split35_northEast_resp.RData")

# Model 3: spline covariate control ------------------------------------------------------------------------------------

# Make splines
require(splines)
require(data.table)
df <- 3
cols <- c("tmmx_lag01_case", "tmmx_lag01_control", "rmax_lag01_case", "rmax_lag01_control")
for (col in cols) {
  dt[ , paste0(col, "_spl", 1:df) := data.table(ns(get(col), df = df))]
}
for (i in 1:df) {
  dt[ , paste0("tmmx_spl", i) := get(paste0("tmmx_lag01_case_spl", i)) - get(paste0("tmmx_lag01_control_spl", i))]
  dt[ , paste0("rmax_spl", i) := get(paste0("rmax_lag01_case_spl", i)) - get(paste0("rmax_lag01_control_spl", i))]
}

# Remove unnecessary columns to save memory
dt[ , c("tmmx_lag01_case", "tmmx_lag01_control", "rmax_lag01_case", "rmax_lag01_control",
        sapply(cols, function(col) paste0(col, "_spl", 1:df))) := NULL]


# Run model
W_cols <- c(paste0("tmmx_spl", 1:df), paste0("rmax_spl", 1:df))
mod3 <- moretrees::moretrees(X = as.matrix(dt[, c("pm25_blw35", "pm25_abv35")]), 
                             W = as.matrix(dt[ , W_cols, with = FALSE]),
                             y = rep(1, nrow(dt)),
                             initial_values = mod2$mod,
                             outcomes = dt$ccs_added_zeros,
                             max_iter = 3E5,
                             update_hyper_freq = 20,
                             tr = tr, 
                             method = "tree",
                             nrestarts = 1,
                             W_method = "shared",
                             print_freq = 1,  
                             get_ml = TRUE)

# Delete g
moretrees_results <- mod3
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL

# save
save(moretrees_results, sd_tmmx, sd_rmax, file = "./results/mod3_split35_northEast_resp.RData")


