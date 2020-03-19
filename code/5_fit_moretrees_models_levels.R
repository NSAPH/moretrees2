# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/")

# Check for updates on moretrees master branch
# devtools::install_github("emgthomas/moretrees_pkg", ref = "devel",)
require(moretrees)
# note: for some updates, may have to restart R session
require(fst)

# Key parameters
dataset <- "resp" # "cvd" or "resp"
split <- "25" # "0", "25", or "35"
hyper_method <- "full" # "full" or "EB"

# Controls
tol <- 1E-8 
tol_hyper <- 1E-5
max_iter <- 1E5
update_hyper_freq <- 20

# Load data
dt <- read_fst(paste0("./data/merged_admissions_enviro/admissions_enviro_", dataset, ".fst"),
               as.data.table = T, 
               columns = c("id", "adate", 
                           # "ssa_state_cd",
                           "ccs_added_zeros", "pm25_lag01_case", "pm25_lag01_control",
                           "tmmx_lag01_case", "tmmx_lag01_control",
                           "rmax_lag01_case", "rmax_lag01_control"))

# # Keep only north east region
# states_list <- c(7, 20, 22, 30, 41, 47, 31, 33, 39)
# dt <- dt[ssa_state_cd %in% states_list]

# First admission only
dt <- dt[order(id, adate)]
dt <- dt[ , .SD[1], by = id]

if (split %in% c("25", "35")) {
  # Split pm2.5 into above and below EPA threshold
  split_val <- as.numeric(split)
  dt[ , pm25_blw_case := ifelse(pm25_lag01_case <= split_val, pm25_lag01_case, 0)]
  dt[ , pm25_abv_case := ifelse(pm25_lag01_case > split_val, pm25_lag01_case, 0)]
  dt[ , pm25_blw_control := ifelse(pm25_lag01_control <= split_val, pm25_lag01_control, 0)]
  dt[ , pm25_abv_control := ifelse(pm25_lag01_control > split_val, pm25_lag01_control, 0)]
  
  # Get difference between case and control
  dt[ , pm25_blw := pm25_blw_case - pm25_blw_control]
  dt[ , pm25_abv := pm25_abv_case - pm25_abv_control]
  
  # Which columns indicate exposure
  X_cols <- c("pm25_blw", "pm25_abv")
  
} else {
  # Get difference between case and control
  dt[ , pm25 := pm25_lag01_case - pm25_lag01_control]
  
  # Which columns indicate exposure
  X_cols <- "pm25"
}

# Get difference between case and control for covariates
dt[ , tmmx := tmmx_lag01_case - tmmx_lag01_control]
dt[ , rmax := rmax_lag01_case - rmax_lag01_control]

# Divide covariates by their standard deviation
sd_tmmx <- sd(dt$tmmx, na.rm = T)
sd_rmax <- sd(dt$rmax, na.rm = T)
dt[ , tmmx := tmmx / sd_tmmx]
dt[ , rmax := rmax / sd_rmax]

# Keep only necessary variables
dt <- dt[ , c(X_cols, "ccs_added_zeros",
              "tmmx", "rmax", "tmmx_lag01_case", "tmmx_lag01_control",
              "rmax_lag01_case", "rmax_lag01_control"), with = FALSE]

# Remove NA rows (moretrees doesn't do this automatically)
dt <- na.omit(dt)

# Get tree
require(magrittr)
require(igraph)
if (dataset == "cvd") root <- "7"
if (dataset == "resp") root <- "8"
tr <- moretrees::ccs_tree(root)$tr 

# check all outcome codes are leaves of tree
setequal(unique(dt$ccs_added_zeros), names(V(tr))[V(tr)$leaf])

# Some leaves are not outcomes, so take subtree
vids <- unique(dt$ccs_added_zeros)
vids <- ego(tr, order = 100, nodes = vids, mode = "in")
vids <- Reduce(union, vids)
tr <- induced_subgraph(graph = tr, vids = vids)
# check again
setequal(unique(dt$ccs_added_zeros), names(V(tr))[V(tr)$leaf])

# Create hyperparameter levels
V(tr)$levels <- 1
V(tr)$levels[V(tr)$leaf] <- 2

# Model 1: no covariate control ------------------------------------------------------------------------------------------
mod1 <- moretrees::moretrees(X = as.matrix(dt[ , X_cols, with = FALSE]), 
                             W = NULL,
                             y = rep(1, nrow(dt)),
                             hyper_method = hyper_method,
                             outcomes = dt$ccs_added_zeros,
                             max_iter = max_iter,
                             update_hyper_freq = update_hyper_freq,
                             tol = tol_hyper,
                             tol_hyper = tol_hyper,
                             tr = tr,
                             nrestarts = 1,
                             print_freq = 1,  
                             get_ml = TRUE)

# Delete g
moretrees_results <- mod1
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL

# save
save(moretrees_results, file = paste0("./results/mod1_split", split, "_", dataset, "_", hyper_method, "2.RData"))

# Model 2: linear covariate control ------------------------------------------------------------------------------------
mod2 <- moretrees::moretrees(X = as.matrix(dt[ , X_cols, with = FALSE]), 
                             W = as.matrix(dt[ , c("tmmx", "rmax")]),
                             y = rep(1, nrow(dt)), 
                             hyper_method = hyper_method,
                             outcomes = dt$ccs_added_zeros,
                             max_iter = max_iter,
                             update_hyper_freq = update_hyper_freq,
                             tol = tol_hyper,
                             tol_hyper = tol_hyper,
                             tr = tr,
                             nrestarts = 1,
                             print_freq = 1,  
                             get_ml = TRUE)

# Delete g
moretrees_results <- mod2
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL

# save
save(moretrees_results, file = paste0("./results/mod2_split", split, "_", dataset, "_", hyper_method, "2.RData"))

# Model 3: spline covariate control ------------------------------------------------------------------------------------

# spline parameters
require(splines)
require(data.table)
nknots <- 2
q <- seq(0, 1, length.out = nknots + 2)

# Get splines for temperature
sd_tmmx <- sd(unlist(dt[ , c("tmmx_lag01_case", "tmmx_lag01_control")]))
dt[ , tmmx_lag01_case := tmmx_lag01_case / sd_tmmx]
dt[ , tmmx_lag01_control := tmmx_lag01_control / sd_tmmx]
tmmx_knots <- quantile(unlist(dt[ , c("tmmx_lag01_case", "tmmx_lag01_control")]), probs = q)
tmmx_internal_knots <- tmmx_knots[2:(length(tmmx_knots) - 1)]
tmmx_boundary_knots <- c(tmmx_knots[1], tmmx_knots[length(tmmx_knots)])
dt[ , paste0("tmmx_spl", 1:(nknots + 1)) := as.data.table(ns(tmmx_lag01_case, knots = tmmx_internal_knots, Boundary.knots = tmmx_boundary_knots)
                                                          - ns(tmmx_lag01_control, knots = tmmx_internal_knots, Boundary.knots = tmmx_boundary_knots))]

# Get splines for humidity
sd_rmax <- sd(unlist(dt[ , c("rmax_lag01_case", "rmax_lag01_control")]))
dt[ , rmax_lag01_case := rmax_lag01_case / sd_rmax]
dt[ , rmax_lag01_control := rmax_lag01_control / sd_rmax]
rmax_knots <- quantile(unlist(dt[ , c("rmax_lag01_case", "rmax_lag01_control")]), probs = q)
rmax_internal_knots <- rmax_knots[2:(length(rmax_knots) - 1)]
rmax_boundary_knots <- c(rmax_knots[1], rmax_knots[length(rmax_knots)])
dt[ , paste0("rmax_spl", 1:(nknots + 1)) := as.data.table(ns(rmax_lag01_case, knots = rmax_internal_knots, Boundary.knots = rmax_boundary_knots)
                                                          - ns(rmax_lag01_control, knots = rmax_internal_knots, Boundary.knots = rmax_boundary_knots))]

# Remove unnecessary columns to save memory
dt[ , c("tmmx_lag01_case", "tmmx_lag01_control", "rmax_lag01_case", "rmax_lag01_control") := NULL]

# Run model
W_cols <- c(paste0("tmmx_spl", 1:(nknots + 1)), paste0("rmax_spl", 1:(nknots + 1)))
mod3 <- moretrees::moretrees(X = as.matrix(dt[ , X_cols, with = FALSE]), 
                             W = as.matrix(dt[ , W_cols, with = FALSE]),
                             y = rep(1, nrow(dt)),
                             outcomes = dt$ccs_added_zeros,
                             max_iter = max_iter,
                             hyper_method = hyper_method,
                             update_hyper_freq = update_hyper_freq,
                             tol = tol_hyper,
                             tol_hyper = tol_hyper,
                             tr = tr,
                             nrestarts = 1,
                             print_freq = 1,  
                             get_ml = TRUE)

# Delete g
moretrees_results <- mod3
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL

# save
save(moretrees_results, file = paste0("./results/mod3_split", split, "_", dataset, "_", hyper_method, "2.RData"))


