# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/")

## This code runs the main analysis + sensitivity analysis 1

# Check for updates on moretrees master branch
# devtools::install_github("emgthomas/moretrees_pkg", ref = "simplified",)
require(moretrees)
# note: for some updates, may have to restart R session
require(fst)

# Key parameters
dataset <- "cvd" # "cvd" or "resp"
split <- "25" # "0" (no knot) or "25" (piecewise linear pm2.5 effect with knot at pm2.5 = 25)

# Prior parameters
a <- c(1, 1)
b <- c(1, 1)
hyper_fixed = list(a = a, b = b)

# Controls
tol <- 1E-8 
tol_hyper <- 1E-4
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

# First admission only
dt <- dt[order(id, adate)]
dt <- dt[ , .SD[1], by = id]

if (split == "25") {
  # Split pm2.5 into above and below EPA threshold
  split_val <- as.numeric(split)
  dt[ , pm25_blw_case := ifelse(pm25_lag01_case <= split_val, pm25_lag01_case, 0)]
  dt[ , pm25_abv_case := ifelse(pm25_lag01_case > split_val, pm25_lag01_case, 0)]
  dt[ , pm25_blw_control := ifelse(pm25_lag01_control <= split_val, pm25_lag01_control, 0)]
  dt[ , pm25_abv_control := ifelse(pm25_lag01_control > split_val, pm25_lag01_control, 0)]
  
  # Which columns indicate exposure
  X_cols_case <- c("pm25_blw_case", "pm25_abv_case")
  X_cols_control <- c("pm25_blw_control", "pm25_abv_control")
  
} else {
  # Which columns indicate exposure
  X_cols_case <- "pm25_lag01_case"
  X_cols_control <- "pm25_lag01_control"
}

# Divide covariates by their standard deviation
sd_tmmx <- sd(unlist(dt[ , c("tmmx_lag01_case", "tmmx_lag01_control")]), na.rm = T)
dt[ , tmmx_lag01_case := tmmx_lag01_case / sd_tmmx]
dt[ , tmmx_lag01_control := tmmx_lag01_control / sd_tmmx]
sd_rmax <- sd(unlist(dt[ , c("rmax_lag01_case", "rmax_lag01_control")]), na.rm = T)
dt[ , rmax_lag01_case := rmax_lag01_case / sd_rmax]
dt[ , rmax_lag01_control := rmax_lag01_control / sd_rmax]

# Which columns indicate covariates
W_cols_case <- c("tmmx_lag01_case", "rmax_lag01_case")
W_cols_control <- c("tmmx_lag01_control", "rmax_lag01_control")

# Keep only necessary variables
dt <- dt[ , c(X_cols_case, X_cols_control,
              W_cols_case, W_cols_control,
              "ccs_added_zeros"), with = FALSE]

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
mod1 <- moretrees::moretrees(Xcase = as.matrix(dt[ , X_cols_case, with = FALSE]), 
                             Xcontrol = as.matrix(dt[ , X_cols_control, with = FALSE]), 
                             outcomes = dt$ccs_added_zeros,
                             hyper_fixed = hyper_fixed,
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
save(moretrees_results, file = paste0("./results/mod1_split", split, "_", dataset, ".RData"))

# Model 2: linear covariate control ------------------------------------------------------------------------------------
vi_params_init <- mod1$mod$vi_params[c("prob", "mu", "Sigma",
                                       "Sigma_inv", "Sigma_det",
                                       "tau_t", "a_t", "b_t")]
hyperparams_init <- mod1$mod$hyperparams[c("eta", "tau")]
mod2 <- moretrees::moretrees(Xcase = as.matrix(dt[ , X_cols_case, with = FALSE]), 
                             Xcontrol = as.matrix(dt[ , X_cols_control, with = FALSE]), 
                             Wcase = as.matrix(dt[ , W_cols_case, with = FALSE]),
                             Wcontrol = as.matrix(dt[ , W_cols_control, with = FALSE]),
                             outcomes = dt$ccs_added_zeros,
                             hyper_fixed = hyper_fixed,
                             vi_params_init = vi_params_init,
                             hyperparams_init = hyperparams_init,
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
save(moretrees_results, sd_tmmx, sd_rmax, file = paste0("./results/mod2_split", split, "_", dataset, ".RData"))

# Model 3: spline covariate control ------------------------------------------------------------------------------------

# spline parameters
require(splines)
require(data.table)
nknots <- 2
q <- seq(0, 1, length.out = nknots + 2)

# Get splines for temperature
tmmx_knots <- quantile(unlist(dt[ , c("tmmx_lag01_case", "tmmx_lag01_control")]), probs = q)
tmmx_internal_knots <- tmmx_knots[2:(length(tmmx_knots) - 1)]
tmmx_boundary_knots <- c(tmmx_knots[1], tmmx_knots[length(tmmx_knots)])
dt[ , paste0("tmmx_spl_case", 1:(nknots + 1)) := as.data.table(ns(tmmx_lag01_case, knots = tmmx_internal_knots, Boundary.knots = tmmx_boundary_knots))]
dt[ , paste0("tmmx_spl_control", 1:(nknots + 1)) := as.data.table(ns(tmmx_lag01_control, knots = tmmx_internal_knots, Boundary.knots = tmmx_boundary_knots))]

# Get splines for humidity
rmax_knots <- quantile(unlist(dt[ , c("rmax_lag01_case", "rmax_lag01_control")]), probs = q)
rmax_internal_knots <- rmax_knots[2:(length(rmax_knots) - 1)]
rmax_boundary_knots <- c(rmax_knots[1], rmax_knots[length(rmax_knots)])
dt[ , paste0("rmax_spl_case", 1:(nknots + 1)) := as.data.table(ns(rmax_lag01_case, knots = rmax_internal_knots, Boundary.knots = rmax_boundary_knots))]
dt[ , paste0("rmax_spl_control", 1:(nknots + 1)) := as.data.table(ns(rmax_lag01_control, knots = rmax_internal_knots, Boundary.knots = rmax_boundary_knots))]

# Remove unnecessary columns to save memory
dt[ , c("tmmx_lag01_case", "tmmx_lag01_control", "rmax_lag01_case", "rmax_lag01_control") := NULL]

# Which columns indicated covariate splines
W_cols_case <- c(paste0("tmmx_spl_case", 1:(nknots + 1)), paste0("rmax_spl_case", 1:(nknots + 1)))
W_cols_control <- c(paste0("tmmx_spl_control", 1:(nknots + 1)), paste0("rmax_spl_control", 1:(nknots + 1)))

# initial values
vi_params_init <- mod2$mod$vi_params[c("prob", "mu", "Sigma",
                                       "Sigma_inv", "Sigma_det",
                                       "tau_t", "a_t", "b_t")]
hyperparams_init <- mod2$mod$hyperparams[c("eta", "tau")]

# Run model
mod3 <- moretrees::moretrees(Xcase = as.matrix(dt[ , X_cols_case, with = FALSE]), 
                             Xcontrol = as.matrix(dt[ , X_cols_control, with = FALSE]), 
                             Wcase = as.matrix(dt[ , W_cols_case, with = FALSE]),
                             Wcontrol = as.matrix(dt[ , W_cols_control, with = FALSE]),
                             outcomes = dt$ccs_added_zeros,
                             hyper_fixed = hyper_fixed,
                             max_iter = max_iter,
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
save(moretrees_results, sd_tmmx, sd_rmax, file = paste0("./results/mod3_split", split, "_", dataset, ".RData"))


