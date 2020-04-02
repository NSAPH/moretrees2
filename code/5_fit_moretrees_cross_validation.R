# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/")

# Check for updates on moretrees master branch
# devtools::install_github("emgthomas/moretrees_pkg", ref = "devel",)
require(moretrees)
# note: for some updates, may have to restart R session
require(fst)

# Key parameters
dataset <- "cvd" # "cvd" or "resp"
split <- "25" # "0" or "25"
nfolds <- 10

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

# Keep only necessary variables
dt <- dt[ , c(X_cols_case, X_cols_control,
              W_cols_case, W_cols_control,
              "ccs_added_zeros"), with = FALSE]

# Remove NA rows (moretrees doesn't do this automatically)
dt <- na.omit(dt)

# Get folds
dt$folds <- sample(1:nfolds, size = nrow(dt), replace = T)

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

# Create levels for CLR estimates
root <- names(V(tr))[degree(tr, mode = "in") == 0]
levels <- as.numeric(distances(tr, v = root, to = V(tr), mode = "out") + 1)
outcomes_levels <- as.data.frame(matrix(nrow = sum(V(tr)$leaf), ncol = max(levels)))
names(outcomes_levels) <- paste0('level', 1:4)
leaves <- names(V(tr))[V(tr)$leaf]
for(v in 1:length(leaves)) {
  outcomes_levels[v, ] <- names(rev(ego(tr, order = 3, nodes = leaves[v], mode = "in")[[1]]))
}
dt <- merge(dt, outcomes_levels, 
            by.x = "ccs_added_zeros", by.y = "level4",
            all.x = T, all.y = F)
names(dt)[names(dt) == "ccs_added_zeros"] <- "level4"

# set keys
setkeyv(dt, c("folds", paste0("level", 1:4)))

# Load initial values
load(paste0("./results/mod3_split", split, "_", dataset, ".RData"))
vi_params_init <- moretrees_results$mod$vi_params
hyperparams_init <- moretrees_results$mod$hyperparams
rm(moretrees_results)

# Set up parallelization
require(doParalllel)
registerDoParallel(cores = nfolds)

# Out-of-sample prediction via 10-fold CV ------------------------------------------------------------------

# Function for computing log-likelihood component for each outcome
ll_fun <- function(v, beta, theta, X, W, outcomes, outcomes_unique){
  out <- outcomes_unique[v]
  sum(moretrees:::logexpit(beta[[v]] %*% X[outcomes == out, , drop = F] + 
                          theta[[v]] %*% W[outcomes == out, , drop = F]))
}

cv_out <- foreach(i = 1:nfolds, .combine = cbind) %doPar% {
  
  # Grouped CLR estimates
  ml_group <- list()
  for (l in 1:max(levels)) {
    out <- as.list(unique(outcomes_levels[ , paste0("level", l)]))
    ml_group[[l]] <- moretrees:::ml_by_group(X = dt[folds != i, X.cols.case, with = F] - dt[folds != i, X.cols.control, with = F],
                                        W = dt[folds != i, W.cols.case, with = F] - dt[folds != i, W.cols.control, with = F],
                                        y = rep(1, nrow(X.train)),
                                        outcomes = dt[folds != i , paste0("level", l), with = FALSE],
                                        outcome_groups = out,
                                        return_theta = T)
  }
  
  # Run moretreees on training data
  mod <- moretrees::moretrees(Xcase = dt[folds != i, X.cols.case, with = F], 
                              Xcontrol = dt[folds != i, X.cols.control, with = F], 
                              Wcase = dt[folds != i, W.cols.case, with = F],
                              Wcontrol = dt[folds != i, W.cols.control, with = F],
                              outcomes = dt[folds !=i, level4],
                              hyper_fixed = hyper_fixed,
                              max_iter = max_iter,
                              update_hyper_freq = update_hyper_freq,
                              tol = tol_hyper,
                              tol_hyper = tol_hyper,
                              tr = tr,
                              nrestarts = 1,
                              print_freq = 1,  
                              get_ml = TRUE)
  
  # Get test set log likelihoods
  ll.moretrees <- mean(sapply(1:sum(V(tree)$leaf), 
                              ll_fun, 
                              beta = mod$beta_est,
                              theta = mod$theta_est,
                              X = dt[folds == i, X.cols.case, with = F] - dt[folds == i, X.cols.control, with = F],
                              W = dt[folds == i, W.cols.case, with = F] - dt[folds == i, W.cols.control, with = F],
                              outcomes = dt[folds == i, level4]))
  ll.ml <- numeric(ncol(outcomes_levels))
  for(l in 1:length(ll.ml)){
    outcomes_unique <- unique(outcomes_levels[ , paste0("level", l)])
    ll.ml[l] <- mean(sapply(1:length(outcomes_unique), 
                            ll_fun, 
                            beta = ml_group[[l]]$beta_ml[ , paste0("est", 1:length(X.cols.case))],
                            theta = ml_group[[l]]$theta_ml[ , paste0("est", 1:length(W.cols.case))],
                            X = dt[folds == i, X.cols.case, with = F] - dt[folds == i, X.cols.control, with = F],
                            W = dt[folds == i, W.cols.case, with = F] - dt[folds == i, W.cols.control, with = F],
                            outcomes = dt[folds == i, paste0("level", i), with = FALSE],
                            outcomes_unique = outcomes_unique))
  }
  
  # Result
  ll.cv <- as.data.frame(matrix(c(i, ll.moretrees, ll.ml),nrow=1))
  names(ll.cv) <- c("fold", "ll.moretrees", paste0("ll.ml", 1:length(ll.ml)))
  return(ll.cv)
}

############### Save results ###############

save(cv_out, file = paste0("./results/cv_split", split, "_", dataset, ".RData"))
