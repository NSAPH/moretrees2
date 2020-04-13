# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/")

# Check for updates on moretrees master branch
# devtools::install_github("emgthomas/moretrees_pkg")
# note: for some updates, may have to restart R session
require(fst)
require(data.table)

# Key parameters
dataset <- "cvd" # "cvd" or "resp"
split <- "25" # "0" or "25"
mod <- 3
nfolds <- 10
set.seed(292345)
seed <- sample(1:1E6, 2)
if (dataset == "cvd") {
  seed <- seed[1]
} else {
  seed <- seed[2]
}

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

if (mod == 1) {
  W_cols_case <- NULL
  W_cols_control <- NULL
}

if (mod == 2) {
  W_cols_case <- c("tmmx_lag01_case", "rmax_lag01_case")
  W_cols_control <- c("tmmx_lag01_control", "rmax_lag01_control")
}

if (mod == 3) {
  # spline parameters
  require(splines)
  nknots <- 2
  q <- seq(0, 1, length.out = nknots + 2)
  
  # Get splines for temperature
  tmmx_knots <- quantile(unlist(dt[ , c("tmmx_lag01_case", "tmmx_lag01_control")]), probs = q, na.rm = T)
  tmmx_internal_knots <- tmmx_knots[2:(length(tmmx_knots) - 1)]
  tmmx_boundary_knots <- c(tmmx_knots[1], tmmx_knots[length(tmmx_knots)])
  dt[ , paste0("tmmx_spl_case", 1:(nknots + 1)) := as.data.table(ns(tmmx_lag01_case, knots = tmmx_internal_knots, Boundary.knots = tmmx_boundary_knots))]
  dt[ , paste0("tmmx_spl_control", 1:(nknots + 1)) := as.data.table(ns(tmmx_lag01_control, knots = tmmx_internal_knots, Boundary.knots = tmmx_boundary_knots))]
  
  # Get splines for humidity
  rmax_knots <- quantile(unlist(dt[ , c("rmax_lag01_case", "rmax_lag01_control")]), probs = q, na.rm = T)
  rmax_internal_knots <- rmax_knots[2:(length(rmax_knots) - 1)]
  rmax_boundary_knots <- c(rmax_knots[1], rmax_knots[length(rmax_knots)])
  dt[ , paste0("rmax_spl_case", 1:(nknots + 1)) := as.data.table(ns(rmax_lag01_case, knots = rmax_internal_knots, Boundary.knots = rmax_boundary_knots))]
  dt[ , paste0("rmax_spl_control", 1:(nknots + 1)) := as.data.table(ns(rmax_lag01_control, knots = rmax_internal_knots, Boundary.knots = rmax_boundary_knots))]
  
  # Which columns indicated covariate splines
  W_cols_case <- c(paste0("tmmx_spl_case", 1:(nknots + 1)), paste0("rmax_spl_case", 1:(nknots + 1)))
  W_cols_control <- c(paste0("tmmx_spl_control", 1:(nknots + 1)), paste0("rmax_spl_control", 1:(nknots + 1)))
}

# Keep only necessary variables
dt <- dt[ , c(X_cols_case, X_cols_control,
              W_cols_case, W_cols_control,
              "ccs_added_zeros"), with = FALSE]

# Remove NA rows (moretrees doesn't do this automatically)
dt <- na.omit(dt)

# Get 
set.seed(seed)
folds_function <- function(nfolds, n) sample(rep(sample(1:nfolds), length.out = n))
dt[ , folds := folds_function(nfolds, .N), by = ccs_added_zeros]
# Check this worked
# dcast(dt, folds ~ ccs_added_zeros)

# set keys
setkey(dt, "folds")

# Load initial values
load(paste0("./results/mod", mod, "_split", split, "_", dataset, ".RData"))
outcome_groups <- moretrees_results$beta_moretrees$outcomes
rm(moretrees_results)

# Set up parallelization
require(doParallel)
registerDoParallel(cores = nfolds)

# Out-of-sample prediction via 10-fold CV ------------------------------------------------------------------

# Function for computing log-likelihood component for each outcome
ll_fun <- function(v, beta, theta = NULL, Xdiff, Wdiff = NULL, outcomes, outcomes_unique){
  out <- outcomes_unique[[v]]
  if (!is.null(theta)) {
    ll <- as.numeric(moretrees:::logexpit(Xdiff[outcomes %in% out, , drop = F] %*% beta[v, ] + 
                                            Wdiff[outcomes %in% out, , drop = F] %*% theta[v, ]))
  } else {
    ll <- as.numeric(moretrees:::logexpit(Xdiff[outcomes %in% out, , drop = F] %*% beta[v, ]))
  }
  return(ll)
}

ll.cv <- as.data.frame(matrix(nrow = nfolds, ncol = length(outcome_groups) * 2 + 1))
names(ll.cv) <- c("fold", paste0("ll1.group", 1:length(outcome_groups)), paste0("ll2.group", 1:length(outcome_groups)))

for (i in 1:nfolds) {
  
  cat("\n\nFold", i, "\n\n")
  
  ll.ml <- numeric(2)
  n.train <- sum(dt$folds != i)
  # Grouped CLR estimates
  ml_group1 <- moretrees:::ml_by_group(X = dt[folds != i, X_cols_case, with = F] - dt[folds != i, X_cols_control, with = F],
                                       W = dt[folds != i, W_cols_case, with = F] - dt[folds != i, W_cols_control, with = F],
                                       y = rep(1, n.train),
                                       outcomes = dt[folds != i , ccs_added_zeros],
                                       outcome_groups = outcome_groups,
                                       return_theta = T,
                                       return_ci = F,
                                       family = "binomial")
  # Get test set log likelihoods for CLR
  ll.ml1 <- sapply(sapply(X = 1:length(outcome_groups), 
                         FUN = ll_fun, 
                         beta = as.matrix(ml_group1$beta_ml[ , paste0("est", 1:length(X_cols_case)), drop = F]),
                         theta = as.matrix(ml_group1$theta_ml[ , paste0("est", 1:length(W_cols_case)), drop = F]),
                         Xdiff = as.matrix(dt[folds == i, X_cols_case, with = F] - dt[folds == i, X_cols_control, with = F]),
                         Wdiff = as.matrix(dt[folds == i, W_cols_case, with = F] - dt[folds == i, W_cols_control, with = F]),
                         outcomes = dt[folds == i , ccs_added_zeros],
                         outcomes_unique = outcome_groups),
                  mean)
  
  # Grouped CLR estimates
  ml_group2 <- moretrees:::ml_by_group(X = dt[folds != i, W_cols_case, with = F] - dt[folds != i, W_cols_control, with = F],
                                       W = NULL,
                                       y = rep(1, n.train),
                                       outcomes = dt[folds != i , ccs_added_zeros],
                                       outcome_groups = outcome_groups,
                                       return_theta = T,
                                       return_ci = F,
                                       family = "binomial")
  # Get test set log likelihoods for CLR
  ll.ml2 <- sapply(sapply(X = 1:length(outcome_groups), 
                                 FUN = ll_fun, 
                                 beta = as.matrix(ml_group2$beta_ml[ , paste0("est", 1:length(W_cols_case)), drop = F]),
                                 Xdiff = as.matrix(dt[folds == i, W_cols_case, with = F] - dt[folds == i, W_cols_control, with = F]),
                                 outcomes = dt[folds == i , ccs_added_zeros],
                                 outcomes_unique = outcome_groups), mean)
    
  # Result
  ll.cv[i, ] <- c(i, ll.ml1, ll.ml2)
}

############### Save results ###############

# Model 1 includes PM2.5 + covariates; Model 2 includes only covariates
save(ll.cv, file = paste0("./results/cv2_mod", mod, "_split", split, "_", dataset, ".RData"))
