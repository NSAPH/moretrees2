# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

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
dt <- read_fst("../data/merged_admissions_enviro/admissions_enviro_cvd.fst",
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
tr <- moretrees::ccs_tree("7")$tr # note: we have an error here.

# check all outcome codes are leaves of tree
sum(!(dt$ccs_added_zeros %in% names(V(tr))[V(tr)$leaf])) == 0

# Model 1: no covariate control ------------------------------------------------------------------------------------------
set.seed(34564)
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

# Get obs counts by group
obs_counts <- sapply(moretrees_results$beta_moretrees$outcomes,
                     function(out, dat) sum(dat %in% out),
                     dat = dt$ccs_added_zeros)
moretrees_results$beta_moretrees$n_obs <- obs_counts

# save
save(moretrees_results, file = "../results/mod1_split35_northEast.RData")

# Model 2: linear covariate control ------------------------------------------------------------------------------------
# get starting vals based on mod1
mod2_init <- mod1$mod
dsgn <- moretrees:::moretrees_design_tree(y = rep(1, n), X = as.matrix(dt[, c("pm25_blw35", "pm25_abv35")]),
                                          W = as.matrix(dt[ , c("tmmx", "rmax")]),
                                          outcomes = dt$ccs_added_zeros, tr = tr)
pL <- sum(V(tr)$leaf)
p <- length(V(tr))
m <- ncol(dsgn$W)
n <- nrow(dt)
mod2_init$hyperparams$g_eta <- c(mod2_init$hyperparams$g_eta, mean(mod2_init$hyperparams$g_eta))
mod2_init$hyperparams$omega <- mod2_init$hyperparams$tau
mod2_init$vi_params$delta <-lapply(1:p, function(i) matrix(0, nrow = m, ncol = 1))
xi <- mapply(`*`, mod2_init$vi_params$prob, mod2_init$vi_params$mu, SIMPLIFY = F)
Xbeta <- numeric(n) + 0
wwT <- moretrees:::rowOuterProds(dsgn$W)
for (u in 1:pL) {
  beta_u <- Reduce(`+`, xi[dsgn$ancestors[[u]]])
  Xbeta[dsgn$outcomes_units[[u]]] <- dsgn$X[dsgn$outcomes_units[[u]], 
                                  ] %*% beta_u
}
wwT_g_eta <- lapply(X = dsgn$outcomes_units, FUN = moretrees:::xxT_g_eta_fun, 
                    xxT = wwT, g_eta = mod2_init$hyperparams$g_eta, K = m)
for (v in 1:p) {
  leaf_descendants <- dsgn$outcomes_nodes[[v]]
  mod2_init$vi_params$Omega_inv[[v]] <- 2 * Reduce(`+`, wwT_g_eta[leaf_descendants]) + 
    diag(1/mod2_init$hyperparams$omega, nrow = m)
  mod2_init$vi_params$Omega[[v]] <- solve(mod2_init$vi_params$Omega_inv[[v]])
  mod2_init$vi_params$Omega_det[v] <- det(mod2_init$vi_params$Omega[[v]])
  mod2_init$vi_params$delta[[v]] <- mod2_init$vi_params$delta[[v]] * 0
  for (u in dsgn$leaf_descendants) {
    anc_u_mv <- setdiff(dsgn$ancestors[[u]], v)
    units_u <- dsgn$outcomes_units[[u]]
    theta_u_mv <- Reduce(`+`, mod2_init$vi_params$delta[anc_u_mv])
    mod2_init$vi_params$delta[[v]] <- mod2_init$vi_params$delta[[v]] + crossprod(W[units_u, 
                                           , drop = FALSE], (y[units_u]/2 - 2 * g_eta[units_u] * 
                                                               (dsgn$W[units_u, , drop = FALSE] %*% theta_u_mv + 
                                                                  Xbeta[units_u])))
  }
  mod2_init$vi_params$delta[[v]] <- mod2_init$vi_params$Omega[[v]] %*% mod2_init$vi_params$delta[[v]]
}

mod2_init$hyperparams <- moretrees:::update_hyperparams_logistic_moretrees(X = dsgn$X, W = dsgn$W, y = dsgn$y, outcomes_units = dsgn$outcomes_units, ancestors = dsgn$ancestors, n = n, K = 2, p = p, m = m, prob = mod2_init$vi_params$prob, mu = mod2_init$vi_params$mu, Sigma = mod2_init$vi_params$Sigma, Sigma_det = mod2_init$vi_params$Sigma_det, tau_t = mod2_init$vi_params$tau_t, delta = mod2_init$vi_params$delta, Omega = mod2_init$vi_params$Omega, Omega_det = mod2_init$vi_params$Omega_det, a_rho = mod2_init$vi_params$a_rho, b_rho = mod2_init$vi_params$b_rho, omega = mod2_init$hyperparams$omega, tau = mod2_init$hyperparams$tau, model = "moretrees", update_hyper = T)

require(gdata)
keep(dt, mod2_init, tr, sd_tmmx, sd_rmax)
set.seed(987234)
mod2 <- moretrees::moretrees(X = as.matrix(dt[, c("pm25_blw35", "pm25_abv35")]), 
                             W = as.matrix(dt[ , c("tmmx", "rmax")]),
                             y = rep(1, nrow(dt)), 
                             initial_values = mod2_init,
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

# Get obs counts by group
obs_counts <- sapply(moretrees_results$beta_moretrees$outcomes,
                     function(out, dat) sum(dat %in% out),
                     dat = dt$ccs_added_zeros)
moretrees_results$beta_moretrees$n_obs <- obs_counts

# save
save(moretrees_results, sd_tmmx, sd_rmax, file = "../results/mod2_split35_northEast_attempt2.RData")

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
set.seed(9845630)
W_cols <- c(paste0("tmmx_spl", 1:df), paste0("rmax_spl", 1:df))
mod3 <- moretrees::moretrees(X = as.matrix(dt[, c("pm25_blw35", "pm25_abv35")]), 
                             W = as.matrix(dt[ , W_cols, with = FALSE]),
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
moretrees_results <- mod3
moretrees_results$mod$hyperparams$g_eta <- NULL
moretrees_results$mod$hyperparams$eta <- NULL

# Get obs counts by group
obs_counts <- sapply(moretrees_results$beta_moretrees$outcomes,
                     function(out, dat) sum(dat %in% out),
                     dat = dt$ccs_added_zeros)
moretrees_results$beta_moretrees$n_obs <- obs_counts

# save
save(moretrees_results, sd_tmmx, sd_rmax, file = "../results/mod3_split35_northEast.RData")


