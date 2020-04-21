# Make sure working directory is set to moretrees2
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/")

source("./code/results_functions.R")
require(fst)
require(splines)
dataset <- "resp"
root <- c("7", "8")
mod <- 3
split <- "0"
nknots <- 2
q <- seq(0, 1, length.out = nknots + 2)


# Load data
dt <- read_fst(paste0("./data/merged_admissions_enviro/admissions_enviro_", dataset, ".fst"),
               as.data.table = T,
               columns = c("id", "adate",
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

# Get tree
require(magrittr)
require(igraph)
if (dataset == "cvd") root <- "7"
if (dataset == "resp") root <- "8"
tr <- moretrees::ccs_tree(root)$tr
vids <- unique(dt$ccs_added_zeros)
vids <- ego(tr, order = 100, nodes = vids, mode = "in")
vids <- Reduce(union, vids)
tr <- induced_subgraph(graph = tr, vids = vids)
V(tr)$levels <- as.numeric(igraph::distances(tr, v = root,
                                             to = V(tr), mode = "out") + 1)

# get CLR estimates for every group on tree
d <- igraph::diameter(tr)
descendants <- igraph::ego(tr, order = d + 1, nodes = names(V(tr)),
                           mode = "out")
descendants <- sapply(descendants, names)
leaves <- names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0])
outcomes_nodes <- sapply(descendants, function(d, leaves) leaves[leaves %in%
                                                                   d], leaves = leaves, simplify = F)

clr_res <- moretrees:::ml_by_group(X = dt[, X_cols_case, with = F] - dt[, X_cols_control, with = F],
                                     W = dt[, W_cols_case, with = F] - dt[, W_cols_control, with = F],
                                     y = rep(1, nrow(dt)),
                                     outcomes = dt[, ccs_added_zeros],
                                     outcome_groups = outcomes_nodes,
                                     return_theta = T,
                                     return_ci = T,
                                     ci_level = 0.95,
                                     family = "binomial")
pltdat_ml <- clr_res$beta_ml
pltdat_ml$group <- pltdat_ml$outcomes <- NULL
pltdat_ml$node <- names(V(tr))
pltdat_ml$n <- sapply(outcomes_nodes, function(o) sum(dt$ccs_added_zeros %in% o))
pltdat_ml_theta <- clr_res$theta_ml

# Save
save(pltdat_ml, tr, file = paste0("./results/nest_pltdat_", dataset, ".RData"))
save(tmmx_internal_knots, tmmx_boundary_knots, rmax_internal_knots, rmax_boundary_knots, pltdat_ml_theta,
     file = paste0("./results/theta_pltdat_", dataset, ".RData"))

