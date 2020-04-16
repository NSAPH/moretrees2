# Make sure working directory is set to moretrees2
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/")

source("./code/results_functions.R")
require(fst)
require(splines)
dataset <- c("cvd", "resp")
root <- c("7", "8")
mod <- 3
split <- "0"

nknots <- 2
q <- seq(0, 1, length.out = nknots + 2)

for (i in 1:length(dataset)) { # datasets
  
  # Load results
  ds <- dataset[i]
  load(file = paste0("./results/mod", mod, "_split", split, "_", ds, ".RData"))
  
  # Get tree
  tr <- ccs_tree(root[i])$tr
  vids <- unlist(moretrees_results$beta_moretrees$outcomes)
  vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
  vids <- Reduce(union, vids)
  tr <- induced_subgraph(tr, vids)
  moretrees_results$tr <- tr
  class(moretrees_results) <- "moretrees_result"
  
  # Load data
  dt <- read_fst(paste0("./data/merged_admissions_enviro/admissions_enviro_", ds, ".fst"),
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
  
  #   # Divide covariates by their standard deviation
  #   sd_tmmx <- sd(unlist(dt[ , c("tmmx_lag01_case", "tmmx_lag01_control")]), na.rm = T)
  #   dt[ , tmmx_lag01_case := tmmx_lag01_case / sd_tmmx]
  #   dt[ , tmmx_lag01_control := tmmx_lag01_control / sd_tmmx]
  #   sd_rmax <- sd(unlist(dt[ , c("rmax_lag01_case", "rmax_lag01_control")]), na.rm = T)
  #   dt[ , rmax_lag01_case := rmax_lag01_case / sd_rmax]
  #   dt[ , rmax_lag01_control := rmax_lag01_control / sd_rmax]
  
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
  
  # get moretrees individual estimates
  pltdat <- get_moretrees_indiv(moretrees_results, nsim = 100000)
  V(tr)$levels <- as.numeric(igraph::distances(tr, v = root[i],
                                               to = V(tr), mode = "out") + 1)
  d <- igraph::diameter(tr)
  descendants <- igraph::ego(tr, order = d + 1, nodes = pltdat$node, 
                             mode = "out")
  descendants <- sapply(descendants, names)
  leaves <- names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0])
  outcomes_nodes <- sapply(descendants, function(d, leaves) leaves[leaves %in% 
                                                                    d], leaves = leaves, simplify = F)
  pltdat$n <- sapply(outcomes_nodes, function(o) sum(dt$ccs_added_zeros %in% o))
  
  pltdat_ml <- moretrees:::ml_by_group(X = dt[, X_cols_case, with = F] - dt[, X_cols_control, with = F],
                                       W = dt[, W_cols_case, with = F] - dt[, W_cols_control, with = F],
                                       y = rep(1, nrow(dt)),
                                       outcomes = dt[, ccs_added_zeros],
                                       outcome_groups = outcomes_nodes,
                                       return_theta = F,
                                       return_ci = T,
                                       ci_level = 0.95,
                                       family = "binomial")
  
  # plot
  xlab <- "Excess Rate (%)"
  if (i == 1) {
    lab.widths <- c(0.9, 0.9, 0.9, 0.9)
    lab.txt.width <- c(20, 25, 25, 15)
    pdf(file = "./figures/cvd_nested_plot_results.pdf", height = 15, width = 15)
    beta_indiv_plot_fun(pltdat, tr, xlab = xlab, lab.widths = lab.widths,
                        lab.txt.width = lab.txt.width, axis.height = 1.5,
                        cil_min = -5, ciu_max = 5,
                        lab.txt.size = 4, digits = 0, axis.txt.size = 10,
                        plot_depth = 4, wrap_labs = FALSE)
    dev.off()
    pdf(file = "./figures/cvd_nested_plot_results_ml.pdf", height = 15, width = 15)
    beta_indiv_plot_fun(pltdat_ml, tr, xlab = xlab, lab.widths = lab.widths,
                        lab.txt.width = lab.txt.width, axis.height = 1.5,
                        cil_min = -5, ciu_max = 5,
                        lab.txt.size = 4, digits = 0, axis.txt.size = 10,
                        plot_depth = 4, wrap_labs = FALSE)
    dev.off()
  }
  if (i == 2) {
    lab.widths <- c(0.9, 0.9, 0.9, 0.9)
    lab.txt.width <- c(20, 25, 25, 15)
    pdf(file = "./figures/resp_nested_plot_results.pdf", height = 10, width = 15)
    beta_indiv_plot_fun(pltdat, tr, xlab = xlab, lab.widths = lab.widths,
                        lab.txt.width = lab.txt.width, axis.height = 1.2,
                        cil_min = -7, ciu_max = 7,
                        lab.txt.size = 4, digits = 0, axis.txt.size = 10,
                        plot_depth = 4, wrap_labs = FALSE,
                        force_lims = T)
    dev.off()
    pdf(file = "./figures/resp_nested_plot_results_ml.pdf", height = 15, width = 15)
    beta_indiv_plot_fun(pltdat_ml, tr, xlab = xlab, lab.widths = lab.widths,
                        lab.txt.width = lab.txt.width, axis.height = 1.5,
                        cil_min = -5, ciu_max = 5,
                        lab.txt.size = 4, digits = 0, axis.txt.size = 10,
                        plot_depth = 4, wrap_labs = FALSE)
    dev.off()
  }
}