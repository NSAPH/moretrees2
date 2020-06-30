
source("./code/results_functions.R")
require(fst)
require(splines)
require(reshape2)
require(stringr)
dataset <- c("cvd", "resp")
root <- c("7", "8")
mod <- 3
split <- "0"

nknots <- 2
q <- seq(0, 1, length.out = nknots + 2)
nplt <- 1000

for (i in 1:length(dataset)) { # datasets
  ds <- dataset[i]

  load(paste0("./results/theta_pltdat_", ds, ".RData"), verbose = TRUE)
  load(paste0("./results/mod3_split0_", ds, ".Rdata"))
  theta_est <- moretrees_results$theta_est
  theta_est_clr <- pltdat_ml_theta
  theta_est_clr$group <- NULL
  theta_est_clr$outcomes <- sapply(theta_est_clr$outcomes, paste0, collapse = "")
  theta_est_clr <- theta_est_clr[theta_est_clr$outcomes %in% row.names(theta_est), ]
  theta_est_clr <- theta_est_clr[!duplicated(theta_est_clr), ]
  theta_est_clr <- theta_est_clr[match(row.names(theta_est), theta_est_clr$outcomes), ]
  all.equal(theta_est_clr$outcomes, row.names(theta_est))
  theta_est_clr <- theta_est_clr[ , paste0("est", 1:6)]
  
  # Temp
  tmmx <- seq(tmmx_boundary_knots[1], tmmx_boundary_knots[2], length.out = nplt)
  spl_tmmx <- ns(tmmx, knots = tmmx_internal_knots, Boundary.knots = tmmx_boundary_knots)
  theta_tmmx <- theta_est[ , paste0("est", 1:3)]
  tmmx_df <- as.data.frame(spl_tmmx %*% t(theta_tmmx))
  tmmx_clr <- as.data.frame(spl_tmmx %*% t(theta_est_clr[ , 1:3]))
  names(tmmx_df) <- row.names(theta_est)
  
  # Humidity
  rmax <- seq(rmax_boundary_knots[1], rmax_boundary_knots[2], length.out = nplt)
  spl_rmax <- ns(rmax, knots = rmax_internal_knots, Boundary.knots = rmax_boundary_knots)
  theta_rmax <- theta_est[ , paste0("est", 4:6)]
  rmax_df <- as.data.frame(spl_rmax %*% t(theta_rmax))
  rmax_clr <- as.data.frame(spl_tmmx %*% t(theta_est_clr[ , 4:6]))
  names(rmax_df) <- row.names(theta_est)
  
  # Get tree
  tr <- ccs_tree(root[i])$tr
  vids <- unlist(moretrees_results$beta_moretrees$outcomes)
  vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
  vids <- Reduce(union, vids)
  tr <- induced_subgraph(tr, vids)
  
  # Get ancestors
  leaves <- row.names(theta_est)
  d <- igraph::diameter(tr)
  ancestors <- igraph::ego(tr, order = d + 1, nodes = leaves, mode = "in")
  ancestors <- sapply(ancestors, names, simplify = F)
  ancestors <- sapply(ancestors, function(a, nodes) which(nodes %in% a), 
                      nodes = names(V(tr)),
                      simplify = F)
  names(ancestors) <- leaves
  
  # Get confidence intervals
  Omega <- moretrees_results$mod$vi_params$Omega
  theta_var <- list()
  for (v in 1:ncol(tmmx_df)) {
    theta_var[[v]] <- Reduce(`+`, Omega[ancestors[[v]]])
  }
  tmmx_cil <- tmmx_ciu <- tmmx_df
  rmax_cil <- rmax_ciu <- rmax_df
  for (v in 1:ncol(tmmx_df)) {
    tmmx_sd <- sqrt(diag(spl_tmmx %*% theta_var[[v]][1:3, 1:3] %*% t(spl_tmmx)))
    tmmx_cil[ , v] <- tmmx_df[ , v] - 1.96 * tmmx_sd
    tmmx_ciu[ , v] <- tmmx_df[ , v] + 1.96 * tmmx_sd
    rmax_sd <- sqrt(diag(spl_rmax %*% theta_var[[v]][4:6, 4:6] %*% t(spl_rmax)))
    rmax_cil[ , v] <- rmax_df[ , v] - 1.96 * rmax_sd
    rmax_ciu[ , v] <- rmax_df[ , v] + 1.96 * rmax_sd
  }

  # Melt temp
  tmmx_df$tmmx <- tmmx
  tmmx_df <- melt(tmmx_df, id.vars = "tmmx")
  names(tmmx_df)[2:3] <- c("outcome", "lrr")
  tmmx_df$outcome <- str_remove_all(tmmx_df$outcome, "\\.0")
  tmmx_df$outcome <- factor(tmmx_df$outcome,
                             levels = str_remove_all(row.names(theta_est), "\\.0"))
  tmmx_df$cil <- melt(tmmx_cil)$value
  tmmx_df$ciu <- melt(tmmx_ciu)$value
  tmmx_df$tmmx <- tmmx_df$tmmx - 273.15
  tmmx_df$rr <- exp(tmmx_df$lrr)
  tmmx_df$cil <- exp(tmmx_df$cil)
  tmmx_df$ciu <- exp(tmmx_df$ciu)
  tmmx_df$rr_clr <- exp(melt(tmmx_clr)$value)
  
  # Plot temp
  require(ggplot2)
  plt <- ggplot(tmmx_df, aes(x = tmmx, y = rr)) +
    geom_ribbon(aes(ymin = cil, ymax = ciu),
                alpha = 0.5) +
    geom_line() +
    geom_line(aes(x = tmmx, y = rr_clr), col = "red", lty = 2) +
    facet_wrap( . ~ outcome, nrow = 10) +
    theme_minimal() +
    ylab("Rate Ratio") +
    xlab("Temperature (Celsius)")  +
    ylim(0.8, 3)
  pdf(file = paste0("./figures/temp_RR_", ds, ".pdf"),
                    width = 8, height = 10)
  print(plt)
  dev.off()
  
  # Melt humidity
  rmax_df$rmax <- rmax
  rmax_df <- melt(rmax_df, id.vars = "rmax")
  names(rmax_df)[2:3] <- c("outcome", "lrr")
  rmax_df$outcome <- str_remove_all(rmax_df$outcome, "\\.0")
  rmax_df$outcome <- factor(rmax_df$outcome,
                            levels = str_remove_all(row.names(theta_est), "\\.0"))
  rmax_df$cil <- melt(rmax_cil)$value
  rmax_df$ciu <- melt(rmax_ciu)$value
  rmax_df$rmax <- rmax_df$rmax
  rmax_df$rr <- exp(rmax_df$lrr)
  rmax_df$cil <- exp(rmax_df$cil)
  rmax_df$ciu <- exp(rmax_df$ciu)
  rmax_df$rr_clr <- exp(melt(rmax_clr)$value)
  
  # Plot humidity
  require(ggplot2)
  plt <- ggplot(rmax_df, aes(x = rmax, y = rr)) +
    geom_ribbon(aes(ymin = cil, ymax = ciu),
                alpha = 0.5) +
    geom_line() +
    geom_line(aes(x = rmax, y = rr_clr), col = "red", lty = 2) +
    facet_wrap( . ~ outcome, nrow = 10) +
    theme_minimal() +
    ylab("Rate Ratio") +
    xlab("Humidity (%)") +
    ylim(0.8, 1.5)
  pdf(file = paste0("./figures/humidity_RR_", ds, ".pdf"),
      width = 8, height = 10)
  print(plt)
  dev.off()
  
}
