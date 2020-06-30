source("./code/results_functions.R")

### CVD tree ----------------------------------------------------------------
load(file = paste0("./results/mod1_split0_cvd.Rdata"))

# Get tree
tr <- ccs_tree("7")$tr
vids <- unlist(moretrees_results$beta_moretrees$outcomes)
vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
vids <- Reduce(union, vids)
tr <- induced_subgraph(tr, vids)

# Get labels
ccs_labels <- read.csv("./data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")

# Make plot
root <- names(igraph::V(tr))[igraph::degree(tr, mode = "in") == 0]
V(tr)$levels <- as.numeric(igraph::distances(tr, v = root, to = V(tr), mode = "out") + 1)
labels <- data.frame(node = names(V(tr)),
                     nodename = str_remove_all(names(V(tr)), "\\.0"), 
                     levels = V(tr)$levels,
                     stringsAsFactors = F)
labels <- merge(labels, ccs_labels, by.x = "nodename", by.y = "ccs_code",
                all.x = T, all.y = F, sort = F)
labels$label <- as.character(labels$label)
labels$label[labels$node == "7.2.11"] <- "Heart failure and congestive heart failure"
labels$label <- paste0(labels$nodename, ": ", labels$label)
labels$nodename <- NULL

p <- ggtree(tr, ladderize = F)
labels <- labels[match(p$data$label, labels$node), ]
p$data$levels <- labels$levels
for (i in 1:4) {
  p$data[ , paste0('label.long', i)] <- as.character(NA)
  labs <- labels$label[labels$levels == i]
  if (i < 3) {
    labs <- str_wrap(labs, width = 20)
  }
  if (i == 3) {
    labs <- str_wrap(labs, width = 60)
  }
  p$data[p$data$levels == i , paste0('label.long', i)] <- labs
}
node.pos <- c(-0.4, 0.2, 1.85, 1.9)
node.nudge <- c(-0.55,-0.58, -1.6, 0)
node.angle <- c(90, 90, 0, 0)
p$data$x <- node.pos[p$data$levels]
p <- p + geom_tiplab(aes(label = label.long4)) + xlim(c(-3, 5))
for (i in 1:3) {
  p <- p + geom_nodelab(aes_string(label = paste0("label.long", i)),
                        geom = "label",
                        hjust = 0,
                        nudge_x = node.nudge[i],
                        angle = node.angle[i],
                        fill = "white",
                        label.size = NA)
}
p <- p + theme(plot.margin = unit(c(-0.3, -4.7, -0.3, -6.8), unit = "in"))
p <- p + coord_cartesian(clip = "off")
p <- p + scale_y_reverse()
p <- p + geom_text(data = data.frame(x = node.pos + node.nudge + 0.13,
                                     y = min(p$data$y) - 2,
                                     lab = paste0("Level ", 1:4)),
                   aes(x = x, y = y, label = lab),
                   size = 5)

pdf(file = "./figures/ccs_tree_cvd.pdf", width = 13, height = 15)
p
dev.off()

### RD tree ----------------------------------------------------------------
load(file = paste0("./results/mod1_split0_resp.Rdata"))

# Get tree
tr <- ccs_tree("8")$tr
vids <- unlist(moretrees_results$beta_moretrees$outcomes)
vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
vids <- Reduce(union, vids)
tr <- induced_subgraph(tr, vids)

# Get labels
ccs_labels <- read.csv("./data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")

# Make plot
root <- names(igraph::V(tr))[igraph::degree(tr, mode = "in") == 0]
V(tr)$levels <- as.numeric(igraph::distances(tr, v = root, to = V(tr), mode = "out") + 1)
labels <- data.frame(node = names(V(tr)),
                     nodename = str_remove_all(names(V(tr)), "\\.0"), 
                     levels = V(tr)$levels,
                     stringsAsFactors = F)
labels <- merge(labels, ccs_labels, by.x = "nodename", by.y = "ccs_code",
                all.x = T, all.y = F, sort = F)
labels$label <- paste0(labels$nodename, ": ", labels$label)
labels$nodename <- NULL

p <- ggtree(tr, ladderize = F)
labels <- labels[match(p$data$label, labels$node), ]
p$data$levels <- labels$levels
for (i in 1:4) {
  p$data[ , paste0('label.long', i)] <- as.character(NA)
  labs <- labels$label[labels$levels == i]
  if (i == 1) labs <- str_wrap(labs, width = 20)
  if (i == 2) {
    labs <- str_wrap(labs, width = 30)
  }
  if (i == 3) {
    labs <- str_wrap(labs, width = 60)
  }
  p$data[p$data$levels == i , paste0('label.long', i)] <- labs
}
node.pos <- c(-0.4, 0.5, 2.15, 2.2)
node.nudge <- c(-0.55,-0.85, -1.6, 0)
node.angle <- c(90, 90, 0, 0)
p$data$x <- node.pos[p$data$levels]
p <- p + geom_tiplab(aes(label = label.long4)) + xlim(c(-3.1, 5.6))
for (i in 1:3) {
  p <- p + geom_nodelab(aes_string(label = paste0("label.long", i)),
                        geom = "label",
                        hjust = 0,
                        nudge_x = node.nudge[i],
                        angle = node.angle[i],
                        fill = "white",
                        label.size = NA)
}
p <- p + theme(plot.margin = unit(c(-0.3, -4.7, -0.3, -6.8), unit = "in"))
p + coord_cartesian(clip = "off")
p <- p + scale_y_reverse()
p <- p + geom_text(data = data.frame(x = node.pos + node.nudge + 0.13,
                                     y = min(p$data$y) - 2,
                                     lab = paste0("Level ", 1:4)),
                   aes(x = x, y = y, label = lab),
                   size = 5)


pdf(file = "./figures/ccs_tree_resp.pdf", width = 14, height = 11)
p
dev.off()

# Plot prior ---------------------------------------------------------------------
dataset <- c("cvd", "resp")
root <- c("7", "8")
a_rho <- as.matrix(rbind(c(1, 1), c(0.9, 0.5)))
b_rho <- as.matrix(rbind(c(1, 1), c(3 , 2)))
prior <- c("main", "sens4", "sens5")
for (i in 1:2) {
  for (j in 1:3) {
    ds <- dataset[i]
    load(file = paste0("./results/mod1_split0_", ds, ".Rdata"))
    
    # Get tree
    tr <- ccs_tree(root[i])$tr
    vids <- unlist(moretrees_results$beta_moretrees$outcomes)
    vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
    vids <- Reduce(union, vids)
    tr <- induced_subgraph(tr, vids)
    
    # Get levels
    levels <- rep(1, length(V(tr)))
    if (j %in% c(1,2)) {
      leaves <- names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0])
      levels[names(igraph::V(tr)) %in% leaves] <- 2
    }
    
    # Create matrix plot
    rownames.lab.offset <- 18.8 * (i == 1) + 19.1 * (i == 2)
    pltfile3 <- paste0("./figures/prior_", ds, "_", prior[j],".pdf")
    pdf(file = pltfile3, width = 12, height = 9.5)
    print(equal_betas_plot(prob = NULL,
                           show.groups = F,
                           a_rho = a_rho[2 - j %% 2, ], b_rho = b_rho[2 - j %% 2, ],
                           tr = tr,
                           rownames.lab.offset = rownames.lab.offset,
                           levels = levels))
    dev.off()
  }
}

# prior simulations
Nsims <- 100000
dataset <- c("cvd", "resp")
root <- c("7", "8")
a_rho <- as.matrix(rbind(c(1, 1), c(0.9, 0.5)))
b_rho <- as.matrix(rbind(c(1, 1), c(3 , 2)))
prior <- c("main", "sens", "sens2")
require(doParallel)
registerDoParallel(cores = detectCores())
require(foreach)
for (j in 1:3) {
  for (i in 1:2) {
    ds <- dataset[i]
    load(file = paste0("./results/mod1_split0_", ds, ".Rdata"))
    
    # Get tree
    tr <- ccs_tree(root[i])$tr
    vids <- unlist(moretrees_results$beta_moretrees$outcomes)
    vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
    vids <- Reduce(union, vids)
    tr <- induced_subgraph(tr, vids)
    
    # Get levels
    levels <- rep(1, length(V(tr)))
    if (j %in% c(1,2)) {
      levels[names(igraph::V(tr)) %in% leaves] <- 2
    }
    
    # Get ancestor matrix
    A <- igraph::as_adjacency_matrix(tr, sparse = T)
    A <- Matrix::expm(Matrix::t(A))
    A[A > 0 ] <- 1 
    A <- Matrix::Matrix(A, sparse = T)
    A_leaf <- A[V(tr)$leaf, ]
    
    # Run sims
    simsout <- foreach(i = 1:Nsims, .combine = rbind) %dopar% 
      sim.prior.fun(levels, A_leaf, a_rho = a_rho[2 - j %% 2, ],
                    b_rho = b_rho[2 - j %% 2, ])
    # plot(x, dbeta(x, shape1 = 0.5, shape2 = 2), type = "l")
    
    # Make plot
    if (i == 1) breaks = c(1, 20, 40, nrow(A_leaf))
    if (i == 2) breaks = c(1, 10, 20, nrow(A_leaf))
    plt_n <- ggplot(simsout, aes(x = n.groups, y = (..count..)/sum(..count..))) +
      geom_bar() +
      theme_minimal() +
      scale_x_continuous(breaks = breaks,
                         limits = c(0, nrow(A_leaf) + 1)) +
      xlab("Number of Groups") +
      ylab("Prior Probability")
    
    # Median conditional on number of groups
    if (i == 1) {
      breaks = c(0, 5, 10, 20, 58)
      labels = c("1 to 5", "6 to 10", "11 to 20", "21 to 57")
    }
    if (i == 2) {
      breaks = c(0, 5, 10, 20, 32)
      labels = c("1 to 5", "6 to 10", "11 to 20", "21 to 32")
    }
    simsout$n.groups.cat <- cut(x = simsout$n.groups, breaks = breaks)
    
    # Mean conditional on number of groups
    plt_med <- ggplot(simsout, aes(x = n.groups.cat, y = mean.size / median.size)) +
      geom_boxplot() + scale_y_continuous(trans = "log10") +
      xlab("Number of Groups") +
      ylab("Mean / Median Group Size") +
      scale_x_discrete(labels = labels) + 
      theme_minimal()
    
    pdf(file = paste0("./figures/prior_", dataset[i], "_groups_", prior[j],".pdf"),
        width = 8, height = 3)
    grid.arrange(plt_n, plt_med, nrow = 1)
    dev.off()
  }
}

# Covariate plotes --------------------------------------------------------------------------


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
    ylab("Relative Rate") +
    xlab("Temperature (Celsius)") 
  if (i == 1) {
    plt <- plt + ylim(0.8, 2)
  } else {
    plt <- plt + ylim(0.8, 3)
  }
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
    ylab("Relative Rate") +
    xlab("Humidity (%)")
  if (i == 1) {
    plt <- plt + ylim (0.8, 1.2)
  } else {
    plt <- plt + ylim(0.8, 1.5)
  }
  pdf(file = paste0("./figures/humidity_RR_", ds, ".pdf"),
      width = 8, height = 10)
  print(plt)
  dev.off()
  
}

# Nested CLR plots --------------------------------------------------------------------------

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
  ds <- dataset[i]
  
  load(paste0("./results/nest_pltdat_", ds, ".RData"))
  pltdat_ml[ , c("est1", "cil1", "ciu1")] <- 100 * (exp(pltdat_ml[, c("est1", "cil1", "ciu1")] * 10 ) - 1)
  
  # plot
  xlab <- "Excess Rate (%)"
  if (i == 1) {
    lab.widths <- c(0.75, 0.75, 0.75, 0.85)
    lab.txt.width <- c(20, 25, 25, 15)
    pdf(file = "./figures/cvd_nested_clr.pdf", height = 15, width = 14)
    beta_indiv_plot_fun(pltdat_ml, tr, xlab = xlab, lab.widths = lab.widths,
                        lab.txt.width = lab.txt.width, axis.height = 1.5,
                        cil_min = -8, ciu_max = 8,
                        force_lims = TRUE,
                        lab.txt.size = 4, digits = 0, axis.txt.size = 10,
                        plot_depth = 4, wrap_labs = FALSE)
    dev.off()
  }
  if (i == 2) {
    lab.widths <- c(0.75, 0.75, 0.75, 0.75)
    lab.txt.width <- c(20, 25, 25, 15)
    pdf(file = "./figures/resp_nested_clr.pdf", height = 10, width = 14)
    beta_indiv_plot_fun(pltdat_ml, tr, xlab = xlab, lab.widths = lab.widths,
                        lab.txt.width = lab.txt.width, axis.height = 1.3,
                        cil_min = -8, ciu_max = 8,
                        force_lims = TRUE,
                        lab.txt.size = 4, digits = 0, axis.txt.size = 10,
                        plot_depth = 4, wrap_labs = FALSE)
    dev.off()
  }
}

