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

# Make horrible plot
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

# Make horrible plot
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
prior <- c("main", "sens", "sens2")
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
