# ----------------------------------------------------------------------------
# ----------------------------- CVD results ----------------------------------
# ----------------------------------------------------------------------------

source("./code/results_functions.R")
require(xtable)

dataset <- c("cvd", "resp")
root <- c("7", "8")
splits <- c("0", "25")
mod <- 3

for (i in 1:length(dataset)) { # datasets
  ds <- dataset[i]
  
  # Read in results of moretrees model 
  load(file = paste0("./results/mod", mod, "_split0_", ds, ".Rdata"))
  mod0 <- moretrees_results
  load(file = paste0("./results/mod", mod, "_split25_", ds, ".Rdata"))
  mod25 <- moretrees_results
  
  # Get tree
  tr <- ccs_tree(root[i])$tr
  vids <- unlist(moretrees_results$beta_moretrees$outcomes)
  vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
  vids <- Reduce(union, vids)
  tr <- induced_subgraph(tr, vids)
  
  # Create results table
  est0 <- ccs_table(root = root[i], moretrees_results = mod0,
                    tr = tr, digits = 1, mult = 10)
  est0$long_label <- NULL
  row.names(est0) <- NULL
  est0$n_obs <- formatC(est0$n_obs, format="d", big.mark=",")
  est0 <- cbind("Model" = 1, est0)
  est0$ci_est2 <- est0$ci_est1
  est25 <- ccs_table(root = root[i], moretrees_results = mod25,
                     tr = tr, digits = 1, mult = 10)
  est25$long_label <- NULL
  row.names(est25) <- NULL
  est25$n_obs <- formatC(est25$n_obs, format="d", big.mark=",")
  est25 <- cbind("Model" = 2, est25)
  est <- rbind(est0, est25)
  
  align <- c("l", "l", "l", "p{6.5cm}", "r", "r", rep("r", 2))
  display <- c("d", "d", "d", "s", "d", "d", rep("f", 2))
  tabnames <- c("Model", "Group", "CCS codes",
                "$n_{out}$", "$n_{obs}$",
                paste0("\\% change in rate below 25 \\mu g \\cdot m^{-3}$ (95\\%CI)"),
                paste0("\\% change in rate above 25 \\mu g \\cdot m^{-3}$ (95\\%CI)"))
  
  # make xtable
  est_xtable <- xtable(est, align = align,
                       digits = 1, display = display)
  names(est_xtable) <- tabnames
  
  tabfile <- paste0("./figures/mod", mod, "_", ds, "_table.tex")
  write(print(est_xtable, floating = FALSE, include.rownames = FALSE,
              sanitize.text.function = function(x) x),
        file = tabfile)
  
  # Create results table with ML estimates
  est_ml0 <- ccs_table(root = root[i], moretrees_results = mod0,
                       type = "ml",
                       digits = 1, mult = 10)
  est_ml0$short_label <- NULL
  est_ml0$long_label <- NULL
  est_ml0$n_outcomes <- NULL
  est_ml0$n_obs <- NULL
  row.names(est_ml0) <- NULL
  est_ml0$ci_est2 <- est_ml0$ci_est1
  est_ml0 <- cbind("Model" = 1, est_ml0)
  est_ml25 <- ccs_table(root = root[i], moretrees_results = mod25,
                        type = "ml",
                        digits = 1, mult = 10)
  est_ml25$short_label <- NULL
  est_ml25$long_label <- NULL
  est_ml25$n_outcomes <- NULL
  est_ml25$n_obs <- NULL
  row.names(est_ml0) <- NULL
  est_ml25 <- cbind("Model" = 2, est_ml25)
  est_ml <- rbind(est_ml0, est_ml25)
  align <- c("l", "l", "l", rep("r", 2))
  display <- c("d", "d", "d", rep("f", 2))
  est_xtable_ml <- xtable(est_ml, align = align,
                          digits = 3, display = display)
  names(est_xtable_ml) <- tabnames[-c(3, 4, 5)]
  tabfile <- paste0("./figures/mod", mod, "_", ds, "_ml_table.tex")
  write(print(est_xtable_ml, floating = FALSE, include.rownames = FALSE,
              sanitize.text.function = function(x) x),
        file = tabfile)
  
  # Create tree plots
  leaves <- names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0])
  groups.df25 <- data.frame(leaves = leaves, 
                            Group25 = as.factor(mod25$beta_est$group))
  mod0$tr <- tr
  class(mod0) <- "moretrees_result"
  cols_g <- RColorBrewer::brewer.pal(max(mod0$beta_moretrees$group), "Set3")
  p <- plot(mod0, group.text.size = 5,
            group.text.offset = 0.5,
            legend.text.size = 5) %<+% groups.df25 + 
    geom_tippoint(ggplot2::aes(fill = Group25),
                  shape = 21, size = 5,
                  color = "white",
                  stroke = 0.0001,
                  position = position_nudge(x = -1.1)) +
    geom_tiplab(ggplot2::aes(label = Group25),
                size = 5,
                offset = -1.6,
                vjust = 0.5,
                hjust = 0.5) +
    labs(colour = "Model 1 Groups:", fill = "Model 2 Groups:") + 
    theme(legend.position="bottom",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.margin = margin(-8, 20, 0, 20),
          legend.box.margin = margin(0, 0, 0, 0)) +
    guides(colour = guide_legend(order = 1, nrow = 1), 
           fill = guide_legend(order = 2, nrow = 1))
  if (i == 1) p <- p + scale_fill_manual(values = cols_g[c(2, 1, 3, 5, 6)])
  pltfile <- paste0("./figures/mod", mod, "_", ds, "_tree.pdf")
  pdf(file = pltfile, width = 12, height = 2.1)
  p
  dev.off()
  
  # Create matrix plots
  rownames.lab.offset <- 18.8 * (i == 1) + 19.1 * (i == 2)
  pltfile2 <- paste0("./figures/mod", mod, "_split0_", ds, "_matrix.pdf")
  pdf(file = pltfile2, width = 12, height = 9.5)
  print(equal_betas_plot(prob = mod0$mod$vi_params$prob,
                         groups = mod0$beta_est$group,
                         tr = tr,
                         rownames.lab.offset = rownames.lab.offset))
  dev.off()
  pltfile2 <- paste0("./figures/mod", mod, "_split25_", ds, "_matrix.pdf")
  pdf(file = pltfile2, width = 12, height = 9.5)
  print(equal_betas_plot(prob = mod25$mod$vi_params$prob,
                         groups = mod25$beta_est$group,
                         tr = tr,
                         rownames.lab.offset = rownames.lab.offset))
  dev.off()
}

# Plot cross-validation results
require(reshape2)
nfolds <- 10
colnms <- c("Dataset",
            "Model",
            "fold",
            "MOReTreeS",
            "CLR\n(MOReTrees)",
            "CLR\n(Level 1)",
            "CLR\n(Level 2)",
            "CLR\n(Level 3)",
            "CLR\n(Level 4)")
datasetnms <- c("CVD Dataset", "RD Dataset")
cv.res <- as.data.frame(matrix(nrow = 0, ncol = length(colnms)))
names(cv.res) <- colnms
for(i in 1:length(dataset)){
  for (j in 1:length(splits)) {
    load(paste0("./results/cv_mod3_split", splits[j], "_", dataset[i], ".RData"))
    ll.cv <- cbind(rep(datasetnms[i], nfolds),
                   rep(paste0("Model ", j), nfolds),
                   ll.cv)
    names(ll.cv) <- colnms
    cv.res <- rbind(cv.res, ll.cv)
  }
}
cv.res$`CLR\n(MOReTrees)` <- NULL
cv.df <- reshape(cv.res, direction = "long",
                 varying = list(colnms[c(4,6:9)]),
                 times = colnms[c(4,6:9)])
names(cv.df)[4:5] <- c("Method", "ll")
cv.df$Method <- factor(cv.df$Method, levels = colnms[4:9])
cv.df$Model <- factor(cv.df$Model, levels = c("Model 1", "Model 2"))
cv.df$id <- NULL

# Best model CVD
cv.cvd <- subset(cv.res, Dataset == "CVD Dataset")
cv.cvd <- cbind("mod1" = cv.cvd[1:10, 4:8],
                "mod2" = cv.cvd[11:20, 4:8])
cv.cvd.min <- apply(cv.cvd, 1,
                    function(df) names(df)[which.max(df)])
ll.cv.cvd <- c(colMeans(cv.res[cv.res$Dataset == "CVD Dataset" & cv.res$Model == "Model 1", 4:8]), 
               colMeans(cv.res[cv.res$Dataset == "CVD Dataset" & cv.res$Model == "Model 2", 4:8]))
names(ll.cv.cvd)[1:5] <-paste0(names(ll.cv.cvd)[1:5], "_mod1")
names(ll.cv.cvd)[6:10] <-paste0(names(ll.cv.cvd)[6:10], "_mod2")
names(ll.cv.cvd)[order(ll.cv.cvd, decreasing = T)]

# Best model RD
cv.resp <- subset(cv.res, Dataset == "RD Dataset")
cv.resp <- cbind("mod1" = cv.resp[1:10, 4:8],
                 "mod2" = cv.resp[11:20, 4:8])
cv.resp.min <- apply(cv.resp, 1,
                     function(df) names(df)[which.max(df)])
ll.cv.resp <- c(colMeans(cv.res[cv.res$Dataset == "RD Dataset" & cv.res$Model == "Model 1", 4:8]), 
                colMeans(cv.res[cv.res$Dataset == "RD Dataset" & cv.res$Model == "Model 2", 4:8]))
names(ll.cv.resp)[1:5] <-paste0(names(ll.cv.resp)[1:5], "_mod1")
names(ll.cv.resp)[6:10] <-paste0(names(ll.cv.resp)[6:10], "_mod2")
names(ll.cv.resp)[order(ll.cv.resp, decreasing = T)]

# Plot CV results
cv.plot <- ggplot(cv.df, aes(x = Method, y = ll, fill = Model)) + 
  geom_boxplot() +
  facet_wrap(. ~ Dataset, ncol = 2, scales = "free_y") +
  theme_bw(base_size = 22) + 
  scale_fill_grey(start = 0.5, end = 0.8) +
  xlab("Method") +
  ylab("Mean log likelihood\nin test set")

pdf("./figures/cv_plot.pdf",width = 18, height = 4)
cv.plot
dev.off()

