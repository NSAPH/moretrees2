# ----------------------------------------------------------------------------
# ----------------------------- CVD results ----------------------------------
# ----------------------------------------------------------------------------

source("./code/results_functions.R")
hyper_method <- c("full", "EB")
dataset <- c("cvd", "resp")
root <- c("7", "8")
splits <- c("0", "25", "35")
nmods <- 3
approx_ml <- array(dim = c(length(dataset), length(splits), nmods, length(hyper_method)), 
                   dimnames = list(dataset = c("Cardiovascular Data", "Respiratory Data"), 
                                   split = splits, 
                                   Model = 1:nmods,
                                   Method = c("Fully Bayes", "Empirical Bayes")))

for (i in 1:length(dataset)) { # datasets
   for (j in 1:length(splits)) { # splits
      for (mod in 1:nmods) { # models
         for(k in 1:length(hyper_method)) { # method for selecting hyperparams
            
            ds <- dataset[i]
            spl <- splits[j]
            hm <- hyper_method[k]
            
            # Read in results of moretrees model 
            load(file = paste0("./results/mod", mod, "_split", spl, "_", ds, "_", hm, ".Rdata"))
            
            # Get approximate marginal likelihood
            approx_ml[i, j, mod, k] <- moretrees_results$mod$hyperparams$ELBO
            
            # Get tree
            tr <- ccs_tree(root[i])$tr
            vids <- unlist(moretrees_results$beta_moretrees$outcomes)
            vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
            vids <- Reduce(union, vids)
            tr <- induced_subgraph(tr, vids)
            
            # Create results table
            OR_est <- ccs_table(root = root[i], moretrees_results = moretrees_results,
                                tr = tr, digits = 3, mult = 10)
            OR_est$long_label <- NULL
            require(xtable)
            row.names(OR_est) <- NULL
            OR_est$n_obs <- formatC(OR_est$n_obs, format="d", big.mark=",")
            k <- length(moretrees_results$mod$vi_params$mu[[1]])
            align <- c("l", "l", "p{6.5cm}", "r", "r", rep("p{2.2cm}", k))
            display <- c("d", "d", "s", "d", "d", rep("f", k))
            if (spl == "0") {
               tabnames <- c("Group", "CCS codes",
                             "$n_{out}$", "$n_{obs}$",
                             "RR (96\\%CI)")
            } else {
               tabnames <- c("Group", "CCS codes",
                             "$n_{out}$", "$n_{obs}$",
                             paste0("RR below $", spl, " \\mu g \\cdot m^{-3}$ (95\\%CI)"),
                             paste0("RR above $", spl, " \\mu g \\cdot m^{-3}$ (95\\%CI)"))
            }
            
            # make xtable
            OR_xtable <- xtable(OR_est, align = align,
                                digits = 3, display = display)
            names(OR_xtable) <- tabnames
            
            tabfile <- paste0("./figures/mod", mod, "_split", spl, "_", ds, "_", hm, "_table.tex")
            write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
                        sanitize.text.function = function(x) x),
                  file = tabfile)
            
            # Create results table with ML estimates
            OR_est <- ccs_table(root = root[i], moretrees_results = moretrees_results,
                                type = "ml",
                                digits = 3, mult = 10)
            OR_est$short_label <- NULL
            OR_est$long_label <- NULL
            OR_est$n_outcomes <- NULL
            OR_est$n_obs <- NULL
            row.names(OR_est) <- NULL
            align <- c("l", "l", rep("l", k))
            display <- c("d", "d", rep("f", k))
            OR_xtable <- xtable(OR_est, align = align,
                                digits = 3, display = display)
            names(OR_xtable) <- tabnames[-c(2, 3, 4)]
            tabfile <- paste0("./figures/mod", mod, "_split", spl, "_", ds, "_", hm, "_ml_table.tex")
            write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
                        sanitize.text.function = function(x) x),
                  file = tabfile)
            # 
            #          # Create results tree plot
            #          pltfile <- paste0("./figures/mod", mod, "_split", spl, "_", ds, "_", hm, "_tree.pdf")
            #          pdf(file = pltfile, width = 6, height = 2.8)
            #          ccs_plot(root = "8", moretrees_results = moretrees_results,
            #                   tr = tr,
            #                   asp = 1/10, leaf.height = 30, label.dist = 4.5)
            #          dev.off()
            
            # Create matrix plot
            rownames.lab.offset <- 18.8 * (i == 1) + 19.1 * (i == 2)
            pltfile2 <- paste0("./figures/mod", mod, "_split", spl, "_", ds, "_", hm, "_matrix.pdf")
            pdf(file = pltfile2, width = 12, height = 9.5)
            print(equal_betas_plot(prob = moretrees_results$mod$vi_params$prob,
                                   groups = moretrees_results$beta_est$group,
                                   tr = tr,
                                   rownames.lab.offset = rownames.lab.offset))
            dev.off()
         }
      }
   }
}

ccs_labels <- read.csv("./data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")

root <- names(igraph::V(tr))[igraph::degree(tr, mode = "in") == 0]
V(tr)$levels <- as.numeric(igraph::distances(tr, v = root, to = nodes, mode = "out") + 1)
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
p <- p + theme(plot.margin = unit(c(0, -4.7, 0, -6.8), unit = "in"))
p + coord_cartesian(clip = "off")

pdf(file = "./figures/test_tree.pdf", width = 13, height = 15)
p
dev.off()

# Plot prior
for (i in 1:2) {
   ds <- dataset[i]
   load(file = paste0("./results/mod", mod, "_split", spl, "_", ds, "_", hm, ".Rdata"))
   # Get tree
   tr <- ccs_tree(root[i])$tr
   vids <- unlist(moretrees_results$beta_moretrees$outcomes)
   vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
   vids <- Reduce(union, vids)
   tr <- induced_subgraph(tr, vids)
   # Create matrix plot
   rownames.lab.offset <- 18.8 * (i == 1) + 19.1 * (i == 2)
   pltfile3 <- paste0("./figures/prior_", ds, ".pdf")
   pdf(file = pltfile3, width = 12, height = 9.5)
   print(equal_betas_plot(prob = rep(0.5, length(moretrees_results$mod$vi_params$prob)),
                          show.groups = F,
                          tr = tr,
                          rownames.lab.offset = rownames.lab.offset))
   dev.off()
}

# Plot Bayes' factors
require(reshape2)
approx_ml2 <- melt(approx_ml)
approx_ml2$split <- factor(approx_ml2$split)
approx_ml2$Model <- factor(approx_ml2$Model)
ml_plot <- ggplot(approx_ml2) + 
   geom_point(aes(x = Model, y = value, shape = split), size = 2) +
   facet_grid(dataset ~ Method, scales = "free_y") +
   theme_minimal() + xlab("Model") + 
   ylab(expression("Lower Bound for "*pi(Y))) +
   scale_shape_discrete(name = expression(PM[2.5]*" break"),
                        labels = c("No break",
                                   expression("25"*mu*"g"*m^-3),
                                   expression("35"*mu*"g"*m^-3)),
                        solid = F)
pdf(file = paste0("./figures/vb_bound.pdf"), width = 6, height = 4)
ml_plot
dev.off()

# Best model
approx_ml_cvd <- approx_ml2[approx_ml2$dataset == "Cardiovascular Data", ]
approx_ml_cvd[which.max(approx_ml_cvd$value), ]

approx_ml_resp <- approx_ml2[approx_ml2$dataset == "Respiratory Data", ]
approx_ml_resp[which.max(approx_ml_resp$value), ]



