# ----------------------------------------------------------------------------
# ----------------------------- CVD results ----------------------------------
# ----------------------------------------------------------------------------

source("./code/results_functions.R")
require(xtable)

dataset <- c("cvd", "resp")
root <- c("7", "8")
splits <- c("0", "25")
nmods <- 3
approx_ml <- array(dim = c(length(dataset), length(splits), nmods), 
                   dimnames = list(dataset = c("Cardiovascular Data", "Respiratory Data"), 
                                   split = splits, 
                                   Model = 1:nmods))

for (i in 1:length(dataset)) { # datasets
   for (j in 1:length(splits)) { # splits
      for (mod in 1:nmods) { # models
         
         ds <- dataset[i]
         spl <- splits[j]
         
         # Read in results of moretrees model 
         load(file = paste0("./results/mod", mod, "_split", spl, "_", ds, ".Rdata"))
         
         # Get approximate marginal likelihood
         approx_ml[i, j, mod] <- moretrees_results$mod$hyperparams$ELBO
         
         # Get tree
         tr <- ccs_tree(root[i])$tr
         vids <- unlist(moretrees_results$beta_moretrees$outcomes)
         vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
         vids <- Reduce(union, vids)
         tr <- induced_subgraph(tr, vids)
         
         # Create results table
         est <- ccs_table(root = root[i], moretrees_results = moretrees_results,
                             tr = tr, digits = 1, mult = 10)
         est$long_label <- NULL
         row.names(est) <- NULL
         est$n_obs <- formatC(est$n_obs, format="d", big.mark=",")
         k <- length(moretrees_results$mod$vi_params$mu[[1]])
         align <- c("l", "l", "p{6.5cm}", "r", "r", rep("p{2.2cm}", k))
         display <- c("d", "d", "s", "d", "d", rep("f", k))
         if (spl == "0") {
            tabnames <- c("Group", "CCS codes",
                          "$n_{out}$", "$n_{obs}$",
                          "\\% change in rate (95\\%CI)")
         } else {
            tabnames <- c("Group", "CCS codes",
                          "$n_{out}$", "$n_{obs}$",
                          paste0("\\% change in rate below $", spl, " \\mu g \\cdot m^{-3}$ (95\\%CI)"),
                          paste0("\\% change in rate above $", spl, " \\mu g \\cdot m^{-3}$ (95\\%CI)"))
         }
         
         # make xtable
         est_xtable <- xtable(est, align = align,
                             digits = 1, display = display)
         names(est_xtable) <- tabnames
         
         tabfile <- paste0("./figures/mod", mod, "_split", spl, "_", ds, "_table.tex")
         write(print(est_xtable, floating = FALSE, include.rownames = FALSE,
                     sanitize.text.function = function(x) x),
               file = tabfile)
         
         # Create results table with ML estimates
         est_ml <- ccs_table(root = root[i], moretrees_results = moretrees_results,
                             type = "ml",
                             digits = 3, mult = 10)
         est_ml$short_label <- NULL
         est_ml$long_label <- NULL
         est_ml$n_outcomes <- NULL
         est_ml$n_obs <- NULL
         row.names(est_ml) <- NULL
         align <- c("l", "l", rep("l", k))
         display <- c("d", "d", rep("f", k))
         est_xtable_ml <- xtable(est_ml, align = align,
                             digits = 3, display = display)
         names(est_xtable_ml) <- tabnames[-c(2, 3, 4)]
         tabfile <- paste0("./figures/mod", mod, "_split", spl, "_", ds, "_ml_table.tex")
         write(print(est_xtable_ml, floating = FALSE, include.rownames = FALSE,
                     sanitize.text.function = function(x) x),
               file = tabfile)
         
         # Create matrix plot
         rownames.lab.offset <- 18.8 * (i == 1) + 19.1 * (i == 2)
         pltfile2 <- paste0("./figures/mod", mod, "_split", spl, "_", ds, "_matrix.pdf")
         pdf(file = pltfile2, width = 12, height = 9.5)
         print(equal_betas_plot(prob = moretrees_results$mod$vi_params$prob,
                                groups = moretrees_results$beta_est$group,
                                tr = tr,
                                rownames.lab.offset = rownames.lab.offset))
         dev.off()
      }
   }
}

# Plot Bayes' factors
require(reshape2)
approx_ml2 <- melt(approx_ml)
approx_ml2$split <- factor(approx_ml2$split)
approx_ml2$Model <- factor(approx_ml2$Model)
ml_plot <- ggplot(approx_ml2) + 
   geom_point(aes(x = Model, y = value, shape = split), size = 2) +
   facet_wrap(dataset ~ ., nrow = 1, scales = "free_y") +
   theme_minimal() + xlab("Model") + 
   ylab(expression("Lower Bound for "*pi(Y))) +
   scale_shape_discrete(name = expression(PM[2.5]*" break"),
                        labels = c("No break",
                                   expression("25"*mu*"g"*m^-3)),
                        solid = F)
pdf(file = paste0("./figures/vb_bound.pdf"), width = 6, height = 4)
ml_plot
dev.off()

# Best model
approx_ml_cvd <- approx_ml2[approx_ml2$dataset == "Cardiovascular Data", ]
approx_ml_cvd[which.max(approx_ml_cvd$value), ]

approx_ml_resp <- approx_ml2[approx_ml2$dataset == "Respiratory Data", ]
approx_ml_resp[which.max(approx_ml_resp$value), ]



