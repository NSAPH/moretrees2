# ----------------------------------------------------------------------------
# ----------------------------- CVD results ----------------------------------
# ----------------------------------------------------------------------------

source("./code/results_functions.R")
require(xtable)

dataset <- c("cvd", "resp")
root <- c("7", "8")
splits <- c("0", "25")
nmods <- 3
approx_ml <- array(dim = c(length(dataset), length(splits), 2), 
                   dimnames = list(dataset = c("Cardiovascular Data", "Respiratory Data"), 
                                   split = splits, 
                                   Model = 2:nmods))

for (i in 1:length(dataset)) { # datasets
   ds <- dataset[i]
   for (mod in 1:nmods) { # models
      
      # Read in results of moretrees model 
      load(file = paste0("./results/mod", mod, "_split0_", ds, "_sens.Rdata"))
      approx_ml[i, 1, mod - 1] <- moretrees_results$mod$hyperparams$ELBO
      mod0 <- moretrees_results
      load(file = paste0("./results/mod", mod, "_split25_", ds, "_sens.Rdata"))
      approx_ml[i, 2, mod - 1] <- moretrees_results$mod$hyperparams$ELBO
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
      
      tabfile <- paste0("./figures/mod", mod, "_", ds, "_table_sens.tex")
      write(print(est_xtable, floating = FALSE, include.rownames = FALSE,
                  sanitize.text.function = function(x) x),
            file = tabfile)
      
      # # Create results table with ML estimates
      # est_ml0 <- ccs_table(root = root[i], moretrees_results = mod0,
      #                      type = "ml",
      #                      digits = 1, mult = 10)
      # est_ml0$short_label <- NULL
      # est_ml0$long_label <- NULL
      # est_ml0$n_outcomes <- NULL
      # est_ml0$n_obs <- NULL
      # row.names(est_ml0) <- NULL
      # est_ml0$ci_est2 <- est_ml0$ci_est1
      # est_ml0 <- cbind("Model" = 1, est_ml0)
      # est_ml25 <- ccs_table(root = root[i], moretrees_results = mod25,
      #                       type = "ml",
      #                       digits = 1, mult = 10)
      # est_ml25$short_label <- NULL
      # est_ml25$long_label <- NULL
      # est_ml25$n_outcomes <- NULL
      # est_ml25$n_obs <- NULL
      # row.names(est_ml0) <- NULL
      # est_ml25 <- cbind("Model" = 2, est_ml25)
      # est_ml <- rbind(est_ml0, est_ml25)
      # align <- c("l", "l", "l", rep("r", 2))
      # display <- c("d", "d", "d", rep("f", 2))
      # est_xtable_ml <- xtable(est_ml, align = align,
      #                         digits = 3, display = display)
      # names(est_xtable_ml) <- tabnames[-c(3, 4, 5)]
      # tabfile <- paste0("./figures/mod", mod, "_", ds, "_ml_table.tex")
      # write(print(est_xtable_ml, floating = FALSE, include.rownames = FALSE,
      #             sanitize.text.function = function(x) x),
      #       file = tabfile)
      # 
      # # Create matrix plots
      # rownames.lab.offset <- 18.8 * (i == 1) + 19.1 * (i == 2)
      # pltfile2 <- paste0("./figures/mod", mod, "_split0_", ds, "_matrix.pdf")
      # pdf(file = pltfile2, width = 12, height = 9.5)
      # print(equal_betas_plot(prob = mod0$mod$vi_params$prob,
      #                        groups = mod0$beta_est$group,
      #                        tr = tr,
      #                        rownames.lab.offset = rownames.lab.offset))
      # dev.off()
      # pltfile2 <- paste0("./figures/mod", mod, "_split25_", ds, "_matrix.pdf")
      # pdf(file = pltfile2, width = 12, height = 9.5)
      # print(equal_betas_plot(prob = mod25$mod$vi_params$prob,
      #                        groups = mod25$beta_est$group,
      #                        tr = tr,
      #                        rownames.lab.offset = rownames.lab.offset))
      # dev.off()
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
pdf(file = paste0("./figures/vb_bound_sens.pdf"), width = 6, height = 4)
ml_plot
dev.off()

# Best model
approx_ml_cvd <- approx_ml2[approx_ml2$dataset == "Cardiovascular Data", ]
approx_ml_cvd[which.max(approx_ml_cvd$value), ]

approx_ml_resp <- approx_ml2[approx_ml2$dataset == "Respiratory Data", ]
approx_ml_resp[which.max(approx_ml_resp$value), ]



