# ----------------------------------------------------------------------------
# ----------------------------- CVD results ----------------------------------
# ----------------------------------------------------------------------------

source("./code/results_functions.R")
dataset <- c("cvd", "resp")
root <- c("7", "8")
splits <- c("0", "25", "35")
nmods <- 3
bf <- array(dim = c(length(dataset), length(splits), nmods), 
            dimnames = list("dataset" = c("Cardiovascular Data", "Respiratory Data"), 
                            split = splits, "Model" = 1:nmods))

for (i in 1:length(dataset)) { # datasets
   for (j in 1:length(splits)) { # splits
      for (mod in 1:nmods) { # models
         spl <- splits[j]
         # Read in results of moretrees model 
         load(file = paste0("./results/mod", mod, "_split", spl, "_northEast_", dataset[i], ".Rdata"))
         
         bf[i, j, mod] <- moretrees_results$mod$hyperparams$ELBO
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
         
         tabfile <- paste0("./figures/mod", mod, "_split", spl, "_northEast_", dataset[i], "_table.tex")
         write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
                     sanitize.text.function = function(x) x),
               file = tabfile)
         # add Bayes factor to table
         # write(x = paste0("Approximate Bayes Factor = ", round(moretrees_results$mod$hyperparams$ELBO)),
               # file = tabfile, append = T)
         
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
         tabfile <- paste0("./figures/mod", mod, "_split", spl, "_northEast_", dataset[i], "_ml_table.tex")
         write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
                     sanitize.text.function = function(x) x),
               file = tabfile)
         
         # Create results tree plot 
         pltfile <- paste0("./figures/mod", mod, "_split", spl, "_northEast_", dataset[i], "_tree.pdf")
         pdf(file = pltfile, width = 6, height = 2.8)
         ccs_plot(root = "8", moretrees_results = moretrees_results,
                  tr = tr,
                  asp = 1/10, leaf.height = 30, label.dist = 4.5)
         dev.off()

      }
   }
}

# Plot Bayes' factors
require(reshape2)
bf2 <- melt(bf)
bf2$split <- factor(bf2$split)
bf2$Model <- factor(bf2$Model)
bfplot <- ggplot(bf2) + 
   geom_point(aes(x = Model, y = value, shape = split), size = 2) +
   facet_wrap(. ~ dataset, nrow = 1, scales = "free_y") +
   theme_minimal() + xlab("Model") + ylab("Bayes Factor") +
   scale_shape_discrete(name = expression(PM[2.5]*" break"),
       labels = c("No break",
                  expression("25"*mu*"g"*m^-3),
                  expression("35"*mu*"g"*m^-3)),
       solid = F)
pdf(file = "./figures/bf_northEast.pdf", width = 6, height = 2)
bfplot
dev.off()
