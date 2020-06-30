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
   load(file = paste0("./results/mod", mod, "_split0_", ds, "_sens5.Rdata"))
   mod0 <- moretrees_results
   load(file = paste0("./results/mod", mod, "_split25_", ds, "_sens5.Rdata"))
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
   
   tabfile <- paste0("./figures/mod", mod, "_", ds, "_table_sens5.tex")
   write(print(est_xtable, floating = FALSE, include.rownames = FALSE,
               sanitize.text.function = function(x) x),
         file = tabfile)
}


