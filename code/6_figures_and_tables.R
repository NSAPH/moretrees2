
# CVD mod1 -------------------------------------------------------------------
source("./code/results_functions.R")

# Read in results of moretrees model 
load(file = "./results/mod1_split35_northEast.Rdata")

# Create results table 
OR_est <- ccs_table(root = "7", moretrees_results = moretrees_results,
                    digits = 3, mult = 10)
OR_est$long_label <- NULL
require(xtable)
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est,align = c("l", "l", "p{7.8cm}", "c", "p{2.2cm}", "p{2.2cm}"), 
                    digits = 3, display = c("d", "d", "s", "d", "f", "f"))
names(OR_xtable) <- c("Group", "CCS codes", "$n_{out}$", 
                      "RR below $35 \\mu g m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table1.tex")

# Create results table with ML estimates 
OR_est <- ccs_table(root = "7", moretrees_results = moretrees_results,
                    type = "ml",
                    digits = 3, mult = 10)
OR_est$short_label <- NULL
OR_est$long_label <- NULL
OR_est$n_outcomes <- NULL
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "l", "l"), 
                    digits = 3, display = c("d", "d", "f", "f"))
names(OR_xtable) <- c("Group", 
                      "RR below $35 \\mu g m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table1_ml.tex")

# Create results tree plot 
pdf(file = "./figures/mod1_cvd_tree.pdf", width = 6, height = 2.8)
ccs_plot(root = "7", moretrees_results = moretrees_results,
         asp = 1/10, leaf.height = 30, label.dist = 4.5)
dev.off()

# Respiratory disease mod1 ------------------------------------------------------------
rm(list = ls())
source("./code/results_functions.R")

# Read in results of moretrees model 
load(file = "./results/mod1_split35_northEast_resp_run2.Rdata")

# Get tree
tr <- ccs_tree("8")$tr
vids <- unlist(moretrees_results$beta_moretrees$outcomes)
vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
vids <- Reduce(union, vids)
tr <- induced_subgraph(tr, vids)

# Create results table 
OR_est <- ccs_table(root = "8", moretrees_results = moretrees_results,
                    tr = tr, digits = 3, mult = 10)
OR_est$long_label <- NULL
require(xtable)
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est,align = c("l", "l", "p{7.8cm}", "c", "p{2.2cm}", "p{2.2cm}"), 
                    digits = 3, display = c("d", "d", "s", "d", "f", "f"))
names(OR_xtable) <- c("Group", "CCS codes", "$n_{out}$", 
                      "RR below $35 \\mu g m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table2.tex")

# Create results table with ML estimates 
OR_est <- ccs_table(root = "8", moretrees_results = moretrees_results,
                    type = "ml",
                    digits = 3, mult = 10)
OR_est$short_label <- NULL
OR_est$long_label <- NULL
OR_est$n_outcomes <- NULL
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "l", "l"), 
                    digits = 3, display = c("d", "d", "f", "f"))
names(OR_xtable) <- c("Group", 
                      "RR below $35 \\mu g m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table2_ml.tex")

# Create results tree plot 
pdf(file = "./figures/mod1_resp_tree.pdf", width = 6, height = 2.8)
ccs_plot(root = "8", moretrees_results = moretrees_results,
         tr = tr,
         asp = 1/10, leaf.height = 30, label.dist = 4.5)
dev.off()
