# ----------------------------------------------------------------------------
# ----------------------------- CVD results ----------------------------------
# ----------------------------------------------------------------------------

# collecting approximate marginal likelihoods
bayes_factor_cvd <- numeric(3)

# CVD mod1 -------------------------------------------------------------------
source("./code/results_functions.R")

# Read in results of moretrees model 
load(file = "./results/mod1_split35_northEast_cvd.Rdata")

# collect bayes factor
bayes_factor_cvd[1] <- moretrees_results$mod$hyperparams$ELBO

# Create results table 
OR_est <- ccs_table(root = "7", moretrees_results = moretrees_results,
                    digits = 3, mult = 10)
OR_est$long_label <- NULL
OR_est$n_obs <- formatC(OR_est$n_obs, format="d", big.mark=",")
require(xtable)
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "p{7.2cm}", "r", "r", "p{2.2cm}", "p{2.2cm}"), 
                    digits = 3, display = c("d", "d", "s", "d", "f", "f", "f"))
names(OR_xtable) <- c("Group", "CCS codes", 
                      "$n_{out}$", "$n_{obs}$",
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_cvd_mod1.tex")

# Create results table with ML estimates 
OR_est <- ccs_table(root = "7", moretrees_results = moretrees_results,
                    type = "ml",
                    digits = 3, mult = 10)
OR_est$short_label <- NULL
OR_est$long_label <- NULL
OR_est$n_outcomes <- NULL
OR_est$n_obs <- NULL
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "l", "l"), 
                    digits = 3, display = c("d", "d", "f", "f"))
names(OR_xtable) <- c("Group", 
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_cvd_mod1_ml.tex")

# Create results tree plot 
pdf(file = "./figures/mod1_cvd_tree.pdf", width = 6, height = 2.8)
ccs_plot(root = "7", moretrees_results = moretrees_results,
         asp = 1/10, leaf.height = 30, label.dist = 4.5)
dev.off()

# CVD mod2 -------------------------------------------------------------------
source("./code/results_functions.R")

# Read in results of moretrees model 
load(file = "./results/mod2_split35_northEast_cvd_attempt2.RData")

# collect bayes factor
bayes_factor_cvd[2] <- moretrees_results$mod$hyperparams$ELBO

# Create results table 
OR_est <- ccs_table(root = "7", moretrees_results = moretrees_results,
                    digits = 3, mult = 10)
OR_est$long_label <- NULL
OR_est$n_obs <- formatC(OR_est$n_obs, format="d", big.mark=",")
require(xtable)
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "p{7.2cm}", "r", "r", "p{2.2cm}", "p{2.2cm}"), 
                    digits = 3, display = c("d", "d", "s", "d", "f", "f", "f"))
names(OR_xtable) <- c("Group", "CCS codes", 
                      "$n_{out}$", "$n_{obs}$",
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_cvd_mod2.tex")

# Create results table with ML estimates 
OR_est <- ccs_table(root = "7", moretrees_results = moretrees_results,
                    type = "ml",
                    digits = 3, mult = 10)
OR_est$short_label <- NULL
OR_est$long_label <- NULL
OR_est$n_outcomes <- NULL
OR_est$n_obs <- NULL
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "l", "l"), 
                    digits = 3, display = c("d", "d", "f", "f"))
names(OR_xtable) <- c("Group", 
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_cvd_mod2_ml.tex")

# Create results tree plot 
pdf(file = "./figures/mod2_cvd_tree.pdf", width = 6, height = 2.8)
ccs_plot(root = "7", moretrees_results = moretrees_results,
         asp = 1/10, leaf.height = 30, label.dist = 4.5)
dev.off()

# CVD mod3 -------------------------------------------------------------------
source("./code/results_functions.R")

# Read in results of moretrees model 
load(file = "./results/mod3_split35_northEast_cvd.RData")

# collect bayes factor
bayes_factor_cvd[3] <- moretrees_results$mod$hyperparams$ELBO

# Create results table 
OR_est <- ccs_table(root = "7", moretrees_results = moretrees_results,
                    digits = 3, mult = 10)
OR_est$long_label <- NULL
OR_est$n_obs <- formatC(OR_est$n_obs, format="d", big.mark=",")
require(xtable)
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "p{7.2cm}", "r", "r", "p{2.2cm}", "p{2.2cm}"), 
                    digits = 3, display = c("d", "d", "s", "d", "f", "f", "f"))
names(OR_xtable) <- c("Group", "CCS codes", 
                      "$n_{out}$", "$n_{obs}$",
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_cvd_mod3.tex")

# Create results table with ML estimates 
OR_est <- ccs_table(root = "7", moretrees_results = moretrees_results,
                    type = "ml",
                    digits = 3, mult = 10)
OR_est$short_label <- NULL
OR_est$long_label <- NULL
OR_est$n_outcomes <- NULL
OR_est$n_obs <- NULL
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "l", "l"), 
                    digits = 3, display = c("d", "d", "f", "f"))
names(OR_xtable) <- c("Group", 
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_cvd_mod3_ml.tex")

# Create results tree plot 
pdf(file = "./figures/mod3_cvd_tree.pdf", width = 6, height = 2.8)
ccs_plot(root = "7", moretrees_results = moretrees_results,
         asp = 1/10, leaf.height = 30, label.dist = 4.5)
dev.off()

# ----------------------------------------------------------------------------
# ------------------------ Respiratory results -------------------------------
# ----------------------------------------------------------------------------

# collecting approximate marginal likelihoods
rm(list = ls())
bayes_factor_resp <- numeric(3)

# Respiratory disease mod1 ------------------------------------------------------------
source("./code/results_functions.R")

# Read in results of moretrees model 
load(file = "./results/mod1_split35_northEast_resp_run2.Rdata")

# collect bayes factor
bayes_factor_resp[1] <- moretrees_results$mod$hyperparams$ELBO

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
OR_est$n_obs <- formatC(OR_est$n_obs, format="d", big.mark=",")
OR_xtable <- xtable(OR_est,align = c("l", "l", "p{7.6cm}", "r", "r", "p{2.2cm}", "p{2.2cm}"), 
                    digits = 3, display = c("d", "d", "s", "d", "d", "f", "f"))
names(OR_xtable) <- c("Group", "CCS codes", 
                      "$n_{out}$", "$n_{obs}$",
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_resp_mod1.tex")

# Create results table with ML estimates 
OR_est <- ccs_table(root = "8", moretrees_results = moretrees_results,
                    type = "ml",
                    digits = 3, mult = 10)
OR_est$short_label <- NULL
OR_est$long_label <- NULL
OR_est$n_outcomes <- NULL
OR_est$n_obs <- NULL
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "l", "l"), 
                    digits = 3, display = c("d", "d", "f", "f"))
names(OR_xtable) <- c("Group", 
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_resp_mod1_ml.tex")

# Create results tree plot 
pdf(file = "./figures/mod1_resp_tree.pdf", width = 6, height = 2.8)
ccs_plot(root = "8", moretrees_results = moretrees_results,
         tr = tr,
         asp = 1/10, leaf.height = 30, label.dist = 4.5)
dev.off()

# Respiratory disease mod2 ------------------------------------------------------------
source("./code/results_functions.R")

# Read in results of moretrees model 
load(file = "./results/mod2_split35_northEast_resp.Rdata")

# collect bayes factor
bayes_factor_resp[2] <- moretrees_results$mod$hyperparams$ELBO

# Create results table 
OR_est <- ccs_table(root = "8", moretrees_results = moretrees_results,
                    tr = tr, digits = 3, mult = 10)
OR_est$long_label <- NULL
require(xtable)
row.names(OR_est) <- NULL
OR_est$n_obs <- formatC(OR_est$n_obs, format="d", big.mark=",")
OR_xtable <- xtable(OR_est,align = c("l", "l", "p{7.6cm}", "r", "r", "p{2.2cm}", "p{2.2cm}"), 
                    digits = 3, display = c("d", "d", "s", "d", "d", "f", "f"))
names(OR_xtable) <- c("Group", "CCS codes", 
                      "$n_{out}$", "$n_{obs}$",
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_resp_mod2.tex")

# Create results table with ML estimates 
OR_est <- ccs_table(root = "8", moretrees_results = moretrees_results,
                    type = "ml",
                    digits = 3, mult = 10)
OR_est$short_label <- NULL
OR_est$long_label <- NULL
OR_est$n_outcomes <- NULL
OR_est$n_obs <- NULL
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "l", "l"), 
                    digits = 3, display = c("d", "d", "f", "f"))
names(OR_xtable) <- c("Group", 
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_resp_mod2_ml.tex")

# Create results tree plot 
pdf(file = "./figures/mod2_resp_tree.pdf", width = 6, height = 2.8)
ccs_plot(root = "8", moretrees_results = moretrees_results,
         tr = tr,
         asp = 1/10, leaf.height = 30, label.dist = 4.5)
dev.off()

# Respiratory disease mod3 ------------------------------------------------------------
source("./code/results_functions.R")

# Read in results of moretrees model 
load(file = "./results/mod3_split35_northEast_resp.Rdata")

# collect bayes factor
bayes_factor_resp[3] <- moretrees_results$mod$hyperparams$ELBO

# Create results table 
OR_est <- ccs_table(root = "8", moretrees_results = moretrees_results,
                    tr = tr, digits = 3, mult = 10)
OR_est$long_label <- NULL
require(xtable)
row.names(OR_est) <- NULL
OR_est$n_obs <- formatC(OR_est$n_obs, format="d", big.mark=",")
OR_xtable <- xtable(OR_est,align = c("l", "l", "p{7.6cm}", "r", "r", "p{2.2cm}", "p{2.2cm}"), 
                    digits = 3, display = c("d", "d", "s", "d", "d", "f", "f"))
names(OR_xtable) <- c("Group", "CCS codes", 
                      "$n_{out}$", "$n_{obs}$",
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_resp_mod3.tex")

# Create results table with ML estimates 
OR_est <- ccs_table(root = "8", moretrees_results = moretrees_results,
                    type = "ml",
                    digits = 3, mult = 10)
OR_est$short_label <- NULL
OR_est$long_label <- NULL
OR_est$n_outcomes <- NULL
OR_est$n_obs <- NULL
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_est, align = c("l", "l", "l", "l"), 
                    digits = 3, display = c("d", "d", "f", "f"))
names(OR_xtable) <- c("Group", 
                      "RR below $35 \\mu g \\cdot m^{-3}$ (95\\%CI)",
                      "RR above $35 \\mu g \\cdot m^{-3}$ (95\\%CI)")

write(print(OR_xtable, floating = FALSE, include.rownames = FALSE,
            sanitize.text.function = function(x) x),
      file = "./figures/table_resp_mod3_ml.tex")

# Create results tree plot 
pdf(file = "./figures/mod3_resp_tree.pdf", width = 6, height = 2.8)
ccs_plot(root = "8", moretrees_results = moretrees_results,
         tr = tr,
         asp = 1/10, leaf.height = 30, label.dist = 4.5)
dev.off()

