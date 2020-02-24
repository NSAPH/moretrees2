# Read in results of moretrees model
load(file = "./results/attempt1.Rdata")
moretrees_results_init <- moretrees_results
load(file = "./results/attempt1_run2.Rdata")

# Plot ELBO
ELBO <- moretrees_results$mod$ELBO_track
plot.start <- 7600
plot.end <- length(ELBO)
# plot.end <- 8000
plot(plot.start:plot.end, ELBO[plot.start:plot.end], type = "l")
ELBO[plot.end] - ELBO[plot.end - 1]

# Create results table
OR_est <- exp(moretrees_results$beta_moretrees[ , 2:4] * 10)
OR_est$outcomes <- moretrees_results$beta_moretrees$outcomes
View(OR_est)
