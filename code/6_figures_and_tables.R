# Read in results of moretrees model
load(file = "./results/attempt1.Rdata")

# Plot ELBO
ELBO <- moretrees_results$mod$ELBO_track
plot.start <- 30
plot.end <- length(ELBO)
plot(plot.start:plot.end, ELBO[plot.start:plot.end], type = "l")

# Create results table
OR_est <- exp(moretrees_results$beta_moretrees[ , 2:4] * 10)
OR_est$outcomes <- moretrees_results$beta_moretrees$outcomes
View(OR_est)
