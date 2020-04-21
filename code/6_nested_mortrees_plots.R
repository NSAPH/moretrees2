
source("./code/results_functions.R")
require(fst)
require(splines)
dataset <- c("cvd", "resp")
root <- c("7", "8")
mod <- 3
split <- "0"

nknots <- 2
q <- seq(0, 1, length.out = nknots + 2)

for (i in 1:length(dataset)) { # datasets
  ds <- dataset[i]

  load(paste0("./results/nest_pltdat_", ds, ".RData"))
  pltdat_ml[ , c("est1", "cil1", "ciu1")] <- 100 * (exp(pltdat_ml[, c("est1", "cil1", "ciu1")] * 10 ) - 1)

  # plot
  xlab <- "Excess Rate (%)"
  if (i == 1) {
    lab.widths <- c(0.75, 0.75, 0.75, 0.85)
    lab.txt.width <- c(20, 25, 25, 15)
    pdf(file = "./figures/cvd_nested_clr.pdf", height = 15, width = 14)
    beta_indiv_plot_fun(pltdat_ml, tr, xlab = xlab, lab.widths = lab.widths,
                        lab.txt.width = lab.txt.width, axis.height = 1.5,
                        cil_min = -8, ciu_max = 8,
                        force_lims = TRUE,
                        lab.txt.size = 4, digits = 0, axis.txt.size = 10,
                        plot_depth = 4, wrap_labs = FALSE)
    dev.off()
  }
  if (i == 2) {
    lab.widths <- c(0.75, 0.75, 0.75, 0.75)
    lab.txt.width <- c(20, 25, 25, 15)
    pdf(file = "./figures/resp_nested_clr.pdf", height = 10, width = 14)
    beta_indiv_plot_fun(pltdat_ml, tr, xlab = xlab, lab.widths = lab.widths,
                        lab.txt.width = lab.txt.width, axis.height = 1.3,
                        cil_min = -8, ciu_max = 8,
                        force_lims = TRUE,
                        lab.txt.size = 4, digits = 0, axis.txt.size = 10,
                        plot_depth = 4, wrap_labs = FALSE)
    dev.off()
  }
}
