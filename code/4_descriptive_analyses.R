# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/code/")

# Check for updates on moretrees master branch
# devtools::install_github("emgthomas/moretrees_pkg", ref = "devel")
require(moretrees)
# note: for some updates, may have to restart R session
require(fst)

# Relevant columns
cols <- c("id", "adate", "ssa_state_cd",
          "ccs_added_zeros", "pm25_lag01_case", "pm25_lag01_control")

# Select states
states_list <- c(7, 20, 22, 30, 41,
                 47, 31, 33, 39)
# states_list <- c("CT", "ME", "MA", "NH", "RI", 
#                 "VT", "NJ", "NY", "PA")

# Load respiratory data
dt_resp <- read_fst("../data/merged_admissions_enviro/admissions_enviro_resp.fst",
               as.data.table = T, 
               columns = cols)
# Keep only north east region
dt_resp <- dt_resp[ssa_state_cd %in% states_list]

# Load CVD data
dt_cvd <- read_fst("../data/merged_admissions_enviro/admissions_enviro_cvd.fst",
                    as.data.table = T, 
                    columns = cols)
# Keep only north east region
dt_cvd <- dt_cvd[ssa_state_cd %in% states_list]

# Data in text -------------------------------------------------------------------------------------------------

sink(file = "../results/descriptives.txt")

# CVD ---------------------------------------------------------
cat("Total CVD admissions in study period = ", nrow(dt_cvd))

# CVD: first admission only
dt_cvd <- dt_cvd[order(id, adate)]
dt_cvd <- dt_cvd[ , .SD[1], by = id]
n_cvd <- nrow(dt_cvd)
cat("\n\nTotal first admissions for CVD = ", n_cvd)

# CVD: keeping only complete cases
dt_cvd <- na.omit(dt_cvd)
cat("\n\nComplete cases for CVD = ", nrow(dt_cvd), "(", 100 * nrow(dt_cvd) / n_cvd, "%)")

# CVD: median PM2.5 exposure by case status
cat("\n\nPM2.5 summary for cases:\n")
summary(dt_cvd$pm25_lag01_case)
cat("\n\nPM2.5 summary for controls:\n")
summary(dt_cvd$pm25_lag01_control)

# Respiratory disease ------------------------------------------
cat("\n\nTotal respiratory admissions in study period = ", nrow(dt_resp))

# Respiratory: first admission only
dt_resp <- dt_resp[order(id, adate)]
dt_resp <- dt_resp[ , .SD[1], by = id]
n_resp <- nrow(dt_resp)
cat("\n\nTotal first admissions for respiratory disease = ", n_resp)

# Respiratory: keeping only complete cases
dt_resp <- na.omit(dt_resp)
cat("\n\nComplete cases for respiratory disease = ", nrow(dt_resp), "(", 100 * nrow(dt_resp) / n_resp, "%)")

# Respiratory: median PM2.5 exposure by case status
cat("\n\nPM2.5 summary for cases:\n")
summary(dt_resp$pm25_lag01_case)
cat("\n\nPM2.5 summary for controls:\n")
summary(dt_resp$pm25_lag01_control)

sink()

# Table 1 ------------------------------------------------------------------------------------------------------


# Figure 1 -----------------------------------------------------------------------------------------------------
require(ggplot2)
require(moretrees)
require(igraph)
require(reshape2)
require(stringr)

# Get labels for CVD
ccs_labels <- read.csv("../data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")
ccs_labels$label <- paste0(ccs_labels$label, " (", ccs_labels$ccs, ")")
tr_cvd <- ccs_tree("7")$tr
lvls_cvd <- ego(tr_cvd, V(tr_cvd)[V(tr_cvd)$leaf], order = 10, mode = "in")
lvls_cvd <- as.data.frame(t(sapply(lvls_cvd, function(x) rev(names(x)))))
names(lvls_cvd) <- paste0("ccs_lvl", 1:4)
lvls_cvd$lvl4_merge <- lvls_cvd$ccs_lvl4
for (i in 1:4) {
  col_i <- paste0("ccs_lvl", i)
  lvls_cvd[ , col_i]  <- str_remove_all(lvls_cvd[ , col_i], "\\.0")
  names(ccs_labels)[2] <- paste0("label", i)
  lvls_cvd <- merge(lvls_cvd, ccs_labels, by.x = col_i, by.y = "ccs_code",
                    all.x = T, all.y = F, sort = F)
  
}
dt_cvd <- merge(dt_cvd, lvls_cvd[ , c("lvl4_merge", paste0("label", 1:4))], 
                by.x = "ccs_added_zeros", by.y = "lvl4_merge",
                all.x = T, all.y = F)

# Get mean difference in PM25 between case and control by disease
mean.diff.log <- function(x , y) {
  ttest <- t.test(log(x), log(y), paired = TRUE)
  list(est = ttest$estimate, cil = ttest$conf.int[1], ciu = ttest$conf.int[2], n = length(x))
}
for (i in 1:3) {
  dt_cvd[ , paste0(c("est_lvl", "cil_lvl", "ciu_lvl", "n"), i) := mean.diff.log(pm25_lag01_case, pm25_lag01_control),
         by = get(paste0("label", i))]
}
cvd_plot <- unique(dt_cvd[ , paste0(c("label", "est_lvl", "cil_lvl", "ciu_lvl", "n"), rep(1:3, each = 5))])
dt_cvd[ , paste0(c("est_lvl", "cil_lvl", "ciu_lvl"), rep(1:3, each = 3)) := NULL]

# Plot PM25 histograms ---------------------------------------------------
cvd_hist <- ggplot(dt_cvd) + geom_density(aes(x = pm25_lag01_case)) + 
  facet_wrap(label2 ~ ., ncol = 1) + theme_minimal() +
  scale_x_continuous(trans = "log10", lim = c(1, max(dt_cvd$pm25_lag01_case)))

# Nested difference plots -------------------------------------------------------------
require(gridExtra)
plot.count <- 0
grobs <- list()
lims <- c(min(cvd_plot[ , paste0("cil_lvl", 1:3)]),
          max(cvd_plot[ , paste0("ciu_lvl", 1:3)]))
layout <- integer()
for (i in 1:3) {
  lab <- paste0("label", i)
  xmin <- paste0("cil_lvl", i)
  xmax <- paste0("ciu_lvl", i)
  x <- paste0("est_lvl", i)
  for (disease in unique(cvd_plot[ , get(lab)])) {
    dat <- as.data.frame(cvd_plot[get(lab) == disease])
    dat <- dat[ , c(x, xmin, xmax)]
    dat <- unique(dat)
    grob <- ggplot(dat) + 
      geom_point(aes(x = get(x), y = 1)) +
      geom_errorbarh(aes(xmin = get(xmin), xmax = get(xmax), y = 1)) +
      ggtitle(disease) + theme_minimal() + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
      scale_x_continuous(limits = lims)
    plot.count <- plot.count + 1
    grobs[[plot.count]] <- grob
    layout <- c(layout, rep(plot.count, times = sum(cvd_plot[ , get(lab)] == disease)))
  }
}





