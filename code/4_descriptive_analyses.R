# Make sure working directory is set to moretrees2/code
setwd("/nfs/home/E/ethomas/shared_space/ci3_analysis/moretrees2/")
source("./code/results_functions.R")

# Check for updates on moretrees master branch
# devtools::install_github("emgthomas/moretrees_pkg", ref = "devel")
require(moretrees)
# note: for some updates, may have to restart R session
require(fst)

# Relevant columns
cols.omit <- c("ccs_added_zeros", "pm25_lag01_case", "pm25_lag01_control",
               "tmmx_lag01_case", "tmmx_lag01_control",
               "rmax_lag01_case", "rmax_lag01_control")
cols.other <- c("id", "adate", "ssa_state_cd",
          "race_gp", "sex_gp", "age_gp", "dual")
cols <- c(cols.other, cols.omit)

# Select states
states_list <- c(7, 20, 22, 30, 41,
                 47, 31, 33, 39)
# states_list <- c("CT", "ME", "MA", "NH", "RI", 
#                 "VT", "NJ", "NY", "PA")

# Load respiratory data
dt_resp <- read_fst("./data/merged_admissions_enviro/admissions_enviro_resp.fst",
                    as.data.table = T, 
                    columns = cols)
# Keep only north east region
dt_resp <- dt_resp[ssa_state_cd %in% states_list]

# Load CVD data
dt_cvd <- read_fst("./data/merged_admissions_enviro/admissions_enviro_cvd.fst",
                   as.data.table = T, 
                   columns = cols)
# Keep only north east region
dt_cvd <- dt_cvd[ssa_state_cd %in% states_list]

# Data in text -------------------------------------------------------------------------------------------------

sink(file = "./results/descriptives.txt")

# CVD ---------------------------------------------------------
cat("Total CVD admissions in study period = ", nrow(dt_cvd))

# CVD: first admission only
dt_cvd <- dt_cvd[order(id, adate)]
dt_cvd <- dt_cvd[ , .SD[1], by = id]
n_cvd <- nrow(dt_cvd)
cat("\n\nTotal first admissions for CVD = ", n_cvd)

# CVD: keeping only complete cases
dt_cvd <- na.omit(dt_cvd, cols = cols.omit)
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
dt_resp <- na.omit(dt_resp, cols = cols.omit)
cat("\n\nComplete cases for respiratory disease = ", nrow(dt_resp), "(", 100 * nrow(dt_resp) / n_resp, "%)")

# Respiratory: median PM2.5 exposure by case status
cat("\n\nPM2.5 summary for cases:\n")
summary(dt_resp$pm25_lag01_case)
cat("\n\nPM2.5 summary for controls:\n")
summary(dt_resp$pm25_lag01_control)

# Table 1 data

cat("\n\nTable 1 data------------------------------------------------------\n\n")

cat("\n\nCVD:\n")

cat("\n\nSex:\n")
tabyl(dt_cvd$sex_gp)

cat("\n\nAge:\n")
tabyl(dt_cvd$age_gp)

cat("\n\nRace:\n")
tabyl(dt_cvd$race_gp)

cat("\n\nMedicaid:\n")
tabyl(dt_cvd$dual)

cat("\n\nRespiratory disease:\n")

cat("\n\nSex:\n")
tabyl(dt_resp$sex_gp)

cat("\n\nAge:\n")
tabyl(dt_resp$age_gp)

cat("\n\nRace:\n")
tabyl(dt_resp$race_gp)

cat("\n\nMedicaid:\n")
tabyl(dt_resp$dual)

sink()

# Figure 1 -----------------------------------------------------------------------------------------------------

# Get labels for CVD
lvls_cvd <- get_labels('7')
plot_depth <- 4

# fake data for testing
n <- 1000
dt_cvd <- data.table(pm25_lag01_case = abs(rnorm(n, sd = 0.01)),
                     pm25_lag01_control = abs(rnorm(n, sd = 0.01)),
                     ccs_added_zeros = sample(lvls_cvd$lvl4_merge, size = n, replace = T))

# Merge in labels
dt_cvd <- merge(dt_cvd, lvls_cvd, 
                by.x = "ccs_added_zeros", by.y = "lvl4_merge",
                all.x = T, all.y = F, sort = FALSE)
dt_cvd <- dt_cvd[order(match(ccs_added_zeros, lvls_cvd$lvl4_merge))]

# Get plot data
dt_cvd_plot <- dt_plot_fun(dt_cvd, plot_depth = plot_depth)

# Make plot
xlab <- expression(mu*"g"*m^-3)
pdf(file = "./figures/cvd_nested_plot.pdf", height = 8, width = 6)
nested_plots(dt_cvd_plot, xlab = xlab, lab.nudge = 2,
             lab.size = 3, digits = 2, axis.txt.size = 5,
             plot_depth = plot_depth)
dev.off()

# Get labels for RD
lvls_resp <- get_labels('8')

# Merge in labels
dt_resp <- merge(dt_resp, lvls_resp, 
                by.x = "ccs_added_zeros", by.y = "lvl4_merge",
                all.x = T, all.y = F, sort = F)
dt_resp <- dt_resp[order(match(ccs_added_zeros, lvls_resp$lvl4_merge))]

# Get plot data
dt_resp_plot <- dt_plot_fun(dt_resp, lab_var = "ccs_lvl")

# Plot difference in PM25 between case and control by disease ----------------------------------
pdf(file = "./figures/resp_nested_plot.pdf", height = 8, width = 6)
nested_plots(dt_resp_plot, xlab = xlab, lab.nudge = 0.15, lab.size = 3, digits = 2, axis.txt.size = 7)
dev.off()



