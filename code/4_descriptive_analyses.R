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
          "race_gp", "sex_gp", "age_gp", "dual",
          "ozone_lag01_case", "ozone_lag01_control", "zip")
cols <- c(cols.other, cols.omit)

# Select states
# excluded_states <- c(2, 12, 40, 48, 53:63, 97:99)

# Load respiratory data
dt_resp <- read_fst("./data/merged_admissions_enviro/admissions_enviro_resp.fst",
                    as.data.table = T, 
                    columns = cols)
# Exclude states outside contiguous US
# dt_resp <- dt_resp[!(ssa_state_cd %in% excluded_states)]

# Load CVD data
dt_cvd <- read_fst("./data/merged_admissions_enviro/admissions_enviro_cvd.fst",
                   as.data.table = T, 
                   columns = cols)
# # Keep only north east region
# dt_cvd <- dt_cvd[ssa_state_cd %in% states_list]

# Data in text -------------------------------------------------------------------------------------------------
require(janitor)
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

# Respiratory: number of cases per outcome
n_out <- as.integer(table(dt_cvd$ccs_added_zeros))
cat("\n\nNumber of cases per respiratory outcome:\n")
summary(n_out)

# CVD: days above and below threshold
cat("\n\nCase days with PM2.5 > 25:\n")
mean(dt_cvd$pm25_lag01_case > 25)
cat("\n\nControl days with PM2.5 > 25:\n")
mean(dt_cvd$pm25_lag01_control > 25)

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

# Respiratory: number of cases per outcome
n_out <- as.integer(table(dt_resp$ccs_added_zeros))
cat("\n\nNumber of cases per respiratory outcome:\n")
summary(n_out)

# Respiratory: days above and below threshold
cat("\n\nCase days with PM2.5 > 25:\n")
mean(dt_resp$pm25_lag01_case > 25)
cat("\n\nControl days with PM2.5 > 25:\n")
mean(dt_resp$pm25_lag01_control > 25)

# Table 1 data

cat("\n\nTable 1 data------------------------------------------------------\n\n")

cat("\n\nCVD------------------------------------------------------\n\n")

# CVD: summary of PM2.5 exposure by case status
cat("\n\nPM2.5 summary for cases:\n")
summary(dt_cvd$pm25_lag01_case)
cat("\n\nPM2.5 SD for cases:\n")
sd(dt_cvd$pm25_lag01_case)
cat("\n\nPM2.5 summary for controls:\n")
summary(dt_cvd$pm25_lag01_control)
cat("\n\nPM2.5 SD for controls:\n")
sd(dt_cvd$pm25_lag01_control)

# CVD: summary of temperature by case status
cat("\n\nTemperature summary for cases:\n")
summary(dt_cvd$tmmx_lag01_case - 273.15)
cat("\n\nTemperature SD for cases:\n")
sd(dt_cvd$tmmx_lag01_case - 273.15)
cat("\n\nTemperature summary for controls:\n")
summary(dt_cvd$tmmx_lag01_control - 273.15)
cat("\n\nTemperature SD for controls:\n")
sd(dt_cvd$tmmx_lag01_control - 273.15)

# CVD: summary of humidity by case status
cat("\n\nHumidity summary for cases:\n")
summary(dt_cvd$rmax_lag01_case)
cat("\n\nHumidity SD for cases:\n")
sd(dt_cvd$rmax_lag01_case)
cat("\n\nHumidity summary for controls:\n")
summary(dt_cvd$rmax_lag01_control)
cat("\n\nHumidity SD for controls:\n")
sd(dt_cvd$rmax_lag01_control)

# CVD: summary of ozone by case status
cat("\n\nOzone summary for cases:\n")
summary(dt_cvd$ozone_lag01_case)
cat("\n\nOzone SD for cases:\n")
sd(dt_cvd$ozone_lag01_case)
cat("\n\nOzone summary for controls:\n")
summary(dt_cvd$ozone_lag01_control)
cat("\n\nOzone SD for controls:\n")
sd(dt_cvd$ozone_lag01_control)

cat("\n\nRespiratory------------------------------------------------------\n\n")

# Respiratory: summary of PM2.5 exposure by case status
cat("\n\nPM2.5 summary for cases:\n")
summary(dt_resp$pm25_lag01_case)
cat("\n\nPM2.5 SD for cases:\n")
sd(dt_resp$pm25_lag01_case)
cat("\n\nPM2.5 summary for controls:\n")
summary(dt_resp$pm25_lag01_control)
cat("\n\nPM2.5 SD for controls:\n")
sd(dt_resp$pm25_lag01_control)

# resp: summary of temperature by case status
cat("\n\nTemperature summary for cases:\n")
summary(dt_resp$tmmx_lag01_case - 273.15)
cat("\n\nTemperature SD for cases:\n")
sd(dt_resp$tmmx_lag01_case - 273.15)
cat("\n\nTemperature summary for controls:\n")
summary(dt_resp$tmmx_lag01_control - 273.15)
cat("\n\nTemperature SD for controls:\n")
sd(dt_resp$tmmx_lag01_control - 273.15)

# resp: summary of humidity by case status
cat("\n\nHumidity summary for cases:\n")
summary(dt_resp$rmax_lag01_case)
cat("\n\nHumidity SD for cases:\n")
sd(dt_resp$rmax_lag01_case)
cat("\n\nHumidity summary for controls:\n")
summary(dt_resp$rmax_lag01_control)
cat("\n\nHumidity SD for controls:\n")
sd(dt_resp$rmax_lag01_control)

# resp: summary of ozone by case status
cat("\n\nOzone summary for cases:\n")
summary(dt_resp$ozone_lag01_case)
cat("\n\nOzone SD for cases:\n")
sd(dt_resp$ozone_lag01_case)
cat("\n\nOzone summary for controls:\n")
summary(dt_resp$ozone_lag01_control)
cat("\n\nOzone SD for controls:\n")
sd(dt_resp$ozone_lag01_control)

# Table 2 data

cat("\n\nTable 2 data------------------------------------------------------\n\n")

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

# # fake data for testing
# n <- 1000
# dt_cvd <- data.table(pm25_lag01_case = abs(rnorm(n, sd = 0.01)),
#                      pm25_lag01_control = abs(rnorm(n, sd = 0.01)),
#                      ccs_added_zeros = sample(lvls_cvd$lvl4_merge, size = n, replace = T))

# Merge in labels
dt_cvd <- merge(dt_cvd, lvls_cvd, 
                by.x = "ccs_added_zeros", by.y = "lvl4_merge",
                all.x = T, all.y = F, sort = FALSE)
dt_cvd <- dt_cvd[order(match(ccs_added_zeros, lvls_cvd$lvl4_merge))]

# Get plot data
dt_cvd_plot <- dt_plot_fun(dt_cvd, plot_depth = 3)

# Make plot
xlab <- expression(mu*"g"*m^-3)
lab.widths <- c(0.9, 0.9, 0.9)
lab.txt.width <- c(20, 20, 25)
pdf(file = "./figures/cvd_nested_plot.pdf", height = 8, width = 10)
nested_plots(dt_cvd_plot, xlab = xlab, lab.widths = lab.widths,
             lab.txt.width = lab.txt.width, axis.height = 1.5,
             lab.txt.size = 4, digits = 2, axis.txt.size = 10,
             plot_depth = 3)
dev.off()

# Get plot data
dt_cvd_plot2 <- dt_plot_fun(dt_cvd, plot_depth = 4)

# Make plot
lab.widths <- c(0.75,0.9, 0.8, 0.82)
lab.txt.width <- c(20, 25, 25, 15)
pdf(file = "./figures/cvd_nested_plot2.pdf", height = 15, width = 15)
nested_plots(dt_cvd_plot2, xlab = xlab, lab.widths = lab.widths,
             lab.txt.width = lab.txt.width, axis.height = 1.65,
             cil_min = -0.015, ciu_max = 0.025,
             lab.txt.size = 4, digits = 2, axis.txt.size = 10,
             plot_depth = 4)
dev.off()

# Get labels for RD
lvls_resp <- get_labels('8')

# # fake data for testing
# n <- 1000
# dt_resp <- data.table(pm25_lag01_case = abs(rnorm(n, sd = 0.01)),
#                      pm25_lag01_control = abs(rnorm(n, sd = 0.01)),
#                      ccs_added_zeros = sample(lvls_resp$lvl4_merge, size = n, replace = T))

# Merge in labels
dt_resp <- merge(dt_resp, lvls_resp, 
                by.x = "ccs_added_zeros", by.y = "lvl4_merge",
                all.x = T, all.y = F, sort = F)
dt_resp <- dt_resp[order(match(ccs_added_zeros, lvls_resp$lvl4_merge))]

# Get plot data
dt_resp_plot <- dt_plot_fun(dt_resp, plot_depth = 3)

# Plot difference in PM25 between case and control by disease ----------------------------------
xlab <- expression(mu*"g"*m^-3)
lab.widths <- c(0.9, 3.3, 1.2)
lab.txt.width <- c(17, 50, 25)
pdf(file = "./figures/resp_nested_plot.pdf", height = 8, width = 10)
nested_plots(dt_resp_plot, xlab = xlab, lab.widths = lab.widths,
             lab.txt.width = lab.txt.width, axis.height = 1.35,
             lab.txt.size = 4, digits = 2, axis.txt.size = 10,
             plot_depth = 3)
dev.off()

# Get plot data
dt_resp_plot2 <- dt_plot_fun(dt_resp, plot_depth = 4)

# Make plot
lab.widths <- c(0.9, 2.5, 0.9, 0.9)
lab.txt.width <- c(20, 60, 25, 25)
pdf(file = "./figures/resp_nested_plot2.pdf", height = 15, width = 15)
nested_plots(dt_resp_plot2, xlab = xlab, lab.widths = lab.widths,
             lab.txt.width = lab.txt.width, axis.height = 1,
             lab.txt.size = 4, digits = 2, axis.txt.size = 10,
             cil_min = -0.015, ciu_max = 0.025,
             plot_depth = 4)
dev.off()


