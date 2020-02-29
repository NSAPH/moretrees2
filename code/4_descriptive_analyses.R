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

sink()

# Table 1 ------------------------------------------------------------------------------------------------------


# Figure 1 -----------------------------------------------------------------------------------------------------
require(ggplot2)
require(moretrees)
require(igraph)
require(reshape2)

# Get level 2 labels for CVD
ccs_labels <- read.csv("../data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")
tr_cvd <- ccs_tree("7")$tr
lvls_cvd <- ego(tr_cvd, V(tr_cvd)[V(tr_cvd)$leaf], order = 10, mode = "in")
lvls_cvd <- as.data.frame(t(sapply(lvls_cvd, function(x) rev(names(x)))))
names(lvls_cvd) <- paste0("ccs_lvl", 1:4)
lvls_cvd <- merge(lvls_cvd, ccs_labels, by.x = "ccs_lvl2", by.y = "ccs_code",
                  all.x = T, all.y = F)
lvls_cvd$label <- as.character(lvls_cvd$label)
lvls_cvd$label <- paste0(lvls_cvd$label, " (", lvls_cvd$ccs_lvl2, ")")
lvls_cvd$label <- factor(lvls_cvd$label, levels = unique(lvls_cvd$label))
dt_cvd <- merge(dt_cvd, lvls_cvd[ , c("ccs_lvl4", "label")], by.x = "ccs_added_zeros", by.y = "ccs_lvl4",
                all.x = T, all.y = F)

# Reshape
dt_cvd_plot <- dt_cvd[, c("pm25_lag01_case", "pm25_lag01_control", "label")]
dt_cvd_plot <- reshape(dt_cvd_plot, direction = "long", varying = c("pm25_lag01_case", "pm25_lag01_control"),
                       times = c("case", "control"), v.names = "pm25_lag01")

# Plot
cvd_plot <- ggplot(dt_cvd_plot) + geom_density(aes(x = pm25_lag01, fill = time), alpha = 0.2) + 
  facet_wrap(label ~ ., ncol = 1) + theme_minimal() + scale_x_continuous(trans = "log10", limits = c(1, 140))

