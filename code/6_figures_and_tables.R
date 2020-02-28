# Read in results of moretrees model -----------------------------------------
load(file = "./results/attempt4_northEast.Rdata")
sd_pm25 <- 7.48621
sd_tmmx <- 6.029194
sd_rmax <- 12.50401
n <- 3354364

# Create results summary -----------------------------------------------------
OR_est <- exp(moretrees_results$beta_moretrees[ , 2:4] / sd_pm25 * 10)
OR_est$outcomes <- moretrees_results$beta_moretrees$outcomes

# Map CCS codes to diseases --------------------------------------------------
ccs_lvls <- read.csv("./data/Multi_Level_CCS_2015_cleaned/dxmlabel-13.csv")
require(stringr)
ccs_lvls$label <- str_remove(ccs_lvls$label, "\\s\\[.*")
ccs_lvls$label <- str_trim(ccs_lvls$label, side = "right")
ccs_lvls$ccs_code <- str_trim(ccs_lvls$ccs_code, side = "right")
write.csv(ccs_lvls, "./data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv", row.names = F)
require(moretrees)
tr <- ccs_tree("7")
ccs_icd9_mapping <- tr$ccs_icd_mapping
tr <- tr$tr
ccs_zeros <- ccs_icd9_mapping[, c("ccs_original", "ccs_added_zeros")]
ccs_zeros <- ccs_zeros[!duplicated(ccs_zeros), ]
ccs <- merge(ccs_zeros, ccs_lvls, by.x = "ccs_original", by.y = "ccs_code")

# Map outcome groups to disease names -----------------------------------------
digits <- 3
frmt <- paste0("%.",digits,"f")
OR_est$long_label <- character(length = nrow(OR_est))
OR_est$short_label <- character(length = nrow(OR_est))
OR_est$est_ci <- character(length = nrow(OR_est))
firstlower <- function(x) {
  substr(x, 1, 1) <- tolower(substr(x, 1, 1))
  x
}
collapse_labels <- function(ccs_g) {
  ccs_labels <- mapply(function(lab, code) paste0(lab, " (",code,")"),
                       ccs_g$label, ccs_g$ccs_original, USE.NAMES = F)
  if (length(ccs_labels) > 1) {
    ccs_labels[2:length(ccs_labels)] <- sapply(ccs_labels[2:length(ccs_labels)],
                                               firstlower)
  }
  ccs_labels <- paste0(ccs_labels, collapse = "; ")
  return(ccs_labels)
}
require(igraph)
for (g in 1:nrow(OR_est)) {
  ccs_g <- ccs[ccs$ccs_added_zeros %in% OR_est$outcomes[[g]], ]
  OR_est$long_label[g] <- collapse_labels(ccs_g)
  # Collapse siblings
  i <- 1
  while (i <= nrow(ccs_g)) {
    code <- ccs_g$ccs_added_zeros[i]
    # Get siblings from tree
    parent <- names(ego(tr, order = 1, mindist = 1,
                  nodes = code, mode = "in")[[1]])
    sibs <- names(ego(tr, order = 1, mindist = 1,
                nodes = parent, mode = "out")[[1]])
    if (sum(sibs %in% ccs_g$ccs_added_zeros) == length(sibs)) {
      ccs_g$ccs_added_zeros[i] <- parent
      ccs_g$ccs_original[i] <- str_remove_all(parent, "\\.0")
      ccs_g$label[i] <- ccs_lvls$label[ccs_lvls$ccs_code == parent]
      ccs_g <- ccs_g[!(ccs_g$ccs_added_zeros %in% setdiff(sibs, code)), ]
      i <- 1
    } else {
      i <- i + 1
    }
  }
  # Collapse simplified labels into one string
  OR_est$short_label[g] <- collapse_labels(ccs_g)
  # Collapse OR and CI into one string
  est1_frmt <- sprintf(frmt, OR_est$est1[g])
  cil1_frmt <- sprintf(frmt, OR_est$cil1[g])
  ciu1_frmt <- sprintf(frmt, OR_est$ciu1[g])
  OR_est$est_ci1[g] <- paste0(est1_frmt," (",cil1_frmt,", ",ciu1_frmt,")")
}
OR_est$n_outcomes <- sapply(OR_est$outcomes, length)
OR_est$group <- 1:nrow(OR_est)

# Create results table ------------------------------------------------------
OR_tab <- OR_est[ , c("group", "short_label", "n_outcomes", "est_ci1")]
require(xtable)
row.names(OR_est) <- NULL
OR_xtable <- xtable(OR_tab,align=c("l","l","p{8cm}","c","c"),digits=3,
                       display=c("d","d","s","d","f"))
names(OR_xtable) <- c("Group","CCS codes","#Outcomes","OR (95%CI)")

write(print(OR_xtable,floating=FALSE,include.rownames = FALSE),file="./figures/table1.tex")




