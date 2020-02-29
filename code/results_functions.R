require(moretrees)
require(igraph)
require(RColorBrewer)

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

ccs_table <- function(root, moretrees_results, digits = 3, mult = 10) {
  
  ccs_lvls <- read.csv("./data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")
  
  # Get ORs
  k <- ncol(moretrees_results$mod$vi_params$mu[[1]])
  cols <- unlist(lapply(1:k, function(k) paste0(c("est", "cil", "ciu"), k)))
  OR_est <- exp(moretrees_results$beta_moretrees[ , cols] * mult)
  OR_est$outcomes <- moretrees_results$beta_moretrees$outcomes
  
  # Map CCS codes to diseases --------------------------------------------------
  tr <- ccs_tree(root)
  ccs_icd9_mapping <- tr$ccs_icd_mapping
  tr <- tr$tr
  ccs_zeros <- ccs_icd9_mapping[, c("ccs_original", "ccs_added_zeros")]
  ccs_zeros <- ccs_zeros[!duplicated(ccs_zeros), ]
  ccs <- merge(ccs_zeros, ccs_lvls, by.x = "ccs_original", by.y = "ccs_code")
  
  # Map outcome groups to disease names -----------------------------------------
  frmt <- paste0("%.", digits, "f")
  OR_est$long_label <- character(length = nrow(OR_est))
  OR_est$short_label <- character(length = nrow(OR_est))
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
    for (j in 1:k) {
      est_frmt <- sprintf(frmt, OR_est$est1[g])
      cil_frmt <- sprintf(frmt, OR_est$cil1[g])
      ciu_frmt <- sprintf(frmt, OR_est$ciu1[g])
      OR_est[g, paste0("ci_est", j)] <- paste0(est_frmt," (",cil_frmt,", ",ciu_frmt,")")
    }
  }
  OR_est$n_outcomes <- sapply(OR_est$outcomes, length)
  OR_est$root <- 1:nrow(OR_est)
  
  return(OR_est[ , c("root", "short_label", "long_label", "n_outcomes", paste0("ci_est", 1:k))])
}
  
ccs_plot <- function(root, moretrees_results,
                     asp = 1/5, internal.size = 1,
                     leaf.size = 2, ...) {
  
  # Get tree
  tr <- ccs_tree(root)$tr
  
  # Assign groups to tree
  outcomes <- moretrees_results$beta_moretrees$outcomes
  cols_g <- brewer.pal(length(outcomes), "Set1")
  cols <- character(length(V(tr)))
  V(tr)$group <- 0
  for (g in 1:length(outcomes)) {
    V(tr)[outcomes[[g]]]$group <- g
    cols[names(V(tr)) %in% outcomes[[g]]] <- cols_g[g]
  }
  groups <- V(tr)$group
  cols[groups == 0] <- "grey"
  groups[groups == 0] <- ""
  vsize <- rep(internal.size, length(V(tr)))
  vsize[V(tr)$leaf] <- leaf.size
  
  # Create legend
  cols2 <- data.frame(cols, groups)
  cols2 <- cols2[!duplicated(cols2), ]
  cols2 <- cols2[cols2$cols != "grey", ]
  cols2 <- cols2[order(cols2$groups) , ]
  cols2$names <- paste0("Group ", cols2$groups)
  
  # layout
  l <- layout.reingold.tilford(tr)
  l <- layout_as_tree(tr)
  
  # Plot
  plot.igraph(tr, layout = l, vertex.label = NA,
              edge.arrow.size = 0,
              vertex.size = vsize, vertex.color = cols, 
              vertex.frame.color = cols,
              asp = asp)
  legend('bottom', legend = cols2$names,
         pch = 19, pt.cex = 1, col = as.character(cols2$cols),
         bty = "n", horiz = T, text.width = 0.2, cex = 0.5)

}


