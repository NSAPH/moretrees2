require(moretrees)
require(igraph)
require(RColorBrewer)
require(stringr)
require(gridExtra)
require(data.table)
require(ggplot2)

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

ccs_table <- function(root, moretrees_results, digits = 3, mult = 10,
                      type = "moretrees", tr = NULL) {
  
  ccs_lvls <- read.csv("./data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")
  
  # Get ORs
  k <- length(moretrees_results$mod$vi_params$mu[[1]])
  cols <- unlist(lapply(1:k, function(k) paste0(c("est", "cil", "ciu"), k)))
  beta_type <- paste0("beta_", type)
  OR_est <- exp(moretrees_results[beta_type][[1]][ , cols] * mult)
  OR_est$outcomes <- moretrees_results$beta_moretrees$outcomes
  OR_est$n_obs <- moretrees_results$beta_moretrees$n_obs
  vids <- unlist(OR_est$outcomes)
  
  # Map CCS codes to diseases --------------------------------------------------
  ccs_tr <- ccs_tree(root)
  ccs_icd9_mapping <- ccs_tr$ccs_icd_mapping
  ccs_zeros <- ccs_icd9_mapping[, c("ccs_original", "ccs_added_zeros")]
  ccs_zeros <- ccs_zeros[!duplicated(ccs_zeros), ]
  ccs <- merge(ccs_zeros, ccs_lvls, by.x = "ccs_original", by.y = "ccs_code", sort = F)
  
  # Get tree if necessary ------------------------------------------------------
  if (is.null(tr)) tr <- ccs_tr$tr
  ccs <- ccs[ccs$ccs_added_zeros %in% names(V(tr)), ]
  ccs <- ccs[order(match(ccs$ccs_added_zeros, names(V(tr)))), ]
  
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
        parent <- str_remove_all(parent, "\\.0")
        ccs_g$ccs_original[i] <- parent
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
      est_frmt <- sprintf(frmt, OR_est[g, paste0("est", j)])
      cil_frmt <- sprintf(frmt, OR_est[g, paste0("cil", j)])
      ciu_frmt <- sprintf(frmt, OR_est[g, paste0("ciu", j)])
      OR_est[g, paste0("ci_est", j)] <- paste0(est_frmt," (",cil_frmt,", ",ciu_frmt,")")
    }
  }
  OR_est$n_outcomes <- sapply(OR_est$outcomes, length)
  OR_est$group <- 1:nrow(OR_est)
  
  return(OR_est[ , c("group", "short_label", "long_label",
                     "n_outcomes", "n_obs", paste0("ci_est", 1:k))])
}
  
ccs_plot <- function(root, moretrees_results, tr = NULL,
                     asp = 1/5, internal.width = 0.1,
                     leaf.width = 2, internal.height = 0.1,
                     leaf.height = 16,
                     label.dist = 3,
                     label.cex = 0.3,...) {
  
  # Get tree
  if(is.null(tr)) tr <- ccs_tree(root)$tr
  
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
  vwidth <- rep(internal.width, length(V(tr)))
  vwidth[V(tr)$leaf] <- leaf.width
  vheight <- rep(internal.height, length(V(tr)))
  vheight[V(tr)$leaf] <- leaf.height
  
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
  plot.igraph(tr, layout = l,
              edge.arrow.size = 0, vertex.shape = "rectangle",
              vertex.size = vwidth, vertex.size2 = vheight,
              vertex.color = cols, 
              vertex.frame.color = cols,
              vertex.label = groups,
              vertex.label.dist = label.dist,
              vertex.label.degree = pi/2,
              vertex.label.family = "sans",
              vertex.label.cex = label.cex,
              vertex.label.color = "black",
              asp = asp)
  legend('bottom', legend = cols2$names,
         pch = 15, pt.cex = 1, col = as.character(cols2$cols),
         bty = "n", horiz = T, text.width = 0.2, cex = 0.5)

}

get_labels <- function(root) {
  ccs_labels <- read.csv("./data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")
  ccs_labels$label <- paste0(ccs_labels$label, " (", ccs_labels$ccs, ")")
  tr_cvd <- ccs_tree(root)$tr
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
  return(lvls_cvd)
}

mean.diff.log <- function(x , y) {
  ttest <- t.test(log(x), log(y), paired = TRUE)
  list(est = ttest$estimate, cil = ttest$conf.int[1], ciu = ttest$conf.int[2], n = length(x))
}

dt_plot_fun <- function(dt, lab_var = "ccs_lvl", plot_depth = 3) {
  for (i in 1:plot_depth) {
    dt[ , paste0(c("est_lvl", "cil_lvl", "ciu_lvl", "n"), i) := mean.diff.log(pm25_lag01_case, pm25_lag01_control),
        by = get(paste0(lab_var, i))]
    dt[ , paste0("pltlab", i) := paste0(get(paste0(lab_var, i)), " (n = ", get(paste0("n", i)), ")")]
  }
  dt_plot <- unique(dt[ , paste0(c("pltlab", "est_lvl", "cil_lvl", "ciu_lvl", "n"), rep(1:plot_depth, each = 5)), with = FALSE])
  dt[ , paste0(c("pltlab", "est_lvl", "cil_lvl", "ciu_lvl", "n"), rep(1:plot_depth, each = 5)) := NULL]
  return(dt_plot)
}

nested_plots <- function(dt_plot, plot_depth = 3,
                         digits = 1, lab.nudge = 1.5, 
                         errorbar.height = 0.5,
                         lab.size = 1,
                         axis.txt.size = 0.7,
                         xlab = "PM2.5") {
  plot.count <- 0
  grobs <- list()
  lims <- c(min(dt_plot[ , paste0("cil_lvl", 1:plot_depth), with = FALSE]),
            max(dt_plot[ , paste0("ciu_lvl", 1:plot_depth), with = FALSE]))
  x.grid <- round(max(abs(lims)) / 2, digits = digits)
  layout <- integer()
  xaxis.plt <- local({
    dat <- data.table(x = c(-x.grid, x.grid), y = c(0, 0))
    ggplot(dat) + geom_blank() +
      theme(axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.line.x = element_line(colour = "black", size = 0.5),
            axis.text.x = element_text(colour = "black", size = axis.txt.size),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.background = element_blank(),
            plot.margin = margin(t = 0, r = 0, b = 0.1,l = 0, "cm")) +
      xlab(xlab) +
      scale_x_continuous(limits = c(lims[1] - lab.nudge, 
                                    lims[2]),
                         breaks = c(-x.grid, 0, x.grid))

  })
  for (i in 1:plot_depth) {
    lab <- paste0("pltlab", i)
    for (disease in unique(dt_plot[ , get(lab)])) {
      plot.count <- plot.count + 1
      y.times <- sum(dt_plot[ , get(lab)] == disease)
      layout <- c(layout, rep(plot.count, times = y.times))
      grobs[[plot.count]] <- local({
        y.min <- 0
        y.max <- y.times
        y.height <- y.times / 2
        i <- i
        lab <- lab
        xmin <- paste0("cil_lvl", i)
        xmax <- paste0("ciu_lvl", i)
        x <- paste0("est_lvl", i)
        dat <- as.data.frame(dt_plot[get(lab) == disease])
        dat <- dat[ , c(x, xmin, xmax)]
        dat <- unique(dat)
        dat$disease <- disease
        grob <- ggplot(dat) + 
          geom_vline(xintercept = 0, color = "red") +
          geom_vline(xintercept = -x.grid, color = "grey80",
                     lwd = 0.2, lty = 2) +
          geom_vline(xintercept = x.grid, color = "grey80",
                     lwd = 0.2, lty = 2) +
          geom_point(aes(x = get(x), y = y.height)) +
          geom_errorbarh(aes(xmin = get(xmin), xmax = get(xmax), 
                             y = y.height), height = errorbar.height) +
          geom_text(aes(label = disease, 
                        x = lims[1] - lab.nudge, 
                        y = y.height),
                    hjust = 0, size = lab.size) +
          theme_void() +
          theme(panel.border = 
                  element_rect(colour = "black", fill = NA, size = 0.1)) +
          scale_x_continuous(limits = 
                               c(lims[1] - lab.nudge, 
                                 lims[2])) + 
          scale_y_continuous(limits = c(y.min, y.max))
        grob
      })
    }
    plot.count <- plot.count + 1
    grobs[[plot.count]] <- xaxis.plt
    layout <- c(layout, plot.count)
  }
  layout_matrix <- matrix(layout, ncol = plot_depth)
  grid.arrange(grobs = grobs, layout_matrix = layout_matrix, 
               padding = unit(0, "line"))
}
