require(moretrees)
require(igraph)
require(RColorBrewer)
require(stringr)
require(gridExtra)
require(data.table)
require(ggplot2)
require(ggtree)

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
  est <- exp(moretrees_results[beta_type][[1]][ , cols] * mult)
  est <- (est - 1) * 100
  est$outcomes <- moretrees_results$beta_moretrees$outcomes
  est$n_obs <- moretrees_results$beta_moretrees$n_obs
  vids <- unlist(est$outcomes)
  
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
  est$long_label <- character(length = nrow(est))
  est$short_label <- character(length = nrow(est))
  for (g in 1:nrow(est)) {
    ccs_g <- ccs[ccs$ccs_added_zeros %in% est$outcomes[[g]], ]
    ccs_g$label <- as.character(ccs_g$label)
    est$long_label[g] <- collapse_labels(ccs_g)
    # Collapse siblings
    i <- 1
    while (i <= nrow(ccs_g)) {
      code <- ccs_g$ccs_added_zeros[i]
      if (code == root) { # in this case, we have one group with all outcomes
        ccs_g$label[i] <- root 
        ccs_g$label[i] <- as.character(ccs_lvls$label[ccs_lvls$ccs_code == root])
        ccs_g$ccs_original[i] <- root
        i <- nrow(ccs_g) + 1
      } else {
        # Get siblings from tree
        parent <- names(ego(tr, order = 1, mindist = 1,
                            nodes = code, mode = "in")[[1]])
        sibs <- names(ego(tr, order = 1, mindist = 1,
                          nodes = parent, mode = "out")[[1]])
        if (sum(sibs %in% ccs_g$ccs_added_zeros) == length(sibs)) {
          ccs_g$ccs_added_zeros[i] <- parent
          parent <- str_remove_all(parent, "\\.0")
          ccs_g$ccs_original[i] <- parent
          ccs_g$label[i] <- as.character(ccs_lvls$label[ccs_lvls$ccs_code == parent])
          ccs_g <- ccs_g[!(ccs_g$ccs_added_zeros %in% setdiff(sibs, code)), ]
          i <- 1
        } else {
          i <- i + 1
        }
      }
    }
    # Collapse simplified labels into one string
    est$short_label[g] <- collapse_labels(ccs_g)
    # Collapse OR and CI into one string
    for (j in 1:k) {
      est_frmt <- sprintf(frmt, est[g, paste0("est", j)])
      cil_frmt <- sprintf(frmt, est[g, paste0("cil", j)])
      ciu_frmt <- sprintf(frmt, est[g, paste0("ciu", j)])
      est[g, paste0("ci_est", j)] <- paste0(est_frmt," (",cil_frmt,", ",ciu_frmt,")")
    }
  }
  est$n_outcomes <- sapply(est$outcomes, length)
  est$group <- 1:nrow(est)
  
  return(est[ , c("group", "short_label", "long_label",
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
  cols_g <- brewer.pal(length(outcomes), "Set3")
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

dt_plot_fun <- function(dt, plot_depth = 3) {
  for (i in 1:plot_depth) {
    dt[ , paste0(c("est_lvl", "cil_lvl", "ciu_lvl", "n"), i) := mean.diff.log(pm25_lag01_case, pm25_lag01_control),
        by = get(paste0("ccs_lvl", i))]
    if (i < 3) {
      dt[ , paste0("pltlab", i) := paste0(get(paste0("ccs_lvl", i)), ": ", 
                                          str_remove(get(paste0("label", i)), "\\s\\(.*\\)"), " (n = ", get(paste0("n", i)), ")")]
    }
    if (i >= 3) {
      dt[ , paste0("pltlab", i) := paste0(get(paste0("ccs_lvl", i)), " (n = ", get(paste0("n", i)), ")")]
    }
  }
  dt_plot <- unique(dt[ , paste0(c("pltlab", "est_lvl", "cil_lvl", "ciu_lvl", "n"), rep(1:plot_depth, each = 5)), with = FALSE])
  dt[ , paste0(c("pltlab", "est_lvl", "cil_lvl", "ciu_lvl", "n"), rep(1:plot_depth, each = 5)) := NULL]
  return(dt_plot)
}

nested_plots <- function(dt_plot, plot_depth = 3,
                         digits = 1, 
                         lab.widths = rep(1.5, plot_depth),
                         errorbar.height = 0.5,
                         lab.txt.size = 10,
                         lab.txt.width = rep(10, plot_depth),
                         axis.height = 1.5,
                         axis.txt.size = 2,
                         xlab = "PM2.5",
                         cil_min = -0.1,
                         ciu_max = 0.1,
                         wrap_labs = TRUE,
                         force_lims = FALSE) {
  # Remove unnecessary columns
  dt_plot[ , lab_col_num := as.integer(factor(ccs_lvl2, levels = unique(ccs_lvl2)))]
  dt_plot[ , lab_col_num2 := as.integer(factor(ccs_lvl3, levels = unique(ccs_lvl3)))]
  dt_plot <- dt_plot[ , 
                      c(paste0(c("pltlab", "est_lvl", "cil_lvl", "ciu_lvl"), 
                             rep(1:plot_depth, each = 5)),
                        "lab_col_num", "lab_col_num2"), with = FALSE]
  dt_plot <- dt_plot[!duplicated(dt_plot), ]
  # Get some plotting parameters
  dt_plot[ , lab_col := ifelse(lab_col_num %% 2 == 1, "grey70", "grey90")]
  dt_plot[ , lab_col2 := ifelse(lab_col_num %% 2 == 1, 
                                ifelse(lab_col_num2 %% 2 == 1, "grey65", "grey75"), 
                                ifelse(lab_col_num2 %% 2 == 1, "grey85", "grey95"))]
  lims <- c(max(min(dt_plot[ , paste0("cil_lvl", 1:plot_depth), with = FALSE], na.rm = T), cil_min),
            min(max(dt_plot[ , paste0("ciu_lvl", 1:plot_depth), with = FALSE], na.rm = T), ciu_max))
  if (force_lims) lims <- c(cil_min, ciu_max)
  x.ticks <- round(max(abs(lims)) * 3 / 4, digits = digits)
  x.grid <- x.ticks * c(-3/2, -1, -1/2, 1/2, 1, 3/2)
  layout <- integer()
  # Make x axis plot
  xaxis.plt <- local({
    dat <- data.table(x = c(-x.ticks, x.ticks), y = c(0, 0))
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
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, "cm")) +
      xlab(xlab) +
      scale_x_continuous(limits = lims,
                         breaks = c(-x.ticks, 0, x.ticks)) +
      coord_cartesian(expand = FALSE)
  })
  # Make plots
  plot.count <- 0
  grobs <- list()
  for (i in 1:plot_depth) {
    lab <- paste0("pltlab", i)
    diseases <- unique(dt_plot[ , get(lab)])
    # Make label plots
    for (disease in diseases) {
      plot.count <- plot.count + 1
      y.times <- sum(dt_plot[ , get(lab)] == disease)
      layout <- c(layout, rep(plot.count, times = y.times))
      if (i <= 2) lab.col <- as.character(dt_plot[get(lab) == disease]$"lab_col"[1])
      if (i >= 3) lab.col <- as.character(dt_plot[get(lab) == disease]$"lab_col2"[1])
      grobs[[plot.count]] <- local({
        y.height <- y.times / 2
        if (y.times > 1 & wrap_labs) {
          disease_nm <- str_remove(disease, "\\s\\(.*\\)")
          n_nm <- str_extract(disease, "\\(.*\\)")
          disease_wrap <- paste0(str_wrap(disease_nm, lab.txt.width[i]),
                                 "\n", n_nm)
        } else {
          disease_wrap <- disease
        }
        lab.col <- as.character(lab.col)
        grob <- ggplot(data.frame(disease = disease_wrap, x = 0, y = y.height), 
                       aes(x = x, y = y, label = disease)) + 
          geom_text(hjust = 0, size = lab.txt.size) +
          scale_x_continuous(lim = c(0, 1)) +
          theme_void() +
          theme(panel.border = 
                  element_rect(colour = "black", fill = NA, size = 0.3),
                panel.background = element_rect(fill = lab.col))
        grob
      })
    }
    plot.count <- plot.count + 1
    grobs[[plot.count]] <- ggplot() + geom_blank() + theme_void()
    layout <- c(layout, plot.count)
    # Make point range plots
    for (disease in diseases) {
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
        dat[, xmin] <- max(dat[, xmin], lims[1])
        dat[, xmax] <- min(dat[, xmax], lims[2])
        dat <- unique(dat)
        dat$disease <- disease
        grob <- ggplot(dat) + 
          geom_vline(xintercept = 0, color = "red", lwd = 0.4) +
          geom_vline(xintercept = x.grid[1], color = "grey70",
                     lwd = 0.2, lty = 2) +
          geom_vline(xintercept = x.grid[2], color = "grey70",
                     lwd = 0.2, lty = 2) +
          geom_vline(xintercept = x.grid[3], color = "grey70",
                     lwd = 0.2, lty = 2) +
          geom_vline(xintercept = x.grid[4], color = "grey70",
                     lwd = 0.2, lty = 2) +
          geom_vline(xintercept = x.grid[5], color = "grey70",
                     lwd = 0.2, lty = 2) +
          geom_vline(xintercept = x.grid[6], color = "grey70",
                     lwd = 0.2, lty = 2) +
          geom_errorbarh(aes(xmin = get(xmin), xmax = get(xmax), 
                             y = y.height), height = errorbar.height) +
          geom_point(aes(x = get(x), y = y.height)) +
          theme_void() +
          theme(panel.border = 
                  element_rect(colour = "black", fill = NA, size = 0.3)) +
          scale_x_continuous(limits = lims) + 
          scale_y_continuous(limits = c(y.min, y.max)) +
          coord_cartesian(expand = FALSE)
        grob
      })
    }
    plot.count <- plot.count + 1
    grobs[[plot.count]] <- xaxis.plt
    layout <- c(layout, plot.count)
  }
  layout_matrix <- matrix(layout, ncol = plot_depth * 2)
  widths <- c()
  for (i in 1:plot_depth) {
    widths <- c(widths, lab.widths[i], 1)
  }
  heights <- rep(1, nrow(layout_matrix))
  heights[length(heights)] <- axis.height
  grid.arrange(grobs = grobs, layout_matrix = layout_matrix,
               widths = widths, heights = heights,
               padding = unit(0, "line"))
}

equal_betas <- function(v1, v2, prob, levels, a_rho, b_rho, ancestors) {
  anc1 <- ancestors[[v1]]
  anc2 <- ancestors[[v2]]
  non_common_anc <- setdiff(union(anc1, anc2),
                            intersect(anc1, anc2))
  if (!is.null(prob)){
    p_return <- prod(1 - prob[non_common_anc])
  } else {
    n_l <- sapply(1:max(levels), function(l) sum(levels[non_common_anc] == l))
    p_return <- 1
    for (l in 1:length(n_l)) {
      if (n_l[l] > 0) {
        r <- 0:n_l[l]
        p_return <- p_return * prod((b_rho[l] + r) / (a_rho[l] + b_rho[l] + r))
      }
    }
  }
  return(p_return)
}

equal_betas_mat <- function(leaves, prob, levels, a_rho, b_rho, tr) {
  # reorder nodes
  nodes <- names(igraph::V(tr))
  leaves <- names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0])
  # nodes <- c(nodes[!(nodes %in% leaves)], leaves)
  # names(prob) <- nodes
  # levels <- levels[match(names(igraph::V(tr)), nodes)]
  pL <- length(leaves)
  
  # get ancestors
  d <- igraph::diameter(tr)
  ancestors <- igraph::ego(tr, order = d + 1, nodes = leaves, mode = "in")
  ancestors <- lapply(ancestors, names)
  ancestors <- sapply(ancestors, function(a, nodes) which(nodes %in% a), nodes = nodes,
                      simplify = F)
  names(ancestors) <- leaves
  
  # get equality matrix
  pmat <- matrix(nrow = pL, ncol = pL)
  diag(pmat) <- 1
  rownames(pmat) <- colnames(pmat) <- leaves
  for (v1 in 1:(length(leaves) - 1)) {
    for (v2 in (v1 + 1):length(leaves)) {
      pmat[v1, v2] <- equal_betas(v1 = leaves[v1], v2 = leaves[v2],
                                  prob = prob, levels = levels,
                                  a_rho = a_rho, b_rho = b_rho,
                                  ancestors = ancestors)
      pmat[v2, v1] <- pmat[v1, v2]
    }
  }
  
  return(pmat)
}

equal_betas_plot <- function(prob = NULL,
                             groups = NULL,
                             levels = NULL,
                             a_rho = NULL,
                             b_rho = NULL,
                             tr,
                             ccs.text.size = 4,
                             group.text.size = 4,
                             legend.text.size = 14,
                             heatmap.offset = 0.5,
                             heatmap.width = 6,
                             rownames.lab.offset = 0.8,
                             group.lab.offset = 0.3,
                             show.groups = T
) {
  
  leaves <- names(igraph::V(tr)[igraph::degree(tr, mode = "out") == 0])
  if (!(is.null(prob)) & is.null(levels)) {
    # Set default levels if null
    levels <- rep(1, length(V(tr)))
    levels[names(igraph::V(tr)) %in% leaves] <- 2
  }
  pmat <- equal_betas_mat(leaves = leaves, prob = prob, 
                          levels = levels, a_rho = a_rho, b_rho = b_rho, 
                          tr = tr)
  colnames(pmat) <- str_remove_all(colnames(pmat), "\\.0")
  groups.df <- data.frame(leaves = leaves, 
                          leafnames = str_remove_all(leaves, "\\.0"))
  if (is.null(groups)) {
    groups.df$Group <-  as.factor(rep(1, nrow(groups.df)))
  } else {
    groups.df$Group <-  as.factor(groups)
  }
  
  p <- ggtree(tr, ladderize = F)
  # Add tip labels
  p <- p  %<+% groups.df + geom_tiplab(aes(label = leafnames),
                                       size = ccs.text.size,
                                       offset = rownames.lab.offset)
  # Add group colours
  if (show.groups) {
    cols_g <- brewer.pal(length(levels(groups.df$Group)), "Set3")
    p <- p + geom_tippoint(aes(color = Group),
                           shape = 15, size = 4) +
             geom_tiplab(aes(label = Group), size = group.text.size,
                  offset = group.lab.offset)
    if (length(cols_g) == length(levels(groups.df$Group))) {
      p <- p + scale_color_manual(values = cols_g)
    }
  }
  # Add pmat heatmap
  p <- gheatmap(p, pmat, 
                offset = heatmap.offset,
                colnames_angle = 90,
                # colnames_offset_y = cols.names.offset,
                colnames_position = "bottom",
                hjust = 0,
                width = heatmap.width,
                font.size = ccs.text.size)
  # heatmap colors
  p <- p + scale_colour_viridis_c(aesthetics = "fill", 
                                  na.value = "grey90",
                                  limits = c(0, 1)) +
    labs(fill = expression(P("same "*beta[v])))   +
    theme(legend.position = "left",
          legend.text = element_text(size = legend.text.size),
          legend.title = element_text(size = legend.text.size),
          legend.margin = margin(0, -30, 0, 0),
          legend.box.margin = margin(0, 0, 0, 0),
          plot.margin = unit(c(0.5, 0.5, -0.5, 0), "cm")) +
    coord_cartesian(clip = "off")
  p <- p + scale_y_reverse()
  p <- p + guides(colour = guide_legend(order = 1), 
                fill = guide_colourbar(order = 2))
  return(p)
}

## Simulate groups from prior
sim.prior.fun <- function(levels, A_leaf, 
                          a_rho = rep(1, max(levels)), 
                          b_rho = rep(1, max(levels))) {
  rho <- numeric(max(levels))
  for (l in 1:max(levels)) {
    rho[l] <- rbeta(1, shape1 = a_rho[l], shape2 = b_rho[l])
  }
  s <- runif(length(levels)) <= rho[levels]
  gamma <- numeric(length = length(levels)) + 0
  gamma[s] <- rnorm(sum(s), sd = 10)
  beta <- as.numeric(A_leaf %*% gamma)
  groups <- as.integer(as.factor(beta))
  groups.size <- table(groups)
  return(data.frame(n.groups = max(groups), 
           mean.size = mean(groups.size),
           median.size = median(groups.size)))
}

beta_indiv_plot_fun <- function(pltdat, tr, ...) {
  L <- max(V(tr)$levels)
  dt <- t(sapply(igraph::ego(tr, order = L, nodes = V(tr)[V(tr)$leaf], mode = "in"), names))
  dt <- data.table(dt)
  names(dt) <- paste0("ccs_lvl", 4:1)
  for (l in 1:L) {
    dt <- merge(dt, pltdat, by.x = paste0("ccs_lvl", l), by.y = "node", sort = FALSE)
    setnames(dt, c("est1", "cil1", "ciu1", "n"), paste0(c("est_lvl", "cil_lvl", "ciu_lvl", "n"), l))
    dt[ , paste0("pltlab", l) := paste0(str_remove_all(get(paste0("ccs_lvl", l)), "\\.0"),
                                        " (n = ", get(paste0("n", l)), ")")]
    # if (l < L) {
    #   for (code in unique(dt[ , paste0("ccs_lvl", l), with = FALSE][[1]])) {
    #     if(sum(dt[ , paste0("ccs_lvl", l), with = FALSE] == code) == 1) {
    #       dt[dt[ , paste0("ccs_lvl", l), with = FALSE][[1]] == code, 
    #          paste0(c("est_lvl", "cil_lvl", "ciu_lvl"), l) := NA]
    #        dt[dt[ , paste0("ccs_lvl", l), with = FALSE][[1]] == code, 
    #           paste0("pltlab", l) := paste0("blank", get(paste0("pltlab", l)))]
    #     }
    #   }
    # }
  }
  
  p <- nested_plots(dt, ...)
  return(p)
}
