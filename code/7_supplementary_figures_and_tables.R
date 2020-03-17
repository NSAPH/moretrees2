require(ggtree)
require(igraph)

### CVD tree ----------------------------------------------------------------
load(file = paste0("./results/mod1_split0_cvd_full.Rdata"))

# Get tree
tr <- ccs_tree("7")$tr
vids <- unlist(moretrees_results$beta_moretrees$outcomes)
vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
vids <- Reduce(union, vids)
tr <- induced_subgraph(tr, vids)

# Get labels
ccs_labels <- read.csv("./data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")

# Make horrible plot
root <- names(igraph::V(tr))[igraph::degree(tr, mode = "in") == 0]
V(tr)$levels <- as.numeric(igraph::distances(tr, v = root, to = V(tr), mode = "out") + 1)
labels <- data.frame(node = names(V(tr)),
                     nodename = str_remove_all(names(V(tr)), "\\.0"), 
                     levels = V(tr)$levels,
                     stringsAsFactors = F)
labels <- merge(labels, ccs_labels, by.x = "nodename", by.y = "ccs_code",
                all.x = T, all.y = F, sort = F)
labels$label <- paste0(labels$nodename, ": ", labels$label)
labels$nodename <- NULL

p <- ggtree(tr, ladderize = F)
labels <- labels[match(p$data$label, labels$node), ]
p$data$levels <- labels$levels
for (i in 1:4) {
  p$data[ , paste0('label.long', i)] <- as.character(NA)
  labs <- labels$label[labels$levels == i]
  if (i < 3) {
    labs <- str_wrap(labs, width = 20)
  }
  if (i == 3) {
    labs <- str_wrap(labs, width = 60)
  }
  p$data[p$data$levels == i , paste0('label.long', i)] <- labs
}
node.pos <- c(-0.4, 0.2, 1.85, 1.9)
node.nudge <- c(-0.55,-0.58, -1.6, 0)
node.angle <- c(90, 90, 0, 0)
p$data$x <- node.pos[p$data$levels]
p <- p + geom_tiplab(aes(label = label.long4)) + xlim(c(-3, 5))
for (i in 1:3) {
  p <- p + geom_nodelab(aes_string(label = paste0("label.long", i)),
                        geom = "label",
                        hjust = 0,
                        nudge_x = node.nudge[i],
                        angle = node.angle[i],
                        fill = "white",
                        label.size = NA)
}
p <- p + theme(plot.margin = unit(c(0, -4.7, 0, -6.8), unit = "in"))
p + coord_cartesian(clip = "off")

pdf(file = "./figures/ccs_tree_cvd.pdf", width = 13, height = 15)
p
dev.off()

### RD tree ----------------------------------------------------------------
load(file = paste0("./results/mod1_split0_resp_full.Rdata"))

# Get tree
tr <- ccs_tree("8")$tr
vids <- unlist(moretrees_results$beta_moretrees$outcomes)
vids <- ego(graph = tr, order = diameter(tr) + 10, nodes = vids, mode = "in")
vids <- Reduce(union, vids)
tr <- induced_subgraph(tr, vids)

# Get labels
ccs_labels <- read.csv("./data/Multi_Level_CCS_2015_cleaned/dxm_label_clean.csv")

# Make horrible plot
V(tr)$levels <- as.numeric(igraph::distances(tr, v = root, to = V(tr), mode = "out") + 1)
labels <- data.frame(node = names(V(tr)),
                     nodename = str_remove_all(names(V(tr)), "\\.0"), 
                     levels = V(tr)$levels,
                     stringsAsFactors = F)
labels <- merge(labels, ccs_labels, by.x = "nodename", by.y = "ccs_code",
                all.x = T, all.y = F, sort = F)
labels$label <- paste0(labels$nodename, ": ", labels$label)
labels$nodename <- NULL

p <- ggtree(tr, ladderize = F)
labels <- labels[match(p$data$label, labels$node), ]
p$data$levels <- labels$levels
for (i in 1:4) {
  p$data[ , paste0('label.long', i)] <- as.character(NA)
  labs <- labels$label[labels$levels == i]
  if (i == 1) labs <- str_wrap(labs, width = 20)
  if (i == 2) {
    labs <- str_wrap(labs, width = 30)
  }
  if (i == 3) {
    labs <- str_wrap(labs, width = 60)
  }
  p$data[p$data$levels == i , paste0('label.long', i)] <- labs
}
node.pos <- c(-0.4, 0.5, 2.15, 2.2)
node.nudge <- c(-0.55,-0.85, -1.6, 0)
node.angle <- c(90, 90, 0, 0)
p$data$x <- node.pos[p$data$levels]
p <- p + geom_tiplab(aes(label = label.long4)) + xlim(c(-3.1, 5.6))
for (i in 1:3) {
  p <- p + geom_nodelab(aes_string(label = paste0("label.long", i)),
                        geom = "label",
                        hjust = 0,
                        vjust = 1,
                        nudge_x = node.nudge[i],
                        angle = node.angle[i],
                        fill = "white",
                        label.size = NA)
}
p <- p + theme(plot.margin = unit(c(0, -4.7, 0, -6.8), unit = "in"))
p + coord_cartesian(clip = "off")

pdf(file = "./figures/ccs_tree_resp.pdf", width = 14, height = 11)
p
dev.off()
