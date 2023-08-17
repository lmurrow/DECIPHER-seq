library(rlist)
library(ape)
library(matrixStats)

plot_phylo_trees_by_rank <- function(phylo_trees, phylo_partitions){
  require(ape)
  for (i in names(phylo_trees)){
    tree = phylo_trees[[i]]$phylo_tree_pruned
    clusts = phylo_partitions[[i]]
    phy_sub = drop.tip(tree, names(clusts)[which(clusts%in%c("0", "outlier"))])
    plot(phy_sub, cex = 0.7, show.tip.label = F, type = "phylo", no.margin = F, main = i)
    ranks_run = as.numeric(gsub("R", "", gsub("_.*", "", phy_sub$tip.label)))
    tiplabels(pch = 16, cex = 1, 
              col = viridis_pal()(max(ranks_run)-min(ranks_run)+1)[ranks_run])
    tiplabels(pch = 1, lwd = 0.1, cex = rep(1, nclust)[labels], 
              col = rep(alpha("black", 0.3), nclust)[labels])
  }
}

plot_phylo_trees_for_rank <- function(phylo_trees, phylo_partitions, rank_plot = 20){
  for (i in names(phylo_trees)){
    require(ape)
    tree = phylo_trees[[i]]$phylo_tree_pruned
    clusts = phylo_partitions[[i]]
    phy_sub = drop.tip(tree, names(clusts)[which(clusts%in%c("0", "outlier"))])
    phy_sub$tip.label[-grep(paste0("R", rank_plot, "_"), phy_sub$tip.label)] = NA
    plot(phy_sub, cex = 0.7, show.tip.label = T, type = "phylo", no.margin = F, main = i)
    ranks_run = as.numeric(gsub("R", "", gsub("_.*", "", phy_sub$tip.label)))
    ranks_run[ranks_run!=rank_plot] = "grey"
    ranks_run[ranks_run==rank_plot] = "red"
    tiplabels(pch = 16, cex = c(0.2, 1)[as.numeric(as.numeric(gsub("R", "", gsub("_.*", "", orig.labs)))==rank_plot)+1], 
              col = ranks_run)
    tiplabels(pch = 1, lwd = 0.1,cex = c(0.2, 1)[as.numeric(as.numeric(gsub("R", "", gsub("_.*", "", orig.labs)))==rank_plot)+1], 
              col = rep(alpha("black", 0.3), nclust)[labels])
  }
}

