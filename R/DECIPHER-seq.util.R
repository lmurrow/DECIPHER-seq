# Utilities

cell_type_cols=c('#33A02B','#875890','#B2DF8A','#A6CEE3','#F8766D','#1F78B4','#E68D1A','#CEA6A6' ,'#ffcc12')
network_module_cols = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#E41A1C", "#377EB8")

my_theme = theme_bw() + theme(panel.border=element_rect(color="black", fill=NA, size=0.5),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# plot histograms ####
plot_hists = function(phylo_trees, suggested_thresholds = NULL, thresh.use = NULL){
  if (is.null(suggested_thresholds)){
    suggested_thresholds = suggest_dist_thresh(phylo_trees)
  }
  if (suggested_thresholds[1]>suggested_thresholds[2]){
    warning("no common threshold found")
  }
  if (is.null(thresh.use)){
    thresh.use = round(max(suggested_thresholds), 1)
  }
  x.max = max(unlist(lapply(phylo_trees, `[[`, "distance_matrix")))
  p = list()
  for (i in names(phylo_trees)){
    x = phylo_trees[[i]]
    df = data.frame(dist = x[["distance_matrix"]][upper.tri(x[["distance_matrix"]])])
    bins = 100
    p[[i]] = ggplot(df, aes(x = dist)) + 
      geom_histogram(aes(y = stat(count)/sum(count)), bins = bins, fill = "grey") + 
      geom_vline(xintercept = thresh.use, linetype = "dashed", color = "red") +
      my_theme + labs(x = "Pairwise patristic distance", y = "Frequency", title = i) +
      scale_x_continuous(limits = c(0, x.max), expand = c(0,0))
  }
  return(p)
}


# plot phylogenetic trees ####

color_set = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#1B9E77", "#D95F02", "#7570B3",
              "#E7298A", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#6C9A91", "#DF55BD",
              "#DEE14C", "#D19A83", "#DFBBD6", "#89E055", "#859C59", "#E2A550", "#6AE195", "#8664D5", "#D0E588",
              "#D894D6", "#6CE3DE", "#DDD9AA", "#9DE2BA", "#D1586D", "#7EC3E4", "#A33BE3")

plot_phylo_trees <- function(phylo_trees, phylo_partitions){
  for (i in names(phylo_trees)){
    tree = phylo_trees[[i]]$phylo_tree_pruned
    clusts = phylo_partitions[[i]]
    phy_sub = drop.tip(tree, names(clusts)[which(clusts%in%c("0", "outlier"))])
    plot(phy_sub, cex = 0.7, show.tip.label = F, type = "phylo", no.margin = F, main = i)
    labels = clusts[which(clusts%ni%c("0", "outlier"))]
    labels = factor(labels)
    nclust = sum(levels(labels)%ni%c("0", "outlier"))
    tiplabels(pch = 16, cex = rep(1, nclust)[labels], 
              col = color_set[1:(nclust)][labels])
    tiplabels(pch = 1, lwd = 0.1, cex = rep(1, nclust)[labels], 
              col = rep(alpha("black", 0.3), nclust)[labels])
  }
}

# plot module heatmap ####
Heatmap_wrapper <- function(Expression_score_cor, Network){
  mat = Expression_score_cor$cor[names(V(Network$filtered_network)),names(V(Network$filtered_network))]
  diag(mat) = 1
  
  cluster_order = names(table(Network$filtered_modules))[c(2,4,3,1,5,6,7,8)]
  order = NULL
  for (n in cluster_order){
    order = c(order, names(which(Network$filtered_modules==n)))
  }
  
  
  celltypes=factor(gsub("\\.R[1-9].*", "", V(Network$filtered_network)[order]$name), 
                   levels=c("HRpos_Luminal","Basal","Secretory_Luminal","Fibroblast", "Vascular_Endothelial"))
  par(xpd = T)
  heatmap.2(mat[order,order], trace="none", density.info="none", scale="none", 
            breaks=seq(-1,1,length.out=101), 
            ColSideColors = network_module_cols[c(2,4,3,1,5,6,7,8)][factor(Network$filtered_modules[order], levels = cluster_order)],
            RowSideColors = cell_type_cols[celltypes],
            col=colorRampPalette(rev(brewer_pal(palette="RdBu")(10)))(100), margins=c(2,2), 
            Rowv = F, Colv = "Rowv", dendrogram='none', na.color="grey80",
            labRow = NA, labCol = NA)
  legend("left",  inset = -0.05,    
         legend = 1:length(unique(celltypes)),
         col = cell_type_cols, 
         pch = 15, bty = "n", title = "Cell Type")
  legend("top",  inset = -0.2,    
         legend = 1:length(unique(Network$filtered_modules)),
         col = network_module_cols[c(2,4,3,1,5,6,7,8)], 
         pch = 15, bty = "n", title = "Module", 
         ncol = length(unique(Network$filtered_modules)))
}

# Plot fgsea enrichment on network ####
plot_fgsea <- function(Network, fgsea_res, gene_set, fdr_max = 1e-2, fdr_min = 1e-10){
  if (gene_set %ni% fgsea_res$pathway){
    stop("gene_set not in fgsea_res")
  }
  fgsea_res = subset(fgsea_res, pathway==gene_set)
  fdr_vect = fgsea_res$fdr
  names(fdr_vect) = paste0(fgsea_res$Type, ".", fgsea_res$Program)
  fdr_vect[subset(fgsea_res)$NES<0]=1
  fdr_vect[is.na(fdr_vect)] = 1
  fdr_vect = fdr_vect[names(Network$filtered_modules)]
  
  plot(Network$filtered_network, layout = Network$filtered_network_coords, 
       vertex.color=colorRampPalette(c("grey", "red"))(100)[cut(-log10(fdr_vect), breaks = c(0, seq(-log10(fdr_max), -log10(fdr_min), length.out = 99), max(-log10(fdr_vect))), include.lowest = T, right = F)],
       edge.width=0.25, vertex.label.cex=0.5, vertex.size = 7, vertex.label.color="black", 
       vertex.label.family="Helvetica", vertex.frame.color=NA, vertex.label.font=2, vertex.label = NA,
       edge.color = c(NA, "grey20")[factor(E(Network$filtered_network)$sign>0)], main = gene_set)
  legend_image <- as.raster(matrix(rev(colorRampPalette(c("grey", "red"))(100)), ncol=1))
  text(x=1.7, y = c(-0.3,0.3), labels = c(fdr_max, fdr_min), cex = 0.7)
  rasterImage(legend_image, 1.35, -0.35, 1.5,0.35)
}

# Plot association of metadata features with nodes in network ####
Plot_metadata_association <- function(Network, feature_effect_size, plot.title = NULL){
  feature_effect_size = feature_effect_size[,names(V(Network$filtered_network))]
  node_val = as.numeric(feature_effect_size["effect_size",])*as.numeric(feature_effect_size["effect_size_sign",])
  node_size = rep(3, length(V(Network$filtered_network)))
  node_size[feature_effect_size["p_value",]<0.05] = 5
  node_size[feature_effect_size["p_value",]<0.01] = 7
  node_size[feature_effect_size["p_value",]<0.001] = 10
  
  plot(Network$filtered_network, layout = Network$filtered_network_coords, 
       vertex.color=colorRampPalette(c("blue","grey", "red"))(100)[cut(node_val, breaks = seq(-1, 1, length.out = 101), include.lowest = T)],
       edge.width=0.25, vertex.label.cex=0.5, vertex.size = node_size, vertex.label.color="black", 
       vertex.label.family="Helvetica", vertex.frame.color=NA, vertex.label.font=2,
       edge.color = c(NA, "grey20")[factor(E(Network$filtered_network)$sign>0)], vertex.label = NA, main = plot.title)
  legend_image <- as.raster(matrix(rev(colorRampPalette(c("blue","grey", "red"))(100)), ncol=1))
  text(x=1.7, y = c(-0.4,0.2), labels = c("min", "max"), cex = 0.7)
  text(x=1.6, y = 0.5, labels = c("Effect size"), cex = 1)
  rasterImage(legend_image, 1.35, -0.45, 1.5,0.25)
  legend("left",  inset = 0,    
         legend = c("n.s.", "0.05", "0.01", "0.001"),
         col = "grey", pt.cex = c(0.3, 0.5, 0.7, 1)*2,
         pch = 19, bty = "n", title = "p-value")
}
