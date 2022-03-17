##### Set up LIGER objects and run consensus iNMF across range of K values #################
iNMF_ksweep <- function(seurat.object, Type.thresh = 100, Sample.thresh = 10, Batch = TRUE, 
                        scale.factor = 10000, k.max = 40, n.reps = 20, nn.pt = 0.3,
                        max.cores = NULL, output.reps = FALSE, return.scale.data = T){
  if (!"Sample"%in%colnames(seurat.object@meta.data) | 
      !"Type"%in%colnames(seurat.object@meta.data)){
    stop("metadata slot must contain colnames 'Sample' and 'Type'")
  }
  if (Batch){
    if (!"Sample"%in%colnames(seurat.object@meta.data)){
      stop("metadata slot must contain colname 'Batch'")
    }
  }
  if (round(n.reps*nn.pt)<=1){
    stop("n.reps or nn.pt too low to find consensus results")
  }
  library(rliger); library(Matrix); library(matrixStats); library(parallel); 
  library(cluster); library(stats); library(rlist); library(liger); library(RANN)
  
  # Process Seurat object and determine which cell types to analyze
  DefaultAssay(seurat.object)="RNA"
  if ("integrated" %in% Assays(seurat.object)){
    seurat.object[['integrated']] <- NULL
  }
  # Only run analysis on cell types that represented across at least 10 samples with at least 100 cells per sample (or adjust inputs above)
  cell.types <- names(which(colSums(table(seurat.object@meta.data[,c("Sample","Type")])>=Type.thresh)>=Sample.thresh))
  if (length(cell.types)==0){
    stop("Not enough cells to meet input Sample/Type thresholds")
  } else if (length(cell.types)==1){
    warning("Only one cell type met input Sample/Type thresholds")
  }
  print(paste0("Running consensus iNMF on the following ", length(cell.types), " cell types: ", paste(cell.types, collapse = ", ")))
  
  # SET UP LIGER OBJECTS AND RUN K SWEEP
  print(paste0('Creating Liger objects and running k sweep...'))
  names(cell.types) = cell.types
  res <- lapply(cell.types, function(x){
    print(paste0("Cell type ", x))
    seu_sub <- subset(seurat.object, Type==x)
    if (Batch){
      if (length(table(seurat.object$Batch))>1){
        Batch.list <- SplitObject(seu_sub, split.by="Batch")
        Liger.setup=list()
        for (j in 1:length(Batch.list)){
          Liger.setup[[j]]=Batch.list[[j]]@assays$RNA@counts
        }
        names(Liger.setup)=names(Batch.list)
        Liger <- createLiger(Liger.setup)
      } 
    } else {
      Liger <- createLiger(list(Batch1 = seu_sub@assays$RNA@counts))
    }
    # normalize so all cells have same total counts
    Liger <- rliger::normalize(Liger)
    Liger <- selectGenes(Liger)
    # log normalize (this combined with the normalize step above is the same as LogNormalize in Seurat)
    for (k in 1:length(Liger@norm.data)){
      Liger@norm.data[[k]]=log1p(Liger@norm.data[[k]]*scale.factor)
    }
    # scale without centering
    Liger <- rliger::scaleNotCenter(Liger)
    
    # adjust online_iNMF parameters based on dataset size
    minibatch_size <- min(table(seu_sub@meta.data$Batch))
    if (minibatch_size > 5000) {
      minibatch_size = 5000} else if (minibatch_size > 1000) { 
        minibatch_size = floor(minibatch_size/1000)*1000} else if (minibatch_size > 500) { 
          minibatch_size = floor(minibatch_size/500)*500} else if (minibatch_size > 100) { 
            minibatch_size = floor(minibatch_size/100)*100} else minibatch_size = minibatch_size 
    h5_chunk_size <- min(1000, minibatch_size)
    
    Liger_list = list()
    for (k in c(2:k.max)){
      print(paste0('Running K = ', k))
      reps.k = rep(k, n.reps)
      names(reps.k) = paste0("rep", 1:length(reps.k))
      Liger_list[[paste0("R", k)]] = mcmapply(LIGER.run, K = reps.k, Liger = c(Liger), 
                                            minibatch_size = minibatch_size, h5_chunk_size = h5_chunk_size,
                                            SIMPLIFY = F)
      for (i in 1:length(Liger_list[[paste0("R", k)]])){
        colnames(Liger_list[[paste0("R", k)]][[i]][["W"]]) =  paste0("Rep", i, "_", colnames(Liger_list[[paste0("R", k)]][[i]][["W"]])) 
        colnames(Liger_list[[paste0("R", k)]][[i]][["H"]]) = paste0("Rep", i, "_", colnames(Liger_list[[paste0("R", k)]][[i]][["H"]]))
        Liger_list[[paste0("R", k)]][[i]][["V"]] = lapply(Liger_list[[paste0("R", k)]][[i]][["V"]], function(y) {
          colnames(y) = paste0("Rep", i, "_", colnames(y))
          return(y)})
      }
      print(paste0('Done with ', x, ' K = ', k))
    }
    
    # FIND CONSENSUS RESULTS
    print(paste0('Finding consensus results'))
    Consensus.results = mclapply(Liger_list, function(x){
      W_list = lapply(x, '[[', 'W')
      W_2norm = rapply(W_list, f = function(y) {apply(y, MARGIN = 2, FUN = function(y){norm(y, type ="2")})},
                       how = "replace")
      W_list = rapply(W_list, f = function(y) {t(apply(y, MARGIN = 2, FUN = function(y) {y/norm(y, type = "2")}))},
                      how = "replace")
      W.mat = list.rbind(W_list)
      kmeans_clusts = dim(W.mat)[1]/n.reps
      
      # find mean distance to nearest neighbors
      nn.dist = rowMeans(RANN::nn2(as.matrix(W.mat), k = n.reps)$nn.dists[,c(1:round(n.reps*nn.pt))+1])
      names(nn.dist) = rownames(W.mat)
      
      # filter outliers
      min = diff(quantile(nn.dist, probs = c(0.25, 0.75)))*0.5 + quantile(nn.dist, probs = c(0.75))
      if (sum(nn.dist>min)>0) {
        W.mat.filt = W.mat[-which(nn.dist>min),]
      } else {
        W.mat.filt = W.mat
      }
      
      # perform kmeans clustering
      km.res = kmeans(W.mat.filt, centers = kmeans_clusts, nstart = 100, iter.max = 100)
      
      # find consensus W (shared gene loadings)
      W_consensus = matrix(nrow = kmeans_clusts, ncol = ncol(W.mat.filt))
      for (i in seq(kmeans_clusts)){
        row.ind = which(km.res$cluster==i)
        if (length(row.ind) > 1){
          W_consensus[i,] = colMedians(W.mat.filt[row.ind,])
        } else W_consensus[i,] = W.mat.filt[row.ind,]
      }
      rownames(W_consensus) = paste0("R", kmeans_clusts, "_Program", seq(kmeans_clusts))
      colnames(W_consensus) = colnames(W.mat.filt)
      
      # find consensus V (batch-specific gene loadings)
      V_list = lapply(x, `[[`, "V")
      for (i in names(V_list)){
        V_list[[i]] = lapply(V_list[[i]], function(x){t(t(x)*(1/W_2norm[[i]]))})
      }
      V.mat = t(list.cbind(lapply(V_list, list.rbind)))
      if (sum(nn.dist>min)>0) {
        V.mat.filt = V.mat[-which(nn.dist>min),]
      } else {
        V.mat.filt = V.mat
      }
      V_consensus = matrix(nrow = kmeans_clusts, ncol = ncol(V.mat.filt))
      for (i in seq(kmeans_clusts)){
        row.ind = which(km.res$cluster==i)
        if (length(row.ind) > 1){
          V_consensus[i,] = colMedians(V.mat.filt[row.ind,])
        } else V_consensus[i,] = V.mat.filt[row.ind,]
      }
      rownames(V_consensus) = paste0("R", kmeans_clusts, "_Program", seq(kmeans_clusts))
      colnames(V_consensus) = colnames(V.mat.filt)
      Batch.info = 1:length(V_list$rep1)
      names(Batch.info) = sort(names(V_list$rep1))
      V_consensus = lapply(Batch.info, function(x){
        V = V_consensus[,((x-1)*length(colnames(W_consensus))+1):(x*length(colnames(W_consensus)))]
      })
      
      # Solve for H (activity program expression scores) using consensus W and V initializations
      H_consensus_list = lapply(Batch.info, function(x){
        H = solveNNLS(rbind(t(W_consensus) + t(V_consensus[[x]]), 
                            sqrt(5) * t(V_consensus[[x]])), rbind(t(Liger@scale.data[[x]]), matrix(0, 
                                                                                                   dim(W_consensus)[2], dim(Liger@raw.data[[x]])[2])))
      })
      H_consensus_list = lapply(H_consensus_list, t)
      H_consensus = list.rbind(H_consensus_list)
      rownames(H_consensus) = rownames((lapply(x, `[[`, "H"))[[1]])
      colnames(H_consensus) = paste0("R", kmeans_clusts, "_Program", seq(kmeans_clusts))
      consensus_res = list(H = H_consensus, W = W_consensus, V = V_consensus)
      return(consensus_res)
    })
    print(paste0('Done with ', x))
    
    res = list()
    if (return.scale.data){
      res[["scale.data"]] = Liger@scale.data
    }
    if (output.reps) {
      res[["Liger.list"]] = Liger_list
    } 
    res[["consensus.results"]] = Consensus.results
    return(res)
  })
  return(res)
}


# build phylogenetic trees based on activity program gene loadings ####
build_phylo_tree <- function(x){
  
  # add in dummy program as outgroup to root tree (for depth-first search)
  W_list = lapply(x$consensus.results, '[[', "W")
  W_list = lapply(W_list, t)
  #W_list[["root"]] <- as.matrix(data.frame(root = sample(rowMedians(W_list$R2))))
  W_list[["root"]] <- as.matrix(data.frame(root = rowMedians(W_list$R2)))
  W_list = W_list[c("root", names(x$consensus.results))]
  
  # build phylogenetic tree
  cor_mat = cor(list.cbind(W_list))
  phylo_tree = fastme.bal(1-cor_mat)
  phylo_tree = root(phylo_tree, outgroup = "root", resolve.root = T)
  phylo_tree = drop.tip(phylo_tree, "root")
  # convert negative branches to zero and filter 
  phylo_tree_pruned = phylo_tree
  phylo_tree_pruned$edge.length[phylo_tree_pruned$edge.length<0]=0
  dist_mat = cophenetic.phylo(phylo_tree_pruned)
  
  return(list(phylo_tree = phylo_tree, phylo_tree_pruned = phylo_tree_pruned, distance_matrix = dist_mat))
}

# suggest distance threshold for phylogenetic trees ####
suggest_dist_thresh <- function(phylo_trees){
  res_list = lapply(phylo_trees, function(x){
    mat = x$distance_matrix
    diag(mat) = NA
    res = list(min = quantile(colMins(mat, na.rm = T), 0.95), max = max(colMins(mat, na.rm = T)))
  })
  common_min_thresh = max(unlist(lapply(res_list, `[[`, "min")))
  common_max_thresh = min(unlist(lapply(res_list, `[[`, "max")))
  res = c(common_min_thresh, common_max_thresh)
  names(res) = c("min", "max")
  return(res)
}

# Identify "outlier" activity programs representing rare contaminating cells ####
identify_outlier_programs <- function(x, ncells = 50){
  test = lapply(x$consensus.results, `[[`, "H")
  res = lapply(test, function(y) {colMaxs(y)/apply(y, 2, function(z) {mean(sort(z, decreasing = T)[(1:ncells)+1])})})
  return(res)
}

# Partition phylogenetic trees ####
partition_phylo_tree <- function(x, y, dist.thresh = NULL, outlier.thresh = 5){
  if (is.null(dist.thresh)){
    stop("run suggest_dist_thresh")
  }
  tree = x$phylo_tree_pruned
  nodes = 1:tree$Nnode+length(tree$tip.label)
  names(nodes) = nodes
  dist_mat = x$distance_matrix
  res = lapply(nodes, function(x){
    tiplist = tips(tree, x)
    dist = median(dist_mat[tiplist,tiplist][upper.tri(dist_mat[tiplist,tiplist], diag = F)])
  })
  dist = unlist(res)
  ntips<-Ntip(tree)
  cnum <- 0 ## cluster number
  assign <- rep(0,ntips+tree$Nnode) ## cluster assignment
  igraph.tree <- graph.edgelist(tree$edge) ## tree in igraph form
  dfs <- graph.dfs(igraph.tree,root=ntips+1,neimode='out',
                   order=TRUE,dist=TRUE)
  distvec = dist
  ## transverse the tree in depth first order
  for(i in 1:length(dfs$order)){
    node <- dfs$order[i]
    ## skip leaves
    if(node < ntips+1){ next }
    ## If the node's subtree is below the threshold, mark it and
    ## its subtree as members of a new cluster
    if(distvec[node-ntips]<=dist.thresh && assign[node]<=0){
      cnum <- cnum+1
      subtree <- graph.dfs(igraph.tree,node,
                           neimode='out',unreachable=FALSE)$order
      subtree <- subtree[! is.na(subtree)]
      assign[subtree] <- cnum
    }}
  ans <- as.character(assign)
  ans <- ans[1:ntips]
  names(ans) <- tree$tip.label
  
  # identify outlier activity programs
  outliers = unlist(lapply(y, function(z) names(which(z>outlier.thresh))))
  ans[outliers] = "outlier"
  
  # set minimum subtree size to at least 2
  ans = plyr::mapvalues(ans, names(which(table(ans)<=2)), rep("0", length(names(which(table(ans)<=2)))))
  return(ans)
  
}

# Calculate the weighted number of subtrees identified at each K ####
calculate_K_metric <- function(clusters, K.max = NULL){
  if (is.null(K.max)){
    K.max = max(as.numeric(gsub("R", "", gsub("_.*", "", names(clusters)))))
  }
  df = data.frame(program = names(clusters), clust = clusters)
  df$rank = gsub("_.*", "", df$program)
  # filter subtrees containing outlier programs or fewer than 5 activity programs
  df = data.frame(program = names(clusters), clust = clusters)
  df$rank = gsub("_.*", "", df$program)
  df = subset(df, clust%ni%c("0", "outlier"))
  df = df[df$clust%in%names(which(table(df$clust)>=5)),]
  # weight subtrees by the total number of programs in that subtree
  df$clust_weight = as.numeric(plyr::mapvalues(df$clust, from = names(table(df[,c("clust")])), to = table(df[,c("clust")])/40))
  df$rank_clust = paste(df$rank, df$clust, sep = "_")
  # only count subtrees once
  df = df[-which(duplicated(df$rank_clust)),]
  res = df %>%
    group_by(rank) %>%
    summarise(weighted_n_subtrees = sum(clust_weight))
  res = data.frame(res)
  rownames(res) = res$rank
  res = res[paste0("R", 2:K.max),]
}

# suggest K based on phylogenetic partitioning ####
suggest.k = function(K_metrics){
  ind = dim(K_metrics)[1]
  k.suggest = which(round(K_metrics$weighted_n_subtrees)>=round(K_metrics$weighted_n_subtrees[ind]))+1
  k.suggest = intersect(k.suggest, which(diff(K_metrics$weighted_n_subtrees)<0)+1)[1]
}

# Collect results at optimized K and filter outliers ####
NMF_results_opt_k <- function(NMF_results, k, out_score, outlier.thresh = 5) {
  NMF_results = NMF_results$consensus.results[[paste0("R", k)]][c("W", "H")]
  out_score = out_score[[paste0("R", k)]]
  NMF_results$W = NMF_results$W[out_score<=outlier.thresh,]
  NMF_results$H = NMF_results$H[,out_score<=outlier.thresh]
  return(NMF_results)
}

# Calculate each sample's average expression score for each program ###
calc.H.score <- function(NMF_results, metadata){
  H_matrix <- NMF_results$H
  res = data.frame(matrix(nrow = length(unique(metadata$Sample)), ncol= ncol(H_matrix)))
  rownames(res) = unique(metadata$Sample)
  colnames(res) = colnames(H_matrix)
  for (i in unique(metadata$Sample)){
    res[i,] = colMeans(H_matrix[intersect(rownames(H_matrix), rownames(subset(metadata, Sample==i))),])
  }
  return(res)
}

# Build DECIPHER-seq network

Construct_network <- function(Expression_score_cor){
  mat <- Expression_score_cor$sig.cor
  diag(mat) <- 0
  mat[is.na(mat)] <- 0
  mat_pos = mat
  mat_pos[mat_pos<0] = 0
  
  network <- graph_from_adjacency_matrix(mat, weighted=T, mode="undirected", diag=F)
  E(network)$sign = E(network)$weight
  E(network)$weight = abs(E(network)$weight)
  E(network)$sign[E(network)$sign<0] = -1
  E(network)$sign[E(network)$sign>0] = 1
  
  # For visualization purposes, we create the layout from positive edges only, but do community detection (and wTO filtering in the next step) on the signed weighted graph
  network_pos <- graph_from_adjacency_matrix(mat_pos, weighted=T, mode="undirected", diag=F)
  coords = layout_with_fr(network_pos, niter = 10000)
  res = list(network = network, network_coords = coords, mat = mat, network_pos = network_pos)
}

# Perform community detection 

CPM <- function(adjacency_matrix){
  py_run_string("import leidenalg as la; import igraph as ig; import numpy as np")
  py_run_string("G = ig.Graph.Weighted_Adjacency(r.adjacency_matrix.tolist())")
  
  # sweep across a range of resolutions
  py_run_string("optimiser = la.Optimiser()")
  py_run_string("profile = optimiser.resolution_profile(G, la.CPMVertexPartition, 
          weights = 'weight', resolution_range=(0.001, 0.4), number_iterations = 0)")
  sweep = py$profile
  modularity = lapply(sweep, function(x){x$modularity})
  # Use "resolution" that gives max modularity
  partition_use = sweep[[which.max(unlist(modularity))]]
  py_run_string("partition = r.partition_use")
  py_run_string("diff = optimiser.optimise_partition(partition)")
  
  # Optimise this partition
  while(py$diff!=0){
    py_run_string("diff = optimiser.optimise_partition(partition)")
    py_run_string("print(diff)")
  }
  clustering_res = py$partition
  modules = clustering_res$membership + 1
  names(modules) = colnames(adjacency_matrix)
  return(modules)
}

# Filter isolated nodes and modules using weighted topological overlap
Filter_network = function(Network){
  mat = Network$mat
  modules = Network$modules
  topo_overlap = wTO(mat, sign = 'sign')
  diag(topo_overlap) = NA
  clust_topo_overlap = NULL
  for (clust in sort(unique(modules))){
    clust_topo_overlap[clust] = mean(topo_overlap[names(which(modules==clust)),names(which(modules==clust))][upper.tri(topo_overlap[names(which(modules==clust)),names(which(modules==clust))], diag = F)])
  }
  clust_size = table(Network$modules)
  
  # Permuation trial 1: expected mean topo overlap of all nodes within a module (by chance)
  # Use this to filter out modules with poor topological overlap
  reps = 10000
  boot_topo_overlap = matrix(nrow = length(unique(modules)), ncol = reps)
  rownames(boot_topo_overlap) = paste0("Module_", sort(unique(modules)))
  for (rep in 1:reps){
    boot_topo_overlap[,rep] = apply(clust_size, 1, function(x) {
      ind = sample(1:ncol(topo_overlap), x)
      res = mean(topo_overlap[ind,ind][upper.tri(topo_overlap[ind,ind], diag = F)])
    })
  }
  modules_keep = intersect(which(clust_topo_overlap > 
                                   apply(boot_topo_overlap, 1, function(x) {quantile(x, 0.99, na.rm = T)})),
                           which(table(modules)>=4))
  
  # Permutation trial 2: expected mean topo overlap of a node with other nodes in the same module (by chance)
  # Use this to filter out isolated nodes within a module
  reps = 10000
  node_boot_topo_overlap = matrix(nrow = length(unique(modules)), ncol = reps)
  rownames(node_boot_topo_overlap) = paste0("Module_", sort(unique(modules)))
  for (rep in 1:reps){
    node_boot_topo_overlap[,rep] = apply(clust_size, 1, function(x) {
      ind = sample(1:ncol(topo_overlap), x)
      if (x > 1) {
        res = sample(colMeans(topo_overlap[ind,ind], na.rm = T), 1)
      } else res = NA
      return(res)})
  }
  nodes_keep = NULL
  for (clust in names(table(modules)[table(modules)>=4])){
    nodes_keep = c(nodes_keep, names(which(colMeans(topo_overlap[names(which(modules==clust)),names(which(modules==clust))], na.rm = T)>
                                             quantile(node_boot_topo_overlap[paste0("Module_", clust),], 0.99, na.rm = T)))
    )
  }
  keep = intersect(which(modules%in%modules_keep), which(names(modules)%in%nodes_keep))
  keep = which(names(modules)%in%nodes_keep)
  if (min(table(modules[keep]))<4){
    keep = intersect(keep, which(modules%ni%which(table(modules[keep])<4)))
  }
  network_sub=induced_subgraph(Network$network, keep)
  network_sub_pos=induced_subgraph(Network$network_pos, keep)
  coords_sub=layout_with_fr(network_sub_pos, niter = 10000)
  modules_sub=modules[keep]
  names(modules_sub) = V(network_sub)$name
  Network$filtered_network = network_sub
  Network$filtered_network_coords = coords_sub
  Network$filtered_modules = modules_sub
  Network$wTO = topo_overlap
  return(Network)
}

# Gene loading similarity between programs ####
Gene_similarity <- function(NMF_results_atK){
  W_list <- lapply(NMF_results_atK, `[[`, "W")
  genes = unique(unlist(lapply(W_list, colnames)))
  W_list <- lapply(W_list, function(x){
    W_mat = as.data.frame(matrix(NA, nrow = length(genes), ncol = nrow(x)))
    rownames(W_mat) = genes
    colnames(W_mat) = rownames(x)
    W_mat[colnames(x),] = t(x)
    return(W_mat)
  })
  W_comb = list.cbind(W_list)
  W_cor = cor(W_comb, use = "pairwise.complete.obs")
  return(W_cor)
}

# Node p-values for gene loading similarity ####

Permutation_test_gene_cor <- function(Network, gene_correlation_matrix, reps = 10000){
  test_nodes = names(Network$filtered_modules)
  gene_correlation_matrix = gene_correlation_matrix[test_nodes,test_nodes]
  clust_size = table(Network$filtered_modules)
  permutation_res = matrix(nrow = length(clust_size), ncol = reps)
  rownames(permutation_res) = paste0("Module_", names(clust_size))
  for (rep in 1:reps){
    permutation_res[,rep] = apply(clust_size, 1, function(x) {
      ind = sample(1:ncol(gene_correlation_matrix), x)
      if (x > 1) {
        res = sample(colMeans(gene_correlation_matrix[ind,ind], na.rm = T), 1)
      } else res = NA
      return(res)})
  }
  
  node_pval <- rep(NA, ncol(gene_correlation_matrix))
  names(node_pval) <- colnames(gene_correlation_matrix)
  for (node in colnames(gene_correlation_matrix)){
    module = Network$filtered_modules[node]
    temp = mean(gene_correlation_matrix[names(which(Network$filtered_modules == module)), node], na.rm = T)
    node_pval[node] = sum(temp < permutation_res[paste0("Module_", module),])/length(permutation_res[paste0("Module_", module),])
  }
  node_pval[node_pval==0]=1/reps
  return(node_pval)
}

# Infer direct cell-cell interations ####

Infer_direct_interactions <- function(Expression_score, Network, metadata, celltypes.test = NULL, sort.gate = NULL,
                                      adjRsquared.thresh = 0.6){
  if (is.null(celltypes.test)) {
    celltypes.test <- levels(metadata$Type)
    warning("No celltypes.test input. Testing all cells in dataset.")
  } else message(paste0("Testing direct interactions between cell types: ", toString(celltypes.test)))
  if (is.null(sort.gate)) {
    warning(paste0("No sort gate input. Testing all cells of type: ", toString(celltypes.test)))
  } else message(paste0("Testing cells in sort gate: ", toString(sort.gate)))
  mat = Network$mat
  keep = names(Network$filtered_modules)
  mat_sub = mat[keep,keep]
  cell.props = as.data.frame.matrix(prop.table(table(subset(metadata, Sort%in%sort.gate)[,c("Sample","Type")])[,celltypes.test], 1))
  colnames(cell.props) = paste0("freq.", colnames(cell.props))
  Expression_score_comb = list.cbind(Expression_score)
  Expression_score_comb = cbind(Expression_score_comb, cell.props[rownames(Expression_score_comb),])
  res = list()
  for (i in rownames(mat_sub)){
    for (j in colnames(mat_sub)){
      i_type = gsub("\\.R[0-9].*", "", i)
      j_type = gsub("\\.R[0-9].*", "", j)
      if (mat_sub[i,j] > 0 & i_type %in% celltypes.test){
        Factor1 = Expression_score_comb[,i]
        Factor2 = Expression_score_comb[,j]
        freq = Expression_score_comb[,paste0("freq.", i_type)]
        res[["summary"]][[i]][[j]] = summary(lm(Factor2 ~ Factor1 * freq))
        res[["summary_noFreq"]][[i]][[j]] = summary(lm(Factor2 ~ Factor1))
        res[["pval_Fact1"]][[i]][[j]] = summary(lm(Factor2 ~ Factor1 * freq))$coefficients[2,4]
        res[["pval_freq"]][[i]][[j]] = summary(lm(Factor2 ~ Factor1 * freq))$coefficients[3,4]
        res[["pval_interaction"]][[i]][[j]] = summary(lm(Factor2 ~ Factor1 * freq))$coefficients[4,4]
        res[["coef_interaction"]][[i]][[j]] = summary(lm(Factor2 ~ Factor1 * freq))$coefficients[4,1]
        res[["pval_overall"]][[i]][[j]] = glance(summary(lm(Factor2 ~ Factor1 * freq)))$p.value
        res[["adjR"]][[i]][[j]] = summary(lm(Factor2 ~ Factor1 * freq))$adj.r.squared
        res[["adjR_noFreq"]][[i]][[j]] = summary(lm(Factor2 ~ Factor1))$adj.r.squared
      } else {
        res[["summary"]][[i]][[j]] = res[["summary_noFreq"]][[i]][[j]] = res[["pval_Fact1"]][[i]][[j]] = res[["pval_freq"]][[i]][[j]] = 
          res[["pval_interaction"]][[i]][[j]] =  res[["pval_overall"]][[i]][[j]] = res[["coef_interaction"]][[i]][[j]] =
          res[["adjR"]][[i]][[j]] = res[["adjR_noFreq"]][[i]][[j]] = NA
      }
    }
  }
  
  ind = which(unlist(res$pval_interaction)<0.01 & 
                unlist(res$coef_interaction)>0 & 
                unlist(res$pval_freq)>0.05 & 
                unlist(res$pval_Fact1)>0.05 &
                unlist(res$adjR) > adjRsquared.thresh &
                p.adjust(unlist(res$pval_overall), method = "fdr")<0.01)
  
  res_mat = matrix(NA, ncol = ncol(mat_sub), nrow = nrow(mat_sub))
  res_mat[ind] = unlist(res$adjR)[ind]
  
  res = list(results = res, inds = ind, adjacency_matrix = res_mat)
  return(res)
}

# Marker genes  - Identify genes associated with each activity program (also see Kotliar et al.) ####

Marker_gene_analysis <- function(NMF_results_atK, NMF_results){
  W_consensus_k_use = lapply(NMF_results_atK, "[[", "W")
  H_consensus_k_use = lapply(NMF_results_atK, "[[", "H")
  scale.data <- lapply(NMF_results, `[[`, "scale.data")
  W_mat = mcmapply(function(W, H, data) {
    Liger.data = t(list.rbind(data))
    mat=data.frame(matrix(nrow=nrow(Liger.data), ncol=ncol(H)))
    rownames(mat)=rownames(Liger.data)
    colnames(mat)=colnames(H)
    for(gene_ind in rownames(mat)){
      cat('.')
      res = unlist(lapply(1:ncol(mat), function(j){
        summary(lm(Liger.data[gene_ind,] ~ H[colnames(Liger.data),j]))$coef[2,1]
      }))
      mat[gene_ind,] = res
    }
    return(mat)
  }, W = W_consensus_k_use, H = H_consensus_k_use, data = scale.data, SIMPLIFY = F, mc.cores = detectCores())
}

combine_marker_matrices <- function(marker_gene_matrix){
  W_mat = marker_gene_matrix
  genes = unique(unlist(lapply(W_mat, rownames)))
  W_mat_k_use = lapply(W_mat, function(x) {
    res = as.data.frame(matrix(NA, nrow = length(genes), ncol = ncol(x)))
    rownames(res) = genes
    colnames(res) = colnames(x)
    res[rownames(x),] = x
    return(res)
  })
  W_mat_k_use = list.cbind(W_mat_k_use)
}

# Run fgsea, combine results, and calculate fdr ####

fgsea_test <- function(marker_gene_list, Network, path_list){
  res = mclapply(marker_gene_list, function(marker_gene_matrix){
    res_list = apply(marker_gene_matrix, 2, function(x){
      Program_test = as.numeric(x)
      names(Program_test) = rownames(marker_gene_matrix)
      res = fgsea(pathways = path_list, stats = Program_test, 
                  minSize=10, maxSize=Inf, eps=1e-50, nPermSimple=100000)
    })
    res_mat = list.rbind(res_list)
    n_pathways = length(path_list)
    res_mat$Program = rep(colnames(marker_gene_matrix), each = n_pathways)
    return(res_mat)
  }, mc.cores = min(detectCores(), length(marker_gene_list)))
  for (i in names(res)){
    res[[i]]$Type = i
  }
  res = list.rbind(res)
  res$module = paste0(res$Type, "." , res$Program)
  res$module = plyr::mapvalues(res$module, from = names(Network$modules), to = Network$modules)
  res = subset(res, module%in%unique(Network$filtered_modules))
  res$fdr = p.adjust(res$pval, method = "fdr")
  return(res)
}

# Permutation test to calculate fdr of gene set enrichment within a module ####
Get_enrichment_pvals <- function(sets_to_test, Network, nreps = 10000){
  mat = Network$mat[names(Network$filtered_modules), names(Network$filtered_modules)]
  
  # only test for enrichment on nodes within filtered network
  fgsea_test = sets_to_test[paste0(sets_to_test$Type, ".", sets_to_test$Program)%in%names(Network$filtered_modules),]
  
  # only test for enrichment for gene sets with at least five significantly enriched nodes (fdr < 0.01)
  fgsea_test = subset(fgsea_test, pathway %in% names(which(table(subset(fgsea_test, fdr < 0.01&NES>0)$pathway)>=5))) 

  message(paste0("Testing ", length(unique(fgsea_test$pathway)), " gene sets"))
  
  # only test for positive enrichment
  fgsea_test$fdr[fgsea_test$NES<0] = 1
  fgsea_test$fdr[is.na(fgsea_test$fdr)] = 1
  
  Test = unique(fgsea_test$pathway)
  names(Test) = Test
  clust_size = table(Network$filtered_modules)
  names(clust_size) = paste0("Module_", names(clust_size))
  node_boot_fgsea_fdr = lapply(Test, function(x){
    Test_fgsea_test = subset(fgsea_test, pathway==x)
    Test_fdr = Test_fgsea_test$fdr
    names(Test_fdr) = Test_fgsea_test$node_id
    node_boot_fgsea_fdr = mclapply(1:nreps, function(y) {
      apply(clust_size, 1, function(z) {
        ind = sample(1:ncol(mat), z)
        res = sum(Test_fdr[ind] < 0.01, na.rm = T)
        return(res)})
    })
    node_boot_fgsea_fdr = list.cbind(node_boot_fgsea_fdr)
    return(node_boot_fgsea_fdr)
  })
  
  enrichment_pval = list()
  for (test in names(node_boot_fgsea_fdr)){
    Test_fgsea_test = subset(fgsea_test, pathway==test)
    Test_fdr = Test_fgsea_test$fdr
    names(Test_fdr) = paste0(Test_fgsea_test$Type, ".", Test_fgsea_test$Program)
    for (clusters in unique(Network$filtered_modules)){
      enrichment_pval[[test]][[paste0("Module_", clusters)]] =  1- sum(sum(Test_fdr[names(which(Network$filtered_modules==clusters))]<0.01, na.rm = T) > node_boot_fgsea_fdr[[test]][paste0("Module_", clusters),], na.rm = T)/length(node_boot_fgsea_fdr[[test]][paste0("Module_", clusters),])
    }
  }
  enrichment_pval = rapply(enrichment_pval, f = function(x){
    if (x == 0){x = 1/(nreps+1)}
    return(x) }, how = "replace")
}

# Measure effect size and p values for association with metadata features ####
Calculate_metadata_associations <- function(Network, Expression_score, metadata, feature.to.test = NULL, type = c("binary", "continuous")){
  if (is.null(feature.to.test)){
    stop("Choose metadata feature to test")
  } else if (feature.to.test %ni% colnames(metadata)){
    stop("feature.to.test must be a column of the input metadata file")
  }
  H.score.comb = list.cbind(Expression_score)
  
  H.score.comb$feature.to.test = plyr::mapvalues(rownames(H.score.comb), from = metadata$Sample, to = metadata[,feature.to.test], warn_missing = F)
  
  res = data.frame(matrix(nrow = 3, ncol = length(Network$filtered_modules)))
  colnames(res) = names(Network$filtered_modules)
  rownames(res) = c("effect_size", "effect_size_sign", "p_value")
  
  # Binary variable
  if (type == "binary"){
    if (length(levels(factor(H.score.comb$feature.to.test)))>2){
      stop("Create new metadata column with two levels")
    }
    
    test_conditions = levels(factor(metadata[,feature.to.test]))
    group1 = unique(metadata[which(metadata[,feature.to.test]==test_conditions[1]),]$Sample)
    group2 = unique(metadata[which(metadata[,feature.to.test]==test_conditions[2]),]$Sample)
    
    for (i in names(Network$filtered_modules)){
      res[c("effect_size"),i] = as.numeric(wilcox_effsize(H.score.comb, as.formula(paste0(i, "~", "feature.to.test")))$effsize)
      res[c("effect_size_sign", "p_value"),i] = as.numeric(wilcox.test(as.formula(paste0(i, "~", "feature.to.test")), data = H.score.comb,alternative = "two.sided", conf.int = T)[c("estimate", "p.value")])
    }
    res["effect_size_sign",][res["effect_size_sign",]>0] = 1
    res["effect_size_sign",][res["effect_size_sign",]<0] = -1
  } else if (type == "continuous"){
    
    # continuous variable
    if (class(H.score.comb$feature.to.test)!="numeric"){
      warning("Converting metadata feature.to.test to continuous variable")
      H.score.comb$feature.to.test = as.numeric(H.score.comb$feature.to.test)
    } 
    for (i in names(Network$filtered_modules)){
      res[c("effect_size"),i] = cor(H.score.comb[,i], H.score.comb$feature.to.test, use = "pairwise.complete.obs")
      res[c("p_value"),i] = summary(lm(H.score.comb[,i] ~ H.score.comb$feature.to.test))$coefficients[2,4]
    }
    res["effect_size_sign",][res["effect_size",]>0] = 1
    res["effect_size_sign",][res["effect_size",]<0] = -1
    res[c("effect_size"),] = abs(res[c("effect_size"),])
  }
  
  return(res)
}

##### Basic functions #################

solveNNLS <- function(C, B) {
  .Call('_rliger_solveNNLS', PACKAGE = 'rliger', C, B)
}

# wrapper function for running online_iNMF
LIGER.run <- function(K, Liger, minibatch_size, h5_chunk_size, seed = NULL){
  # sample seeds and report for reproducibility
  if (is.null(seed)) {seed = sample(1:1000000, 1)}
  Liger <- online_iNMF(Liger, k = K, max.epochs=10, seed = seed, miniBatch_size = minibatch_size, h5_chunk_size = h5_chunk_size, verbose = F)
  cat(".")
  W_res = t(Liger@W)
  colnames(W_res) = paste0("R", K, "_Program", 1:K)
  H_res = list.rbind(Liger@H)
  colnames(H_res) = paste0("R", K, "_Program", 1:K)
  V_res = Liger@V
  V_res = lapply(V_res, function(y){
    rownames(y) = paste0("R", K, "_Program", 1:K)
    y = t(y)
    return(y)
  })
  params = c(minibatch_size = minibatch_size, h5_chunk_size = h5_chunk_size)
  return(list(W = W_res, H = H_res, V = V_res))
}

'%ni%' <- Negate('%in%')
