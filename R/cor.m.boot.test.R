cor.m.boot.test <- function (mat, ...) {
  require(parallel)
  mat <- as.matrix(mat)
  n_pairs <- ncol(mat) * (ncol(mat)-1) / 2
  i_vec <- rep(0, n_pairs)
  j_vec <- rep(0, n_pairs)
  idx <- 1
  for(i in seq_len(ncol(mat))) {
    for(j in seq_len(i-1)) {
      i_vec[idx] <- i
      j_vec[idx] <- j
      idx <- idx + 1
    }
  }
  res = mcmapply(i_vec, j_vec, FUN = function(i,j) {
    tmp = boot.cor.bca.complete(x = mat[, i], y = mat[, j], null.hyp = 0, alternative = c("two.sided"), 
                                conf.level = 0.95, type = NULL, R = 10000, ncpus = 1)
    res = list(cor = tmp$Observed, p = tmp$p.value, lowCI = tmp$Confidence.limits[1], uppCI = tmp$Confidence.limits[2])
    }, mc.cores = detectCores())
  cor.mat = p.mat = lowCI.mat = uppCI.mat = matrix(rep(NA, ncol(mat)^2), ncol = ncol(mat))
  dimnames(cor.mat) <- dimnames(p.mat) <- dimnames(lowCI.mat) <- dimnames(uppCI.mat) <- list(colnames(mat), colnames(mat))
  for(n in seq_along(i_vec)) {
    cor.mat[i_vec[n], j_vec[n]] <- cor.mat[j_vec[n], i_vec[n]] <- unlist(res["cor",n])
    p.mat[i_vec[n], j_vec[n]] <- p.mat[j_vec[n], i_vec[n]] <- unlist(res["p",n])
    lowCI.mat[i_vec[n], j_vec[n]] <- lowCI.mat[j_vec[n], i_vec[n]] <- unlist(res["lowCI",n])
    uppCI.mat[i_vec[n], j_vec[n]] <- uppCI.mat[j_vec[n], i_vec[n]] <- unlist(res["uppCI",n])
  }
  results = list(cor = cor.mat, p = p.mat, lowCI = lowCI.mat, uppCI = uppCI.mat)
  results
}