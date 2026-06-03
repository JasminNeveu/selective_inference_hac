get.merged.clusters <- function(hcl,k){
  n <- length(hcl$order)
  merge_step <- n-(k-1)

  c1_internal <- hcl$merge[merge_step, 1]
  c2_internal <- hcl$merge[merge_step, 2]

  clusters_k <- cutree(hcl, k = k)

  get_leaf_index <- function(id) {
    if (id < 0) return(-id)
    get_leaf_index(hcl$merge[id, 1])
  }

  leaf1 <- get_leaf_index(c1_internal)
  leaf2 <- get_leaf_index(c2_internal)

  c1 <- clusters_k[leaf1]
  c2 <- clusters_k[leaf2]
  c(unname(c1),unname(c2))
}


compute.pvals <- function(hcl, X, U, Sigma, Y, UY, precUY, linkage, ndraws, sample_split, nY, return_Sigma, return_X_clus, dismat, kmax,parallel_config) {
  k_seq <- 2:kmax

  run_one_k <- function(k) {
    clusters <- get.merged.clusters(hcl, k)
    hc_test  <- PCIdep::test.clusters.hc(
      X = X, U = U, Sigma = Sigma, Y = Y,
      UY = UY, precUY = precUY, linkage = linkage,
      NC = k, clusters = clusters, ndraws = ndraws,
      sample_split = sample_split, nY = nY,
      return_Sigma = return_Sigma, return_X_clus = return_X_clus,
      hcl = hcl, dismat = dismat
    )

    hc_test$pval <- ifelse(hc_test$pval < 2.2e-16, 2.2e-16, hc_test$pval)
    print(paste("k =",k,"clusters:",clusters[1],"-",clusters[2],"---", hc_test$pval))
    hc_test$pval

  }

  if (is.null(parallel_config) || !parallel_config$enabled) {
    return(setNames(sapply(k_seq, run_one_k), k_seq))
  }

  workers <- parallel_config$workers %||% (parallel::detectCores() - 1)
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  pvals <- unlist(furrr::future_map(k_seq, run_one_k))
  stats::setNames(pvals, k_seq)
}

get.data <- function(hcl,X, U,Sigma , Y , UY , precUY  , linkage,ndraws, sample_split, nY, return_Sigma, return_X_clus, dismat,kmax,parallel_config){
  hcldata <- ggdendro::dendro_data(hcl, type = "rectangle")
  hcldata$segments <- hcldata$segments[order(-hcldata$segments$y), ] # ordonne les segments pour qu'ils soient lister de haut en bas de l'arbre
  n <- nrow(hcldata$segments)
  pvec <- compute.pvals(hcl,X, U,Sigma , Y , UY , precUY  , linkage,ndraws, sample_split, nY, return_Sigma, return_X_clus, dismat,kmax,parallel_config)
  pvec_rep <- rep(pvec, each = 4)
  hcldata$segments$pval <- NA
  hcldata$segments$pval[1:length(pvec_rep)] <- pvec_rep

  # ------Add labels pvalues ------
  merge_points <- hcldata$segments %>% #
    filter(y == yend) %>%
    group_by(y) %>%
    summarise(x_min = min(x), x_max = max(x)) %>%
    arrange(desc(y)) %>%
    mutate(
      x_mid = (x_min + x_max) / 2,
      label = c(signif(pvec,3), rep(NA, n() - length(pvec)))
    )
  hcldata$merge_points <- merge_points
  return(hcldata)
}




