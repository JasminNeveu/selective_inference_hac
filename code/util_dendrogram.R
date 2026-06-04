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

# TODO: nb minimal pts dans chaque clusters genre 30/50 pts.
compute.pvals <- function(hcl, X, U, Sigma, Y, UY, precUY, linkage, ndraws, sample_split, nY, return_Sigma, return_X_clus, dismat, kmax,parallel_config,min_pts = 30) {
  k_seq <- 2:kmax

  run_one_k <- function(k) {
    clusters <- get.merged.clusters(hcl, k)
    #if(){
    #warning(....)
    #return (NA)
    #}
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
  merge_points <- hcldata$segments %>%
    filter(y == yend) %>%
    group_by(y) %>%
    #.groups ?
    summarise(x_min = min(x), x_max = max(x)) %>%
    arrange(desc(y)) %>%
    mutate(
      x_mid = (x_min + x_max) / 2,
      k = row_number()+1,
      label = c(signif(pvec,3), rep(NA, n() - length(pvec)))
    )
  hcldata$merge_points <- merge_points
  return(hcldata)
}

get.individuals.merged.clusters <- function(hcl, k) {
  n <- length(hcl$order) 
  if (length(k) != 1 || is.na(k)) {
    stop("k must be a single integer.")
  }
  if (k < 1 || k > n) {
    stop("k must be an integer between 1 and n - 1.")
  }
  merge_step <- n - (k-1)
  c1_internal <- hcl$merge[merge_step, 1] 
  c2_internal <- hcl$merge[merge_step, 2]
  get_individuals <- function(id) {
    if (id < 0) { # cas d'un indiv
      return(-id)
    } else { # cas d'un cluster
      left  <- get_individuals(hcl$merge[id, 1])
      right <- get_individuals(hcl$merge[id, 2])
      return(c(left, right))
    }
  }
  
  individuals_1 <- get_individuals(c1_internal)
  individuals_2 <- get_individuals(c2_internal)
  
  return(list(unname(individuals_1), unname(individuals_2)))
}

get.individuals <- function(hcl,k){
  indiv <- get.individuals.merged.clusters(hcl,k+1)
  c(indiv[[1]],indiv[[2]])
}

labels.groups <- function(data_pvalues,hcl,groups, groups_labels) {

  clusters_list <- lapply(data_pvalues$merge_points$k, function(k) {
    cl <- get.individuals.merged.clusters(hcl, k)
    unique(c(cl[[1]], cl[[2]]))
  })
  cluster_sizes <- vapply(clusters_list, length, numeric(1))

  for (i in seq_along(groups)) {

    group <- groups[[i]]
    label_grp <- groups_labels[i]

    tmp <- lapply(clusters_list, function(indiv) {

      tab <- table(group[indiv])

      list(
        class = names(tab)[which.max(tab)],
        prop  = max(tab) / sum(tab)
      )
    })

    data_pvalues$merge_points[[paste0(label_grp, "_class")]] <-
      vapply(tmp, `[[`, character(1), "class")

    data_pvalues$merge_points[[paste0(label_grp, "_prop")]] <-
      vapply(tmp, `[[`, numeric(1), "prop")
  }
  data_pvalues$merge_points$cluster_size <- cluster_sizes
  return(data_pvalues)
}
