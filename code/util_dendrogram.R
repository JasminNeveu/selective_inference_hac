# TODO: finir documentation

#' @title Find label of clusters merged
#' @description
#' Determines the two clusters merged at the step \code{n-(k-1)} of the HAC algorithm.
#'
#' @param  hierarchical clustering object of class \code{"hclust"}, as returned by
#'   \code{fastcluster::hclust()} or \code{stats::hclust()} on the squared Euclidean distance matrix of \code{X}.
#'
#' @param k Integer. Number of clusters.
#'
#' @return A list of two integers
#'
#' @example

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


#' @title Compute p-values
#' @description
#'  Compute the p-values until the \code{(n-(kmax -1))} step of the HAC algorithm.
#'
#'
#' @param X An \eqn{n \times p} data matrix assumed to arise from a matrix normal distribution
#'   \eqn{\mathcal{MN}(M, U, \Sigma)}.
#' @param kmax Integer. Number of p-values to compute.
#' @param U An \eqn{n \times n} positive-definite matrix describing the dependence structure between the rows of \code{X}. If \code{NULL},
#'   observations are assumed to be independent and \code{U} is set to the \eqn{n \times n} identity matrix.
#' @param Sigma A \eqn{p \times p} positive-definite matrix describing the dependence structure between the columns of \code{X}. If \code{NULL},
#'   \code{Sigma} is over-estimated from an auxiliary independent sample \code{Y} (in the sense of the Loewner partial order).
#' @param Y If \code{Sigma} is \code{NULL}, an independent copy of \code{X} used to estimate \code{Sigma}. It must have the same number of columns as
#'   \code{X}.
#' @param UY If \code{Sigma} is \code{NULL}, an \eqn{n_Y \times n_Y} positive-definite matrix describing the dependence structure between the
#'   rows of \code{Y}. If \code{NULL} and \code{precUY} is not provided, the identity matrix is used by default.
#' @param precUY The inverse of \code{UY}. Supplying \code{precUY} may improve computational efficiency. If \code{UY} is provided but \code{precUY} is
#'   \code{NULL}, \code{precUY} is computed internally by matrix inversion.
#' @param linkage Character string specifying the linkage criterion used in hierarchical clustering. Must be one of \code{"single"},
#'   \code{"average"}, \code{"centroid"}, \code{"ward.D"}, \code{"median"}, \code{"mcquitty"}, or \code{"complete"}.
#'   When \code{hcl} is provided, the linkage is inferred from \code{hcl$method} and this argument is overridden.
#'   If \code{linkage} is supplied explicitly alongside \code{hcl} but does not match \code{hcl$method}, a warning
#'   is emitted and the method stored in \code{hcl} takes precedence.
#' @param hcl An optional precomputed hierarchical clustering object of class \code{"hclust"}, as returned by
#'   \code{fastcluster::hclust()} or \code{stats::hclust()} on the squared Euclidean distance matrix of \code{X}.
#'   When supplied, the clustering step is skipped and the linkage criterion is inferred from \code{hcl$method},
#'   overriding the \code{linkage} argument. Passing a precomputed \code{hcl} is useful when testing several
#'   cluster pairs from the same clustering, as it avoids recomputing the dendrogram each time.
#'   The number of leaves in \code{hcl} must equal \code{nrow(X)}; an error is raised otherwise.
#'   Ignored (with a warning) when \code{sample_split = TRUE}, because sample splitting changes \code{X}
#'   after the clustering would have been computed.
#' @param dismat An optional precomputed squared Euclidean distance object of class \code{"dist"}, as returned by
#'   \code{stats::dist(X, method = "euclidean")^2}. When supplied alongside \code{hcl}, both the distance object
#'   and the dendrogram computations are skipped, which is useful when testing several cluster pairs from the same
#'   clustering. When supplied without \code{hcl}, a warning is emitted and the dendrogram is computed from
#'   \code{dismat}. When \code{hcl} is supplied without \code{dismat} and the linkage is not \code{"complete"},
#'   a warning is emitted and \code{dismat} is recomputed internally. Ignored (with a warning) when
#'   \code{sample_split = TRUE}, because sample splitting changes \code{X} after the distance object would
#'   have been computed.
#' @param ndraws Integer. Number of Monte Carlo samples used to approximate the p-value when \code{linkage = "complete"}. Ignored otherwise.
#' @param sample_split Logical. Whether to use sample splitting to estimate \code{Sigma} when \code{Sigma = NULL}. Ignored when \code{Sigma} is provided by the user.
#' @param nY Integer. If \code{Y} is not provided and \code{sample_split = TRUE}, the number of rows of the auxiliary sample \code{Y} used to estimate \code{Sigma}. If \code{nY} is \code{NULL}, half of the rows of \code{X} are used for estimation. Ignored when \code{Sigma} is provided by the user.
#' @param return_Sigma Logical. Whether to include the column covariance matrix used in the test in the returned list. Ignored when \code{Sigma} is provided by the user. Default is \code{FALSE}.
#' @param return_X_clus Logical. If sample splitting is performed to estimate \code{Sigma}, whether to include the data matrix used for clustering in the returned list. Ignored when \code{sample_split = FALSE} (as the same data matrix is used for clustering and testing). If further analysis of the retrieved clusters is desired, we recommend setting \code{return_X_clus = TRUE} when \code{sample_split = TRUE} to avoid confusion. Default is \code{FALSE}.
#'
#'
#' @return A list of size \code{nrow(X) -1 } of p-values.


compute.pvals <- function(hcl,X, U,Sigma ,Y ,UY , precUY ,linkage,ndraws, sample_split,nY,return_Sigma,return_X_clus,dismat,kmax) {
  pvec <- numeric(kmax-1)
  for (k in 2:kmax) {
    i <- k-1
    clusters <- get.merged.clusters(hcl, k)
    hc_test <- PCIdep::test.clusters.hc(
      X = X,
      U = U,
      Sigma = Sigma,
      Y = Y,
      UY = UY,
      precUY = precUY,
      linkage = linkage,
      NC = k,
      clusters =clusters,
      ndraws = ndraws,
      sample_split = sample_split,
      nY = nY,
      return_Sigma = return_Sigma,
      return_X_clus = return_X_clus,
      hcl= hcl,
      dismat= dismat
    )
    pvec[i] <- hc_test$pval
  }
  names(pvec) <- 2:kmax
  pvec
}


#' @title Build dataframe
#' @description
#' Build dataframe for \code{plot.pdendrogram} function.
#'
#' @param X An \eqn{n \times p} data matrix assumed to arise from a matrix normal distribution
#'   \eqn{\mathcal{MN}(M, U, \Sigma)}.
#' @param U An \eqn{n \times n} positive-definite matrix describing the dependence structure between the rows of \code{X}. If \code{NULL},
#'   observations are assumed to be independent and \code{U} is set to the \eqn{n \times n} identity matrix.
#' @param Sigma A \eqn{p \times p} positive-definite matrix describing the dependence structure between the columns of \code{X}. If \code{NULL},
#'   \code{Sigma} is over-estimated from an auxiliary independent sample \code{Y} (in the sense of the Loewner partial order).
#' @param Y If \code{Sigma} is \code{NULL}, an independent copy of \code{X} used to estimate \code{Sigma}. It must have the same number of columns as
#'   \code{X}.
#' @param UY If \code{Sigma} is \code{NULL}, an \eqn{n_Y \times n_Y} positive-definite matrix describing the dependence structure between the
#'   rows of \code{Y}. If \code{NULL} and \code{precUY} is not provided, the identity matrix is used by default.
#' @param precUY The inverse of \code{UY}. Supplying \code{precUY} may improve computational efficiency. If \code{UY} is provided but \code{precUY} is
#'   \code{NULL}, \code{precUY} is computed internally by matrix inversion.
#' @param linkage Character string specifying the linkage criterion used in hierarchical clustering. Must be one of \code{"single"},
#'   \code{"average"}, \code{"centroid"}, \code{"ward.D"}, \code{"median"}, \code{"mcquitty"}, or \code{"complete"}.
#'   When \code{hcl} is provided, the linkage is inferred from \code{hcl$method} and this argument is overridden.
#'   If \code{linkage} is supplied explicitly alongside \code{hcl} but does not match \code{hcl$method}, a warning
#'   is emitted and the method stored in \code{hcl} takes precedence.
#' @param hcl An optional precomputed hierarchical clustering object of class \code{"hclust"}, as returned by
#'   \code{fastcluster::hclust()} or \code{stats::hclust()} on the squared Euclidean distance matrix of \code{X}.
#'   When supplied, the clustering step is skipped and the linkage criterion is inferred from \code{hcl$method},
#'   overriding the \code{linkage} argument. Passing a precomputed \code{hcl} is useful when testing several
#'   cluster pairs from the same clustering, as it avoids recomputing the dendrogram each time.
#'   The number of leaves in \code{hcl} must equal \code{nrow(X)}; an error is raised otherwise.
#'   Ignored (with a warning) when \code{sample_split = TRUE}, because sample splitting changes \code{X}
#'   after the clustering would have been computed.
#' @param dismat An optional precomputed squared Euclidean distance object of class \code{"dist"}, as returned by
#'   \code{stats::dist(X, method = "euclidean")^2}. When supplied alongside \code{hcl}, both the distance object
#'   and the dendrogram computations are skipped, which is useful when testing several cluster pairs from the same
#'   clustering. When supplied without \code{hcl}, a warning is emitted and the dendrogram is computed from
#'   \code{dismat}. When \code{hcl} is supplied without \code{dismat} and the linkage is not \code{"complete"},
#'   a warning is emitted and \code{dismat} is recomputed internally. Ignored (with a warning) when
#'   \code{sample_split = TRUE}, because sample splitting changes \code{X} after the distance object would
#'   have been computed.
#' @param ndraws Integer. Number of Monte Carlo samples used to approximate the p-value when \code{linkage = "complete"}. Ignored otherwise.
#' @param sample_split Logical. Whether to use sample splitting to estimate \code{Sigma} when \code{Sigma = NULL}. Ignored when \code{Sigma} is provided by the user.
#' @param nY Integer. If \code{Y} is not provided and \code{sample_split = TRUE}, the number of rows of the auxiliary sample \code{Y} used to estimate \code{Sigma}. If \code{nY} is \code{NULL}, half of the rows of \code{X} are used for estimation. Ignored when \code{Sigma} is provided by the user.
#' @param return_Sigma Logical. Whether to include the column covariance matrix used in the test in the returned list. Ignored when \code{Sigma} is provided by the user. Default is \code{FALSE}.
#' @param return_X_clus Logical. If sample splitting is performed to estimate \code{Sigma}, whether to include the data matrix used for clustering in the returned list. Ignored when \code{sample_split = FALSE} (as the same data matrix is used for clustering and testing). If further analysis of the retrieved clusters is desired, we recommend setting \code{return_X_clus = TRUE} when \code{sample_split = TRUE} to avoid confusion. Default is \code{FALSE}.
#'
#'
#' @return A list with components
#' \describe{
#' \item{labels} Label data
#' \item{segments} Line segment data
#' \describe{
#' \item{pval} p-values for each step of the HAC algorithm
#' \item{merge_points} location to label p-values
#'  }
#' }
#'
#' @example

get.data <- function(hcl,X, U,Sigma , Y , UY , precUY  , linkage,ndraws, sample_split, nY, return_Sigma, return_X_clus, dismat,kmax){
  hcldata <- ggdendro::dendro_data(hcl, type = "rectangle")
  hcldata$segments <- hcldata$segments[order(-hcldata$segments$y), ] # ordonne les segments pour qu'ils soient lister de haut en bas de l'arbre
  n <- nrow(hcldata$segments)
  pvec <- compute.pvals(hcl,X, U,Sigma , Y , UY , precUY  , linkage,ndraws, sample_split, nY, return_Sigma, return_X_clus, dismat,kmax)
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
      label = c(round(pvec, 6), rep(NA, n() - length(pvec)))
    )
  hcldata$merge_points <- merge_points
  return(hcldata)
}

