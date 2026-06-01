library(ggplot2)
library(MASS)
library(dplyr)
library(ggdendro)
library(fastcluster)
library(PCIdep)
library(scales)

#' @title p-Dendrogram for hierarchical clustering
#' @description
#' Generates a dendrogram with p-values derived from test of equality of means between two clusters.
#' The p-values are computed using the `test.clusters.hc` function from the library `PCIdep`.
#'
#'
#' @param X An \eqn{n \times p} data matrix assumed to arise from a matrix normal distribution
#'   \eqn{\mathcal{MN}(M, U, \Sigma)}.
#' @param kmax Integer. Number of computed p-values. Default is \code{nrow(X) - 1}.
#' @param labels_pvalues Logical. Whether to label p-values in the p-dendrogram. Default is \code{FALSE}.
#' @param treshold Integer. Threshold for signifiance of the p-value. Default is \code{0.05}.
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
#' @param low,high Colours for low and high ends of the gradient.
#' @param mid Colour for mid point
#' @param return_pvalues Logical. Wether to return the list of p-values computed. Default is \code{FALSE}.
#'
#'
#'
#' @return A ggplot object representing the p-dendrogram.
#'
#' @examples
#' n <- 100
#' n1 <- floor(n / 4)
#' n2 <- floor(n / 4)
#' n3 <- floor(n/4)
#' n4 <- n - n1 - n2 -n3 
#'
#' mu1 <- c(0, 0)
#' mu2 <- c(7.5, 2)
#' mu3 <- c(0, 6)
#' mu4 <- c(6,8)
#'
#' sigma <- diag(2)
#'
#' gen_data <- function(){
#'   X1 <- mvrnorm(n1, mu = mu1, Sigma = sigma)
#'   X2 <- mvrnorm(n2, mu = mu2, Sigma = sigma)
#'   X3 <- mvrnorm(n3, mu = mu3, Sigma = sigma)
#'   X4 <- mvrnorm(n4, mu = mu4, Sigma = sigma)
#'
#'   X <- rbind(X1, X2, X3,X4)
#'   X
#' }
#'
#' X <- gen_data()
#' Y <- gen_data()
#' labels <- factor(c(rep(1, n1), rep(2, n2), rep(3, n3),rep(4,n4)))
#' df <- data.frame(x = X[,1], y = X[,2], cluster = labels)
#'
#'
#' dismat <- dist(df %>% select(-cluster),method="euclidean")^2
#' hcl <- hclust(dismat,method="ward.D")
#' p_dendro <- plot.pdendrogram(hcl=hcl,X=X,Y=Y,linkage="ward.D",dismat=dismat,kmax=10,labels_pvalues = TRUE,)
#' p_dendro

# TODO: gerer return_pvalues argument...
plot.pdendrogram <- function(X, U = NULL, Sigma = NULL, Y = NULL, UY = NULL, precUY = NULL, linkage = 'ward.D', hcl = NULL, ndraws = 2000, sample_split = FALSE, nY = NULL, return_Sigma = FALSE, return_X_clus = FALSE, dismat = NULL,kmax=nrow(X)-1, labels_pvalues=FALSE,treshold = 0.05,low=NULL,mid=NULL,high=NULL,return_pvalues = FALSE){
  n <- nrow(X)
  if(kmax < 1 || kmax > n - 1){
    stop("kmax must be an integer between 1 and n - 1.")
  }
  if(treshold < 0 || treshold > 1){
    warning("treshold must be an integer between 0 and 1.")
  }
  filled <- c(!is.null(low), !is.null(mid), !is.null(high))
  partial_defined <- sum(filled) > 0 && sum(filled) < 3
  if(partial_defined){ # si low ou mid ou high est renseigné, alors tous doivent l'être sinon c'est mal défini
    stop("Arguments 'low', 'mid' and 'high' must be provided together (either all three or none).") 
  }

  if(!is.null(low) && !is.null(mid) && !is.null(high)){
    colours = c(low,mid,high)
  }else{
    colours = c("#B2182B","#D9C27A","#2166AC") #option 3
  }

  data <- get.data(hcl,X, U,Sigma , Y , UY , precUY  , linkage,ndraws, sample_split, nY, return_Sigma, return_X_clus, dismat,kmax)

  p_floor <- 2.2e-16
  breaks <- c(p_floor,treshold, 1)
  labels <- c(expression("<10"^-16), as.character(treshold), "1")
  vals <- rescale(log10(c(p_floor,treshold, 1)))

  p <- ggplot() +
    geom_segment(
      data = data$segments,
      aes(x = x, y = y, xend = xend, yend = yend, colour = pval) 
    ) + 
    labs(title="p-dendrogram",x = "Observations",y="Merging distances") + 
    scale_color_gradientn(
      # colours = c("#B2182B","#F4D35E","#2166AC"), #option 1
      # colours = c("#C51B7D","#E9D8A6","#008B8B"), #option 2
      # colours = c("#B2182B","#D9C27A","#2166AC"), #option 3
      # colours = c("#D55E00","#F0E442","#0072B2"), #option 4,
      colours = colours,
      values  = vals,
      trans   = "log10",
      limits  = c(p_floor, 1),
      breaks  = breaks,
      labels  = labels,
      name    = "p-value",
      na.value = "grey"
    ) +
    theme_minimal()

  if(labels_pvalues){
    p <- p + 
    geom_text(
      data  = data$merge_points %>% filter(label < treshold),
      aes(x = x_mid, y = y, label = label),
      vjust = 1.5,
      size  = 3,
    ) 
  }
  # if(return_pvalues){return_list$pvalues <- data$segments$pval <- NA} faut créer un objet return_list mais c'est chiant après dans l'usage...
  p
}


