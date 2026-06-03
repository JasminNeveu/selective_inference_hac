library(ggplot2)
library(MASS)
library(dplyr)
library(ggdendro)
library(fastcluster)
library(PCIdep)
library(scales)

build.data <- function(X,kmax = nrow(X), U = NULL, Sigma = NULL, Y = NULL, UY = NULL, precUY = NULL, linkage = 'ward.D', hcl = NULL, ndraws = 2000, sample_split = FALSE, nY = NULL, return_Sigma = FALSE, return_X_clus = FALSE, dismat = NULL,parallel_config=NULL){
  n <- nrow(X)
  if(kmax < 1 || kmax > n - 1){
    stop("kmax must be an integer between 1 and n - 1.")
  }
  data_pvalues <- get.data(hcl,X, U,Sigma , Y , UY , precUY  , linkage,ndraws, sample_split, nY, return_Sigma, return_X_clus, dismat,kmax,parallel_config)
  data_pvalues
}


plot.pdendrogram <- function(data_pvalues,labels_pvalues=FALSE,treshold = 0.05,low=NULL,mid=NULL,high=NULL){
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


  p_floor <- 2.2e-16
  breaks <- c(p_floor,treshold, 1)
  labels <- c(expression("<10"^-16), as.character(treshold), "1")
  vals <- rescale(log10(c(p_floor,treshold, 1)))

  p <- ggplot() +
    geom_segment(
      data= data_pvalues$segments,
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
      #TODO: gerer pbs pvalues < treshold plotly
      data  = data_pvalues$merge_points %>% filter(label < treshold),
      aes(x = x_mid, y = y, label = label),
      vjust = 1.5,
      size  = 3,
    ) 
  }
  p
}

# rm(list=ls())
# library(fastcluster)
# library(ggplot2)
# data <- data.frame(
#   x = c(1,2,5,5.6,6),
#   y = c(2,2.5,4.7,5,4)
# )
#
# ggplot(data=data,aes(x=x,y=y))+geom_point()
# hcl <- hclust(dist(data,method="euclidean"),method="ward.D")
# plot(hcl)
#

get.individuals.merged.clusters <- function(hcl, k) {
  n <- length(hcl$order)
  if (k < 1 || k > n - 1) {
    stop("k must be an integer between 1 and n - 1.")
  }
  merge_step <- n - (k - 1)
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

