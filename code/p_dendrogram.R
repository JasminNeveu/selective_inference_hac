library(ggplot2)
library(MASS)
library(dplyr)
library(ggdendro)
library(fastcluster)
library(PCIdep)
library(scales)

build.data <- function(X,kmax = nrow(X)-5, U = NULL, Sigma = NULL, Y = NULL, UY = NULL, precUY = NULL, linkage = 'ward.D', hcl = NULL, ndraws = 2000, sample_split = FALSE, nY = NULL, return_Sigma = FALSE, return_X_clus = FALSE, dismat = NULL,parallel_config=NULL){
  n <- nrow(X)
  if(kmax < 1 || kmax > n - 1){
    stop("kmax must be an integer between 1 and n - 1.")
  }
  data_pvalues <- get.data(hcl,X, U,Sigma , Y , UY , precUY  , linkage,ndraws, sample_split, nY, return_Sigma, return_X_clus, dismat,kmax,parallel_config)
  data_pvalues
}


plot.pdendrogram <- function(hcl,data_pvalues,labels_pvalues=FALSE,treshold = 0.05,low=NULL,mid=NULL,high=NULL,groups=NULL,groups_labels=NULL){
  # check paremeters
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

  # plot branches
  p <- ggplot() +
    geom_segment(
      data= data_pvalues$segments,
      aes(x = x, y = y, xend = xend, yend = yend, colour = pval,text="") 
    ) + 
    labs(title="p-dendrogram",x = "Observations",y="Merging distances") + 
    scale_color_gradientn(
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

  # plot labels



  if(labels_pvalues){
    # add groups labels
    if(!is.null(groups)){
      if (is.null(groups_labels)) {
        groups_labels <- paste0("group", seq_along(groups))
      }
      if (length(groups_labels) != length(groups)) {
        stop("groups_labels must have the same length as groups.")
      }
    }

    data_pvalues <- labels.groups(data_pvalues,hcl,groups,groups_labels)
    df <- data_pvalues$merge_points 
    data_pvalues$merge_points$hover_text <- vapply(
      seq_len(nrow(df
      )),
      function(i) {
        row <- df[i, , drop = TRUE]
        txt <- paste0("\n size = ",row[["cluster_size"]])

        for(label in groups_labels) {
          txt <- paste0(
            txt,
            "\n",
            label, ": ",
            row[[paste0(label, "_class")]],
            " (",
            round(100 * as.numeric(row[[paste0(label, "_prop")]]), 1),
            "%)"
          )
        }
        txt
      },
      character(1)
    )

    p <- p +
      geom_text(
        data = data_pvalues$merge_points %>% filter(!is.null(label)),
        aes(
          x = x_mid,
          y = y,
          label = label,
          text=  paste0("p-value = ",label, "\n k = ", k,hover_text)
        ),
        vjust = 1.5,
        size  = 3
      )
  }
  p
}
