library(PCIdep)
library(ggplot2)
library(reshape2)
library(dplyr)
library(rsample)
library(fastcluster)
pak::pak("lucylgao/clusterpval")
source("utils.R")

CONFIG <- list(
  seed = 42,
  split_prop = 0.2,
  n_pcs = 11,
  pc_search_seq = seq(1, 30, by = 2),
  dist_method = "euclidean",
  hclust_method = "ward.D",
  output_dir = "output",
  nb_cluster = 26,
  lib = "PCIdep"

)

main <- function() {
  raw <- load_data()
  data <- raw$data
  super_pop <- raw$super_pop
  X <- data %>% select(-ID, -pop)
  hcl_full <- fit_hclust(X)

  clusters <- cutree(hcl_full, k = CONFIG$nb_cluster)
  pop_list <- build_pop_list(clusters, data$pop, super_pop)
  cluster_mapping <- build_cluster_mapping(pop_list,CONFIG$nb_cluster)

  # ground_truth <- build_absolute(pop_list, CONFIG$nb_cluster)

  # optim_pcs_df <- optim_pcs(X, CONFIG$nb_cluster, ground_truth)
  # save_plot(plot_pc_optim(optim_pcs_df), "distances_pcs.png")
  print(paste("Computing p-values with",CONFIG$n_pcs,"PCs", "for", CONFIG$nb_cluster,"clusters", "with",CONFIG$lib, "package"))
  pvec <- compute_pvals(X, CONFIG$nb_cluster, CONFIG$n_pcs,CONFIG$lib)
 
  save_plot(
    plot_pval_heatmap(pvec, CONFIG$nb_cluster, cluster_mapping),
    sprintf("heatmap_%s_pval_pc%02d.png", CONFIG$lib,CONFIG$n_pcs)
  )
  }



main()



