param_check <- function() {
  if (CONFIG$split_prop > 1 || CONFIG$split_prop < 0) {
    stop("split_prop should be in [0,1]")
  }
  if (!CONFIG$dist_method %in% c("euclidean")) {
    stop("dist_method should be one of: euclidean")
  }
  if (!CONFIG$hclust_method %in% c(
    "single", "average", "centroid",
    "ward.D", "median", "mcquitty", "complete"
  )) {
    stop('hclust_method must be one of: single, average, centroid, ward.D, median, mcquitty, complete')
  }
  if (CONFIG$nb_cluster < 2 || CONFIG$nb_cluster > 100) {
    stop("nb_cluster should be in [2,100]")
  }
  if (!CONFIG$lib %in% c("PCIdep", "clusterpval")) {
    stop("lib must be one of: PCIdep, clusterpval")
  }
  print("Params checked")
}

load_data <- function() {
  data      <- readRDS("data/1KGP_100PC.Rda")
  super_pop <- read.table(
    "data/20131219.populations.tsv",
    sep    = "\t",
    header = TRUE
  )
  list(data = data, super_pop = super_pop)
}

pair_to_idx <- function(i, j, nb_cluster) {
  (i - 1) * (nb_cluster - i / 2) + (j - i)
}

n_pairs <- function(k){
  return(k * (k - 1) / 2)
}

fit_hclust <- function(X) {
  hclust(
    dist(X, method = CONFIG$dist_method),
    method = CONFIG$hclust_method
  )
}

build_pop_list <- function(clusters, pop_col, super_pop) {
  nb_cluster <- max(clusters)

  majority_super <- vapply(seq_len(nb_cluster), function(k) {
    pops_in_k <- pop_col[clusters == k]
    super_pop %>%
      filter(Population.Code %in% pops_in_k) %>%
      count(Super.Population) %>%
      slice_max(n, n = 1, with_ties = FALSE) %>%
      pull(Super.Population)
  }, character(1))

  super_pop_levels <- unique(super_pop$Super.Population)
  setNames(
    lapply(super_pop_levels, function(sp) which(majority_super == sp)),
    super_pop_levels
  )
}

compute_pvals <- function(X, nb_cluster, pc, lib) {

  set.seed(CONFIG$seed)

  split <- initial_split(X, prop = CONFIG$split_prop)
  Y <- training(split)[, seq_len(pc)]
  data_clust <- testing(split)

  hcl <- fit_hclust(data_clust)

  pairs <- which(upper.tri(matrix(NA, nb_cluster, nb_cluster)), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1], pairs[, 2]), ]

  pvec <- numeric(n_pairs(nb_cluster))

  for (r in seq_len(nrow(pairs))) {

    i <- pairs[r, 1]
    j <- pairs[r, 2]
    idx <- pair_to_idx(i, j, nb_cluster)

    p <- switch(
      lib,

      "PCIdep" = {
        PCIdep::test.clusters.hc(
          hcl = hcl,
          X = data_clust[, 1:pc],
          cluster = c(i, j),
          NC = nb_cluster,
          Y = Y
        )$pvalue
      },

      "clusterpval" = {
        clusterpval::test_hier_clusters_exact(
          as.matrix(data_clust[, 1:pc]),
          link = "ward.D",
          K = nb_cluster,
          k1 = i,
          k2 = j,
          hcl = hcl
        )$pval
      },

      stop("lib must be one of: PCIdep, clusterpval")
    )

    pvec[idx] <- p

    print(paste("Pair (", i, ",", j, ") - p =", p))
  }

  pvec
}


build_absolute <- function(pop_list, nb_cluster) {
  vapply(
    seq_len(n_pairs(nb_cluster)),
    function(idx) {
      for (i in seq_len(nb_cluster - 1)) {
        for (j in (i + 1):nb_cluster) {
          if (pair_to_idx(i, j, nb_cluster) == idx) {
            same <- any(vapply(pop_list, function(pop) i %in% pop && j %in% pop, logical(1)))
            return(as.integer(same))
          }
        }
      }
      NA_integer_
    },
    integer(1)
  )
}

pval_l1_distance <- function(pvec, ground_truth){ 
  return(sum(abs(pvec - ground_truth)))}



optim_pcs <- function(X, nb_cluster, ground_truth, pc_list = CONFIG$pc_search_seq) {
  results <- lapply(pc_list, function(pc) {
    message(sprintf("Testing pc = %d ...", pc))
    pvec <- compute_pvals(X, nb_cluster, pc)
    data.frame(pc = pc, distance = pval_l1_distance(pvec, ground_truth))
  })
  do.call(rbind, results)
}

build_cluster_mapping <- function(pop_list,nb_cluster) {
  data.frame(
    position  = 1:nb_cluster,
    cluster   = unname(unlist(lapply(pop_list, sort))),
    super_pop = rep(names(pop_list), lengths(pop_list))
  )
}

plot_pval_heatmap <- function(pvals_flat, nb_cluster, cluster_mapping) {
  pmat <- matrix(NA, nrow = nb_cluster, ncol = nb_cluster)
  for (i in seq_len(nb_cluster - 1)) {
    for (j in (i + 1):nb_cluster) {
      idx        <- pair_to_idx(i, j, nb_cluster)
      pmat[i, j] <- pvals_flat[idx]
    }
  }

  ordered_clusters <- cluster_mapping$cluster
  id_to_pos        <- setNames(cluster_mapping$position, cluster_mapping$cluster)

  df <- melt(pmat, varnames = c("k1", "k2"), value.name = "pval") %>%
    filter(!is.na(pval)) %>%
    mutate(
      k1 = id_to_pos[as.character(k1)],
      k2 = id_to_pos[as.character(k2)]
    )

  df_sym <- bind_rows(df, rename(df, k1 = k2, k2 = k1)) %>%
    distinct(k1, k2, .keep_all = TRUE) %>%
    filter(k1 >= k2)

  pop_bounds <- cluster_mapping %>%
    group_by(super_pop) %>%
    summarise(
      start = min(position) - 0.5,
      end   = max(position) + 0.5,
      mid   = mean(position),
      .groups = "drop"
    )

  p <- ggplot(df_sym, aes(x = k1, y = k2, fill = pval)) +
    geom_tile(color = "black", linewidth = 0.3) +
    scale_x_continuous(breaks = seq_len(nb_cluster), labels = ordered_clusters) +
    scale_y_continuous(breaks = seq_len(nb_cluster), labels = ordered_clusters) +
    scale_fill_gradient(
      low      = "#ffffff",
      high     = "#1b263b",
      name     = "p-value",
      na.value = "white"
    ) +
    coord_equal(clip = "off") +
    theme_minimal() +
    theme(
      panel.grid  = element_blank(),
      axis.title  = element_blank(),
      axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
      plot.margin = margin(t = 10, r = 10, b = 50, l = 50)
    )

  for (i in seq_len(nrow(pop_bounds))) {
    p <- p +
      annotate("rect",
        xmin = pop_bounds$start[i], xmax = pop_bounds$end[i],
        ymin = -2,                  ymax = -0.7,
        fill = "#fff", color = "black"
      ) +
      annotate("text",
        x     = pop_bounds$mid[i], y = -1.35,
        label = pop_bounds$super_pop[i], size = 3
      ) +
      annotate("rect",
        ymin = pop_bounds$start[i], ymax = pop_bounds$end[i],
        xmin = -2,                  xmax = -0.7,
        fill = "#fff", color = "black"
      ) +
      annotate("text",
        y     = pop_bounds$mid[i], x = -1.35,
        label = pop_bounds$super_pop[i], size = 3, angle = 90
      )
  }

  p
}

plot_pc_optim <- function(optim_pcs_df) {
  ggplot(optim_pcs_df, aes(x = pc, y = distance)) +
    geom_point() +
    scale_x_continuous(breaks = seq(min(optim_pcs_df$pc), max(optim_pcs_df$pc), by = 5)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank())
}

save_plot <- function(plot, filename, width = 8, height = 8, dpi = 300) {
  dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(CONFIG$output_dir, filename)
  ggsave(path, plot = plot, dpi = dpi, width = width, height = height)
  message("Saved: ", path)
  invisible(path)
}



