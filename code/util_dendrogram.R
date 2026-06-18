get.merged.clusters <- function(hcl, k) {
  n <- length(hcl$order)
  merge_step <- n - (k - 1)

  c1_internal <- hcl$merge[merge_step, 1]
  c2_internal <- hcl$merge[merge_step, 2]

  clusters_k <- cutree(hcl, k = k)

  get_leaf_index <- function(id) {
    if (id < 0) {
      return(-id)
    }
    get_leaf_index(hcl$merge[id, 1])
  }

  leaf1 <- get_leaf_index(c1_internal)
  leaf2 <- get_leaf_index(c2_internal)

  c1 <- clusters_k[leaf1]
  c2 <- clusters_k[leaf2]
  c(unname(c1), unname(c2))
}

compute.pvals <- function(
  hcl,
  X,
  U,
  Sigma,
  Y,
  UY,
  precUY,
  linkage,
  ndraws,
  sample_split,
  nY,
  return_Sigma,
  return_X_clus,
  dismat,
  kmax,
  parallel_config
) {
  k_seq <- 2:kmax

  run_one_k <- function(k) {
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
      clusters = clusters,
      ndraws = ndraws,
      sample_split = sample_split,
      nY = nY,
      return_Sigma = return_Sigma,
      return_X_clus = return_X_clus,
      hcl = hcl,
      dismat = dismat
    )

    hc_test$pval <- ifelse(hc_test$pval < 2.2e-16, 2.2e-16, hc_test$pval)
    print(paste(
      "k =",
      k,
      "clusters:",
      clusters[1],
      "-",
      clusters[2],
      "---",
      hc_test$pval
    ))
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

get.data <- function(
  hcl,
  X,
  U,
  Sigma,
  Y,
  UY,
  precUY,
  linkage,
  ndraws,
  sample_split,
  nY,
  return_Sigma,
  return_X_clus,
  dismat,
  kmax,
  parallel_config
) {
  hcldata <- ggdendro::dendro_data(hcl, type = "rectangle")
  hcldata$segments <- hcldata$segments[order(-hcldata$segments$y), ] # ordonne les segments pour qu'ils soient lister de haut en bas de l'arbre
  n <- nrow(hcldata$segments)
  pvec <- compute.pvals(
    hcl,
    X,
    U,
    Sigma,
    Y,
    UY,
    precUY,
    linkage,
    ndraws,
    sample_split,
    nY,
    return_Sigma,
    return_X_clus,
    dismat,
    kmax,
    parallel_config
  )
  pvec_rep <- rep(pvec, each = 4)
  hcldata$segments$pval <- NA
  hcldata$segments$pval[1:length(pvec_rep)] <- pvec_rep

  # ------Add labels pvalues ------
  merge_points <- hcldata$segments %>%
    filter(y == yend) %>%
    group_by(y) %>%
    summarise(x_min = min(x), x_max = max(x)) %>%
    arrange(desc(y)) %>%
    mutate(
      x_mid = (x_min + x_max) / 2,
      k = row_number() + 1,
      label = c(signif(pvec, 3), rep(NA, n() - length(pvec)))
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
  merge_step <- n - (k - 1)
  c1_internal <- hcl$merge[merge_step, 1]
  c2_internal <- hcl$merge[merge_step, 2]
  get_individuals <- function(id) {
    if (id < 0) {
      # cas d'un indiv
      return(-id)
    } else {
      # cas d'un cluster
      left <- get_individuals(hcl$merge[id, 1])
      right <- get_individuals(hcl$merge[id, 2])
      return(c(left, right))
    }
  }

  individuals_1 <- get_individuals(c1_internal)
  individuals_2 <- get_individuals(c2_internal)

  return(list(unname(individuals_1), unname(individuals_2)))
}

get.individuals <- function(hcl, k) {
  indiv <- get.individuals.merged.clusters(hcl, k + 1)
  c(indiv[[1]], indiv[[2]])
}

clustering.type <- function(x, tol = 2.1e-16) {
  if (is.factor(x) || is.character(x) || is.integer(x)) {
    return("hard")
  }

  if (is.list(x)) {
    same_length <- length(unique(lengths(x))) == 1
    valid_memberships <- all(vapply(
      x,
      function(v) {
        is.numeric(v) &&
          all(v >= -tol) &&
          all(v <= 1 + tol) &&
          abs(sum(v) - 1) < tol
      },
      logical(1)
    ))

    if (same_length && valid_memberships) {
      return("soft")
    }
  }

  if (is.matrix(x) || is.data.frame(x)) {
    x <- as.matrix(x)
    if (
      is.numeric(x) &&
        all(x >= -tol) &&
        all(x <= 1 + tol) &&
        all(abs(rowSums(x) - 1) < tol)
    ) {
      return("soft")
    }
  }

  return("unknown")
}

validate_groups_labels <- function(
  groups,
  groups_labels,
  groups_subclasses = NULL
) {
  if (is.null(groups_labels)) {
    return(invisible(NULL))
  }

  if (length(groups_labels) != length(groups)) {
    stop("'groups_labels' must have the same length as 'groups'.")
  }

  if (is.null(groups_subclasses)) {
    groups_subclasses <- vector("list", length(groups))
  }

  if (length(groups_subclasses) != length(groups)) {
    stop("'groups_subclasses' must have the same length as 'groups'.")
  }

  for (i in seq_along(groups)) {
    if (
      !is.character(groups_labels[[i]]) ||
        length(groups_labels[[i]]) != 1
    ) {
      stop(
        paste0(
          "groups_labels[[",
          i,
          "]] must be a character(1)."
        )
      )
    }

    type <- clustering.type(groups[[i]])

    if (type == "soft") {
      n_sub <- if (is.list(groups[[i]])) {
        length(groups[[i]])
      } else {
        ncol(groups[[i]])
      }

      if (is.null(groups_subclasses[[i]])) {
        stop(
          paste0(
            "groups_subclasses[[",
            i,
            "]] must be provided for soft clusters."
          )
        )
      }

      if (length(groups_subclasses[[i]]) != n_sub) {
        stop(
          paste0(
            "groups_subclasses[[",
            i,
            "]] has ",
            length(groups_subclasses[[i]]),
            " elements but group has ",
            n_sub,
            " subclasses."
          )
        )
      }
    }
  }

  invisible(NULL)
}

labels.groups <- function(
  data_pvalues,
  hcl,
  groups,
  groups_labels,
  groups_subclasses = NULL
) {
  if (is.null(groups_subclasses)) {
    groups_subclasses <- vector("list", length(groups))
  }

  clusters_list <- lapply(data_pvalues$merge_points$k, function(k) {
    cl <- get.individuals.merged.clusters(hcl, k)
    unique(c(cl[[1]], cl[[2]]))
  })

  cluster_sizes <- vapply(clusters_list, length, numeric(1))

  for (i in seq_along(groups)) {
    group <- groups[[i]]
    group_name <- groups_labels[[i]]
    subclasses <- groups_subclasses[[i]]

    type_group <- clustering.type(group)

    tmp <- lapply(clusters_list, function(indiv) {
      if (type_group == "hard") {
        tab <- table(group[indiv])

        idx <- which.max(tab)

        list(
          class = names(tab)[idx],
          prop = as.numeric(tab[idx] / sum(tab))
        )
      } else {
        mean_vals <- colMeans(group[indiv, , drop = FALSE])

        idx <- which.max(mean_vals)

        list(
          class = subclasses[idx],
          prop = as.numeric(mean_vals[idx])
        )
      }
    })

    data_pvalues$merge_points[[paste0(group_name, "_class")]] <-
      vapply(tmp, `[[`, character(1), "class")

    data_pvalues$merge_points[[paste0(group_name, "_prop")]] <-
      vapply(tmp, `[[`, numeric(1), "prop")
  }

  data_pvalues$merge_points$cluster_size <- cluster_sizes

  data_pvalues
}


normalize_hard <- function(group, ord) {
  data.frame(
    x = seq_along(ord),
    group_value = as.character(group)[ord],
    subclass_label = NA_character_,
    group_type = "hard",
    stringsAsFactors = FALSE
  )
}


normalize_soft <- function(group, subclasses, ord) {
  if (is.list(group)) {
    mat <- do.call(cbind, group)
  } else {
    mat <- as.matrix(group)
  }
  mat <- mat[ord, , drop = FALSE]

  n_indiv <- nrow(mat)
  n_subclass <- ncol(mat)

  data.frame(
    x = rep(seq_len(n_indiv), times = n_subclass),
    group_value = as.numeric(mat),
    subclass_label = rep(subclasses, each = n_indiv),
    group_type = "soft",
    stringsAsFactors = FALSE
  )
}

build_leaves <- function(
  groups,
  ord,
  labels = NULL,
  subclasses = NULL,
  tol = 2.2e-15
) {
  n <- length(groups)

  if (is.null(names(groups))) {
    names(groups) <- paste0("group_", seq_len(n))
  }
  if (is.null(labels)) {
    labels <- as.list(names(groups))
  }
  if (is.null(subclasses)) {
    subclasses <- vector("list", n)
  }

  rows <- lapply(seq_len(n), function(i) {
    group <- groups[[i]]
    group_name <- labels[[i]]
    type <- clustering.type(group, tol)

    if (type == "unknown") {
      stop(sprintf("Groupe '%s' : type inconnu (ni hard ni soft).", group_name))
    }

    if (type == "hard") {
      df <- normalize_hard(group, ord)
    } else {
      sub <- subclasses[[i]]
      if (is.null(sub) || (length(sub) == 1 && is.na(sub))) {
        k <- if (is.matrix(group) || is.data.frame(group)) {
          ncol(group)
        } else {
          length(group[[1]])
        }
        sub <- paste0("sub_", seq_len(k))
      }
      df <- normalize_soft(group, sub, ord)
    }

    df$group_id <- i
    df$group_name <- group_name
    df
  })

  result <- do.call(rbind, rows)
  result[, c(
    "x",
    "group_id",
    "group_name",
    "subclass_label",
    "group_value",
    "group_type"
  )]
}

assign_strip_positions <- function(
  leaves_long,
  strip_height_factor,
  base_y = 0,
  gap = 0.1
) {
  n_indiv <- max(leaves_long$x, na.rm = TRUE)
  strip_height <- strip_height_factor * n_indiv
  gap_size <- gap * n_indiv

  leaves_long$strip_y <- -((leaves_long$group_id - 0.5) *
    strip_height +
    gap_size)

  list(data = leaves_long, strip_height = strip_height)
}

prepare_soft_stack <- function(soft_df) {
  soft_df <- soft_df[order(soft_df$x, soft_df$subclass_label), ]
  soft_df$x_min <- soft_df$x - 0.5
  soft_df$x_max <- soft_df$x + 0.5

  chunks <- split(soft_df, list(soft_df$group_id, soft_df$x), drop = TRUE)
  chunks <- lapply(chunks, function(g) {
    g$y_min <- c(0, head(cumsum(g$group_value), -1))
    g$y_max <- cumsum(g$group_value)
    g
  })
  do.call(rbind, chunks)
}


add_cluster_strips <- function(p, leaves_long, strip_height) {
  soft_ids <- sort(unique(leaves_long$group_id[
    leaves_long$group_type == "soft"
  ]))
  hard_ids <- sort(unique(leaves_long$group_id[
    leaves_long$group_type == "hard"
  ]))

  if (length(soft_ids) > 0) {
    soft_all <- leaves_long[leaves_long$group_type == "soft", ]
    all_subclasses <- unique(soft_all$subclass_label)
    colors <- setNames(
      scales::hue_pal()(length(all_subclasses)),
      all_subclasses
    )
    stacked_all <- prepare_soft_stack(soft_all)
    legend_title <- paste(unique(soft_all$group_name), collapse = ", ")

    p <- p +
      ggplot2::geom_rect(
        data = stacked_all,
        ggplot2::aes(
          xmin = x_min,
          xmax = x_max,
          ymin = strip_y - (0.5 - y_min) * strip_height,
          ymax = strip_y - (0.5 - y_max) * strip_height,
          fill = subclass_label,
          text = subclass_label
        )
      ) +
      ggplot2::scale_fill_manual(values = colors, name = legend_title)

    if (length(hard_ids) > 0) p <- p + ggnewscale::new_scale_fill()
  }

  for (i in seq_along(hard_ids)) {
    gid <- hard_ids[[i]]
    hard_df <- leaves_long[leaves_long$group_id == gid, ]
    name <- unique(hard_df$group_name)

    p <- p +
      ggplot2::geom_tile(
        data = hard_df,
        ggplot2::aes(x = x, y = strip_y, fill = group_value),
        width = 1,
        height = strip_height
      ) +
      ggplot2::scale_fill_discrete(name = name)

    if (i < length(hard_ids)) p <- p + ggnewscale::new_scale_fill()
  }

  p
}
