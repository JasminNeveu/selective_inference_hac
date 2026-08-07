library(ggplot2)
library(clue)
library(progressr)
library(scales)
library(patchwork)
library(MASS)
library(dplyr)
library(ggdendro)
library(fastcluster)
library(PCIdep)
library(ggnewscale)

plot.pdendrogram <- function(
  #TODO: checker parametres
  pvals,
  hcl,
  labels_pvalues = FALSE,
  alpha = 0.05,
  strip_height_factor = 1,
  low = NULL,
  mid = NULL,
  high = NULL,
  groups = NULL,
  groups_labels = NULL,
  groups_subclasses = NULL,
  log_height = FALSE
) {
  # check paremeters
  if (alpha < 0 || alpha > 1) {
    warning("alpha must be an integer between 0 and 1.")
  }
  filled <- c(!is.null(low), !is.null(mid), !is.null(high))
  partial_defined <- sum(filled) > 0 && sum(filled) < 3
  if (partial_defined) {
    # si low ou mid ou high est renseigné, alors tous doivent l'être sinon c'est mal défini
    stop(
      "Arguments 'low', 'mid' and 'high' must be provided together (either all three or none)."
    )
  }
  if (!is.null(low) && !is.null(mid) && !is.null(high)) {
    colours <- c(low, mid, high)
  } else {
    colours <- c("#B2182B", "#D9C27A", "#2166AC") #option 3
  }

  if (log_height) {
    hcl$height <- log10(hcl$height + 1)
  }

  p_floor <- 2.2e-16
  breaks <- c(p_floor, alpha, 1)
  # TODO: voir si on ne peut pas faire mieux que plain texte comme ça...
  labels <- c("< 2.2⁻¹⁶", as.character(alpha), "1")
  vals <- rescale(log10(c(p_floor, alpha, 1)))

  data_pvalues <- data.pdendrogram(pvals, hcl)

  # plot branches
  p <- ggplot() +
    geom_segment(
      data = data_pvalues$segments,
      aes(x = x, y = y, xend = xend, yend = yend, colour = pval)
    ) +
    labs(title = "p-dendrogram", x = "Observations", y = "Merging distances") +
    scale_color_gradientn(
      colours = colours,
      values = vals,
      trans = "log10",
      limits = c(p_floor, 1),
      breaks = breaks,
      labels = labels,
      na.value = "grey",
      name = "p-value"
    ) +
    theme_minimal()

  # plot groups individuals
  if (!is.null(groups)) {
    validate_groups_labels(
      groups,
      groups_labels,
      groups_subclasses
    )

    leaves_long <- build_leaves(
      groups = groups,
      ord = hcl$order,
      labels = groups_labels,
      subclasses = groups_subclasses
    )
    tmp <- assign_strip_positions(leaves_long, strip_height_factor)

    leaves_long <- tmp$data
    strip_height <- tmp$strip_height
    p <- add_cluster_strips(p, leaves_long, strip_height)
  }

  # add pvalues labels (and plotly hover text)
  if (labels_pvalues) {
    data_pvalues <- labels.groups(
      data_pvalues,
      hcl,
      groups,
      groups_labels,
      groups_subclasses
    )
    df <- data_pvalues$merge_points
    data_pvalues$merge_points$hover_text <- vapply(
      seq_len(nrow(df)),
      function(i) {
        row <- df[i, , drop = TRUE]

        txt <- paste0("\nsize = ", row[["cluster_size"]])
        for (group_name in groups_labels) {
          txt <- paste0(
            txt,
            "\n",
            group_name,
            ": ",
            row[[paste0(group_name, "_class")]],
            " (",
            round(
              100 * as.numeric(row[[paste0(group_name, "_prop")]]),
              1
            ),
            "%)"
          )
        }

        txt
      },
      character(1)
    )
    ymax <- max(data_pvalues$merge_points$y)
    p <- p +
      geom_text(
        data = data_pvalues$merge_points %>%
          filter(!is.na(label)),
        aes(
          x = x_mid,
          y = y,
          label = label,
          text = paste0("p-value = ", label, "\nk = ", k, hover_text)
        ),

        vjust = "top",
        size = 3
      )
  }
  p
}


plotly.pdendrogram <- function(p, slider = FALSE) {
  if (!is_ggplot(p)) {
    stop("p should be an ggplot object.")
  }
  if (!is.logical(slider)) {
    stop("slider should be a boolean.")
  }

  if (slider) {
    p$layers$geom_text$mapping <- modifyList(
      p$layers$geom_text$mapping,
      aes(frame = k)
    )
  }

  ggplotly(p, tooltip = "text") %>% style(textposition = "bottom")
}


selective.cutree <- function(
  hcl,
  pvals,
  min_pts = 10,
  alpha = 0.05
) {
  node_sizes <- get.node.sizes(hcl)
  n_obs <- length(hcl$order)
  labels <- rep(0L, n_obs)
  cluster_id <- 1L

  # Fonction récursive de visite de l'arbre
  # On garde 'peeling_root' en mémoire pour annuler le bruit si cul-de-sac
  visit <- function(row, peeling_root) {
    if (row < 0) {
      return()
    }

    c1 <- hcl$merge[row, 1]
    c2 <- hcl$merge[row, 2]

    s1 <- if (c1 < 0) 1 else node_sizes[c1]
    s2 <- if (c2 < 0) 1 else node_sizes[c2]

    # Cas 1 : Cul-de-sac (les deux enfants sont trop petits)
    if (s1 < min_pts && s2 < min_pts) {
      leafs <- get.node.leafs(hcl, peeling_root, node_sizes)
      labels[leafs] <<- cluster_id
      cluster_id <<- cluster_id + 1L
      return()
    }

    # Cas 2 : Grignotage asymétrique (gauche trop petit, on continue à droite)
    if (s1 < min_pts && s2 >= min_pts) {
      if (c2 > 0) {
        visit(c2, peeling_root)
      }
      return()
    }

    # Cas 3 : Grignotage asymétrique (droite trop petit, on continue à gauche)
    if (s1 >= min_pts && s2 < min_pts) {
      if (c1 > 0) {
        visit(c1, peeling_root)
      }
      return()
    }

    # Cas 4 : Les deux enfants sont assez grands (test de validité)
    if (s1 >= min_pts && s2 >= min_pts) {
      # FIX 2 : Traduction de l'indice de ligne (row) vers le nombre de clusters (k)
      k <- n_obs - row + 1
      p_val <- pvals[as.character(k)]

      if (!is.na(p_val) && p_val < alpha) {
        # Split significatif : les enfants deviennent de nouvelles racines de pelage
        if (c1 > 0) {
          visit(c1, peeling_root = c1)
        }
        if (c2 > 0) visit(c2, peeling_root = c2)
      } else {
        # Noyau dur trouvé : on valide le noyau courant
        leafs <- get.node.leafs(hcl, row, node_sizes)
        labels[leafs] <<- cluster_id
        cluster_id <<- cluster_id + 1L
      }
      return()
    }
  }

  root_node <- nrow(hcl$merge)
  visit(root_node, peeling_root = root_node)

  return(labels)
}

plot.heatmap <- function(
  true_labels,
  predicted_clusters,
  super_population = NULL,
  prop = FALSE,
  predicted_clusters_name = "Clusters"
) {
  M <- build.heatmap.matrix(true_labels, predicted_clusters, prop)
  df <- as.data.frame(as.table(M))
  names(df) <- c("Cluster", "Label", "Value")

  legend_text <- if (prop) "Proportion" else "Count"

  if (!is.null(super_population)) {
    col_labels <- colnames(M)

    # Matching sécurisé s'adaptant au type (entier, caractère, facteur)
    matched_superpop <- super_population[as.character(col_labels)]

    if (any(is.na(matched_superpop))) {
      warning("Certains labels n'ont pas de super_population associée.")
    }

    unique_superpops <- unique(na.omit(matched_superpop))
    n_superpops <- length(unique_superpops)

    if (n_superpops <= 8) {
      palette_colors <- RColorBrewer::brewer.pal(max(3, n_superpops), "Set1")[
        1:n_superpops
      ]
    } else {
      palette_colors <- rainbow(n_superpops)
    }

    superpop_colors <- setNames(palette_colors, unique_superpops)
    label_colors <- unname(superpop_colors[matched_superpop])
  }

  p <- ggplot(
    df,
    aes(x = Label, y = factor(Cluster, levels = rownames(M)), fill = Value)
  ) +
    geom_tile(color = "black", linewidth = 0.3) +
    coord_fixed() +
    scale_fill_gradient(low = "white", high = "#1b263b") +
    labs(
      x = "True population",
      y = predicted_clusters_name,
      fill = legend_text
    ) +
    theme_minimal() +
    theme(panel.grid = element_blank())

  if (!is.null(super_population)) {
    col_labels_orig <- type.convert(colnames(M), as.is = TRUE)

    legend_df <- data.frame(
      Label = col_labels_orig,
      Cluster = rownames(M)[1],
      super_population = matched_superpop
    )

    angle_text <- if (is.character(true_labels)) 90 else 0
    p <- p +
      geom_point(
        data = legend_df,
        aes(x = Label, y = Cluster, color = super_population),
        alpha = 0,
        inherit.aes = FALSE
      ) +
      scale_color_manual(
        values = superpop_colors,
        name = "Super-population",
        guide = guide_legend(override.aes = list(alpha = 1, size = 3))
      ) +
      theme(
        axis.text.x = element_text(
          angle = angle_text,
          hjust = 1,
          vjust = 0.5,
          colour = label_colors
        )
      )
  } else {
    p <- p +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        )
      )
  }

  p
}
