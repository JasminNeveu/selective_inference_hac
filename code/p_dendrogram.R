library(ggplot2)
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
        ,
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

selective.cutree <- function(hcl, k_list) {
  n <- length(hcl$order)
  labels <- numeric(n)
  for (i in seq_len(length(k_list))) {
    k <- k_list[i]
    indiv <- get.individuals(hcl, k)
    labels[indiv] <- i
  }
  labels
}

plot.heatmap <- function(
  clusters,
  population,
  super_population = NULL,
  cluster_name = "Cluster",
  population_name = "Population"
) {
  stopifnot(length(clusters) == length(population))

  tab <- table(population, clusters)
  tab_prop <- prop.table(tab, margin = 1)

  best_cluster <- apply(tab_prop, 1, which.max)
  pop_order <- names(sort(best_cluster))

  df <- as.data.frame(tab_prop)
  colnames(df) <- c("population", "cluster", "prop")

  df$population <- factor(df$population, levels = pop_order)
  df$cluster <- factor(df$cluster, levels = sort(unique(clusters)))

  p <- ggplot(df, aes(population, cluster, fill = prop)) +
    geom_tile(color = "black", linewidth = 0.3) +
    coord_fixed() +
    labs(
      x = population_name,
      y = cluster_name,
      fill = "Proportion"
    ) +
    scale_fill_gradient(
      low = "white",
      high = "#1b263b"
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      )
    )

  if (!is.null(super_population)) {
    stopifnot(length(super_population) == length(population))

    label_df <- data.frame(
      population = population,
      super_population = super_population
    ) |>
      dplyr::distinct()

    label_df <- label_df[match(pop_order, label_df$population), ]

    superpop_colors <- c(
      AFR = "#ff9f1c",
      AMR = "#e71d36",
      EAS = "#2ec4b6",
      EUR = "#011627",
      SAS = "#984EA3"
    )

    label_colors <- superpop_colors[label_df$super_population]

    p <- p +
      geom_point(
        data = transform(
          df,
          super_population = label_df$super_population[
            match(population, label_df$population)
          ]
        ),
        aes(color = super_population),
        alpha = 0,
        inherit.aes = TRUE
      ) +
      scale_color_manual(
        values = superpop_colors,
        name = "Super-population",
        guide = guide_legend(
          override.aes = list(alpha = 1, size = 3)
        )
      ) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5,
          colour = label_colors
        )
      )
  }

  p
}
