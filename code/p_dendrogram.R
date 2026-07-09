library(ggplot2)
library(scales)
library(patchwork)
library(MASS)
library(dplyr)
library(ggdendro)
library(fastcluster)
library(PCIdep)
library(ggnewscale)

plot.pdendrogram <- function(
  pvals,
  hcl,
  labels_pvalues = FALSE,
  threshold = 0.05,
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
  if (threshold < 0 || threshold > 1) {
    warning("threshold must be an integer between 0 and 1.")
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
  breaks <- c(p_floor, threshold, 1)
  # TODO: voir si on ne peut pas faire mieux que plain texte comme ça...
  labels <- c("< 2.2⁻¹⁶", as.character(threshold), "1")
  vals <- rescale(log10(c(p_floor, threshold, 1)))

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


plot.ground.truth <- function(
  pvals,
  hcl,
  threshold = 0.05,
  population,
  X,
  log_axis = FALSE,
  epsilon = 0
) {
  mean_hat <- get.mean_hat(X, population)
  kmax <- length(pvals)
  df <- data.frame(
    x = numeric(0),
    y = numeric(0),
    status = character(0)
  )
  for (k in 2:kmax) {
    individuals <- get.individuals.merged.clusters(hcl, k)
    indiv1 <- individuals[[1]]
    indiv2 <- individuals[[2]]

    mu_1 <- colMeans(mean_hat[indiv1, , drop = FALSE])
    mu_2 <- colMeans(mean_hat[indiv2, , drop = FALSE])
    nu_t_mu <- norm(mu_1 - mu_2)

    pval <- pvals[k - 1]

    status <- dplyr::case_when(
      pval < threshold & nu_t_mu > epsilon ~ "TP",
      pval >= threshold & nu_t_mu < epsilon ~ "TN",
      pval < threshold & nu_t_mu < epsilon ~ "FP",
      pval >= threshold & nu_t_mu > epsilon ~ "FN"
    )

    df <- rbind(
      df,
      data.frame(
        x = nu_t_mu,
        y = pval,
        status = status
      )
    )
  }

  p <- ggplot(df, aes(x = x, y = y, color = status)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      x = expression(nu^T * mu),
      y = "p-value",
      title = "Ground truth"
    ) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    annotate(
      "text",
      x = Inf,
      y = threshold,
      label = paste0("alpha = ", threshold),
      hjust = 1.1,
      vjust = -0.5,
      color = "red"
    )
  if (log_axis) {
    p <- p + scale_y_continuous(trans = 'log10')
  }
  p
}


plotly.pdendrogram <- function(p, slider = FALSE) {
  # TODO: test if sliders == boolean sinon error
  if (slider) {
    p$layers$geom_text$mapping <- modifyList(
      p$layers$geom_text$mapping,
      aes(frame = k)
    )
  }

  ggplotly(p, tooltip = "text") %>% style(textposition = "bottom")
}
