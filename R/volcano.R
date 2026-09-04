.pf_volcano_measure_labels <- c(
  mean_difference = "Mean difference",
  odds_ratio = "Log odds ratio",
  incidence_rate_ratio = "Log incidence-rate ratio",
  common_odds_ratio = "Log common odds ratio"
)

.pf_volcano_threshold <- function(data, significance, alpha) {
  if (identical(significance, "none")) {
    return(data.frame(
      method = character(), y = numeric(), threshold_label = character(),
      stringsAsFactors = FALSE
    ))
  }

  if (identical(significance, "bonferroni")) {
    family_size <- .pf_plot_family_size(data)
    if (family_size < 1L) {
      stop("No successful tests are available for Bonferroni correction.",
           call. = FALSE)
    }
    y <- -log10(alpha / family_size)
    return(data.frame(
      method = "Bonferroni",
      y = y,
      threshold_label = paste0(
        "Bonferroni: ", format(signif(y, 3L), trim = TRUE)
      ),
      stringsAsFactors = FALSE
    ))
  }

  passing <- data$status == "ok" & is.finite(data$neg_log10_p) &
    is.finite(data$q_value) & data$q_value <= alpha
  if (!any(passing)) {
    return(data.frame(
      method = character(), y = numeric(), threshold_label = character(),
      stringsAsFactors = FALSE
    ))
  }
  y <- min(data$neg_log10_p[passing])
  data.frame(
    method = "BH FDR",
    y = y,
    threshold_label = paste0(
      "BH FDR: ", format(signif(y, 3L), trim = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

.pf_volcano_select_labels <- function(data, label, top_label_n, significance,
                                      threshold, alpha, missing) {
  if (!is.null(label)) {
    return(.pf_select_ids(data, label, "label", missing))
  }
  if (!top_label_n) return(character())

  significant <- switch(
    significance,
    bonferroni = if (nrow(threshold)) data$y_plot >= threshold$y[[1L]] else
      rep(FALSE, nrow(data)),
    bh = is.finite(data$q_value) & data$q_value <= alpha,
    none = rep(FALSE, nrow(data))
  )
  ordering <- order(
    !significant, -data$y_plot, -abs(data$x_plot), data$phenotype,
    method = "radix"
  )
  data$phenotype[utils::head(ordering, top_label_n)]
}

#' Plot a PheWAS volcano result
#'
#' Plots fitted coefficients against the authoritative stored `-log10(p)`
#' values using the same group legend, colors, effect-direction symbols, theme,
#' and annotation style as [plot_phewas_manhattan()]. Coefficients remain on
#' their fitted model or link scale. Results containing different effect
#' measures are separated into free-scale facets, each centered at zero,
#' because their magnitudes are not directly comparable. Mean differences also
#' require common outcome units and transformations before their magnitudes can
#' be compared within a panel.
#'
#' @details The default Bonferroni threshold uses the recorded original
#'   successful testing family, including successful rows omitted from the plot
#'   because their estimate is not finite. Threshold text is placed inside the
#'   plotting panel rather than added as an axis tick.
#'
#' @param x A `phewas_result`.
#' @param significance Testing threshold to draw: `"bonferroni"`, `"bh"`, or
#'   `"none"`.
#' @param label,highlight Explicit phenotype values to label or ring.
#' @param category_order Optional character vector giving every plotted group
#'   in display order.
#' @param x_limit,y_limit Optional reviewed axis limits. Limits that clip data
#'   or the testing threshold error. `x_limit` must be symmetric around zero.
#' @param missing What to do when selected phenotypes are unavailable.
#' @param title Optional in-plot title.
#' @param label_mode Label text: phenotype `description`, `phenotype`, or
#'   `"both"`.
#' @param top_label_n Number of strongest associations to label automatically
#'   when `label` is `NULL`. Use zero to disable automatic labels.
#' @param label_wrap_width Approximate characters per label line.
#' @param group_legend_columns Number of columns in the bottom group legend.
#' @param group_wrap_width Approximate characters per group-label line.
#' @param point_size Point size in millimetres.
#' @param base_size,label_size Base and association-label font sizes in points.
#'
#' @return A `ggplot` object with a `phewas_audit` attribute.
#' @export
plot_phewas_volcano <- function(
    x,
    significance = "bonferroni",
    label = NULL,
    highlight = NULL,
    category_order = NULL,
    x_limit = NULL,
    y_limit = NULL,
    missing = c("error", "drop"),
    title = NULL,
    label_mode = c("description", "phenotype", "both"),
    top_label_n = 0L,
    label_wrap_width = 26L,
    group_legend_columns = 5L,
    group_wrap_width = 24L,
    point_size = 2.2,
    base_size = 12,
    label_size = base_size) {
  significance <- match.arg(significance, c("bonferroni", "bh", "none"))
  missing <- match.arg(missing)
  label_mode <- match.arg(label_mode)
  if (!is.null(title) &&
      (length(title) != 1L || is.na(title) || !is.character(title))) {
    stop("title must be NULL or one nonmissing character value.", call. = FALSE)
  }
  integer_settings <- list(
    top_label_n = top_label_n,
    label_wrap_width = label_wrap_width,
    group_legend_columns = group_legend_columns,
    group_wrap_width = group_wrap_width
  )
  valid_integer <- vapply(integer_settings, function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value == floor(value)
  }, logical(1))
  if (any(!valid_integer)) {
    stop("Label and group count settings must be finite whole numbers.",
         call. = FALSE)
  }
  if (top_label_n < 0L || top_label_n > 50L) {
    stop("top_label_n must be between 0 and 50.", call. = FALSE)
  }
  if (label_wrap_width < 8L || group_wrap_width < 8L ||
      group_legend_columns < 1L) {
    stop("Wrap widths must be at least 8 and group_legend_columns at least 1.",
         call. = FALSE)
  }
  numeric_settings <- list(
    point_size = point_size, base_size = base_size, label_size = label_size
  )
  valid_numeric <- vapply(numeric_settings, function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value > 0
  }, logical(1))
  if (any(!valid_numeric)) {
    stop("Plot sizes must be finite and greater than zero.", call. = FALSE)
  }

  all_data <- .pf_plot_result(x)
  needed <- c("estimate", "effect_measure")
  absent <- setdiff(needed, names(all_data))
  if (length(absent)) {
    stop("Result is missing volcano column(s): ", paste(absent, collapse = ", "),
         ".", call. = FALSE)
  }
  all_data$estimate <- suppressWarnings(as.numeric(all_data$estimate))
  all_data$effect_measure <- trimws(as.character(all_data$effect_measure))
  all_data$neg_log10_p <- suppressWarnings(as.numeric(all_data$neg_log10_p))
  all_data$q_value <- suppressWarnings(as.numeric(all_data$q_value))
  group_data <- all_data[
    all_data$status == "ok" & is.finite(all_data$neg_log10_p), , drop = FALSE
  ]
  category_order <- .pf_plot_group_order(group_data, category_order)
  color_by_group <- .pf_group_colors(group_data, category_order)
  data <- group_data[is.finite(group_data$estimate), , drop = FALSE]
  if (!nrow(data)) {
    stop("No successful finite associations are available to plot.",
         call. = FALSE)
  }

  measures <- unique(data$effect_measure)
  invalid_measures <- setdiff(measures, names(.pf_volcano_measure_labels))
  if (anyNA(measures) || any(!nzchar(measures)) || length(invalid_measures)) {
    shown <- unique(c(
      invalid_measures,
      if (anyNA(measures) || any(!nzchar(measures))) "missing" else character()
    ))
    stop(
      "Volcano plots require supported effect_measure values; found: ",
      paste(shown, collapse = ", "), ".", call. = FALSE
    )
  }
  measure_order <- names(.pf_volcano_measure_labels)
  measure_order <- measure_order[measure_order %in% measures]
  panel_order <- unname(.pf_volcano_measure_labels[measure_order])
  data$effect_panel <- factor(
    unname(.pf_volcano_measure_labels[data$effect_measure]),
    levels = panel_order
  )
  data$x_plot <- data$estimate
  data$y_plot <- data$neg_log10_p
  data$effect_direction <- .pf_effect_direction(data)

  alpha <- if (identical(significance, "none")) {
    tryCatch(.pf_fdr_threshold(x, x, NULL), error = function(error) NA_real_)
  } else {
    .pf_fdr_threshold(x, x, NULL)
  }
  thresholds <- .pf_volcano_threshold(all_data, significance, alpha)
  thresholds$y_plot <- thresholds$y

  x_limit <- .pf_check_limits(data$x_plot, x_limit, "x-axis")
  if (!is.null(x_limit)) {
    tolerance <- sqrt(.Machine$double.eps) * max(1, abs(x_limit))
    if (x_limit[[1L]] >= 0 || x_limit[[2L]] <= 0 ||
        abs(abs(x_limit[[1L]]) - abs(x_limit[[2L]])) > max(tolerance)) {
      stop("x-axis limits must be symmetric around zero.", call. = FALSE)
    }
  }
  y_limit <- .pf_check_limits(
    c(data$y_plot, thresholds$y), y_limit, "y-axis"
  )

  highlight <- .pf_select_ids(data, highlight, "highlight", missing)
  label <- .pf_volcano_select_labels(
    data, label, as.integer(top_label_n), significance, thresholds, alpha,
    missing
  )
  data$plot_label <- NA_character_
  labeled <- data$phenotype %in% label
  data$plot_label[labeled] <- .pf_manhattan_label_text(
    data[labeled, , drop = FALSE], label_mode, as.integer(label_wrap_width)
  )

  data$group <- factor(data$group, levels = category_order)
  data <- data[order(
    data$y_plot, data$effect_panel, data$group, data$phenotype,
    method = "radix"
  ), , drop = FALSE]

  if (is.null(x_limit)) {
    symmetric_limits <- do.call(rbind, lapply(
      split(data, data$effect_panel, drop = TRUE),
      function(panel_data) {
        bound <- max(abs(panel_data$x_plot), na.rm = TRUE)
        if (!is.finite(bound) || bound == 0) bound <- 1
        data.frame(
          effect_panel = panel_data$effect_panel[[1L]],
          x_plot = c(-bound, bound), y_plot = 0,
          stringsAsFactors = FALSE
        )
      }
    ))
    symmetric_limits$effect_panel <- factor(
      symmetric_limits$effect_panel, levels = panel_order
    )
  } else {
    symmetric_limits <- data.frame(
      effect_panel = factor(rep(panel_order, each = 2L), levels = panel_order),
      x_plot = rep(x_limit, times = length(panel_order)),
      y_plot = 0,
      stringsAsFactors = FALSE
    )
  }

  y_max <- max(c(group_data$neg_log10_p, thresholds$y), na.rm = TRUE)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  display_y_limit <- if (is.null(y_limit)) c(0, y_max * 1.12) else y_limit
  y_axis <- .pf_plot_y_axis(
    display_y_limit, snap_upper = is.null(y_limit)
  )
  display_y_limit <- y_axis$limits
  y_breaks <- y_axis$breaks
  y_labels <- format(
    signif(y_breaks, digits = 5L), trim = TRUE, scientific = FALSE
  )

  effect_levels <- levels(droplevels(data$effect_direction))
  effect_shapes <- c(
    "Positive effect" = 24,
    "Negative effect" = 25,
    "Unknown effect" = 21
  )

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = x_plot, y = y_plot)) +
    ggplot2::geom_blank(
      data = symmetric_limits,
      ggplot2::aes(x = x_plot, y = y_plot), inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(
      yintercept = 0, color = "#B8B8B8", linewidth = 0.4
    ) +
    ggplot2::geom_vline(
      xintercept = 0, color = "#8A8A8A", linewidth = 0.55
    )
  if (nrow(thresholds)) {
    line_types <- c("BH FDR" = "dotdash", "Bonferroni" = "dashed")
    line_colors <- c("BH FDR" = "#536B82", "Bonferroni" = "#222222")
    plot <- plot +
      ggplot2::geom_hline(
        data = thresholds,
        ggplot2::aes(
          yintercept = y_plot, linetype = method, color = method
        ),
        inherit.aes = FALSE, linewidth = 0.75, show.legend = FALSE
      ) +
      ggplot2::scale_linetype_manual(values = line_types, guide = "none") +
      ggplot2::scale_color_manual(values = line_colors, guide = "none")
  }
  plot <- plot +
    ggplot2::geom_point(
      ggplot2::aes(fill = group, shape = effect_direction),
      color = "#222222", size = point_size, stroke = 0.35, alpha = 0.88
    ) +
    .pf_phewas_group_scale(
      color_by_group,
      category_order,
      group_wrap_width,
      group_legend_columns,
      show_legend = TRUE
    ) +
    ggplot2::scale_shape_manual(
      values = effect_shapes,
      breaks = effect_levels,
      name = NULL,
      drop = TRUE
    ) +
    ggplot2::scale_x_continuous(
      limits = x_limit,
      expand = ggplot2::expansion(mult = c(0.06, 0.06))
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = display_y_limit,
      expand = ggplot2::expansion(mult = c(0, 0.015))
    ) +
    ggplot2::labs(
      title = title,
      x = if (length(panel_order) == 1L) panel_order[[1L]] else
        "Estimate on fitted model scale",
      y = expression(-log[10](italic(p)))
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        order = 2,
        nrow = 1,
        override.aes = list(
          fill = "#808080", color = "#222222", size = 3.8,
          stroke = 0.4, alpha = 1
        )
      )
    ) +
    .pf_phewas_theme(base_size, show_x_text = TRUE) +
    ggplot2::coord_cartesian(clip = "off")

  if (length(panel_order) > 1L) {
    plot <- plot + ggplot2::facet_wrap(
      ggplot2::vars(effect_panel), scales = "free_x"
    )
  }
  if (nrow(thresholds)) {
    plot <- plot + ggplot2::geom_label(
      data = thresholds,
      ggplot2::aes(x = Inf, y = y_plot, label = threshold_label),
      inherit.aes = FALSE,
      hjust = 1.03, vjust = -0.3,
      size = base_size * 0.78 / 2.845276,
      color = "#333333", fill = "#FFFFFFE6", linewidth = 0.25,
      label.padding = grid::unit(0.12, "lines"),
      label.r = grid::unit(0.08, "lines"),
      show.legend = FALSE
    )
  }
  plot <- .pf_phewas_add_highlights(plot, data, highlight, point_size)
  plot <- .pf_phewas_add_labels(plot, data, label, label_size)

  attr(plot, "phewas_audit") <- list(
    plot = "volcano",
    phenotypes = data$phenotype,
    labels = label,
    label_mode = label_mode,
    highlights = highlight,
    significance = significance,
    threshold = thresholds[c("method", "y")],
    testing_family_size = .pf_plot_family_size(all_data),
    effect_measures = measure_order,
    effect_panels = panel_order,
    group_order = category_order,
    x_limit = x_limit,
    y_limit = y_limit
  )
  plot
}
