.plot_global_variables <- c(
  "x", "y", "neg_log10_p", "group", "yintercept", "method",
  "selected_description", "x_plot", "y_plot", "group_x",
  "plot_description", "count", "set", "description", "run",
  "native_effect", "native_conf_low", "native_conf_high", "effect_panel",
  "xintercept", "effect_direction", "plot_label", "threshold_label"
)
utils::globalVariables(.plot_global_variables)

.pf_select_ids <- function(data, ids, argument, missing = c("error", "drop")) {
  missing <- match.arg(missing)
  if (is.null(ids)) return(character())
  ids <- unique(as.character(ids))
  absent <- setdiff(ids, data$phenotype)
  if (length(absent) && missing == "error") {
    stop(argument, " contains phenotype value(s) unavailable for this plot: ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }
  intersect(ids, data$phenotype)
}

.pf_check_limits <- function(values, limits, axis = "axis") {
  if (is.null(limits)) return(NULL)
  if (!is.numeric(limits) || length(limits) != 2L || any(!is.finite(limits)) ||
      limits[[1L]] >= limits[[2L]]) {
    stop(axis, " limits must be two increasing finite numbers.", call. = FALSE)
  }
  values <- values[is.finite(values)]
  if (length(values) && (min(values) < limits[[1L]] || max(values) > limits[[2L]])) {
    stop(axis, " limits would clip plotted data; review or remove the limits.",
         call. = FALSE)
  }
  as.numeric(limits)
}

.pf_natural_order <- function(x) {
  x <- as.character(x)
  first_number <- regexpr("[0-9]+(?:[.][0-9]+)?", x, perl = TRUE)
  prefix <- ifelse(
    first_number > 0L,
    tolower(substr(x, 1L, pmax(0L, first_number - 1L))),
    tolower(x)
  )
  number <- rep(NA_real_, length(x))
  has_number <- first_number > 0L
  match_length <- attr(first_number, "match.length")
  number[has_number] <- suppressWarnings(as.numeric(substr(
    x[has_number], first_number[has_number],
    first_number[has_number] + match_length[has_number] - 1L
  )))
  order(prefix, is.na(number), number, tolower(x), x, method = "radix")
}

.pf_plot_wrap <- function(x, width) {
  vapply(as.character(x), function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    paragraphs <- strsplit(value, "\n", fixed = TRUE)[[1L]]
    wrapped <- unlist(lapply(paragraphs, function(paragraph) {
      lines <- strwrap(paragraph, width = width)
      if (length(lines)) lines else ""
    }), use.names = FALSE)
    paste(wrapped, collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

.pf_valid_plot_color <- function(x) {
  if (is.na(x) || !nzchar(x)) return(FALSE)
  isTRUE(tryCatch({
    grDevices::col2rgb(x)
    TRUE
  }, error = function(e) FALSE))
}

.pf_group_colors <- function(data, group_order) {
  colors <- stats::setNames(rep(NA_character_, length(group_order)), group_order)
  for (group_name in group_order) {
    declared <- unique(data$color[
      data$group == group_name & !is.na(data$color) & nzchar(data$color)
    ])
    if (length(declared) > 1L) {
      stop("Each group must have at most one declared color.", call. = FALSE)
    }
    if (length(declared)) {
      if (!.pf_valid_plot_color(declared[[1L]])) {
        stop("Group `", group_name, "` has an invalid color: ",
             declared[[1L]], ".", call. = FALSE)
      }
      colors[[group_name]] <- declared[[1L]]
    }
  }
  missing <- is.na(colors)
  if (any(missing)) {
    fallback <- grDevices::hcl.colors(
      max(3L, sum(missing)), palette = "Dark 3"
    )
    colors[missing] <- fallback[seq_len(sum(missing))]
  }
  colors
}

.pf_plot_group_order <- function(data, category_order = NULL) {
  observed <- unique(as.character(data$group))
  if (!is.null(category_order)) {
    category_order <- as.character(category_order)
    if (anyNA(category_order) || any(!nzchar(category_order)) ||
        anyDuplicated(category_order) || !setequal(category_order, observed)) {
      stop(
        "category_order must contain every plotted group exactly once.",
        call. = FALSE
      )
    }
    return(category_order)
  }

  group_numbers <- split(data$groupnum, data$group)
  has_one_number <- vapply(group_numbers, function(value) {
    value <- unique(value[is.finite(value)])
    length(value) == 1L
  }, logical(1))
  if (all(has_one_number)) {
    order_map <- vapply(group_numbers, function(value) {
      unique(value[is.finite(value)])[[1L]]
    }, numeric(1))
    return(names(sort(order_map, method = "radix")))
  }
  sort(observed, method = "radix")
}

.pf_phewas_group_guide <- function(columns, order = 1L) {
  ggplot2::guide_legend(
    order = order,
    ncol = as.integer(columns),
    byrow = TRUE,
    override.aes = list(
      shape = 21, color = "#242424", size = 3.8,
      stroke = 0.35, alpha = 1
    )
  )
}

.pf_plot_y_axis <- function(limits, n = 5L, snap_upper = FALSE) {
  limits <- suppressWarnings(as.numeric(limits))
  if (length(limits) != 2L || any(!is.finite(limits)) ||
      limits[[1L]] >= limits[[2L]]) {
    stop("Plot limits must be two increasing finite numbers.", call. = FALSE)
  }
  limits <- signif(limits, digits = 12L)
  breaks <- pretty(limits, n = as.integer(n))
  spacing <- diff(sort(unique(breaks)))
  spacing <- spacing[is.finite(spacing) & spacing > 0]
  if (isTRUE(snap_upper) && length(spacing)) {
    candidates <- breaks[breaks >= limits[[2L]]]
    if (length(candidates) &&
        candidates[[1L]] - limits[[2L]] <= min(spacing) * 0.01) {
      limits[[2L]] <- candidates[[1L]]
    }
  }
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(limits))
  breaks <- breaks[
    breaks >= limits[[1L]] - max(tolerance) &
      breaks <= limits[[2L]] + max(tolerance)
  ]
  list(limits = limits, breaks = breaks)
}

.pf_plot_index_breaks <- function(n, target = 5L) {
  n <- as.integer(n)
  target <- as.integer(target)
  if (n < 1L || target < 1L) return(integer())
  unique(as.integer(floor(seq.int(1L, n, length.out = min(n, target)))))
}

.pf_phewas_group_scale <- function(colors, group_order, wrap_width, columns,
                                   show_legend = TRUE) {
  ggplot2::scale_fill_manual(
    values = colors,
    breaks = group_order,
    labels = .pf_plot_wrap(group_order, as.integer(wrap_width)),
    name = "Phenotype group",
    guide = if (isTRUE(show_legend)) {
      .pf_phewas_group_guide(columns)
    } else {
      "none"
    },
    drop = FALSE
  )
}

.pf_phewas_theme <- function(base_size, show_x_text = TRUE) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(
        color = "#E8E8E8", linewidth = 0.4
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.x = grid::unit(1.25, "lines"),
      axis.line = ggplot2::element_line(color = "#222222", linewidth = 0.7),
      axis.ticks = ggplot2::element_line(color = "#222222", linewidth = 0.65),
      axis.ticks.length = grid::unit(5, "pt"),
      axis.text.x = if (isTRUE(show_x_text)) {
        ggplot2::element_text(
          size = base_size, hjust = 0.5, vjust = 0.5,
          lineheight = 0.92, color = "#222222"
        )
      } else {
        ggplot2::element_blank()
      },
      axis.ticks.x = if (isTRUE(show_x_text)) {
        ggplot2::element_line(color = "#222222", linewidth = 0.65)
      } else {
        ggplot2::element_blank()
      },
      axis.text.y = ggplot2::element_text(size = base_size, color = "#222222"),
      axis.title = ggplot2::element_text(
        color = "#111111", face = "bold", size = base_size * 1.05
      ),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 9)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      strip.background = ggplot2::element_rect(
        fill = "#F2F2F2", color = "#D6D6D6", linewidth = 0.45
      ),
      strip.text = ggplot2::element_text(
        face = "bold", color = "#222222", size = base_size
      ),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.direction = "horizontal",
      legend.text = ggplot2::element_text(size = base_size * 0.88),
      legend.title = ggplot2::element_text(
        face = "bold", size = base_size * 0.92, hjust = 0.5
      ),
      legend.title.position = "top",
      legend.key.width = grid::unit(1.15, "lines"),
      legend.spacing.x = grid::unit(0.5, "lines"),
      plot.title = ggplot2::element_text(
        face = "bold", color = "#111111", size = base_size * 1.3,
        margin = ggplot2::margin(b = 10)
      ),
      plot.margin = ggplot2::margin(t = 8, r = 24, b = 8, l = 10)
    )
}

.pf_phewas_add_highlights <- function(plot, data, highlights, point_size) {
  if (!length(highlights)) return(plot)
  highlight_data <- data[data$phenotype %in% highlights, , drop = FALSE]
  ring_size <- max(4.2, point_size * 2)
  plot +
    ggplot2::geom_point(
      data = highlight_data,
      ggplot2::aes(x = x_plot, y = y_plot),
      inherit.aes = FALSE, shape = 21, fill = NA, color = "white",
      size = ring_size, stroke = 1.5, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = highlight_data,
      ggplot2::aes(x = x_plot, y = y_plot),
      inherit.aes = FALSE, shape = 21, fill = NA, color = "#111111",
      size = ring_size, stroke = 0.85, show.legend = FALSE
    )
}

.pf_phewas_add_labels <- function(plot, data, labels, label_size) {
  if (!length(labels)) return(plot)
  label_data <- data[data$phenotype %in% labels, , drop = FALSE]
  plot + ggrepel::geom_label_repel(
    data = label_data,
    ggplot2::aes(x = x_plot, y = y_plot, label = plot_label),
    inherit.aes = FALSE,
    seed = 20260819,
    size = label_size / 2.845276,
    color = "#111111",
    fill = "#FFFFFFF2",
    linewidth = 0.28,
    label.padding = grid::unit(0.18, "lines"),
    label.r = grid::unit(0.12, "lines"),
    box.padding = 0.65,
    point.padding = grid::unit(0.6, "lines"),
    min.segment.length = 0,
    segment.color = "#444444",
    segment.size = 0.45,
    force = 1.4,
    force_pull = 0.45,
    direction = "both",
    max.overlaps = Inf,
    max.iter = 50000,
    max.time = 10,
    show.legend = FALSE
  )
}

.pf_effect_direction <- function(data) {
  estimate <- suppressWarnings(as.numeric(.pf_column(data, "estimate", NA_real_)))
  direction <- rep("Unknown effect", nrow(data))
  direction[is.finite(estimate) & estimate > 0] <- "Positive effect"
  direction[is.finite(estimate) & estimate < 0] <- "Negative effect"

  unresolved <- direction == "Unknown effect"
  if (any(unresolved)) {
    native <- suppressWarnings(as.numeric(.pf_column(
      data, "native_effect", NA_real_
    )))
    measure <- as.character(.pf_column(data, "effect_measure", NA_character_))
    difference <- measure %in% c("mean_difference", "difference", "beta")
    use_difference <- unresolved & is.finite(native) & difference
    direction[use_difference & native > 0] <- "Positive effect"
    direction[use_difference & native < 0] <- "Negative effect"
    use_ratio <- unresolved & is.finite(native) & !difference
    direction[use_ratio & native > 1] <- "Positive effect"
    direction[use_ratio & native < 1] <- "Negative effect"
  }

  factor(
    direction,
    levels = c("Positive effect", "Negative effect", "Unknown effect")
  )
}

.pf_manhattan_axis_scale <- function(values, mode, break_top, top_size) {
  values <- suppressWarnings(as.numeric(values))
  finite_values <- values[is.finite(values)]
  max_y <- if (length(finite_values)) max(finite_values) else 1
  max_y <- max(1, max_y)
  request_break <- identical(mode, "on") ||
    (identical(mode, "auto") && max_y > break_top / (1 - top_size))
  enabled <- request_break && break_top / (1 - top_size) <= max_y

  rescale <- function(value) suppressWarnings(as.numeric(value))
  if (enabled) {
    top_plot_range <- break_top / (1 - top_size) - break_top
    top_data_range <- max_y - break_top
    rescale <- function(value) {
      value <- suppressWarnings(as.numeric(value))
      output <- value
      above <- is.finite(value) & value > break_top
      output[above] <- break_top +
        (value[above] - break_top) / (top_data_range / top_plot_range)
      output
    }
  }

  lower_top <- if (enabled) break_top else max_y
  tick_values <- pretty(c(0, lower_top), n = 5)
  tick_values <- tick_values[tick_values >= 0 & tick_values <= lower_top]
  if (enabled && !any(abs(tick_values - break_top) < 1e-8)) {
    tick_values <- c(tick_values, break_top)
  }
  if (enabled) tick_values <- c(tick_values, max_y)
  tick_positions <- rescale(tick_values)
  ordering <- order(tick_positions)
  tick_values <- tick_values[ordering]
  tick_positions <- tick_positions[ordering]
  keep <- !duplicated(round(tick_positions, digits = 8L))

  list(
    enabled = enabled,
    max_y = max_y,
    break_top = break_top,
    top_size = top_size,
    rescale = rescale,
    tick_values = tick_values[keep],
    tick_positions = tick_positions[keep]
  )
}

.pf_manhattan_labels <- function(data, label, top_label_n, cutoff, missing) {
  manual <- !is.null(label)
  if (manual) return(.pf_select_ids(data, label, "label", missing))
  if (!top_label_n) return(character())

  significant <- if (is.finite(cutoff)) {
    is.finite(data$q_value) & data$q_value <= cutoff
  } else {
    rep(FALSE, nrow(data))
  }
  neg_log10_q <- suppressWarnings(as.numeric(.pf_column(
    data, "neg_log10_q", NA_real_
  )))
  neg_log10_q[!is.finite(neg_log10_q)] <- -Inf
  order_index <- order(
    !significant, -neg_log10_q, -data$neg_log10_p, data$phenotype,
    method = "radix"
  )
  data$phenotype[utils::head(order_index, top_label_n)]
}

.pf_manhattan_label_text <- function(data, label_mode, wrap_width) {
  description <- as.character(data$description)
  missing_description <- is.na(description) | !nzchar(description)
  description[missing_description] <- data$phenotype[missing_description]
  text <- switch(
    label_mode,
    description = description,
    phenotype = data$phenotype,
    both = ifelse(
      description == data$phenotype,
      data$phenotype,
      paste0(description, "\n(", data$phenotype, ")")
    )
  )
  .pf_plot_wrap(text, wrap_width)
}

.pf_plot_result <- function(x) {
  data <- .pf_associations(x, "x")
  needed <- c("phenotype", "status", "p_value", "q_value", "neg_log10_p")
  absent <- setdiff(needed, names(data))
  if (length(absent)) {
    stop("Result is missing plot column(s): ", paste(absent, collapse = ", "),
         ".", call. = FALSE)
  }
  data
}

.pf_plot_family_size <- function(data) {
  recorded <- suppressWarnings(as.numeric(.pf_column(
    data, "testing_family_size", NA_real_
  )))
  present <- !is.na(recorded)
  if (any(present & (!is.finite(recorded) | recorded < 1 |
                     recorded != floor(recorded)))) {
    stop("The result records an invalid testing_family_size value.",
         call. = FALSE)
  }
  recorded <- unique(recorded[present])
  if (length(recorded) > 1L) {
    stop("The result records conflicting testing_family_size values.",
         call. = FALSE)
  }
  observed <- sum(data$status == "ok", na.rm = TRUE)
  if (length(recorded)) {
    if (recorded[[1L]] < observed) {
      stop("testing_family_size is smaller than the successful result family.",
           call. = FALSE)
    }
    return(as.integer(recorded[[1L]]))
  }
  observed
}

#' Plot a PheWAS Manhattan-style result
#'
#' Uses canonical group colors, effect-direction triangles, group separators,
#' directly annotated testing thresholds, boxed labels, and bottom legends.
#' Stored `neg_log10_p` values remain authoritative, so probabilities that
#' underflowed to zero are still plottable. Group and phenotype order are
#' deterministic.
#'
#' @param x A `phewas_result`.
#' @param significance Any of `"bh"`, `"bonferroni"`, and `"none"`.
#' @param label,highlight Explicit phenotype values to label or ring.
#' @param category_order Optional character vector giving every group in
#'   display order. The argument name is retained for API compatibility;
#'   plotted metadata use the canonical `group` field.
#' @param y_limit Optional reviewed limits on the original `-log10(p)` scale.
#'   Limits that clip data or testing thresholds error.
#' @param missing What to do when selected phenotypes are unavailable.
#' @param title Optional in-plot title. `NULL` leaves room for an external
#'   report or slide title.
#' @param label_mode Label text: phenotype `description`, `phenotype`, or
#'   `"both"`.
#' @param top_label_n Number of automatically selected labels when `label` is
#'   `NULL`. Significant associations rank first; use zero to disable.
#' @param label_wrap_width Approximate characters per label line.
#' @param group_display Where to show phenotype groups. The default `"legend"`
#'   matches the volcano plot and leaves a numeric phenotype index on the
#'   Manhattan x-axis. `"auto"` uses centered group labels for compact group
#'   sets and a legend for crowded group sets; `"x_axis"` always places group
#'   labels on the Manhattan x-axis.
#' @param group_legend_columns Number of columns in the group legend.
#' @param group_label_rows Rows used to dodge group labels on the x-axis.
#' @param group_wrap_width Approximate characters per group-label line.
#' @param point_size Point size in millimetres.
#' @param base_size,label_size Base and association-label font sizes in points.
#' @param axis_break_mode Optional upper-tail compression: `"off"`, `"auto"`,
#'   or `"on"`.
#' @param axis_break_top Compression threshold on the original `-log10(p)`
#'   scale.
#' @param axis_break_top_size Fraction of the plotting range allotted above
#'   the compression threshold.
#'
#' @return A `ggplot` object with a `phewas_audit` attribute.
#' @export
plot_phewas_manhattan <- function(
    x,
    significance = "bonferroni",
    label = NULL,
    highlight = NULL,
    category_order = NULL,
    y_limit = NULL,
    missing = c("error", "drop"),
    title = NULL,
    label_mode = c("description", "phenotype", "both"),
    top_label_n = 0L,
    label_wrap_width = 26L,
    group_display = c("legend", "auto", "x_axis"),
    group_legend_columns = 5L,
    group_label_rows = 3L,
    group_wrap_width = 24L,
    point_size = 2.2,
    base_size = 12,
    label_size = base_size,
    axis_break_mode = c("off", "auto", "on"),
    axis_break_top = 55,
    axis_break_top_size = 0.125) {
  missing <- match.arg(missing)
  label_mode <- match.arg(label_mode)
  group_display_requested <- match.arg(group_display)
  axis_break_mode <- match.arg(axis_break_mode)
  significance <- unique(match.arg(
    significance, c("bh", "bonferroni", "none"), several.ok = TRUE
  ))
  if ("none" %in% significance && length(significance) > 1L) {
    stop("`none` cannot be combined with significance thresholds.", call. = FALSE)
  }
  if (!is.null(title) &&
      (length(title) != 1L || is.na(title) || !is.character(title))) {
    stop("title must be NULL or one nonmissing character value.", call. = FALSE)
  }
  integer_settings <- list(
    top_label_n = top_label_n,
    label_wrap_width = label_wrap_width,
    group_legend_columns = group_legend_columns,
    group_label_rows = group_label_rows,
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
      group_legend_columns < 1L || group_label_rows < 1L) {
    stop("Wrap widths must be at least 8 and group layout counts at least 1.",
         call. = FALSE)
  }
  numeric_settings <- list(
    point_size = point_size, base_size = base_size, label_size = label_size,
    axis_break_top = axis_break_top
  )
  valid_numeric <- vapply(numeric_settings, function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value > 0
  }, logical(1))
  if (any(!valid_numeric)) {
    stop("Plot sizes and axis_break_top must be finite and greater than zero.",
         call. = FALSE)
  }
  if (!is.numeric(axis_break_top_size) || length(axis_break_top_size) != 1L ||
      is.na(axis_break_top_size) || !is.finite(axis_break_top_size) ||
      axis_break_top_size <= 0 ||
      axis_break_top_size >= 1) {
    stop("axis_break_top_size must be strictly between zero and one.",
         call. = FALSE)
  }

  data <- .pf_plot_result(x)
  successful <- data$status == "ok"
  testing_family_size <- .pf_plot_family_size(data)
  data <- data[
    successful & is.finite(data$neg_log10_p), , drop = FALSE
  ]
  if (!nrow(data)) {
    stop("No successful finite associations are available to plot.",
         call. = FALSE)
  }
  highlight <- .pf_select_ids(data, highlight, "highlight", missing)

  category_order <- .pf_plot_group_order(data, category_order)
  color_by_group <- .pf_group_colors(data, category_order)
  data$group <- factor(data$group, levels = category_order)
  phenotype_order <- suppressWarnings(as.numeric(.pf_column(
    data, "phenotype_order", NA_real_
  )))
  phenotype_order[!is.finite(phenotype_order)] <- Inf
  natural_rank <- integer(nrow(data))
  natural_rank[.pf_natural_order(data$phenotype)] <- seq_len(nrow(data))
  data <- data[order(
    data$group, phenotype_order, natural_rank, data$phenotype,
    method = "radix"
  ), , drop = FALSE]
  data$x <- seq_len(nrow(data))
  data$x_plot <- data$x
  data$effect_direction <- .pf_effect_direction(data)

  cutoff <- tryCatch(
    .pf_fdr_threshold(x, x, NULL), error = function(e) NA_real_
  )
  label <- .pf_manhattan_labels(
    data, label, as.integer(top_label_n), cutoff, missing
  )
  data$plot_label <- NA_character_
  labeled <- data$phenotype %in% label
  data$plot_label[labeled] <- .pf_manhattan_label_text(
    data[labeled, , drop = FALSE], label_mode, as.integer(label_wrap_width)
  )

  thresholds <- data.frame(
    method = character(), y = numeric(), stringsAsFactors = FALSE
  )
  if ("bh" %in% significance) {
    if (!is.finite(cutoff)) {
      stop("The result does not record an FDR threshold.", call. = FALSE)
    }
    passing <- is.finite(data$q_value) & data$q_value <= cutoff
    if (any(passing)) {
      thresholds <- rbind(thresholds, data.frame(
        method = "BH FDR", y = min(data$neg_log10_p[passing]),
        stringsAsFactors = FALSE
      ))
    }
  }
  if ("bonferroni" %in% significance) {
    if (!is.finite(cutoff)) {
      stop("The result does not record a significance threshold.", call. = FALSE)
    }
    thresholds <- rbind(thresholds, data.frame(
      method = "Bonferroni", y = -log10(cutoff / testing_family_size),
      stringsAsFactors = FALSE
    ))
  }

  y_limit <- .pf_check_limits(
    c(data$neg_log10_p, thresholds$y), y_limit, "y-axis"
  )
  axis_scale <- .pf_manhattan_axis_scale(
    c(data$neg_log10_p, thresholds$y, y_limit),
    mode = axis_break_mode,
    break_top = axis_break_top,
    top_size = axis_break_top_size
  )
  data$y_plot <- axis_scale$rescale(data$neg_log10_p)
  thresholds$y_plot <- axis_scale$rescale(thresholds$y)

  spans <- do.call(rbind, lapply(
    split(data, data$group, drop = TRUE),
    function(z) data.frame(
      group = as.character(z$group[[1L]]),
      start = min(z$x),
      end = max(z$x),
      midpoint = mean(range(z$x)),
      stringsAsFactors = FALSE
    )
  ))
  separators <- if (nrow(spans) > 1L) spans$end[-nrow(spans)] + 0.5 else numeric()
  group_labels <- .pf_plot_wrap(spans$group, as.integer(group_wrap_width))
  group_display <- if (identical(group_display_requested, "auto")) {
    if (nrow(spans) > 12L) "legend" else "x_axis"
  } else {
    group_display_requested
  }

  y_plot_max <- max(
    c(data$y_plot, thresholds$y_plot, axis_scale$rescale(y_limit)),
    na.rm = TRUE
  )
  if (!is.finite(y_plot_max) || y_plot_max <= 0) y_plot_max <- 1
  threshold_annotations <- data.frame(
    y_plot = numeric(), threshold_label = character(), stringsAsFactors = FALSE
  )
  if (nrow(thresholds)) {
    remaining <- seq_len(nrow(thresholds))
    tolerance <- 0.035 * y_plot_max
    while (length(remaining)) {
      seed <- remaining[[1L]]
      cluster <- remaining[
        abs(thresholds$y_plot[remaining] - thresholds$y_plot[[seed]]) <=
          tolerance
      ]
      threshold_labels <- paste0(
        thresholds$method[cluster], ": ",
        format(signif(thresholds$y[cluster], digits = 3L), trim = TRUE)
      )
      threshold_annotations <- rbind(
        threshold_annotations,
        data.frame(
          y_plot = max(thresholds$y_plot[cluster]),
          threshold_label = paste(threshold_labels, collapse = "  |  "),
          stringsAsFactors = FALSE
        )
      )
      remaining <- setdiff(remaining, cluster)
    }
  }
  display_limits <- if (is.null(y_limit)) {
    c(0, y_plot_max * 1.12)
  } else {
    axis_scale$rescale(y_limit)
  }
  if (axis_scale$enabled) {
    axis_positions <- axis_scale$tick_positions
    axis_values <- axis_scale$tick_values
  } else {
    y_axis <- .pf_plot_y_axis(
      display_limits, snap_upper = is.null(y_limit)
    )
    display_limits <- y_axis$limits
    axis_positions <- y_axis$breaks
    axis_values <- axis_positions
  }
  axis_labels <- format(
    signif(axis_values, digits = 5L), trim = TRUE, scientific = FALSE
  )
  if (axis_scale$enabled) {
    at_break <- abs(axis_values - axis_scale$break_top) < 1e-8
    axis_labels[at_break] <- paste0(axis_labels[at_break], "  //")
  }
  axis_order <- order(axis_positions)
  axis_positions <- axis_positions[axis_order]
  axis_labels <- axis_labels[axis_order]

  effect_levels <- levels(droplevels(data$effect_direction))
  effect_shapes <- c(
    "Positive effect" = 24,
    "Negative effect" = 25,
    "Unknown effect" = 21
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = x_plot, y = y_plot))
  if (length(separators)) {
    plot <- plot + ggplot2::geom_vline(
      xintercept = separators, color = "#D8D8D8", linewidth = 0.45
    )
  }
  if (nrow(thresholds)) {
    line_types <- c("BH FDR" = "dotdash", "Bonferroni" = "dashed")
    line_colors <- c("BH FDR" = "#536B82", "Bonferroni" = "#222222")
    plot <- plot +
      ggplot2::geom_hline(
        data = thresholds,
        ggplot2::aes(yintercept = y_plot, linetype = method, color = method),
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
      show_legend = identical(group_display, "legend")
    ) +
    ggplot2::scale_shape_manual(
      values = effect_shapes,
      breaks = effect_levels,
      name = NULL,
      drop = TRUE
    ) +
    ggplot2::scale_x_continuous(
      breaks = if (identical(group_display, "x_axis")) {
        spans$midpoint
      } else {
        .pf_plot_index_breaks(nrow(data))
      },
      labels = if (identical(group_display, "x_axis")) {
        group_labels
      } else {
        as.character(.pf_plot_index_breaks(nrow(data)))
      },
      expand = ggplot2::expansion(mult = c(0.012, 0.012)),
      guide = ggplot2::guide_axis(n.dodge = if (
        identical(group_display, "x_axis")
      ) as.integer(group_label_rows) else 1L)
    ) +
    ggplot2::scale_y_continuous(
      breaks = axis_positions,
      labels = axis_labels,
      limits = display_limits,
      expand = ggplot2::expansion(mult = c(0, 0.015))
    ) +
    ggplot2::labs(
      title = title,
      x = if (identical(group_display, "legend")) {
        "Phenotypes ordered by group"
      } else {
        NULL
      },
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
    .pf_phewas_theme(
      base_size,
      show_x_text = TRUE
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if (nrow(threshold_annotations)) {
    plot <- plot + ggplot2::geom_label(
      data = threshold_annotations,
      ggplot2::aes(x = Inf, y = y_plot, label = threshold_label),
      inherit.aes = FALSE,
      hjust = 1.03,
      vjust = -0.3,
      size = base_size * 0.78 / 2.845276,
      color = "#333333",
      fill = "#FFFFFFE6",
      linewidth = 0.25,
      label.padding = grid::unit(0.12, "lines"),
      label.r = grid::unit(0.08, "lines"),
      show.legend = FALSE
    )
  }

  plot <- .pf_phewas_add_highlights(plot, data, highlight, point_size)
  plot <- .pf_phewas_add_labels(plot, data, label, label_size)

  attr(plot, "phewas_audit") <- list(
    plot = "manhattan",
    phenotypes = data$phenotype,
    labels = label,
    label_mode = label_mode,
    highlights = highlight,
    thresholds = thresholds[c("method", "y")],
    testing_family_size = testing_family_size,
    group_order = category_order,
    group_display_requested = group_display_requested,
    group_display = group_display,
    axis_break = list(
      mode = axis_break_mode,
      enabled = axis_scale$enabled,
      top = axis_break_top,
      top_size = axis_break_top_size
    )
  )
  plot
}

#' Plot concordance between two PheWAS results
#'
#' @param comparison A `phewas_comparison` from [compare_phewas()].
#' @param scale Plot standardized coefficients (`"standardized"`) or compatible
#'   native effects (`"native"`).
#' @param label Phenotype values to label.
#' @param x_limit,y_limit Optional reviewed limits that error if they clip data.
#' @param missing What to do when selected IDs are unavailable.
#'
#' @return A `ggplot` object.
#' @export
plot_phewas_concordance <- function(
    comparison, scale = c("standardized", "native"), label = NULL,
    x_limit = NULL, y_limit = NULL, missing = c("error", "drop")) {
  if (!inherits(comparison, "phewas_comparison")) {
    stop("comparison must be a phewas_comparison.", call. = FALSE)
  }
  scale <- match.arg(scale)
  missing <- match.arg(missing)
  data <- as.data.frame(
    data.table::copy(comparison$matched), stringsAsFactors = FALSE
  )
  if (!nrow(data)) stop("The comparison has no mutually successful phenotypes.", call. = FALSE)
  label <- .pf_select_ids(data, label, "label", missing)
  if (scale == "native") {
    if (!isTRUE(comparison$native_compatible)) {
      stop("Native effects are incompatible: ",
           paste(comparison$native_incompatibility, collapse = ", "), ".",
           call. = FALSE)
    }
    data$x_plot <- data$native_effect_x
    data$y_plot <- data$native_effect_y
    axis_name <- unique(data$effect_measure_x)[[1L]]
    neutral <- if (identical(axis_name, "mean_difference")) 0 else 1
  } else {
    data$x_plot <- data$z_x
    data$y_plot <- data$z_y
    axis_name <- "standardized estimate (estimate / SE)"
    neutral <- 0
  }
  data$plot_description <- ifelse(
    data$phenotype %in% label, data$description_x, NA_character_
  )
  x_limit <- .pf_check_limits(data$x_plot, x_limit, "x-axis")
  y_limit <- .pf_check_limits(data$y_plot, y_limit, "y-axis")
  run_labels <- unname(comparison$run_labels)

  plot <- ggplot2::ggplot(data, ggplot2::aes(x = x_plot, y = y_plot)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey65", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = neutral, color = "grey85", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = neutral, color = "grey85", linewidth = 0.3) +
    ggplot2::geom_point(ggplot2::aes(color = group_x), alpha = 0.8, size = 2) +
    ggplot2::labs(
      x = paste(run_labels[[1L]], axis_name),
      y = paste(run_labels[[2L]], axis_name), color = "Group"
    ) +
    ggplot2::theme_bw()
  if (length(label)) {
    plot <- plot + ggrepel::geom_text_repel(
      data = data[data$phenotype %in% label, , drop = FALSE],
      ggplot2::aes(label = plot_description), size = 3, max.overlaps = Inf,
      show.legend = FALSE
    )
  }
  if (!is.null(x_limit) || !is.null(y_limit)) {
    plot <- plot + ggplot2::coord_cartesian(xlim = x_limit, ylim = y_limit)
  }
  attr(plot, "phewas_audit") <- list(
    plot = paste0("concordance_", scale), phenotypes = data$phenotype,
    labels = label, run_labels = run_labels
  )
  plot
}

#' Plot compatible native effects from a PheWAS comparison
#'
#' @inheritParams plot_phewas_concordance
#' @return A `ggplot` object.
#' @export
plot_phewas_effects <- function(comparison, label = NULL, x_limit = NULL,
                                y_limit = NULL, missing = c("error", "drop")) {
  plot_phewas_concordance(
    comparison, scale = "native", label = label, x_limit = x_limit,
    y_limit = y_limit, missing = missing
  )
}

#' Plot FDR signal overlap by phenotype group
#'
#' @param comparison A `phewas_comparison` from [compare_phewas()].
#'
#' @return A `ggplot` object.
#' @export
plot_phewas_category_overlap <- function(comparison) {
  if (!inherits(comparison, "phewas_comparison")) {
    stop("comparison must be a phewas_comparison.", call. = FALSE)
  }
  source <- comparison$group_overlap
  if (!nrow(source)) stop("The comparison has no group information.", call. = FALSE)
  data <- rbind(
    data.frame(group = source$group, set = "both", count = source$both),
    data.frame(group = source$group, set = "only_x", count = source$only_x),
    data.frame(group = source$group, set = "only_y", count = source$only_y)
  )
  data$group <- factor(data$group, levels = rev(sort(unique(data$group))))
  data$set <- factor(data$set, levels = c("both", "only_x", "only_y"))
  set_labels <- c(
    both = "Both",
    only_x = paste0(comparison$run_labels[["x"]], " only"),
    only_y = paste0(comparison$run_labels[["y"]], " only")
  )
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = count, y = group, fill = set)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_discrete(labels = set_labels) +
    ggplot2::labs(x = "FDR-significant phenotypes", y = NULL, fill = "Signal set") +
    ggplot2::theme_bw()
  attr(plot, "phewas_audit") <- list(
    plot = "group_overlap", fdr_threshold = comparison$fdr_threshold,
    data = data
  )
  plot
}

.pf_forest_table <- function(result, phenotype_ids, run_label, missing) {
  data <- .pf_associations(result, run_label)
  needed <- c("phenotype", "status", "effect_measure", "native_effect",
              "native_conf_low", "native_conf_high", "direction", "reference")
  absent <- setdiff(needed, names(data))
  if (length(absent)) {
    stop("Result `", run_label, "` is missing forest column(s): ",
         paste(absent, collapse = ", "), ".", call. = FALSE)
  }
  available <- data$phenotype[data$status == "ok" &
    is.finite(data$native_effect) & is.finite(data$native_conf_low) &
    is.finite(data$native_conf_high)]
  absent_ids <- setdiff(phenotype_ids, available)
  if (length(absent_ids) && missing == "error") {
    stop("Result `", run_label, "` lacks successful native estimates for: ",
         paste(absent_ids, collapse = ", "), ".", call. = FALSE)
  }
  data <- data[data$phenotype %in% intersect(phenotype_ids, available), , drop = FALSE]
  data$run <- run_label
  data
}

#' Draw a multi-run forest plot
#'
#' Effect measures are placed in free-scale facets so differences, odds
#' ratios, incidence-rate ratios, and ordinal odds ratios are never put on a
#' shared axis.
#'
#' @param results A `phewas_result` or a list of results.
#' @param phenotype_ids Explicit phenotype values to include, in display order.
#' @param run_labels Optional labels, one per result.
#' @param missing Whether an absent or unsuccessful selection is an error.
#' @param limits Optional reviewed numeric x limits; clipping is an error.
#'
#' @return A `ggplot` object.
#' @export
plot_phewas_forest <- function(results, phenotype_ids, run_labels = NULL,
                               missing = c("error", "drop"), limits = NULL) {
  missing <- match.arg(missing)
  phenotype_ids <- unique(as.character(phenotype_ids))
  if (!length(phenotype_ids) || anyNA(phenotype_ids) || any(!nzchar(phenotype_ids))) {
    stop("phenotype_ids must explicitly name at least one phenotype.", call. = FALSE)
  }
  if (is.data.frame(results) || inherits(results, "phewas_result")) results <- list(results)
  if (!is.list(results) || !length(results)) stop("results must contain at least one result.", call. = FALSE)
  if (is.null(run_labels)) {
    run_labels <- vapply(seq_along(results), function(i) {
      .pf_analysis_label(results[[i]], paste0("run_", i))
    }, character(1))
  }
  if (length(run_labels) != length(results) || anyNA(run_labels) ||
      any(!nzchar(run_labels)) || anyDuplicated(run_labels)) {
    stop("run_labels must provide one unique, nonempty label per result.", call. = FALSE)
  }
  tables <- Map(function(result, run_label) {
    .pf_forest_table(result, phenotype_ids, run_label, missing)
  }, results, as.character(run_labels))
  data <- do.call(rbind, tables)
  rownames(data) <- NULL
  if (!nrow(data)) stop("No selected successful associations remain.", call. = FALSE)

  limits <- .pf_check_limits(c(data$native_conf_low, data$native_conf_high), limits, "x-axis")
  description_map <- vapply(phenotype_ids, function(id) {
    values <- data$description[data$phenotype == id]
    if (length(values)) values[[1L]] else id
  }, character(1))
  names(description_map) <- phenotype_ids
  if (anyDuplicated(description_map)) {
    duplicated_descriptions <- unique(description_map[duplicated(description_map) |
      duplicated(description_map, fromLast = TRUE)])
    replace <- description_map %in% duplicated_descriptions
    description_map[replace] <- paste0(
      description_map[replace], " [", names(description_map)[replace], "]"
    )
  }
  data$description <- factor(
    description_map[data$phenotype], levels = rev(unname(description_map))
  )
  data$run <- factor(data$run, levels = run_labels)
  measure_labels <- c(
    mean_difference = "Difference",
    odds_ratio = "Odds ratio",
    incidence_rate_ratio = "Incidence-rate ratio",
    common_odds_ratio = "Common odds ratio"
  )
  direction_labels <- c(
    phenotypes_as_outcomes = "phenotypes as outcomes",
    phenotypes_as_predictors = "phenotypes as predictors"
  )
  readable_measure <- unname(measure_labels[data$effect_measure])
  readable_measure[is.na(readable_measure)] <- data$effect_measure[is.na(readable_measure)]
  readable_direction <- unname(direction_labels[data$direction])
  readable_direction[is.na(readable_direction)] <- data$direction[is.na(readable_direction)]
  data$effect_panel <- paste(readable_measure, readable_direction, sep = " - ")
  ratio_measures <- c("odds_ratio", "incidence_rate_ratio", "common_odds_ratio",
                      "OR", "IRR", "ordinal_OR")
  neutral <- data.frame(
    effect_panel = unique(data$effect_panel), stringsAsFactors = FALSE
  )
  measure_by_panel <- setNames(data$effect_measure, data$effect_panel)
  neutral$xintercept <- ifelse(
    unname(measure_by_panel[neutral$effect_panel]) %in% ratio_measures, 1, 0
  )
  dodge <- ggplot2::position_dodge(width = 0.55)
  plot <- ggplot2::ggplot(
    data, ggplot2::aes(x = native_effect, y = description, color = run)
  ) +
    ggplot2::geom_vline(
      data = neutral, ggplot2::aes(xintercept = xintercept), inherit.aes = FALSE,
      color = "grey70", linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = native_conf_low, xmax = native_conf_high),
      orientation = "y", width = 0.2, position = dodge
    ) +
    ggplot2::geom_point(position = dodge, size = 2) +
    ggplot2::facet_wrap(~effect_panel, scales = "free_x") +
    ggplot2::labs(x = "Effect estimate (95% CI)", y = NULL, color = "Run") +
    ggplot2::theme_bw()
  if (!is.null(limits)) plot <- plot + ggplot2::coord_cartesian(xlim = limits)
  attr(plot, "phewas_audit") <- list(
    plot = "forest", phenotypes = phenotype_ids, run_labels = run_labels,
    data = data[c("phenotype", "description", "run", "effect_measure", "direction",
                  "native_effect", "native_conf_low", "native_conf_high")]
  )
  plot
}

#' Save a PheWAS plot and optional audit sidecar
#'
#' @param plot A ggplot returned by a phewasFlow plotting function.
#' @param filename One or more destination image or PDF paths.
#' @param width,height,units,dpi Passed to [ggplot2::ggsave()].
#' @param audit_file Optional JSON path for the plot's selections and source
#'   data summary.
#' @param ... Additional arguments passed to [ggplot2::ggsave()].
#'
#' @return The normalized output path or paths, invisibly.
#' @export
save_phewas_plot <- function(plot, filename, width = 10, height = 6,
                             units = "in", dpi = 300, audit_file = NULL, ...) {
  if (!inherits(plot, "ggplot")) stop("plot must be a ggplot object.", call. = FALSE)
  if (!is.character(filename) || !length(filename) || anyNA(filename) ||
      any(!nzchar(filename)) || anyDuplicated(filename)) {
    stop("filename must contain one or more unique, non-empty paths.",
         call. = FALSE)
  }
  parents <- unique(dirname(filename))
  for (parent in parents) {
    if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
      stop("Could not create plot directory: ", parent, call. = FALSE)
    }
  }
  for (path in filename) {
    ggplot2::ggsave(
      filename = path, plot = plot, width = width, height = height,
      units = units, dpi = dpi, ...
    )
  }
  if (!is.null(audit_file)) {
    audit <- attr(plot, "phewas_audit", exact = TRUE)
    if (is.null(audit)) stop("The plot has no phewasFlow audit metadata.", call. = FALSE)
    audit_parent <- dirname(audit_file)
    if (!dir.exists(audit_parent) && !dir.create(audit_parent, recursive = TRUE)) {
      stop("Could not create audit directory: ", audit_parent, call. = FALSE)
    }
    jsonlite::write_json(audit, audit_file, pretty = TRUE, auto_unbox = TRUE,
                         dataframe = "rows", na = "null")
  }
  invisible(normalizePath(filename, winslash = "/", mustWork = TRUE))
}
