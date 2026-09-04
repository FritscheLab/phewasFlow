make_plot_result <- function(analysis_id = "plot_run", direction = "phenotypes_as_outcomes",
                             measure = "mean_difference") {
  out <- data.frame(
    analysis_id = analysis_id,
    phenotype = c("p1", "p2", "p3"),
    status = "ok",
    description = c("One", "Two", "Three"),
    group = c("A", "A", "B"),
    groupnum = c(1, 1, 2),
    color = c("#3366AA", "#3366AA", "#AA6633"),
    direction = direction,
    response = if (direction == "phenotypes_as_outcomes") c("p1", "p2", "p3") else "endpoint",
    predictor = if (direction == "phenotypes_as_outcomes") "pgs" else c("p1", "p2", "p3"),
    outcome_type = "continuous",
    predictor_type = "numeric",
    transformation = "none",
    reference = NA_character_,
    outcome_reference = NA_character_,
    predictor_reference = NA_character_,
    effect_measure = measure,
    estimate = c(0.2, -0.1, 0.4),
    std_error = c(0.1, 0.1, 0.2),
    p_value = c(0, 0.02, 0.4),
    q_value = c(0, 0.03, 0.4),
    neg_log10_p = c(400, -log10(0.02), -log10(0.4)),
    native_effect = c(0.2, -0.1, 0.4),
    native_conf_low = c(0.0, -0.3, 0.0),
    native_conf_high = c(0.4, 0.1, 0.8),
    stringsAsFactors = FALSE
  )
  class(out) <- c("phewas_result", "data.frame")
  attr(out, "spec") <- list(analysis_id = analysis_id, fdr_threshold = 0.05)
  out
}

make_crowded_plot_result <- function() {
  groups <- c(
    "Infections", "Endocrine/Metab", "Neoplasms", "Blood/Immune",
    "Mental", "Neurological", "Sense organs", "Cardiovascular",
    "Respiratory", "Gastrointestinal", "Genitourinary", "Musculoskeletal",
    "Dermatological", "Neonatal", "Congenital", "Symptoms", "Pregnancy",
    "Genetic"
  )
  counts <- c(32L, 28L, 24L, 22L, 30L, 26L, 18L, 20L, 18L, 24L,
              20L, 34L, 12L, 3L, 4L, 8L, 3L, 2L)
  group_number <- rep(seq_along(groups), counts)
  n <- length(group_number)
  out <- make_plot_result()[rep(1L, n), , drop = FALSE]
  rownames(out) <- NULL
  out$phenotype <- sprintf("P%04d", seq_len(n))
  out$description <- paste("Phenotype", seq_len(n))
  out$group <- groups[group_number]
  out$groupnum <- group_number
  colors <- grDevices::hcl.colors(length(groups), palette = "Dark 3")
  out$color <- colors[group_number]
  out$response <- out$phenotype
  out$estimate <- rep(c(0.2, -0.2), length.out = n)
  out$p_value <- seq(0.001, 0.2, length.out = n)
  out$q_value <- pmin(1, out$p_value * 2)
  out$neg_log10_p <- -log10(out$p_value)
  out$native_effect <- out$estimate
  out$native_conf_low <- out$estimate - 0.1
  out$native_conf_high <- out$estimate + 0.1
  class(out) <- c("phewas_result", "data.frame")
  attr(out, "spec") <- list(analysis_id = "crowded", fdr_threshold = 0.05)
  out
}

test_that("Manhattan plot uses authoritative underflow-safe values", {
  result <- make_plot_result()
  default_plot <- plot_phewas_manhattan(result)
  expect_identical(
    attr(default_plot, "phewas_audit")$thresholds$method, "Bonferroni"
  )
  plot <- plot_phewas_manhattan(
    result, significance = "none", label = "p1", highlight = "p2"
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(max(plot$data$neg_log10_p), 400)
  built <- ggplot2::ggplot_build(plot)
  expect_lte(built$layout$panel_params[[1L]]$y.range[[1L]], 0)
  expect_equal(attr(plot, "phewas_audit")$labels, "p1")
  expect_identical(attr(plot, "phewas_audit")$phenotypes, c("p1", "p2", "p3"))
  expect_identical(attr(plot, "phewas_audit")$group_display_requested, "legend")
  expect_identical(attr(plot, "phewas_audit")$group_display, "legend")
  expect_identical(
    as.character(plot$data$effect_direction),
    c("Positive effect", "Negative effect", "Positive effect")
  )
  expect_false(any(c(
    "column", "phenotype_id", "phecode", "phenotype_column", "label",
    "category", "category_order"
  ) %in% names(plot$data)))
  expect_error(plot_phewas_manhattan(result, label = "not_present"),
               "unavailable")
  expect_error(plot_phewas_manhattan(result, significance = "none",
                                     y_limit = c(0, 20)), "clip")
})

test_that("crowded group names use an ordered legend without dropping groups", {
  result <- make_crowded_plot_result()
  groups <- unique(result$group)

  automatic <- plot_phewas_manhattan(
    result, significance = "none", group_display = "auto"
  )
  audit <- attr(automatic, "phewas_audit")
  expect_identical(audit$group_display_requested, "auto")
  expect_identical(audit$group_display, "legend")
  expect_identical(audit$group_order, groups)
  expect_length(automatic$scales$get_scales("x")$breaks, 5L)
  expect_identical(
    automatic$labels$x, "Phenotypes ordered by group"
  )
  expect_identical(automatic$scales$get_scales("fill")$breaks, groups)
  expect_s3_class(ggplot2::ggplotGrob(automatic), "gtable")

  forced_axis <- plot_phewas_manhattan(
    result, significance = "none", group_display = "x_axis"
  )
  expect_identical(
    attr(forced_axis, "phewas_audit")$group_display_requested, "x_axis"
  )
  expect_identical(attr(forced_axis, "phewas_audit")$group_display, "x_axis")
  expect_length(forced_axis$scales$get_scales("x")$breaks, length(groups))

  forced_legend <- plot_phewas_manhattan(
    make_plot_result(), significance = "none", group_display = "legend"
  )
  expect_identical(
    attr(forced_legend, "phewas_audit")$group_display_requested, "legend"
  )
  expect_identical(attr(forced_legend, "phewas_audit")$group_display, "legend")
  expect_identical(forced_legend$scales$get_scales("x")$breaks, 1:3)
  expect_identical(forced_legend$labels$x, "Phenotypes ordered by group")
})

test_that("Manhattan and volcano views share one visual grammar", {
  result <- make_plot_result()
  common <- list(
    significance = "bonferroni",
    label = "p1",
    highlight = "p2",
    group_legend_columns = 2L,
    base_size = 11,
    label_size = 9
  )
  manhattan <- do.call(plot_phewas_manhattan, c(list(x = result), common))
  volcano <- do.call(plot_phewas_volcano, c(list(x = result), common))

  manhattan_audit <- attr(manhattan, "phewas_audit")
  volcano_audit <- attr(volcano, "phewas_audit")
  expect_identical(manhattan_audit$group_display, "legend")
  expect_identical(manhattan_audit$group_order, volcano_audit$group_order)

  manhattan_fill <- manhattan$scales$get_scales("fill")
  volcano_fill <- volcano$scales$get_scales("fill")
  expect_identical(manhattan_fill$name, "Phenotype group")
  expect_identical(manhattan_fill$name, volcano_fill$name)
  expect_identical(manhattan_fill$breaks, volcano_fill$breaks)
  expect_identical(
    manhattan_fill$get_labels(manhattan_fill$breaks),
    volcano_fill$get_labels(volcano_fill$breaks)
  )
  expect_identical(
    manhattan_fill$palette(length(manhattan_fill$breaks)),
    volcano_fill$palette(length(volcano_fill$breaks))
  )
  expect_false(manhattan_fill$drop)
  expect_false(volcano_fill$drop)
  expect_identical(manhattan_fill$guide$params, volcano_fill$guide$params)
  expect_identical(manhattan_fill$guide$params$override.aes$shape, 21)

  manhattan_shape <- manhattan$scales$get_scales("shape")
  volcano_shape <- volcano$scales$get_scales("shape")
  expect_identical(manhattan_shape$breaks, volcano_shape$breaks)
  expect_identical(
    manhattan_shape$palette(length(manhattan_shape$breaks)),
    volcano_shape$palette(length(volcano_shape$breaks))
  )

  primary_point <- function(plot) {
    Filter(function(layer) {
      inherits(layer$geom, "GeomPoint") && inherits(layer$data, "waiver")
    }, plot$layers)[[1L]]
  }
  manhattan_point <- primary_point(manhattan)
  volcano_point <- primary_point(volcano)
  expect_identical(manhattan_point$aes_params, volcano_point$aes_params)
  expect_identical(names(manhattan_point$mapping), names(volcano_point$mapping))

  label_layer <- function(plot) {
    Filter(function(layer) inherits(layer$geom, "GeomLabelRepel"),
           plot$layers)[[1L]]
  }
  expect_identical(
    label_layer(manhattan)$aes_params,
    label_layer(volcano)$aes_params
  )
  expect_identical(
    label_layer(manhattan)$geom_params,
    label_layer(volcano)$geom_params
  )

  threshold_layer <- function(plot, geom) {
    Filter(function(layer) {
      inherits(layer$geom, geom) && is.data.frame(layer$data) &&
        "threshold_label" %in% names(layer$data)
    }, plot$layers)[[1L]]
  }
  manhattan_threshold <- threshold_layer(manhattan, "GeomLabel")
  volcano_threshold <- threshold_layer(volcano, "GeomLabel")
  expect_identical(
    manhattan_threshold$aes_params,
    volcano_threshold$aes_params
  )
  expect_identical(
    manhattan_threshold$geom_params,
    volcano_threshold$geom_params
  )
  threshold_line <- function(plot) {
    Filter(function(layer) {
      inherits(layer$geom, "GeomHline") && is.data.frame(layer$data) &&
        "method" %in% names(layer$data)
    }, plot$layers)[[1L]]
  }
  expect_identical(
    threshold_line(manhattan)$aes_params,
    threshold_line(volcano)$aes_params
  )
  expect_identical(
    threshold_line(manhattan)$geom_params,
    threshold_line(volcano)$geom_params
  )

  shared_theme <- c(
    "panel.grid.major.y", "panel.grid.minor", "axis.line", "axis.ticks",
    "axis.ticks.length", "axis.text.y", "axis.title", "axis.title.y",
    "legend.position", "legend.box", "legend.box.just", "legend.direction",
    "legend.text", "legend.title", "legend.title.position",
    "legend.key.width", "legend.spacing.x", "plot.title", "plot.margin"
  )
  for (field in shared_theme) {
    expect_identical(manhattan$theme[[field]], volcano$theme[[field]],
                     info = field)
  }
  expect_identical(manhattan$labels$y, volcano$labels$y)
  expect_equal(
    manhattan$scales$get_scales("y")$breaks,
    volcano$scales$get_scales("y")$breaks
  )
  expect_identical(
    manhattan$scales$get_scales("y")$labels,
    volcano$scales$get_scales("y")$labels
  )
  expect_equal(
    ggplot2::ggplot_build(manhattan)$layout$panel_params[[1L]]$y.range,
    ggplot2::ggplot_build(volcano)$layout$panel_params[[1L]]$y.range
  )
  expect_identical(manhattan$labels$x, "Phenotypes ordered by group")
})

test_that("shared y axes include the highest regular tick", {
  for (case in list(c(10.7, 12), c(11.9, 12), c(12.5, 14))) {
    result <- make_plot_result()
    result$neg_log10_p <- c(case[[1L]], 7.4, 0.4)
    result$p_value <- 10^(-result$neg_log10_p)
    result$q_value <- pmin(1, result$p_value * nrow(result))

    manhattan <- plot_phewas_manhattan(result, significance = "none")
    volcano <- plot_phewas_volcano(result, significance = "none")
    expect_true(case[[2L]] %in% manhattan$scales$get_scales("y")$breaks)
    expect_identical(
      manhattan$scales$get_scales("y")$breaks,
      volcano$scales$get_scales("y")$breaks
    )
  }

  result$estimate[[1L]] <- NA_real_
  manhattan <- plot_phewas_manhattan(result, significance = "none")
  volcano <- plot_phewas_volcano(result, significance = "none")
  expect_identical(
    manhattan$scales$get_scales("y")$breaks,
    volcano$scales$get_scales("y")$breaks
  )
})

test_that("canonical phenotype metadata drives descriptions and natural order", {
  result <- make_plot_result()
  result$phenotype <- c("DX_10", "DX_2", "DX_1")
  result$description <- c("Diagnosis ten", "Diagnosis two", "Diagnosis one")
  result$group <- c("First", "First", "Second")
  result$groupnum <- c(1, 1, 2)
  result$color <- c("#3366AA", "#3366AA", "#AA6633")

  plot <- plot_phewas_manhattan(
    result, significance = "none", label = "DX_10"
  )

  expect_identical(plot$data$phenotype, c("DX_2", "DX_10", "DX_1"))
  selected <- plot$data$plot_label[!is.na(plot$data$plot_label)]
  expect_identical(unname(selected), "Diagnosis ten")
})

test_that("Manhattan plot uses group blocks, labeled thresholds, and boxed labels", {
  result <- make_plot_result()
  result$testing_family_size <- 100L
  plot <- plot_phewas_manhattan(
    result,
    significance = "bonferroni",
    label = "p1",
    highlight = "p2",
    label_mode = "both",
    group_display = "x_axis"
  )

  audit <- attr(plot, "phewas_audit")
  expect_identical(audit$group_order, c("A", "B"))
  expect_identical(audit$highlights, "p2")
  expect_identical(audit$thresholds$method, "Bonferroni")
  expect_identical(audit$testing_family_size, 100L)
  expect_equal(audit$thresholds$y, -log10(0.05 / 100))
  expect_equal(plot$scales$get_scales("x")$breaks, c(1.5, 3))
  expect_false(any(grepl(
    "Bonferroni", plot$scales$get_scales("y")$labels, fixed = TRUE
  )))
  threshold_layers <- Filter(function(layer) {
    is.data.frame(layer$data) && "threshold_label" %in% names(layer$data) &&
      any(grepl("Bonferroni", layer$data$threshold_label, fixed = TRUE))
  }, plot$layers)
  expect_length(threshold_layers, 1L)
  expect_match(
    threshold_layers[[1L]]$data$threshold_label, "Bonferroni: [0-9]"
  )
  expect_identical(
    plot$data$plot_label[plot$data$phenotype == "p1"],
    "One\n(p1)"
  )
  expect_gte(length(ggplot2::ggplot_build(plot)$data), 6L)
})

test_that("automatic labels and optional upper-tail compression are deterministic", {
  result <- make_plot_result()
  plot <- plot_phewas_manhattan(
    result,
    significance = "none",
    top_label_n = 2,
    axis_break_mode = "on",
    axis_break_top = 20
  )

  audit <- attr(plot, "phewas_audit")
  expect_identical(audit$labels, c("p1", "p2"))
  expect_true(audit$axis_break$enabled)
  expect_lt(max(plot$data$y_plot), max(plot$data$neg_log10_p))
  expect_error(plot_phewas_manhattan(result, top_label_n = 51), "between 0 and 50")
  expect_error(
    plot_phewas_manhattan(result, axis_break_top_size = 1),
    "strictly between"
  )

  result$q_value <- c(0.01, 0.2, 0.3)
  filled <- plot_phewas_manhattan(
    result, significance = "none", top_label_n = 3
  )
  expect_identical(
    attr(filled, "phewas_audit")$labels,
    c("p1", "p2", "p3")
  )
})

test_that("plotting accepts data.table phewas results without mutation", {
  result <- data.table::as.data.table(make_plot_result())
  class(result) <- c("phewas_result", "data.table", "data.frame")
  attr(result, "spec") <- list(analysis_id = "dt", fdr_threshold = 0.05)
  before <- data.table::copy(result)

  plot <- plot_phewas_manhattan(result, significance = "none")

  expect_s3_class(plot, "ggplot")
  expect_identical(result, before)
})

test_that("comparison plots enforce native compatibility", {
  x <- make_plot_result("x")
  y <- make_plot_result("y")
  comparison <- compare_phewas(x, y)

  expect_s3_class(plot_phewas_concordance(comparison), "ggplot")
  expect_s3_class(plot_phewas_effects(comparison), "ggplot")
  overlap <- plot_phewas_category_overlap(comparison)
  expect_s3_class(overlap, "ggplot")
  expect_true("group" %in% names(overlap$data))
  expect_false("category" %in% names(overlap$data))
  expect_identical(attr(overlap, "phewas_audit")$plot, "group_overlap")

  reverse <- make_plot_result("reverse", direction = "phenotypes_as_predictors")
  incompatible <- compare_phewas(x, reverse)
  expect_error(plot_phewas_effects(incompatible), "incompatible")
})

test_that("forest plots separate effect and direction semantics", {
  forward <- make_plot_result("forward")
  reverse <- make_plot_result(
    "reverse", direction = "phenotypes_as_predictors", measure = "odds_ratio"
  )
  reverse$native_effect <- exp(reverse$estimate)
  reverse$native_conf_low <- exp(reverse$estimate - 1.96 * reverse$std_error)
  reverse$native_conf_high <- exp(reverse$estimate + 1.96 * reverse$std_error)

  plot <- plot_phewas_forest(
    list(forward, reverse), c("p1", "p2"),
    run_labels = c("Forward", "Reverse")
  )

  expect_s3_class(plot, "ggplot")
  expect_equal(length(unique(plot$data$effect_panel)), 2)
  expect_true(all(c("phenotype", "description") %in% names(plot$data)))
  expect_false(any(c(
    "column", "phenotype_id", "phecode", "phenotype_column", "label",
    "category", "category_order"
  ) %in% names(plot$data)))
  expect_error(
    plot_phewas_forest(forward, c("p1", "absent")), "lacks successful"
  )
  expect_error(
    plot_phewas_forest(forward, "p1", limits = c(-0.05, 0.05)), "clip"
  )
})

test_that("plot writer can create an image and JSON audit", {
  plot <- plot_phewas_manhattan(make_plot_result(), significance = "none")
  image <- tempfile(fileext = ".png")
  pdf <- tempfile(fileext = ".pdf")
  audit <- tempfile(fileext = ".json")
  on.exit(unlink(c(image, pdf, audit)), add = TRUE)

  paths <- NULL
  expect_silent(paths <- save_phewas_plot(
    plot, c(image, pdf), width = 4, height = 3, audit_file = audit
  ))
  expect_length(paths, 2L)
  expect_true(file.exists(image))
  expect_true(file.exists(pdf))
  expect_true(file.exists(audit))
  expect_error(save_phewas_plot(plot, character()), "one or more")
})
