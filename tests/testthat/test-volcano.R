make_volcano_result <- function() {
  result <- data.frame(
    analysis_id = "volcano_test",
    phenotype = paste0("p", seq_len(8L)),
    status = "ok",
    description = c(
      "Strong positive", "Strong negative", "Moderate positive",
      "Moderate negative", "Small positive", "Small negative",
      "Near zero one", "Near zero two"
    ),
    group = rep(c("Cardiovascular", "Musculoskeletal"), each = 4L),
    groupnum = rep(c(1, 2), each = 4L),
    color = rep(c("#2B8CBE", "#7B2CBF"), each = 4L),
    estimate = c(1.3, -1.1, 0.7, -0.6, 0.3, -0.25, 0.08, -0.05),
    effect_measure = "mean_difference",
    p_value = c(0, 1e-12, 1e-8, 1e-6, 0.002, 0.02, 0.3, 0.8),
    q_value = c(0, 8e-12, 8e-8, 8e-6, 0.0032, 0.0267, 0.3429, 0.8),
    neg_log10_p = c(400, 12, 8, 6, -log10(0.002), -log10(0.02),
                    -log10(0.3), -log10(0.8)),
    testing_family_size = 12L,
    stringsAsFactors = FALSE
  )
  class(result) <- c("phewas_result", "data.frame")
  attr(result, "spec") <- list(
    analysis_id = "volcano_test", fdr_threshold = 0.05
  )
  result
}

test_that("volcano defaults to an in-panel Bonferroni threshold", {
  result <- make_volcano_result()
  plot <- plot_phewas_volcano(result)
  audit <- attr(plot, "phewas_audit")

  expect_s3_class(plot, "ggplot")
  expect_identical(audit$plot, "volcano")
  expect_identical(audit$significance, "bonferroni")
  expect_identical(audit$testing_family_size, 12L)
  expect_identical(audit$labels, character())
  expect_equal(audit$threshold$y, -log10(0.05 / 12))
  expect_identical(audit$threshold$method, "Bonferroni")
  expect_equal(plot$data$x_plot, plot$data$estimate)
  expect_equal(plot$data$y_plot, plot$data$neg_log10_p)
  expect_equal(max(plot$data$y_plot), 400)
  expect_identical(plot$labels$x, "Mean difference")
  expect_identical(audit$group_order,
                   c("Cardiovascular", "Musculoskeletal"))
  expect_identical(
    plot$scales$get_scales("fill")$breaks,
    c("Cardiovascular", "Musculoskeletal")
  )

  threshold_labels <- Filter(function(layer) {
    inherits(layer$geom, "GeomLabel") && is.data.frame(layer$data) &&
      "threshold_label" %in% names(layer$data)
  }, plot$layers)
  expect_length(threshold_labels, 1L)
  expect_match(threshold_labels[[1L]]$data$threshold_label, "Bonferroni")
  y_scale <- plot$scales$get_scales("y")
  expect_false(is.character(y_scale$labels) &&
                 any(grepl("Bonferroni", y_scale$labels, fixed = TRUE)))
})

test_that("mixed model scales use clearly named zero-centered facets", {
  result <- make_volcano_result()
  result$effect_measure[5:8] <- "odds_ratio"
  result$estimate[5:8] <- c(-0.12, 0.08, -0.03, 0.02)

  plot <- plot_phewas_volcano(
    result, significance = "none", top_label_n = 0
  )
  audit <- attr(plot, "phewas_audit")

  expect_s3_class(plot$facet, "FacetWrap")
  expect_true(plot$facet$params$free$x)
  expect_identical(
    audit$effect_measures,
    c("mean_difference", "odds_ratio")
  )
  expect_identical(
    audit$effect_panels,
    c("Mean difference", "Log odds ratio")
  )
  expect_identical(plot$labels$x, "Estimate on fitted model scale")

  built <- ggplot2::ggplot_build(plot)
  ranges <- lapply(built$layout$panel_params, function(panel) panel$x.range)
  expect_true(all(vapply(ranges, function(range) {
    isTRUE(all.equal(abs(range[[1L]]), abs(range[[2L]]), tolerance = 1e-8))
  }, logical(1))))
})

test_that("volcano labels and highlights are deterministic and explicit", {
  result <- make_volcano_result()
  plot <- plot_phewas_volcano(
    result,
    significance = "none",
    label = c("p2", "p1"),
    highlight = "p3",
    label_mode = "both"
  )
  audit <- attr(plot, "phewas_audit")

  expect_identical(audit$labels, c("p2", "p1"))
  expect_identical(audit$highlights, "p3")
  expect_identical(
    plot$data$plot_label[plot$data$phenotype == "p1"],
    "Strong positive\n(p1)"
  )
  label_layers <- Filter(function(layer) {
    inherits(layer$geom, "GeomLabelRepel")
  }, plot$layers)
  expect_length(label_layers, 1L)
  expect_identical(label_layers[[1L]]$geom_params$seed, 20260819)
  expect_error(plot_phewas_volcano(result, label = "absent"), "unavailable")
  expect_s3_class(
    plot_phewas_volcano(result, label = "absent", missing = "drop"),
    "ggplot"
  )
})

test_that("volcano supports BH and validates scientific scales", {
  result <- make_volcano_result()
  bh <- plot_phewas_volcano(
    result, significance = "bh", top_label_n = 0
  )
  expect_identical(attr(bh, "phewas_audit")$threshold$method, "BH FDR")
  expect_equal(attr(bh, "phewas_audit")$threshold$y, -log10(0.02))

  result$effect_measure[[1L]] <- "unsupported_ratio"
  expect_error(plot_phewas_volcano(result), "supported effect_measure")

  result <- make_volcano_result()
  expect_error(
    plot_phewas_volcano(result, x_limit = c(-1, 1)), "clip"
  )
  expect_error(
    plot_phewas_volcano(result, x_limit = c(-2, 3)), "symmetric"
  )
  expect_error(
    plot_phewas_volcano(result, y_limit = c(0, 20)), "clip"
  )
  expect_error(
    plot_phewas_volcano(result, top_label_n = 51), "between 0 and 50"
  )
})

test_that("Bonferroni keeps the recorded family after plot filtering", {
  result <- make_volcano_result()
  result$estimate[[8L]] <- NA_real_
  plot <- plot_phewas_volcano(result, top_label_n = 0)
  audit <- attr(plot, "phewas_audit")

  expect_length(audit$phenotypes, 7L)
  expect_identical(audit$testing_family_size, 12L)
  expect_equal(audit$threshold$y, -log10(0.05 / 12))

  result$testing_family_size[[8L]] <- 11L
  expect_error(plot_phewas_volcano(result), "conflicting")
})

test_that("volcano preserves the full plotted group legend after estimate filtering", {
  result <- make_volcano_result()
  result$estimate[result$group == "Musculoskeletal"] <- NA_real_

  plot <- plot_phewas_volcano(result, significance = "none")
  audit <- attr(plot, "phewas_audit")
  fill <- plot$scales$get_scales("fill")

  expect_identical(
    audit$group_order,
    c("Cardiovascular", "Musculoskeletal")
  )
  expect_identical(fill$breaks, audit$group_order)
  expect_false(fill$drop)
  expect_identical(levels(plot$data$group), audit$group_order)
  expect_error(
    plot_phewas_volcano(result, group_wrap_width = 7),
    "Wrap widths"
  )
})
