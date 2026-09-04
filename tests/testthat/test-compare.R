make_comparison_result <- function(
    analysis_id = "run_a", ids = c("p1", "p2", "p3"), shift = 0,
    direction = "phenotypes_as_outcomes") {
  n <- length(ids)
  index <- match(ids, c("p1", "p2", "p3"))
  if (anyNA(index)) index <- seq_len(n)
  out <- data.frame(
    analysis_id = analysis_id,
    phenotype = ids,
    status = rep("ok", n),
    description = paste("Phenotype", ids),
    group = c("A", "B", "A")[index],
    groupnum = c(1, 2, 1)[index],
    color = c("#3366AA", "#AA6633", "#3366AA")[index],
    direction = direction,
    response = if (direction == "phenotypes_as_outcomes") ids else "endpoint",
    predictor = if (direction == "phenotypes_as_outcomes") "pgs" else ids,
    outcome_type = "continuous",
    predictor_type = "numeric",
    transformation = if (direction == "phenotypes_as_outcomes") "zscore" else "none",
    reference = NA_character_,
    outcome_reference = NA_character_,
    predictor_reference = NA_character_,
    effect_measure = "mean_difference",
    estimate = index / 10 + shift,
    std_error = rep(0.1, n),
    p_value = c(0.001, 0.02, 0.3)[index],
    q_value = c(0.003, 0.03, 0.3)[index],
    neg_log10_p = -log10(c(0.001, 0.02, 0.3)[index]),
    native_effect = index / 10 + shift,
    native_conf_low = index / 10 + shift - 0.196,
    native_conf_high = index / 10 + shift + 0.196,
    stringsAsFactors = FALSE
  )
  class(out) <- c("phewas_result", "data.frame")
  attr(out, "spec") <- list(analysis_id = analysis_id, fdr_threshold = 0.05)
  out
}

test_that("compare_phewas reports matched concordance and signal overlap", {
  x <- make_comparison_result("x")
  y <- make_comparison_result("y", shift = 0.02)

  comparison <- compare_phewas(x, y)

  expect_s3_class(comparison, "phewas_comparison")
  expect_s3_class(comparison$matched, "data.table")
  expect_s3_class(comparison$metrics, "data.table")
  expect_s3_class(comparison$group_overlap, "data.table")
  expect_equal(comparison$metrics$n_mutually_successful, 3)
  expect_equal(comparison$metrics$fdr_intersection, 2)
  expect_equal(comparison$metrics$fdr_jaccard, 1)
  expect_equal(comparison$metrics$sign_concordance, 1)
  expect_true(comparison$metrics$native_compatible)
  expect_equal(comparison$metrics$native_slope, 1, tolerance = 1e-12)
  expect_equal(sum(comparison$group_overlap$both), 2)
  expect_true("group" %in% names(comparison$group_overlap))
  expect_identical(comparison$matched$color_x, x$color)
  expect_identical(comparison$matched$color_y, y$color)
  legacy <- c(
    "column", "phenotype_id", "phecode", "phenotype_column", "label",
    "category", "category_order"
  )
  expect_false(any(legacy %in% names(comparison$matched)))
})

test_that("unambiguous legacy results normalize without input mutation", {
  x <- make_comparison_result("x")
  names(x)[names(x) == "phenotype"] <- "phenotype_id"
  x$phenotype_column <- x$phenotype_id
  names(x)[names(x) == "description"] <- "label"
  names(x)[names(x) == "group"] <- "category"
  names(x)[names(x) == "groupnum"] <- "category_order"
  x$phecode <- x$phenotype_id
  x <- data.table::as.data.table(x)
  class(x) <- c("phewas_result", "data.table", "data.frame")
  attr(x, "spec") <- list(analysis_id = "x", fdr_threshold = 0.05)
  before <- data.table::copy(x)

  comparison <- compare_phewas(x, x)

  expect_identical(x, before)
  expect_identical(comparison$matched$phenotype, c("p1", "p2", "p3"))
  expect_identical(
    comparison$matched$description_x,
    paste("Phenotype", c("p1", "p2", "p3"))
  )
  expect_identical(comparison$matched$groupnum_x, c(1, 2, 1))
  expect_false(any(c(
    "column", "phenotype_id", "phecode", "phenotype_column", "label",
    "category", "category_order"
  ) %in% names(comparison$matched)))
})

test_that("partial overlap uses mutual successes only for correlations", {
  x <- make_comparison_result("x")
  y <- make_comparison_result("y", ids = c("p2", "p3"))
  y$status[[2L]] <- "skipped"
  y$q_value[[2L]] <- NA_real_

  comparison <- compare_phewas(x, y, fdr_threshold = 0.05)

  expect_equal(comparison$metrics$n_mutually_successful, 1)
  expect_equal(comparison$matched$phenotype, "p2")
  expect_true(is.na(comparison$metrics$standardized_estimate_correlation))
})

test_that("direction changes disable native metrics but retain standardized comparison", {
  x <- make_comparison_result("forward")
  y <- make_comparison_result(
    "reverse", direction = "phenotypes_as_predictors"
  )

  comparison <- compare_phewas(x, y)

  expect_false(comparison$native_compatible)
  expect_true("direction" %in% comparison$native_incompatibility)
  expect_equal(comparison$metrics$n_mutually_successful, 3)
  expect_true(is.finite(comparison$metrics$standardized_estimate_correlation))
})

test_that("duplicate IDs and inverted phenotype references are rejected", {
  duplicate <- make_comparison_result()
  duplicate$phenotype[[3L]] <- "p1"
  expect_error(compare_phewas(duplicate, make_comparison_result("y")),
               "duplicated phenotype")

  ambiguous <- make_comparison_result()
  ambiguous$phenotype_id <- ambiguous$phenotype
  ambiguous$phenotype_id[[1L]] <- "other"
  expect_error(
    compare_phewas(ambiguous, make_comparison_result("y")),
    "cannot be safely normalized.*phenotype"
  )

  x <- make_comparison_result("x")
  y <- make_comparison_result("y")
  x$outcome_reference[[1L]] <- "0"
  y$outcome_reference[[1L]] <- "1"
  expect_error(compare_phewas(x, y), "reference levels")
})

test_that("mixed native effect scales do not produce aggregate native metrics", {
  x <- make_comparison_result("x")
  y <- make_comparison_result("y")
  x$effect_measure[[3L]] <- y$effect_measure[[3L]] <- "odds_ratio"

  comparison <- compare_phewas(x, y)

  expect_false(comparison$native_compatible)
  expect_true("mixed effect measures" %in% comparison$native_incompatibility)
  expect_true(is.na(comparison$metrics$native_correlation))
})

test_that("different fixed outcomes disable native effect comparison", {
  x <- make_comparison_result("x", direction = "phenotypes_as_predictors")
  y <- make_comparison_result("y", direction = "phenotypes_as_predictors")
  y$response <- "other_endpoint"

  comparison <- compare_phewas(x, y)

  expect_false(comparison$native_compatible)
  expect_true("response" %in% comparison$native_incompatibility)
  expect_true(is.na(comparison$metrics$native_correlation))
  expect_true(is.finite(comparison$metrics$standardized_estimate_correlation))
})
