make_contrast_fixture <- function() {
  set.seed(771)
  n <- 180
  pgs <- stats::rnorm(n, mean = 1, sd = 2)
  age <- stats::rnorm(n, 55, 9)
  added <- stats::rnorm(n)
  y1 <- 0.55 * pgs + 0.9 * added + 0.01 * age + stats::rnorm(n)
  y2 <- -0.30 * pgs + 0.02 * age + stats::rnorm(n)
  y1[c(2, 8, 11)] <- NA_real_
  y2[c(3, 7)] <- NA_real_
  added[c(1, 4, 9)] <- NA_real_
  data <- data.frame(pgs = pgs, y1 = y1, y2 = y2, age = age, added = added)
  metadata <- data.frame(
    phenotype = c("y1", "y2"),
    description = c("Outcome one", "Outcome two"),
    group = c("A", "B"),
    groupnum = c(1, 2),
    color = c("#3366AA", "#AA6633"),
    variable_type = c("numeric", "numeric"),
    outcome_type = c("continuous", "continuous"),
    reference = NA_character_,
    levels = NA_character_,
    scores = NA_character_,
    offset = NA_character_,
    stringsAsFactors = FALSE
  )
  eligibility <- list(
    min_n = 20L, min_cases = 2L, min_controls = 2L,
    min_outcome_levels = 3L
  )
  base <- phewas_spec(
    "base", anchor = "pgs", direction = "phenotypes_as_outcomes",
    phenotypes = metadata, covariates = "age", anchor_type = "numeric",
    anchor_transform = "zscore", eligibility = eligibility,
    fdr_threshold = 0.05
  )
  adjusted <- phewas_spec(
    "adjusted", anchor = "pgs", direction = "phenotypes_as_outcomes",
    phenotypes = metadata, covariates = c("age", "added"),
    anchor_type = "numeric", anchor_transform = "zscore",
    eligibility = eligibility, fdr_threshold = 0.05
  )
  list(data = data, base = base, adjusted = adjusted)
}

test_that("nested contrasts use identical rows and shared corrections", {
  fixture <- make_contrast_fixture()
  contrast <- run_phewas_contrast(
    fixture$data, fixture$base, fixture$adjusted
  )

  expect_s3_class(contrast, "phewas_contrast")
  expect_s3_class(contrast$associations, "data.table")
  expect_s3_class(contrast$diagnostics, "data.table")
  expect_equal(nrow(contrast$associations), 2)
  expect_true(all(contrast$associations$status == "ok"))
  expect_true(all(contrast$associations$n_complete > 100))
  expect_true(all(nzchar(contrast$associations$sample_fingerprint)))
  base_rows <- stats::complete.cases(fixture$data[c("pgs", "age", "added")])
  expect_equal(
    unique(contrast$associations$transformation_mean),
    mean(fixture$data$pgs[base_rows])
  )
  expect_equal(
    unique(contrast$associations$transformation_source_n), sum(base_rows)
  )
  expect_equal(
    contrast$associations$base_q_value,
    stats::p.adjust(contrast$associations$base_p_value, method = "BH"),
    tolerance = 1e-12
  )
  expect_equal(
    contrast$associations$adjusted_q_value,
    stats::p.adjust(contrast$associations$adjusted_p_value, method = "BH"),
    tolerance = 1e-12
  )
  expect_true(contrast$associations$partial_r_squared[
    contrast$associations$phenotype == "y1"
  ] > 0)
  expect_false(any(c(
    "column", "phenotype_id", "phecode", "phenotype_column", "label",
    "category", "category_order"
  ) %in% names(contrast$associations)))
  expect_true("phenotype" %in% names(contrast$diagnostics))
})

test_that("contrast subset order is specification order", {
  fixture <- make_contrast_fixture()
  contrast <- run_phewas_contrast(
    fixture$data, fixture$base, fixture$adjusted,
    phenotype_ids = c("y2", "y1")
  )
  expect_equal(contrast$associations$phenotype, c("y1", "y2"))
  expect_error(
    run_phewas_contrast(fixture$data, fixture$base, fixture$adjusted, "unknown"),
    "Unknown phenotype"
  )
})

test_that("contrast rejects non-nested or scientifically different specs", {
  fixture <- make_contrast_fixture()
  same_covariates <- fixture$adjusted
  same_covariates$covariates <- fixture$base$covariates
  expect_error(
    run_phewas_contrast(fixture$data, fixture$base, same_covariates),
    "must add at least one"
  )

  different_threshold <- fixture$adjusted
  different_threshold$fdr_threshold <- 0.1
  expect_error(
    run_phewas_contrast(fixture$data, fixture$base, different_threshold),
    "differ outside"
  )
})

test_that("contrast preserves canonical annotations and exclusion provenance", {
  fixture <- make_contrast_fixture()
  fixture$data$cohort <- rep("keep", nrow(fixture$data))
  fixture$data$cohort[c(1, 5, 10)] <- "drop"
  fixture$base$exclusion_column <- "cohort"
  fixture$base$exclusion_values <- "drop"
  fixture$adjusted$exclusion_column <- "cohort"
  fixture$adjusted$exclusion_values <- "drop"
  for (name in c("base", "adjusted")) {
    fixture[[name]]$phenotypes$description <- c("Long Y1", "Long Y2")
    fixture[[name]]$phenotypes$group <- c("Group one", "Group two")
    fixture[[name]]$phenotypes$groupnum <- c(11, 22)
    fixture[[name]]$phenotypes$color <- c("#112233", "#445566")
  }

  contrast <- run_phewas_contrast(
    fixture$data, fixture$base, fixture$adjusted
  )

  expect_equal(nrow(contrast$associations), 2L)
  expect_identical(contrast$associations$phenotype, c("y1", "y2"))
  expect_identical(contrast$associations$description, c("Long Y1", "Long Y2"))
  expect_identical(contrast$associations$group, c("Group one", "Group two"))
  expect_equal(contrast$associations$groupnum, c(11, 22))
  expect_identical(contrast$associations$color, c("#112233", "#445566"))
  expect_true(all(contrast$associations$exclusion_active))
  expect_true(all(contrast$associations$exclusion_column == "cohort"))
  expect_true(all(contrast$associations$n_before_exclusion == nrow(fixture$data)))
  expect_true(all(contrast$associations$n_excluded == 3L))
  expect_true(all(contrast$associations$n_after_exclusion == nrow(fixture$data) - 3L))
})

test_that("unexpected fit failures remain contrast diagnostics", {
  fixture <- make_contrast_fixture()
  testthat::local_mocked_bindings(
    .fit_one_association = function(...) stop("Unexpected fit failure"),
    .package = "phewasFlow"
  )

  contrast <- run_phewas_contrast(fixture$data, fixture$base, fixture$adjusted)

  expect_equal(nrow(contrast$associations), 2)
  expect_true(all(contrast$associations$status == "error"))
  expect_true(all(grepl("Unexpected fit failure", contrast$diagnostics$message)))
  expect_true(all(is.na(contrast$associations$base_estimate)))
  expect_equal(contrast$run_metadata$testing_universe, 0)
})
