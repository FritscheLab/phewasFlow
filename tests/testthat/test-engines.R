make_engine_fixture <- function() {
  set.seed(20260818)
  n <- 320
  x <- stats::rnorm(n)
  covariate <- stats::rnorm(n)
  binary <- ifelse(1.1 * x + 0.2 * covariate + stats::rlogis(n) > 0,
                   "case", "control")
  continuous <- 1.4 * x - 0.3 * covariate + stats::rnorm(n, sd = 0.8)
  log_exposure <- stats::runif(n, -0.4, 0.5)
  count <- stats::rnbinom(
    n, mu = exp(log_exposure + 0.35 * x + 0.1 * covariate), size = 1.7
  )
  ordinal <- cut(
    1.0 * x + 0.2 * covariate + stats::rlogis(n),
    breaks = c(-Inf, -0.7, 0.6, Inf),
    labels = c("low", "middle", "high"), ordered_result = TRUE
  )
  data <- data.frame(x, covariate, binary, continuous, count, ordinal,
                     log_exposure)
  metadata <- data.frame(
    phenotype = c("binary", "continuous", "count", "ordinal"),
    description = c("Binary", "Continuous", "Count", "Ordinal"),
    group = "Synthetic", groupnum = 1,
    color = "#336699",
    variable_type = c("binary", "numeric", "count", "ordinal"),
    outcome_type = c("binary", "continuous", "count", "ordinal"),
    reference = c("control", NA, NA, NA),
    offset = c(NA, NA, "log_exposure", NA),
    stringsAsFactors = FALSE
  )
  metadata$levels <- I(list(
    character(), character(), character(), c("low", "middle", "high")
  ))
  metadata$scores <- I(replicate(4, numeric(), simplify = FALSE))
  spec <- phewas_spec(
    "engines", "x", "phenotypes_as_outcomes", metadata,
    covariates = "covariate", anchor_type = "numeric",
    eligibility = list(
      min_n = 50L, min_cases = 10L, min_controls = 10L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  list(data = data, metadata = metadata, spec = spec)
}

test_that("Firth logistic engine matches a direct logistf fit", {
  fixture <- make_engine_fixture()
  result <- run_phewas(
    fixture$data, fixture$spec, phenotype_ids = "binary"
  )
  direct_data <- transform(
    fixture$data,
    binary = as.integer(binary == "case")
  )
  direct <- logistf::logistf(
    binary ~ x + covariate, data = direct_data, pl = TRUE,
    control = logistf::logistf.control(maxit = 1000),
    plcontrol = logistf::logistpl.control(maxit = 1000)
  )

  expect_identical(result$status, "ok")
  expect_identical(result$test_method, "penalized_likelihood_ratio")
  expect_equal(result$estimate, unname(stats::coef(direct)[["x"]]),
               tolerance = 1e-10)
  expect_equal(result$conf_low, unname(direct$ci.lower[["x"]]),
               tolerance = 1e-10)
  expect_equal(result$native_effect, exp(result$estimate), tolerance = 1e-14)
})

test_that("linear engine matches lm and reports difference effects", {
  fixture <- make_engine_fixture()
  result <- run_phewas(
    fixture$data, fixture$spec, phenotype_ids = "continuous"
  )
  direct <- stats::lm(continuous ~ x + covariate, data = fixture$data)

  expect_identical(result$status, "ok")
  expect_identical(result$engine, "lm")
  expect_identical(result$effect_measure, "mean_difference")
  expect_equal(result$estimate, unname(stats::coef(direct)[["x"]]),
               tolerance = 1e-12)
  expect_equal(result$native_effect, result$estimate)
})

test_that("negative-binomial engine honors a declared log-offset", {
  fixture <- make_engine_fixture()
  result <- run_phewas(fixture$data, fixture$spec, phenotype_ids = "count")
  direct <- MASS::glm.nb(
    count ~ x + covariate + offset(log_exposure), data = fixture$data,
    link = log
  )

  expect_identical(result$status, "ok")
  expect_identical(result$offset, "log_exposure")
  expect_identical(result$effect_measure, "incidence_rate_ratio")
  expect_equal(result$estimate, unname(stats::coef(direct)[["x"]]),
               tolerance = 1e-10)
  expect_equal(result$theta, direct$theta, tolerance = 1e-10)
})

test_that("negative-binomial engine also fits without an offset", {
  fixture <- make_engine_fixture()
  row <- fixture$metadata[fixture$metadata$phenotype == "count", ]
  row$offset <- NA_character_
  spec <- phewas_spec(
    "count_without_offset", "x", "phenotypes_as_outcomes", row,
    covariates = "covariate", anchor_type = "numeric",
    eligibility = fixture$spec$eligibility, fdr_threshold = 0.05
  )
  result <- run_phewas(fixture$data, spec)
  direct <- MASS::glm.nb(count ~ x + covariate, data = fixture$data, link = log)

  expect_identical(result$status, "ok")
  expect_true(is.na(result$offset))
  expect_equal(result$estimate, unname(stats::coef(direct)[["x"]]),
               tolerance = 1e-10)
})

test_that("ordinal engine uses declared level order and higher-category odds", {
  fixture <- make_engine_fixture()
  result <- run_phewas(fixture$data, fixture$spec, phenotype_ids = "ordinal")
  direct <- ordinal::clm(ordinal ~ x + covariate, data = fixture$data,
                         link = "logit")

  expect_identical(result$status, "ok")
  expect_identical(result$ordinal_levels, "low|middle|high")
  expect_identical(result$effect_measure, "common_odds_ratio")
  expect_gt(result$estimate, 0)
  expect_equal(result$estimate, unname(direct$beta[["x"]]),
               tolerance = 1e-10)
})

test_that("separation is fitted by the real Firth engine", {
  data <- data.frame(
    x = c(rep(-1, 12), rep(1, 12)),
    separated = c(rep("control", 12), rep("case", 12))
  )
  metadata <- data.frame(
    phenotype = "separated", description = "Separated",
    group = "Synthetic", groupnum = 1, color = "black",
    variable_type = "binary", outcome_type = "binary",
    reference = "control", stringsAsFactors = FALSE
  )
  spec <- phewas_spec(
    "separation", "x", "phenotypes_as_outcomes", metadata,
    anchor_type = "numeric",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  result <- run_phewas(data, spec)

  expect_identical(result$status, "ok")
  expect_true(is.finite(result$estimate))
  expect_true(is.finite(result$conf_high))
  expect_identical(result$engine, "logistf")
})

test_that("failed negative-binomial dispersion estimation is not a success", {
  set.seed(3)
  data <- data.frame(x = stats::rnorm(100))
  data$count <- stats::rpois(100, exp(0.3 * data$x))
  metadata <- data.frame(
    phenotype = "count", description = "Count", group = "Synthetic",
    groupnum = 1, color = "black", variable_type = "count",
    outcome_type = "count"
  )
  spec <- phewas_spec(
    "dispersion", "x", "phenotypes_as_outcomes", metadata,
    anchor_type = "numeric",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  result <- run_phewas(data, spec)

  expect_identical(result$status, "error")
  expect_identical(result$reason_code, "nonconvergence")
  expect_false(result$converged)
  expect_match(result$warnings, "iteration limit reached|alternation limit reached")
  expect_true(is.na(result$p_value))
  expect_true(is.na(result$q_value))
  expect_equal(result$testing_family_size, 0)
})
