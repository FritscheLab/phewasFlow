make_run_fixture <- function() {
  set.seed(101)
  n <- 180
  pgs <- stats::rnorm(n, mean = 2, sd = 3)
  age <- stats::rnorm(n, 55, 8)
  p1 <- 0.7 * pgs + 0.2 * age + stats::rnorm(n)
  p2 <- -0.4 * pgs + stats::rnorm(n)
  binary <- ifelse(0.5 * scale(pgs)[, 1] + stats::rlogis(n) > 0, "yes", "no")
  binary[c(1, 2, 3)] <- NA_character_
  data <- data.frame(pgs, age, p1, p2, binary)
  metadata <- data.frame(
    phenotype = c("p1", "p2", "binary"),
    description = c("P1", "P2", "Binary"),
    group = c("A", "B", "A"), groupnum = c(1, 2, 1),
    color = c("#111111", "#222222", "#111111"),
    variable_type = c("numeric", "numeric", "binary"),
    outcome_type = c("continuous", "continuous", "binary"),
    reference = c(NA, NA, "no"), stringsAsFactors = FALSE
  )
  spec <- phewas_spec(
    "run", "pgs", "phenotypes_as_outcomes", metadata,
    covariates = "age", anchor_type = "numeric", anchor_transform = "zscore",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  list(data = data, metadata = metadata, spec = spec)
}

test_that("run returns canonical rows and one global testing family", {
  fixture <- make_run_fixture()
  result <- run_phewas(fixture$data, fixture$spec)

  expect_s3_class(result, "phewas_result")
  expect_true(data.table::is.data.table(result))
  expect_identical(class(result), c("phewas_result", "data.table", "data.frame"))
  mutable_copy <- data.table::copy(result)
  expect_warning(mutable_copy[, audit_flag := TRUE], NA)
  expect_true(all(mutable_copy$audit_flag))
  expect_identical(result$phenotype, fixture$metadata$phenotype)
  expect_identical(result$description, fixture$metadata$description)
  expect_identical(result$group, fixture$metadata$group)
  expect_identical(result$groupnum, fixture$metadata$groupnum)
  expect_identical(result$color, fixture$metadata$color)
  expect_false(any(c(
    "phenotype_id", "phecode", "phenotype_column", "column", "label",
    "category", "category_order"
  ) %in% names(result)))
  expect_true(all(c(
    "effective_formula", "effective_covariates",
    "dropped_invariant_covariates"
  ) %in% names(result)))
  expect_true(all(result$status == "ok"))
  expect_equal(result$q_value, stats::p.adjust(result$p_value, "BH"),
               tolerance = 1e-14)
  expect_equal(result$bonferroni_p,
               stats::p.adjust(result$p_value, "bonferroni"), tolerance = 1e-14)
  expect_true(all(result$testing_family_size == 3L))
  expect_equal(unique(result$transformation_mean), mean(fixture$data$pgs),
               tolerance = 1e-14)
  expect_equal(unique(result$transformation_sd), stats::sd(fixture$data$pgs),
               tolerance = 1e-14)
  expect_identical(result$n_complete[result$phenotype == "binary"], 177L)
  expect_null(attr(result, "fit", exact = TRUE))
  expect_s3_class(attr(result, "spec"), "phewas_spec")
  expect_true(is.list(attr(result, "run_metadata")))
  printed <- capture.output(print(result))
  expect_true(any(grepl("p_value_text", printed, fixed = TRUE)))
  expect_true(any(grepl("neg_log10_p", printed, fixed = TRUE)))
  selected_print <- capture.output(print(
    result[, .(phenotype, description, group, groupnum)]
  ))
  expect_true(any(grepl("group", selected_print, fixed = TRUE)))
  expect_true(any(grepl("groupnum", selected_print, fixed = TRUE)))
})

test_that("phenotype is both the canonical key and data column", {
  set.seed(102)
  n <- 80L
  data <- data.frame(
    pgs = stats::rnorm(n),
    ph2 = stats::rnorm(n)
  )
  metadata <- data.frame(
    phenotype = "ph2",
    description = "Long phenotype description",
    group = "Phenotype group",
    groupnum = 2,
    color = "#336699",
    variable_type = "numeric",
    outcome_type = "continuous",
    stringsAsFactors = FALSE
  )
  spec <- phewas_spec(
    "metadata", "pgs", "phenotypes_as_outcomes", metadata,
    anchor_type = "numeric",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )

  result <- run_phewas(data, spec)

  expect_identical(result$phenotype, "ph2")
  expect_identical(result$description, "Long phenotype description")
  expect_identical(result$group, "Phenotype group")
  expect_identical(result$groupnum, 2)
})

test_that("run output never retains matching legacy spec aliases", {
  fixture <- make_run_fixture()
  legacy_spec <- fixture$spec
  legacy_spec$phenotypes$phenotype_id <- legacy_spec$phenotypes$phenotype
  legacy_spec$phenotypes$phenotype_column <- legacy_spec$phenotypes$phenotype
  legacy_spec$phenotypes$column <- legacy_spec$phenotypes$phenotype
  legacy_spec$phenotypes$label <- legacy_spec$phenotypes$description
  legacy_spec$phenotypes$category <- legacy_spec$phenotypes$group
  legacy_spec$phenotypes$category_order <- legacy_spec$phenotypes$groupnum

  result <- run_phewas(fixture$data, legacy_spec)
  result_metadata <- attr(result, "spec", exact = TRUE)$phenotypes

  expect_false(any(c(
    "phenotype_id", "phecode", "phenotype_column", "column", "label",
    "category", "category_order"
  ) %in% names(result_metadata)))
  expect_identical(result_metadata$phenotype, fixture$metadata$phenotype)
})

test_that("exclusion markers filter both directions before scaling and completeness", {
  fixture <- make_run_fixture()
  n <- nrow(fixture$data)
  fixture$data$cohort_marker <- "keep"
  fixture$data$cohort_marker[1:7] <- "drop"
  fixture$data$cohort_marker[8:10] <- "withdrawn"
  fixture$data$cohort_marker[11] <- NA_character_
  fixture$data$pgs[1] <- 1e6
  fixture$data$age[12] <- NA_real_
  expected_base <- setdiff(seq_len(n), c(1:10, 12))

  forward_spec <- phewas_spec(
    "excluded_forward", "pgs", "phenotypes_as_outcomes",
    fixture$metadata, covariates = "age", anchor_type = "numeric",
    anchor_transform = "zscore", exclusion_column = "cohort_marker",
    exclusion_values = c("drop", "withdrawn"),
    eligibility = fixture$spec$eligibility, fdr_threshold = 0.05
  )
  forward <- run_phewas(fixture$data, forward_spec)

  expect_true(all(forward$n_before_exclusion == n))
  expect_true(all(forward$n_excluded == 10L))
  expect_true(all(forward$n_after_exclusion == n - 10L))
  expect_equal(
    unique(forward$transformation_mean),
    mean(fixture$data$pgs[expected_base]), tolerance = 1e-14
  )
  forward_metadata <- attr(forward, "run_metadata")
  expect_identical(forward_metadata$exclusion_column, "cohort_marker")
  expect_identical(forward_metadata$exclusion_values, c("drop", "withdrawn"))
  expect_identical(forward_metadata$n_excluded, 10L)

  reverse_metadata <- data.frame(
    phenotype = "p1", description = "P1",
    group = "A", groupnum = 1, color = "black",
    variable_type = "numeric"
  )
  reverse_spec <- phewas_spec(
    "excluded_reverse", "p2", "phenotypes_as_predictors",
    reverse_metadata, outcome_type = "continuous", covariates = "age",
    exclusion_column = "cohort_marker", exclusion_values = c("drop", "withdrawn"),
    eligibility = fixture$spec$eligibility, fdr_threshold = 0.05
  )
  reverse <- run_phewas(fixture$data, reverse_spec)

  expect_identical(reverse$status, "ok")
  expect_identical(reverse$n_excluded, 10L)
  expect_identical(reverse$n_base, as.integer(length(expected_base)))
  expect_identical(reverse$n_complete, as.integer(length(expected_base)))
})

test_that("expected data failures remain explicit rows", {
  fixture <- make_run_fixture()
  fixture$data$constant <- 1
  fixture$data$duplicate_age <- fixture$data$age
  extra <- data.frame(
    phenotype = c("constant", "p1"),
    description = c("Constant", "Rank deficient"),
    group = "C", groupnum = 3, color = "black",
    variable_type = "numeric", outcome_type = "continuous",
    reference = NA_character_, stringsAsFactors = FALSE
  )

  constant_spec <- phewas_spec(
    "constant", "pgs", "phenotypes_as_outcomes", extra[1, ],
    anchor_type = "numeric",
    eligibility = fixture$spec$eligibility, fdr_threshold = 0.05
  )
  constant <- run_phewas(fixture$data, constant_spec)
  expect_identical(constant$status, "skipped")
  expect_identical(constant$reason_code, "constant_outcome")
  expect_true(is.na(constant$q_value))

  rank_spec <- phewas_spec(
    "rank", "pgs", "phenotypes_as_outcomes", extra[2, ],
    covariates = c("age", "duplicate_age"), anchor_type = "numeric",
    eligibility = fixture$spec$eligibility, fdr_threshold = 0.05
  )
  rank <- run_phewas(fixture$data, rank_spec)
  expect_identical(rank$status, "skipped")
  expect_identical(rank$reason_code, "rank_deficient")
  expect_identical(rank$effective_covariates, "age|duplicate_age")
  expect_true(is.na(rank$dropped_invariant_covariates))
  expect_match(rank$effective_formula, "`age` \\+ `duplicate_age`$")

  binary_metadata <- data.frame(
    phenotype = "none", description = "No cases",
    group = "C", groupnum = 3, color = "black",
    variable_type = "binary", outcome_type = "binary", reference = "no",
    stringsAsFactors = FALSE
  )
  fixture$data$none <- "no"
  binary_spec <- phewas_spec(
    "none", "pgs", "phenotypes_as_outcomes", binary_metadata,
    anchor_type = "numeric", eligibility = fixture$spec$eligibility,
    fdr_threshold = 0.05
  )
  no_cases <- run_phewas(fixture$data, binary_spec)
  expect_identical(no_cases$status, "skipped")
  expect_identical(no_cases$reason_code, "insufficient_cases")
})

test_that("outcome-specific invariant covariates are omitted from the fit", {
  skip_if_not_installed("logistf")
  set.seed(204)
  n_per_sex <- 80L
  n <- 2L * n_per_sex
  SNPSEX <- factor(
    rep(c("female", "male"), each = n_per_sex),
    levels = c("female", "male", "unknown")
  )
  x <- stats::rnorm(n)
  age <- stats::rnorm(n, 55, 8)
  array <- factor(
    sample(c("A", "B"), n, replace = TRUE),
    levels = c("A", "B", "unused")
  )
  all_sexes <- rep(c("control", "case"), length.out = n)
  male_only <- rep(NA_character_, n)
  male_only[SNPSEX == "male"] <- rep(
    c("control", "case"), length.out = n_per_sex
  )
  data <- data.frame(x, age, SNPSEX, array, all_sexes, male_only)
  metadata <- data.frame(
    phenotype = c("all_sexes", "male_only"),
    description = c("All sexes", "Male-only outcome"),
    group = "Synthetic", groupnum = 1, color = "#336699",
    variable_type = "binary", outcome_type = "binary",
    reference = "control", stringsAsFactors = FALSE
  )
  spec <- phewas_spec(
    "sex_specific", "x", "phenotypes_as_outcomes", metadata,
    covariates = c("age", "SNPSEX", "array"), anchor_type = "numeric",
    eligibility = list(
      min_n = 50L, min_cases = 20L, min_controls = 20L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )

  result <- run_phewas(data, spec)
  all_sexes_result <- result[result$phenotype == "all_sexes"]
  male_only_result <- result[result$phenotype == "male_only"]

  expect_true(all(result$status == "ok"))
  expect_identical(
    all_sexes_result$effective_covariates, "age|SNPSEX|array"
  )
  expect_true(is.na(all_sexes_result$dropped_invariant_covariates))
  expect_match(
    all_sexes_result$effective_formula,
    "`age` \\+ `SNPSEX` \\+ `array`$"
  )
  expect_identical(male_only_result$effective_covariates, "age|array")
  expect_identical(
    male_only_result$dropped_invariant_covariates, "SNPSEX"
  )
  expect_identical(
    male_only_result$effective_formula,
    "`male_only` ~ `x` + `age` + `array`"
  )
  expect_match(
    male_only_result$formula,
    "`age` \\+ `SNPSEX` \\+ `array`$"
  )

  direct_data <- data[SNPSEX == "male", , drop = FALSE]
  direct_data$male_binary <- as.integer(direct_data$male_only == "case")
  direct <- logistf::logistf(
    male_binary ~ x + age + array,
    data = direct_data,
    pl = TRUE,
    control = logistf::logistf.control(maxit = 1000L),
    plcontrol = logistf::logistpl.control(maxit = 1000L)
  )
  expect_equal(
    male_only_result$estimate,
    unname(stats::coef(direct)[["x"]]),
    tolerance = 1e-10
  )
})

test_that("log-scale correction remains authoritative after underflow", {
  rows <- data.frame(
    status = c("ok", "ok", "skipped"),
    log_p = c(-1000, log(0.02), NA_real_),
    stringsAsFactors = FALSE
  )
  corrected <- phewasFlow:::.apply_multiple_testing(rows, 0.05)

  expect_identical(corrected$p_value, NULL)
  expect_equal(corrected$q_value[[1L]], 0)
  expect_equal(corrected$neg_log10_q[[1L]],
               -( -1000 + log(2)) / log(10), tolerance = 1e-12)
  expect_match(corrected$q_value_text[[1L]], "e-4")
  expect_true(corrected$fdr_significant[[1L]])
  expect_true(is.na(corrected$q_value[[3L]]))
})

test_that("phenotypes can be binary and ordinal trend predictors", {
  set.seed(55)
  n <- 220
  binary <- sample(c("absent", "present"), n, replace = TRUE)
  stage <- sample(c("low", "middle", "high"), n, replace = TRUE)
  stage_score <- c(low = 0, middle = 1, high = 3)[stage]
  endpoint <- 0.6 * (binary == "present") + 0.4 * stage_score + stats::rnorm(n)
  data <- data.frame(binary, stage, endpoint)
  metadata <- data.frame(
    phenotype = c("binary", "stage"),
    description = c("Binary", "Stage"), group = "Predictors",
    groupnum = 1, color = "#445566",
    variable_type = c("binary", "ordinal"),
    reference = c("absent", NA), stringsAsFactors = FALSE
  )
  metadata$levels <- I(list(character(), c("low", "middle", "high")))
  metadata$scores <- I(list(numeric(), c(0, 1, 3)))
  spec <- phewas_spec(
    "reverse", "endpoint", "phenotypes_as_predictors", metadata,
    outcome_type = "continuous",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  result <- run_phewas(data, spec)

  expect_true(all(result$status == "ok"))
  expect_identical(result$response, rep("endpoint", 2))
  expect_identical(result$predictor, c("binary", "stage"))
  expect_identical(result$scores[result$phenotype == "stage"], "0|1|3")
  expect_true(all(result$estimate > 0))
})

test_that("reverse scans apply binary phenotype count thresholds", {
  n <- 120L
  endpoint <- rep(0:31, length.out = n)
  data <- data.frame(
    endpoint = endpoint,
    few_nonreference = c(rep(1L, 49L), rep(0L, n - 49L)),
    few_reference = c(rep(0L, 49L), rep(1L, n - 49L)),
    eligible = rep(c(0L, 1L), length.out = n)
  )
  metadata <- data.frame(
    phenotype = c("few_nonreference", "few_reference", "eligible"),
    description = c("Few non-reference", "Few reference", "Eligible"),
    group = "Diagnoses", groupnum = 1, color = "#445566",
    variable_type = "binary", reference = "0",
    stringsAsFactors = FALSE
  )
  spec <- phewas_spec(
    "reverse_ordinal_thresholds", "endpoint", "phenotypes_as_predictors",
    metadata,
    outcome_type = "ordinal", anchor_type = "ordinal",
    anchor_levels = as.character(0:31),
    eligibility = list(
      min_n = 100L, min_cases = 50L, min_controls = 50L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )

  result <- run_phewas(data, spec)

  expect_identical(
    result$status,
    c("skipped", "skipped", "ok")
  )
  expect_identical(
    result$reason_code,
    c("insufficient_cases", "insufficient_controls", "ok")
  )
  expect_match(result$message[[1L]], "49 non-reference.*min_cases is 50")
  expect_match(result$message[[2L]], "49 reference.*min_controls is 50")
  expect_identical(result$engine, rep("clm", 3L))
  expect_identical(result$outcome_type, rep("ordinal", 3L))
  expect_identical(result$effect_measure, rep("common_odds_ratio", 3L))
  expect_identical(result$ordinal_levels[[3L]], paste(0:31, collapse = "|"))
  expect_equal(result$native_effect[[3L]], exp(result$estimate[[3L]]))
})

test_that("ordinal predictor defaults are deterministic and recorded", {
  set.seed(56)
  stage <- rep(c("low", "middle", "high"), each = 50)
  endpoint <- match(stage, c("low", "middle", "high")) + stats::rnorm(150)
  metadata <- data.frame(
    phenotype = "stage", description = "Stage",
    group = "Predictors", groupnum = 1, color = "#445566",
    variable_type = "ordinal", reference = NA_character_,
    levels = "low|middle|high", scores = NA_character_,
    stringsAsFactors = FALSE
  )
  spec <- phewas_spec(
    "ordinal_default", "endpoint", "phenotypes_as_predictors", metadata,
    outcome_type = "continuous",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  result <- run_phewas(data.frame(stage, endpoint), spec)

  expect_identical(result$status, "ok")
  expect_identical(result$scores, "1|2|3")
  expect_gt(result$estimate, 0)
})

test_that("worker counts reject nonintegers and overflow before fitting", {
  fixture <- make_run_fixture()
  for (workers in list(Inf, NaN, 1.5, .Machine$integer.max + 1, 0, "2")) {
    expect_error(
      run_phewas(fixture$data, fixture$spec, workers = workers),
      "`workers` must be one positive integer"
    )
  }
})

test_that("sequential and PSOCK backends are numerically identical", {
  probe <- tryCatch(parallel::makePSOCKcluster(1L), error = function(error) NULL)
  if (is.null(probe)) {
    skip("This test environment does not permit local server sockets.")
  }
  parallel::stopCluster(probe)

  fixture <- make_run_fixture()
  sequential <- run_phewas(fixture$data, fixture$spec, backend = "sequential")
  parallel <- run_phewas(fixture$data, fixture$spec, backend = "psock", workers = 2L)

  columns <- c(
    "phenotype", "status", "reason_code", "estimate", "std_error",
    "conf_low", "conf_high", "log_p", "q_value", "bonferroni_p",
    "effective_formula", "effective_covariates",
    "dropped_invariant_covariates", "sample_fingerprint"
  )
  expect_equal(as.data.frame(parallel)[columns],
               as.data.frame(sequential)[columns], tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_identical(attr(parallel, "run_metadata")$backend, "psock")
})
