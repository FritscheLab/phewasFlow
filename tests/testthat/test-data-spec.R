test_that("strict assembly preserves the sample spine and reports unmatched IDs", {
  sample <- data.frame(id = c("p2", "p1", "p3"), enrolled = TRUE)
  phenotype <- data.frame(id = c("p1", "p2", "outside"), trait = c(1, 0, 1))
  anchor <- data.frame(id = c("p3", "p1"), pgs = c(0.4, -0.2))
  phenotype_before <- phenotype

  result <- assemble_phewas_data(
    sample, phenotypes = phenotype, anchors = anchor, id = "id"
  )

  expect_s3_class(result, "data.frame")
  expect_true(data.table::is.data.table(result))
  expect_identical(result$id, sample$id)
  expect_equal(result$trait, c(0, 1, NA))
  expect_equal(result$pgs, c(NA, -0.2, 0.4))
  expect_identical(phenotype, phenotype_before)

  report <- attr(result, "join_report")
  expect_true(data.table::is.data.table(report))
  expect_equal(report$unmatched_input_rows[report$table == "phenotypes"], 1)
  expect_equal(report$unmatched_sample_rows[report$table == "phenotypes"], 1)
  expect_equal(report$unmatched_sample_rows[report$table == "anchors"], 1)
})

test_that("data.table inputs are copied and never modified by reference", {
  sample <- data.table::data.table(id = c("p1", "p2"), selected = TRUE)
  phenotypes <- data.table::data.table(id = c("p1", "p2"), trait = c(0L, 1L))
  sample_before <- data.table::copy(sample)
  phenotypes_before <- data.table::copy(phenotypes)

  result <- assemble_phewas_data(sample, phenotypes = phenotypes, id = "id")
  result[, trait := trait + 10L]

  expect_true(data.table::is.data.table(result))
  expect_identical(sample, sample_before)
  expect_identical(phenotypes, phenotypes_before)
})

test_that("assembly rejects unsafe keys and ambiguous columns", {
  expect_error(
    assemble_phewas_data(
      data.frame(id = c(1L, 1L)), data.frame(id = 1L, x = 1), id = "id"
    ),
    "one row per ID"
  )
  expect_error(
    assemble_phewas_data(
      data.frame(id = c("", "p2")), data.frame(id = "p2", x = 1), id = "id"
    ),
    "missing or blank"
  )
  expect_error(
    assemble_phewas_data(
      data.frame(id = 1L, x = 1), data.frame(id = 1L, x = 2), id = "id"
    ),
    "Column collision"
  )
  expect_error(
    assemble_phewas_data(
      data.frame(id = 1L), data.frame(id = "1", x = 2), id = "id"
    ),
    "different type"
  )
})

test_that("spec validation keeps explicit phenotype and ordinal metadata", {
  metadata <- data.frame(
    phenotype = c("binary_trait", "ordered_trait"),
    description = c("Binary", "Ordered"),
    group = c("A", "B"),
    groupnum = c(1, 2),
    color = c("#112233", "#445566"),
    variable_type = c("binary", "ordinal"),
    outcome_type = c("binary", "ordinal"),
    reference = c("0", NA),
    levels = c(NA, "low|middle|high"),
    scores = NA_character_,
    offset = NA_character_,
    stringsAsFactors = FALSE
  )
  eligibility <- list(
    min_n = 20L, min_cases = 5L, min_controls = 5L,
    min_outcome_levels = 3L
  )
  spec <- phewas_spec(
    "explicit", anchor = "pgs", direction = "phenotypes_as_outcomes",
    phenotypes = metadata, covariates = "age", anchor_type = "numeric",
    anchor_transform = "zscore", eligibility = eligibility,
    fdr_threshold = 0.05
  )

  expect_s3_class(spec, "phewas_spec")
  expect_identical(spec$eligibility, eligibility)
  expect_identical(
    phewasFlow:::.metadata_sequence(spec$phenotypes$levels[[2L]]),
    c("low", "middle", "high")
  )
  expect_identical(validate_phewas_spec(spec), spec)
  expect_error(
    validate_phewas_spec(spec, data.frame(pgs = 1, age = 2, binary_trait = 0)),
    "missing declared columns"
  )
})

test_that("binary phenotype metadata defaults type and reference", {
  eligibility <- list(
    min_n = 20L, min_cases = 5L, min_controls = 5L,
    min_outcome_levels = 3L
  )
  metadata <- data.table::data.table(
    phenotype = c("default_absent", "default_blank", "custom", "ordered"),
    description = c("Absent", "Blank", "Custom", "Ordered"),
    group = c("A", "B", "C", "D"),
    groupnum = 1:4,
    color = c("#112233", "#223344", "#334455", "#445566"),
    variable_type = c(NA, " ", "binary", "ordinal"),
    outcome_type = c("binary", "binary", "binary", "ordinal"),
    reference = c(NA, " ", "control", NA),
    levels = c(NA, NA, NA, "low|middle|high")
  )
  before <- data.table::copy(metadata)

  spec <- phewas_spec(
    "binary_defaults", "pgs", "phenotypes_as_outcomes", metadata,
    anchor_type = "numeric", eligibility = eligibility,
    fdr_threshold = 0.05
  )

  expect_identical(metadata, before)
  expect_identical(
    spec$phenotypes$variable_type,
    c("binary", "binary", "binary", "ordinal")
  )
  expect_identical(
    spec$phenotypes$reference,
    c("0", "0", "control", NA_character_)
  )
  expect_identical(validate_phewas_spec(spec), spec)

  core <- data.frame(
    phenotype = "diagnosis", description = "Diagnosis",
    group = "Diagnoses", groupnum = 1, color = "#112233",
    stringsAsFactors = FALSE
  )
  forward <- transform(core, outcome_type = "binary")
  forward_spec <- phewas_spec(
    "forward_defaults", "pgs", "phenotypes_as_outcomes", forward,
    anchor_type = "numeric", eligibility = eligibility,
    fdr_threshold = 0.05
  )
  expect_identical(forward_spec$phenotypes$variable_type, "binary")
  expect_identical(forward_spec$phenotypes$reference, "0")
  expect_error(
    phewas_spec(
      "forward_missing_outcome_type", "pgs", "phenotypes_as_outcomes", core,
      anchor_type = "numeric", eligibility = eligibility,
      fdr_threshold = 0.05
    ),
    "outcome_type"
  )

  reverse_spec <- phewas_spec(
    "reverse_defaults", "endpoint", "phenotypes_as_predictors", core,
    outcome_type = "continuous", anchor_type = "numeric",
    eligibility = eligibility, fdr_threshold = 0.05
  )
  expect_identical(reverse_spec$phenotypes$variable_type, "binary")
  expect_identical(reverse_spec$phenotypes$reference, "0")
})

test_that("spec rejects inference-prone or inconsistent declarations", {
  metadata <- data.frame(
    phenotype = "trait", description = "Trait",
    group = "A", groupnum = 1, color = "black",
    variable_type = "binary", outcome_type = "binary",
    reference = NA_character_, stringsAsFactors = FALSE
  )
  eligibility <- list(
    min_n = 10L, min_cases = 2L, min_controls = 2L,
    min_outcome_levels = 3L
  )
  expect_error(
    phewas_spec(
      "bad", "pgs", "phenotypes_as_outcomes", metadata,
      anchor_type = "numeric", eligibility = eligibility
    ),
    "fdr_threshold.*explicitly supplied"
  )
  defaulted <- phewas_spec(
    "defaulted", "pgs", "phenotypes_as_outcomes", metadata,
    anchor_type = "numeric", eligibility = eligibility,
    fdr_threshold = 0.05
  )
  expect_identical(defaulted$phenotypes$reference, "0")
  metadata$reference <- "0"
  expect_error(
    phewas_spec(
      "bad", "pgs", "phenotypes_as_outcomes", metadata,
      eligibility = eligibility, fdr_threshold = 0.05
    ),
    "anchor_type"
  )
  expect_error(
    phewas_spec(
      "bad", "pgs", "phenotypes_as_outcomes", metadata,
      anchor_type = "binary", anchor_reference = "0",
      anchor_transform = "zscore", eligibility = eligibility,
      fdr_threshold = 0.05
    ),
    "Only a numeric predictor"
  )
  expect_error(
    phewas_spec(
      "bad", "pgs", "phenotypes_as_outcomes", metadata,
      anchor_type = "numeric",
      eligibility = list(min_n = 10L, min_cases = 2L),
      fdr_threshold = 0.05
    ),
    "must explicitly define"
  )
})

test_that("canonical metadata defines the analysis-data contract", {
  metadata <- data.table::data.table(
    phenotype = "250.2",
    description = "Type 2 diabetes",
    group = "Endocrine",
    groupnum = 3,
    color = "#7C3AED",
    variable_type = "binary",
    outcome_type = "binary",
    reference = "0"
  )
  before <- data.table::copy(metadata)
  eligibility <- list(
    min_n = 20L, min_cases = 5L, min_controls = 5L,
    min_outcome_levels = 3L
  )
  spec <- phewas_spec(
    "icon", "pgs", "phenotypes_as_outcomes", metadata,
    anchor_type = "numeric", eligibility = eligibility,
    fdr_threshold = 0.05
  )

  expect_identical(metadata, before)
  expect_identical(spec$phenotypes$phenotype, "250.2")
  expect_identical(spec$phenotypes$description, "Type 2 diabetes")
  expect_identical(spec$phenotypes$group, "Endocrine")
  expect_identical(spec$phenotypes$groupnum, 3)
  expect_false(any(c(
    "phenotype_id", "phecode", "phenotype_column", "column", "label", "category",
    "category_order"
  ) %in% names(spec$phenotypes)))
  expect_identical(validate_phewas_spec(spec), spec)
})

test_that("unambiguous legacy metadata normalizes without retaining aliases", {
  legacy <- data.frame(
    phenotype_id = "250.2",
    phecode = "250.2",
    phenotype_column = "250.2",
    column = "250.2",
    label = "Type 2 diabetes",
    category = "Endocrine",
    category_order = "3",
    color = "#7C3AED",
    variable_type = "binary",
    outcome_type = "binary",
    reference = "0",
    stringsAsFactors = FALSE
  )
  before <- legacy
  spec <- phewas_spec(
    "legacy", "pgs", "phenotypes_as_outcomes", legacy,
    anchor_type = "numeric",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )

  expect_identical(legacy, before)
  expect_identical(spec$phenotypes$phenotype, "250.2")
  expect_identical(spec$phenotypes$description, "Type 2 diabetes")
  expect_identical(spec$phenotypes$group, "Endocrine")
  expect_identical(spec$phenotypes$groupnum, 3)
  expect_false(any(c(
    "phenotype_id", "phecode", "phenotype_column", "column", "label", "category",
    "category_order"
  ) %in% names(spec$phenotypes)))
})

test_that("conflicting canonical and legacy metadata is rejected", {
  metadata <- data.frame(
    phenotype = "250.2",
    phenotype_id = "diabetes",
    description = "Type 2 diabetes",
    group = "Endocrine",
    groupnum = 3,
    color = "#7C3AED",
    variable_type = "binary",
    outcome_type = "binary",
    reference = "0",
    stringsAsFactors = FALSE
  )
  make_spec <- function(x) {
    phewas_spec(
      "ambiguous", "pgs", "phenotypes_as_outcomes", x,
      anchor_type = "numeric",
      eligibility = list(
        min_n = 20L, min_cases = 5L, min_controls = 5L,
        min_outcome_levels = 3L
      ),
      fdr_threshold = 0.05
    )
  }

  expect_error(make_spec(metadata), "`phenotype`.*`phenotype_id`.*disagree")

  aliases_only <- metadata
  aliases_only$phenotype <- NULL
  aliases_only$phecode <- "250.2"
  expect_error(
    make_spec(aliases_only),
    "`phenotype`.*`phecode`.*disagree"
  )

  matching <- metadata
  matching$phenotype_id <- matching$phenotype
  matching$label <- "Different description"
  expect_error(make_spec(matching), "`description`.*`label`.*disagree")

  mismatched_column <- metadata
  mismatched_column$phenotype_id <- mismatched_column$phenotype
  mismatched_column$column <- "phecode_250_2"
  expect_error(
    make_spec(mismatched_column),
    "`phenotype`.*`column`.*disagree"
  )

  mismatched_phenotype_column <- metadata
  mismatched_phenotype_column$phenotype_id <-
    mismatched_phenotype_column$phenotype
  mismatched_phenotype_column$phenotype_column <- "phecode_250_2"
  expect_error(
    make_spec(mismatched_phenotype_column),
    "`phenotype`.*`phenotype_column`.*disagree"
  )

  column_only <- metadata
  column_only$phenotype <- NULL
  column_only$phenotype_id <- NULL
  column_only$column <- "250.2"
  column_only_spec <- make_spec(column_only)
  expect_identical(column_only_spec$phenotypes$phenotype, "250.2")
  expect_false("column" %in% names(column_only_spec$phenotypes))
})

test_that("groups have one unique order and one valid color", {
  metadata <- data.frame(
    phenotype = c("a1", "a2", "b1"),
    description = c("A one", "A two", "B one"),
    group = c("A", "A", "B"),
    groupnum = c(1, 1, 2),
    color = c("#FF0000", "#FF0000", "blue"),
    variable_type = "numeric",
    outcome_type = "continuous",
    stringsAsFactors = FALSE
  )
  make_spec <- function(x) {
    phewas_spec(
      "groups", "pgs", "phenotypes_as_outcomes", x,
      anchor_type = "numeric",
      eligibility = list(
        min_n = 20L, min_cases = 5L, min_controls = 5L,
        min_outcome_levels = 3L
      ),
      fdr_threshold = 0.05
    )
  }

  expect_s3_class(make_spec(metadata), "phewas_spec")

  inconsistent_order <- metadata
  inconsistent_order$groupnum[[2L]] <- 3
  expect_error(make_spec(inconsistent_order), "one `groupnum`.*A")

  inconsistent_color <- metadata
  inconsistent_color$color[[2L]] <- "blue"
  expect_error(make_spec(inconsistent_color), "one `color`.*A")

  invalid_color <- metadata
  invalid_color$color[[1L]] <- "not-a-real-color"
  expect_error(make_spec(invalid_color), "invalid R color")

  tied_order <- metadata
  tied_order$groupnum[tied_order$group == "B"] <- 1
  expect_error(make_spec(tied_order), "uniquely order groups")
})

test_that("exclusion declarations require a paired, type-safe marker", {
  metadata <- data.frame(
    phenotype = "trait", description = "Trait",
    group = "A", groupnum = 1, color = "black",
    variable_type = "numeric", outcome_type = "continuous"
  )
  eligibility <- list(
    min_n = 10L, min_cases = 2L, min_controls = 2L,
    min_outcome_levels = 3L
  )
  expect_error(
    phewas_spec(
      "bad", "pgs", "phenotypes_as_outcomes", metadata,
      anchor_type = "numeric", exclusion_values = "drop",
      eligibility = eligibility, fdr_threshold = 0.05
    ),
    "requires an `exclusion_column`"
  )
  spec <- phewas_spec(
    "excluded", "pgs", "phenotypes_as_outcomes", metadata,
    anchor_type = "numeric", exclusion_column = "cohort_marker",
    exclusion_values = c("drop", "withdrawn"), eligibility = eligibility,
    fdr_threshold = 0.05
  )
  expect_error(
    validate_phewas_spec(
      spec,
      data.frame(pgs = 1:12, trait = 1:12, cohort_marker = 0L)
    ),
    "not type-compatible"
  )
})

test_that("mutable specifications recheck every model declaration", {
  example <- phewas_example_data(n = 100)
  spec <- phewas_spec(
    "validation", "pgs", "phenotypes_as_outcomes",
    example$forward_metadata, anchor_type = "numeric", covariates = "age",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  changes <- list(
    list(field = "anchor", value = "age", message = "also be a covariate"),
    list(field = "direction", value = "invalid", message = "arg"),
    list(field = "anchor_transform", value = "invalid", message = "arg"),
    list(field = "anchor_type", value = "binary", message = "anchor_reference"),
    list(field = "covariates", value = "pheno_binary", message = "scanned phenotypes"),
    list(field = "fdr_threshold", value = 2, message = "fdr_threshold")
  )
  for (change in changes) {
    invalid <- spec
    invalid[[change$field]] <- change$value
    expect_error(validate_phewas_spec(invalid), change$message)
    expect_error(run_phewas(example$data, invalid), change$message)
  }

  for (value in list(Inf, NaN, .Machine$integer.max + 1, 20.5)) {
    invalid <- spec
    invalid$eligibility$min_n <- value
    expect_error(validate_phewas_spec(invalid), "must be one integer")
  }
  invalid <- spec
  invalid$eligibility <- c(invalid$eligibility, list(min_n = 2L))
  expect_error(validate_phewas_spec(invalid), "explicitly named list")
})

test_that("factor metadata keeps its declared numeric values", {
  metadata <- data.frame(
    phenotype = c("a", "b"), description = c("A", "B"),
    group = c("A", "B"), groupnum = factor(c("20", "10")),
    color = c("red", "blue"), outcome_type = "binary"
  )
  spec <- phewas_spec(
    "factor_metadata", "pgs", "phenotypes_as_outcomes", metadata,
    anchor_type = "numeric",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  expect_equal(spec$phenotypes$groupnum, c(20, 10))
  expect_equal(
    phewasFlow:::.metadata_sequence(factor(c("1", "5", "10")), numeric = TRUE),
    c(1, 5, 10)
  )
})

test_that("synthetic sample sizes reject truncation and overflow", {
  for (value in list(100.5, "100", Inf, NaN, .Machine$integer.max + 1)) {
    expect_error(phewas_example_data(n = value), "must be one integer")
  }
})
