test_that("synthetic workflow runs from models through shards and plots", {
  example <- phewas_example_data(n = 220, seed = 401)
  expect_true(all(vapply(example, data.table::is.data.table, logical(1))))
  eligibility <- list(
    min_n = 100L, min_cases = 10L, min_controls = 10L,
    min_outcome_levels = 3L
  )
  direct_spec <- phewas_spec(
    "end_to_end", "pgs", "phenotypes_as_outcomes",
    example$forward_metadata, covariates = "age", anchor_type = "numeric",
    anchor_transform = "zscore", eligibility = eligibility,
    fdr_threshold = 0.05
  )
  direct <- run_phewas(example$data, direct_spec)
  expect_identical(direct$engine, c("logistf", "lm", "glm.nb", "clm"))
  expect_true(all(direct$status == "ok"))
  expect_equal(direct$q_value, stats::p.adjust(direct$p_value, "BH"),
               tolerance = 1e-12)

  reverse_spec <- phewas_spec(
    "end_to_end_reverse", "endpoint_continuous",
    "phenotypes_as_predictors", example$reverse_metadata,
    covariates = "age", outcome_type = "continuous",
    eligibility = eligibility, fdr_threshold = 0.05
  )
  reverse <- run_phewas(example$data, reverse_spec)
  expect_true(all(reverse$status == "ok"))
  expect_identical(names(reverse), names(direct))

  base_spec <- direct_spec
  base_spec$analysis_id <- "end_to_end_base"
  base_spec$covariates <- character()
  contrast <- run_phewas_contrast(example$data, base_spec, direct_spec)
  expect_true(all(contrast$associations$status == "ok"))
  expect_true(any(is.finite(contrast$associations$delta_r_squared)))

  directory <- withr::local_tempdir()
  data_path <- file.path(directory, "analysis.tsv")
  metadata_path <- file.path(directory, "metadata.tsv")
  config_path <- file.path(directory, "analysis.yml")
  output_path <- file.path(directory, "output")
  data.table::fwrite(example$data, data_path, sep = "\t")
  metadata <- example$forward_metadata
  metadata$levels <- vapply(metadata$levels, paste, collapse = "|", character(1))
  metadata$scores <- vapply(metadata$scores, paste, collapse = "|", character(1))
  data.table::fwrite(metadata, metadata_path, sep = "\t", na = "")
  yaml::write_yaml(list(
    schema_version = 1,
    inputs = list(data = data_path, phenotype_metadata = metadata_path),
    analysis = list(
      analysis_id = "end_to_end", direction = "phenotypes_as_outcomes",
      anchor = "pgs", anchor_type = "numeric", anchor_transform = "zscore",
      covariates = "age", fdr_threshold = 0.05,
      eligibility = eligibility
    ),
    execution = list(shards = 2L, backend = "sequential", workers = 1L),
    output = list(directory = output_path)
  ), config_path)

  create_phewas_manifest(config_path)
  run_phewas_shard(config_path, 1L)
  run_phewas_shard(config_path, 2L)
  combined <- combine_phewas_shards(config_path)
  combined <- combined[match(direct$phenotype, combined$phenotype), ]
  expect_equal(combined$estimate, direct$estimate, tolerance = 1e-10)
  expect_equal(combined$q_value, direct$q_value, tolerance = 1e-10)

  comparison <- compare_phewas(direct, combined)
  expect_equal(comparison$metrics$n_mutually_successful, 4)
  expect_s3_class(ggplot2::ggplot_build(plot_phewas_manhattan(direct)),
                  "ggplot_built")
  expect_s3_class(
    ggplot2::ggplot_build(plot_phewas_volcano(direct, top_label_n = 0L)),
    "ggplot_built"
  )
  expect_s3_class(ggplot2::ggplot_build(plot_phewas_concordance(comparison)),
                  "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(plot_phewas_category_overlap(comparison)),
                  "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(plot_phewas_forest(
    direct, phenotype_ids = direct$phenotype
  )), "ggplot_built")
})
