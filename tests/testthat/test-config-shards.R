make_shard_fixture <- function(directory, shards = 2L) {
  example <- phewas_example_data(n = 120, seed = 19)
  metadata <- data.frame(
    phenotype = c("pheno_continuous", "endpoint_continuous"),
    description = c("Continuous phenotype", "Continuous endpoint"),
    group = c("Measurements", "Endpoints"),
    groupnum = c(2L, 1L),
    color = c("#10B981", "#3B82F6"),
    variable_type = c("numeric", "numeric"),
    outcome_type = c("continuous", "continuous"),
    stringsAsFactors = FALSE
  )
  saveRDS(example$data, file.path(directory, "analysis.rds"))
  saveRDS(metadata, file.path(directory, "metadata.rds"))
  config <- list(
    schema_version = 1,
    inputs = list(
      data = "analysis.rds",
      phenotype_metadata = "metadata.rds"
    ),
    analysis = list(
      analysis_id = "config_shard_test",
      direction = "phenotypes_as_outcomes",
      anchor = "pgs",
      anchor_type = "numeric",
      anchor_transform = "none",
      covariates = c("age", "sex"),
      fdr_threshold = 0.05,
      eligibility = list(
        min_n = 20L,
        min_cases = 5L,
        min_controls = 5L,
        min_outcome_levels = 3L
      )
    ),
    execution = list(shards = shards, backend = "sequential", workers = 1L),
    output = list(directory = "output"),
    plot = list(
      output = "output/phewas.png",
      significance = "bh",
      label = list(),
      highlight = list()
    )
  )
  path <- file.path(directory, "phewasflow.yml")
  yaml::write_yaml(config, path)
  path
}

make_bundle_result <- function(seed = 43L) {
  example <- phewas_example_data(n = 120L, seed = seed)
  metadata <- example$forward_metadata[
    example$forward_metadata$outcome_type == "continuous", , drop = FALSE
  ]
  spec <- phewas_spec(
    analysis_id = "bundle_round_trip",
    anchor = "pgs",
    direction = "phenotypes_as_outcomes",
    phenotypes = metadata,
    covariates = "age",
    anchor_type = "numeric",
    eligibility = list(
      min_n = 20L, min_cases = 5L, min_controls = 5L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  run_phewas(example$data, spec)
}

complete_shard_fixture <- function(directory) {
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  path <- make_shard_fixture(directory)
  manifest <- create_phewas_manifest(path)
  shard_paths <- vapply(1:2, function(shard_id) {
    run_phewas_shard(path, shard_id)
  }, character(1L))
  result <- combine_phewas_shards(path)
  list(
    config = path,
    manifest = manifest,
    shards = shard_paths,
    result = result,
    output = file.path(directory, "output")
  )
}

test_that("YAML paths are relative to the config and produce a valid spec", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  config <- read_phewasflow_config(path)

  expect_s3_class(config, "phewasflow_config")
  expect_true(all(file.exists(unlist(config$inputs))))
  expect_identical(config$execution$shards, 2L)
  expect_identical(config$analysis$eligibility$min_outcome_levels, 3L)
  expect_identical(config$plot$group_display, "legend")

  spec <- config_to_phewas_spec(config)
  expect_s3_class(spec, "phewas_spec")
  expect_identical(spec$analysis_id, "config_shard_test")
  expect_identical(
    spec$phenotypes$phenotype,
    c("pheno_continuous", "endpoint_continuous")
  )
})

test_that("numeric-looking phenotype names remain exact strings", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory, shards = 1L)
  metadata <- readRDS(file.path(directory, "metadata.rds"))[1L, , drop = FALSE]
  metadata$phenotype <- "001"
  metadata_path <- file.path(directory, "leading-zero-metadata.tsv")
  data.table::fwrite(metadata, metadata_path, sep = "\t", na = "NA")

  data <- readRDS(file.path(directory, "analysis.rds"))
  data[["001"]] <- data$pheno_continuous
  saveRDS(data, file.path(directory, "analysis.rds"), version = 3L)

  raw <- yaml::read_yaml(path)
  raw$inputs$phenotype_metadata <- basename(metadata_path)
  yaml::write_yaml(raw, path)

  config <- read_phewasflow_config(path)
  spec <- config_to_phewas_spec(config)
  expect_identical(spec$phenotypes$phenotype, "001")
  expect_true("001" %in% names(data))

  manifest <- create_phewas_manifest(config)
  expect_identical(manifest$phenotype, "001")
  run_phewas_shard(config, 1L)
  result <- combine_phewas_shards(config)
  expect_identical(result$phenotype, "001")
  expect_identical(result$status, "ok")
})

test_that("TSV metadata receives binary phenotype defaults", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory, shards = 1L)
  metadata <- data.frame(
    phenotype = "pheno_binary",
    description = "Binary phenotype",
    group = "Diagnoses",
    groupnum = 1L,
    color = "#3B82F6",
    outcome_type = "binary",
    stringsAsFactors = FALSE
  )
  metadata_path <- file.path(directory, "binary-metadata.tsv")
  data.table::fwrite(metadata, metadata_path, sep = "\t", na = "NA")

  raw <- yaml::read_yaml(path)
  raw$inputs$phenotype_metadata <- basename(metadata_path)
  yaml::write_yaml(raw, path)

  spec <- config_to_phewas_spec(read_phewasflow_config(path))
  expect_identical(spec$phenotypes$variable_type, "binary")
  expect_identical(spec$phenotypes$reference, "0")
})

test_that("the scheduler-neutral CLI validates a configuration", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  output <- utils::capture.output(
    status <- phewasflow_cli(c("validate", "--config", path))
  )

  expect_identical(status, 0L)
  output <- paste(output, collapse = "\n")
  expect_match(output, '"status": "valid"')
  expect_match(output, '"phenotype_count": 2')
})

test_that("plot group display is configurable and validated", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  raw <- yaml::read_yaml(path)
  raw$plot$group_display <- "legend"
  yaml::write_yaml(raw, path)

  config <- read_phewasflow_config(path)
  expect_identical(config$plot$group_display, "legend")

  raw$plot$group_display <- "labels"
  yaml::write_yaml(raw, path)
  expect_error(
    read_phewasflow_config(path),
    "`plot.group_display` must be one of: auto, x_axis, legend"
  )
})

test_that("plot output accepts a PNG, PDF, or stem and defaults to Bonferroni", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  raw <- yaml::read_yaml(path)
  raw$plot$output <- "figures/publication"
  raw$plot$significance <- NULL
  yaml::write_yaml(raw, path)

  config <- read_phewasflow_config(path)
  expect_identical(config$plot$significance, "bonferroni")
  expect_identical(basename(config$plot$output), "publication")
  output_paths <- phewasFlow:::.pf_plot_output_paths(config$plot$output)
  expect_identical(
    names(output_paths),
    c("manhattan_png", "manhattan_pdf", "volcano_png", "volcano_pdf")
  )
  expect_identical(
    unname(basename(output_paths)),
    c(
      "publication.png", "publication.pdf",
      "publication-volcano.png", "publication-volcano.pdf"
    )
  )

  expect_identical(
    unname(basename(phewasFlow:::.pf_plot_output_paths(
      file.path(directory, "named.PDF")
    ))),
    c(
      "named.png", "named.pdf", "named-volcano.png", "named-volcano.pdf"
    )
  )

  raw$plot$output <- "figures/publication.svg"
  yaml::write_yaml(raw, path)
  expect_error(
    read_phewasflow_config(path, check_data = FALSE),
    "`plot.output` must be a filename stem or end in `.png` or `.pdf`"
  )
})

test_that("plot CLI writes both plot types as PNG and PDF", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  raw <- yaml::read_yaml(path)
  raw$plot$group_display <- "legend"
  yaml::write_yaml(raw, path)
  seen_groups <- character()
  seen_significance <- character()
  saved <- data.frame(
    plot = character(), path = character(), stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    read_phewas_bundle = function(...) structure(list(), class = "phewas_result"),
    plot_phewas_manhattan = function(..., significance, group_display) {
      seen_groups <<- c(seen_groups, group_display)
      seen_significance <<- c(
        seen_significance, paste0("manhattan:", significance)
      )
      structure(list(kind = "manhattan"), class = c("ggplot", "list"))
    },
    plot_phewas_volcano = function(..., significance) {
      seen_significance <<- c(
        seen_significance, paste0("volcano:", significance)
      )
      structure(list(kind = "volcano"), class = c("ggplot", "list"))
    },
    save_phewas_plot = function(plot, filename, ...) {
      saved <<- rbind(
        saved,
        data.frame(plot = plot$kind, path = filename, stringsAsFactors = FALSE)
      )
      invisible(filename)
    },
    .package = "phewasFlow"
  )

  output <- utils::capture.output(
    status <- phewasflow_cli(c("plot", "--config", path))
  )
  expect_identical(status, 0L)
  response <- jsonlite::fromJSON(paste(output, collapse = "\n"))
  configured_paths <- phewasFlow:::.pf_plot_output_paths(
    read_phewasflow_config(path)$plot$output
  )
  expect_identical(response$plot, unname(configured_paths[["manhattan_png"]]))
  expect_identical(
    unlist(response$plots, use.names = TRUE),
    configured_paths
  )
  expect_identical(saved$path, unname(configured_paths))
  expect_identical(
    saved$plot,
    c("manhattan", "manhattan", "volcano", "volcano")
  )

  override <- file.path(directory, "figures", "custom.pdf")
  invisible(utils::capture.output(
    phewasflow_cli(c(
      "plot", "--config", path, "--output", override,
      "--group-display", "x_axis"
    ))
  ))
  override_paths <- phewasFlow:::.pf_plot_output_paths(override)
  expect_identical(tail(saved$path, 4L), unname(override_paths))
  expect_identical(seen_groups, c("legend", "x_axis"))
  expect_identical(
    seen_significance,
    c("manhattan:bh", "volcano:bh", "manhattan:bh", "volcano:bh")
  )
  expect_error(
    phewasflow_cli(c(
      "plot", "--config", path, "--group-display", "labels"
    )),
    "`--group-display` must be auto, x_axis, or legend"
  )
  expect_error(
    phewasflow_cli(c(
      "plot", "--config", path, "--output", "publication.svg"
    )),
    "`--output` must be a filename stem or end in `.png` or `.pdf`"
  )
})

test_that("configuration errors identify missing files and explicit thresholds", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  raw <- yaml::read_yaml(path)
  raw$inputs$data <- "does-not-exist.rds"
  yaml::write_yaml(raw, path)
  expect_error(
    read_phewasflow_config(path),
    "Configured file does not exist"
  )

  path <- make_shard_fixture(directory)
  raw <- yaml::read_yaml(path)
  raw$analysis$eligibility$min_outcome_levels <- NULL
  yaml::write_yaml(raw, path)
  expect_error(read_phewasflow_config(path), "Missing eligibility field")

  path <- make_shard_fixture(directory)
  raw <- yaml::read_yaml(path)
  raw$analysis$covaraite <- "age"
  yaml::write_yaml(raw, path)
  expect_error(read_phewasflow_config(path), "Unknown `analysis` field")
})

test_that("configuration carries an explicit pre-analysis exclusion", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  data <- readRDS(file.path(directory, "analysis.rds"))
  data$X999 <- c(1L, rep(0L, nrow(data) - 1L))
  saveRDS(data, file.path(directory, "analysis.rds"))

  raw <- yaml::read_yaml(path)
  raw$analysis$exclusion_column <- "X999"
  raw$analysis$exclusion_values <- 1L
  yaml::write_yaml(raw, path)
  config <- read_phewasflow_config(path)
  spec <- config_to_phewas_spec(config)
  expect_identical(spec$exclusion_column, "X999")
  expect_equal(spec$exclusion_values, 1)

  raw$analysis$exclusion_column <- NULL
  yaml::write_yaml(raw, path)
  expect_error(read_phewasflow_config(path), "only allowed with")

  raw$analysis$exclusion_column <- "X999"
  raw$analysis$exclusion_values <- NULL
  yaml::write_yaml(raw, path)
  expect_error(read_phewasflow_config(path), "must be a nonempty vector")
})

test_that("configuration rejects nonfinite and out-of-range numeric settings", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  raw <- yaml::read_yaml(path)
  for (value in c(Inf, NaN, .Machine$integer.max + 1, 1.5)) {
    invalid <- raw
    invalid$execution$shards <- value
    expect_error(
      validate_phewasflow_config(invalid, config_path = path),
      "`execution.shards` must be a positive integer"
    )
  }
  for (field in c("width", "height", "dpi")) {
    invalid <- raw
    invalid$plot[[field]] <- Inf
    expect_error(
      validate_phewasflow_config(invalid, config_path = path),
      paste0("`plot.", field, "` must be a finite positive number")
    )
  }
})

test_that("version 2 manifests are deterministic by phenotype and detect stale inputs", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  first_path <- file.path(directory, "first.tsv")
  second_path <- file.path(directory, "second.tsv")
  first <- create_phewas_manifest(path, first_path)
  second <- create_phewas_manifest(path, second_path)

  columns <- names(first)
  expect_identical(columns, c(
    "manifest_version", "run_fingerprint", "analysis_fingerprint",
    "inputs_fingerprint", "package_version", "shard_id", "target_order",
    "phenotype", "shard_fingerprint"
  ))
  expect_identical(unclass(first[columns]), unclass(second[columns]))
  expect_identical(unique(first$manifest_version), "2")
  expect_identical(first$phenotype, sort(first$phenotype))
  expect_false("phenotype_id" %in% names(first))
  expect_identical(first$shard_id, c(1L, 2L))
  expect_false(file.exists(paste0(first_path, ".sha256")))

  data <- readRDS(file.path(directory, "analysis.rds"))
  data$pgs[[1L]] <- data$pgs[[1L]] + 0.01
  saveRDS(data, file.path(directory, "analysis.rds"))
  expect_error(
    run_phewas_shard(path, 1L, manifest_path = first_path),
    "does not match the current configuration or input files"
  )
})

test_that("stale version 1 manifests and shard artifacts are rejected clearly", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  manifest <- create_phewas_manifest(path)
  manifest_path <- attr(manifest, "manifest_path", exact = TRUE)

  stale_manifest <- data.table::fread(
    manifest_path, sep = "\t", data.table = FALSE
  )
  stale_manifest$manifest_version <- "1"
  names(stale_manifest)[names(stale_manifest) == "phenotype"] <- "phenotype_id"
  data.table::fwrite(stale_manifest, manifest_path, sep = "\t")
  expect_error(
    run_phewas_shard(path, 1L),
    "Manifest version 1 is stale"
  )

  create_phewas_manifest(path)
  shard_path <- run_phewas_shard(path, 1L)
  stale_artifact <- readRDS(shard_path)
  stale_artifact$artifact_version <- "1"
  stale_artifact$phenotype_ids <- stale_artifact$phenotypes
  stale_artifact$phenotypes <- NULL
  saveRDS(stale_artifact, shard_path, version = 3L)
  expect_error(
    run_phewas_shard(path, 1L),
    "Shard artifact version 1 is stale"
  )
})

test_that("shards restart safely and combine with exact coverage", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  data <- readRDS(file.path(directory, "analysis.rds"))
  data$X999 <- c(rep(1L, 3L), rep(0L, nrow(data) - 3L))
  saveRDS(data, file.path(directory, "analysis.rds"))
  raw <- yaml::read_yaml(path)
  raw$analysis$exclusion_column <- "X999"
  raw$analysis$exclusion_values <- 1L
  yaml::write_yaml(raw, path)
  manifest <- create_phewas_manifest(path)
  first <- run_phewas_shard(path, 1L)

  expect_error(combine_phewas_shards(path), "Shard set is not exact")
  expect_identical(run_phewas_shard(path, 1L), first)
  run_phewas_shard(path, 2L)
  result <- combine_phewas_shards(path)

  expect_s3_class(result, "phewas_result")
  expect_identical(as.character(result$phenotype), manifest$phenotype)
  expect_true(all(c(
    "effective_formula", "effective_covariates",
    "dropped_invariant_covariates"
  ) %in% names(result)))
  expect_true(all(result$status == "ok"))
  expect_true(all(is.finite(result$q_value)))
  direct <- run_phewas(
    readRDS(file.path(directory, "analysis.rds")),
    config_to_phewas_spec(path),
    backend = "sequential"
  )
  direct <- as.data.frame(direct)
  direct <- direct[match(result$phenotype, direct$phenotype), , drop = FALSE]
  expect_equal(result$estimate, direct$estimate, tolerance = 1e-12)
  expect_equal(result$log_p, direct$log_p, tolerance = 1e-12)
  expect_equal(result$q_value, direct$q_value, tolerance = 1e-12)
  expect_identical(unique(result$testing_family_size), 2L)
  run_metadata <- attr(result, "run_metadata", exact = TRUE)
  expect_true(run_metadata$exclusion_active)
  expect_identical(run_metadata$exclusion_column, "X999")
  expect_equal(run_metadata$exclusion_values, 1)
  expect_identical(run_metadata$n_before_exclusion, 120L)
  expect_identical(run_metadata$n_excluded, 3L)
  expect_identical(run_metadata$n_after_exclusion, 117L)
  output <- file.path(directory, "output")
  expect_setequal(list.files(output), c("results.rds", "results.tsv"))
  state <- file.path(output, ".phewasflow")
  expect_setequal(
    list.files(state, recursive = TRUE),
    c("manifest.tsv", "shards/shard-0001.rds", "shards/shard-0002.rds")
  )
  expect_equal(read_phewas_bundle(output), result, ignore_attr = TRUE)
  expect_s3_class(read_phewas_bundle(output), "data.table")
  stored <- readRDS(file.path(output, "results.rds"))
  expect_setequal(
    names(attributes(stored)),
    c("names", "row.names", "class", "spec", "run_metadata")
  )
  artifact <- readRDS(first)
  expect_identical(artifact$artifact_version, "2")
  expect_identical(artifact$phenotypes, manifest$phenotype[manifest$shard_id == 1L])
  expect_false("phenotype_ids" %in% names(artifact))
  expect_false(any(c("data", "fit", "fits", "models") %in% names(artifact)))

  # An unreadable shard is rejected.
  writeBin(charToRaw("corrupt"), first)
  expect_error(
    combine_phewas_shards(path, overwrite = TRUE),
    "Could not read shard artifact"
  )
})

test_that("result bundles contain only the RDS and TSV outputs", {
  directory <- withr::local_tempdir()
  result <- make_bundle_result()
  attr(result, "bundle_directory") <- file.path(directory, "private-inputs")
  attr(result, "participant_data") <- data.frame(id = "synthetic-participant")
  bundle <- file.path(directory, "bundle")
  write_phewas_bundle(result, bundle, run_metadata = list(run_fingerprint = "test"))
  restored <- read_phewas_bundle(bundle)

  expect_setequal(list.files(bundle), c("results.rds", "results.tsv"))
  stored <- readRDS(file.path(bundle, "results.rds"))
  expect_setequal(
    names(attributes(stored)),
    c("names", "row.names", "class", "spec", "run_metadata")
  )
  expect_null(attr(restored, "bundle_directory", exact = TRUE))
  expect_null(attr(restored, "participant_data", exact = TRUE))
  expect_s3_class(restored, "phewas_result")
  tsv <- data.table::fread(
    file.path(bundle, "results.tsv"), sep = "\t", data.table = FALSE
  )
  expect_identical(names(tsv), names(as.data.frame(restored)))
  expect_true(all(c(
    "phenotype", "description", "group", "groupnum", "color"
  ) %in% names(tsv)))
  expect_false(any(c(
    "phenotype_id", "phecode", "phenotype_column", "column", "label", "category",
    "category_order"
  ) %in% names(tsv)))
  expect_identical(
    attr(restored, "run_metadata", exact = TRUE)$run_fingerprint,
    "test"
  )
  expect_error(write_phewas_bundle(result, bundle), "Results already exist")
})

test_that("bundle writes require every canonical phenotype annotation", {
  directory <- withr::local_tempdir()
  result <- make_bundle_result()

  for (field in c("phenotype", "description", "group", "groupnum", "color")) {
    incomplete <- data.table::copy(result)
    incomplete[[field]] <- NULL
    expect_error(
      write_phewas_bundle(
        incomplete, file.path(directory, paste0("missing-", field)),
        run_metadata = list(run_fingerprint = paste0("missing-", field))
      ),
      paste0("missing required column.*", field)
    )
  }

  legacy_metadata <- attr(result, "run_metadata", exact = TRUE)
  legacy_metadata$phenotype_ids <- legacy_metadata$phenotypes
  legacy_metadata$phenotypes <- NULL
  expect_error(
    write_phewas_bundle(
      result, file.path(directory, "legacy-run-metadata"),
      run_metadata = legacy_metadata
    ),
    "legacy run metadata field `phenotype_ids`"
  )

  legacy_spec_result <- data.table::copy(result)
  legacy_spec <- attr(legacy_spec_result, "spec", exact = TRUE)
  legacy_spec$phenotypes$column <- legacy_spec$phenotypes$phenotype
  attr(legacy_spec_result, "spec") <- legacy_spec
  expect_error(
    write_phewas_bundle(
      legacy_spec_result, file.path(directory, "legacy-spec-metadata")
    ),
    "noncanonical specification metadata"
  )
})

test_that("legacy bundles are normalized on read only when aliases agree", {
  directory <- withr::local_tempdir()
  result <- make_bundle_result()
  specification <- attr(result, "spec", exact = TRUE)
  run_metadata <- attr(result, "run_metadata", exact = TRUE)
  run_metadata$run_fingerprint <- "legacy-result"
  run_metadata$phenotype_ids <- run_metadata$phenotypes
  run_metadata$phenotypes <- NULL
  legacy <- as.data.frame(data.table::copy(result), stringsAsFactors = FALSE)
  names(legacy)[names(legacy) == "phenotype"] <- "phenotype_id"
  legacy$phecode <- legacy$phenotype_id
  legacy$phenotype_column <- legacy$phenotype_id
  legacy$label <- legacy$description
  legacy$category <- legacy$group
  legacy$category_order <- legacy$groupnum
  class(legacy) <- c("phewas_result", "data.frame")
  attr(legacy, "spec") <- specification
  attr(legacy, "run_metadata") <- run_metadata

  expect_error(
    write_phewas_bundle(legacy, file.path(directory, "new-write")),
    "missing required column.*phenotype"
  )
  canonical_with_alias <- data.table::copy(result)
  canonical_with_alias$phenotype_id <- canonical_with_alias$phenotype
  expect_error(
    write_phewas_bundle(
      canonical_with_alias, file.path(directory, "aliased-write"),
      run_metadata = list(run_fingerprint = "aliased-write")
    ),
    "legacy result column"
  )

  legacy_directory <- file.path(directory, "legacy")
  dir.create(legacy_directory)
  saveRDS(legacy, file.path(legacy_directory, "results.rds"), version = 3L)
  restored <- read_phewas_bundle(legacy_directory)
  expect_equal(
    as.data.frame(restored), as.data.frame(result),
    ignore_attr = TRUE
  )
  expect_false(any(c(
    "phenotype_id", "phecode", "phenotype_column", "column", "label", "category",
    "category_order"
  ) %in% names(restored)))
  expect_identical(
    attr(restored, "run_metadata", exact = TRUE)$phenotypes,
    result$phenotype
  )
  expect_false(
    "phenotype_ids" %in% names(attr(restored, "run_metadata", exact = TRUE))
  )

  conflicting_directory <- file.path(directory, "legacy-conflict")
  dir.create(conflicting_directory)
  conflicting <- legacy
  conflicting$phenotype_column[[1L]] <- paste0(
    conflicting$phenotype_column[[1L]], "_other"
  )
  saveRDS(
    conflicting, file.path(conflicting_directory, "results.rds"), version = 3L
  )
  expect_error(
    read_phewas_bundle(conflicting_directory),
    "cannot be safely normalized.*phenotype"
  )
})

test_that("canonical result rows also migrate legacy attributes", {
  directory <- withr::local_tempdir()
  result <- make_bundle_result()

  legacy_spec <- attr(result, "spec", exact = TRUE)
  metadata <- legacy_spec$phenotypes
  names(metadata)[names(metadata) == "phenotype"] <- "phenotype_id"
  metadata$column <- metadata$phenotype_id
  names(metadata)[names(metadata) == "description"] <- "label"
  names(metadata)[names(metadata) == "group"] <- "category"
  names(metadata)[names(metadata) == "groupnum"] <- "category_order"
  legacy_spec$phenotypes <- metadata

  legacy_run_metadata <- attr(result, "run_metadata", exact = TRUE)
  legacy_run_metadata$run_fingerprint <- "legacy-attributes"
  legacy_run_metadata$phenotype_ids <- legacy_run_metadata$phenotypes
  legacy_run_metadata$phenotypes <- NULL

  portable <- as.data.frame(data.table::copy(result), stringsAsFactors = FALSE)
  class(portable) <- c("phewas_result", "data.frame")
  attr(portable, "spec") <- legacy_spec
  attr(portable, "run_metadata") <- legacy_run_metadata
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  saveRDS(portable, file.path(directory, "results.rds"), version = 3L)

  restored <- read_phewas_bundle(directory)
  restored_spec <- attr(restored, "spec", exact = TRUE)
  restored_metadata <- attr(restored, "run_metadata", exact = TRUE)
  expect_identical(restored_spec$phenotypes$phenotype, result$phenotype)
  expect_false(any(c(
    "phenotype_id", "column", "label", "category", "category_order"
  ) %in% names(restored_spec$phenotypes)))
  expect_identical(restored_metadata$phenotypes, result$phenotype)
  expect_false("phenotype_ids" %in% names(restored_metadata))
})

test_that("combined results are validated and the TSV is repaired on reuse", {
  directory <- withr::local_tempdir()
  fixture <- complete_shard_fixture(directory)
  rds_path <- file.path(fixture$output, "results.rds")
  tsv_path <- file.path(fixture$output, "results.tsv")
  expected_tsv <- readLines(tsv_path, warn = FALSE)

  unlink(tsv_path)
  reused <- suppressMessages(combine_phewas_shards(fixture$config))
  expect_equal(reused, fixture$result, ignore_attr = TRUE)
  expect_identical(readLines(tsv_path, warn = FALSE), expected_tsv)

  writeLines("stale table", tsv_path)
  suppressMessages(combine_phewas_shards(fixture$config))
  expect_identical(readLines(tsv_path, warn = FALSE), expected_tsv)

  tampered <- readRDS(rds_path)
  tampered$estimate[[1L]] <- tampered$estimate[[1L]] + 1
  saveRDS(tampered, rds_path, version = 3L)
  expect_error(
    combine_phewas_shards(fixture$config),
    "do not match the validated shard set"
  )
  restored <- combine_phewas_shards(fixture$config, overwrite = TRUE)
  expect_equal(restored, fixture$result, ignore_attr = TRUE)
  expect_identical(readLines(tsv_path, warn = FALSE), expected_tsv)
})

test_that("readable shard artifacts cannot hide schema or metadata drift", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  create_phewas_manifest(path)
  first <- run_phewas_shard(path, 1L)
  second <- run_phewas_shard(path, 2L)
  original_first <- readRDS(first)
  original_second <- readRDS(second)

  altered <- original_first
  altered$results$log_p <- NULL
  saveRDS(altered, first, version = 3L)
  expect_error(combine_phewas_shards(path), "missing required column")

  saveRDS(original_first, first, version = 3L)
  altered <- original_first
  altered$results$color <- NULL
  saveRDS(altered, first, version = 3L)
  expect_error(combine_phewas_shards(path), "missing required column.*color")

  saveRDS(original_first, first, version = 3L)
  altered <- original_first
  altered$results$status[[1L]] <- "mystery"
  saveRDS(altered, first, version = 3L)
  expect_error(combine_phewas_shards(path), "invalid association status")

  saveRDS(original_first, first, version = 3L)
  altered <- original_first
  altered$results$analysis_id[[1L]] <- "other_analysis"
  saveRDS(altered, first, version = 3L)
  expect_error(combine_phewas_shards(path), "inconsistent specification")

  saveRDS(original_first, first, version = 3L)
  altered <- original_second
  altered$run_metadata$n_base <- altered$run_metadata$n_base + 1L
  saveRDS(altered, second, version = 3L)
  expect_error(combine_phewas_shards(path), "run metadata are inconsistent")
})

test_that("unexpected shard IDs beyond four digits are rejected", {
  directory <- withr::local_tempdir()
  fixture <- complete_shard_fixture(directory)
  extra <- file.path(dirname(fixture$shards[[1L]]), "shard-10000.rds")
  expect_true(file.copy(fixture$shards[[1L]], extra))
  expect_error(
    combine_phewas_shards(fixture$config),
    "unexpected: shard-10000.rds"
  )
})

test_that("configured and overridden output paths cannot overwrite inputs", {
  directory <- withr::local_tempdir()
  path <- make_shard_fixture(directory)
  input_path <- file.path(directory, "analysis.rds")
  input_hash <- phewasFlow:::.pf_file_sha256(input_path)

  expect_error(
    create_phewas_manifest(path, path = input_path),
    "collides with a protected analysis path"
  )
  expect_identical(phewasFlow:::.pf_file_sha256(input_path), input_hash)

  raw <- yaml::read_yaml(path)
  raw$plot$output <- "analysis.rds"
  yaml::write_yaml(raw, path)
  expect_error(
    read_phewasflow_config(path, check_data = FALSE),
    "`plot.output` must be a filename stem or end in `.png` or `.pdf`"
  )
  expect_identical(phewasFlow:::.pf_file_sha256(input_path), input_hash)

  # Every derived companion path is protected, not only the configured name.
  path <- make_shard_fixture(directory)
  input_pdf <- file.path(directory, "publication.pdf")
  expect_true(file.copy(input_path, input_pdf))
  raw <- yaml::read_yaml(path)
  raw$inputs$data <- basename(input_pdf)
  raw$plot$output <- "publication.png"
  yaml::write_yaml(raw, path)
  expect_error(
    read_phewasflow_config(path, check_data = FALSE),
    "would overwrite an input.*publication.pdf"
  )

  fixture <- complete_shard_fixture(file.path(directory, "complete"))
  result_rds <- file.path(fixture$output, "results.rds")
  result_hash <- phewasFlow:::.pf_file_sha256(result_rds)
  expect_error(
    phewasflow_cli(c(
      "plot", "--config", fixture$config, "--output", result_rds
    )),
    "`--output` must be a filename stem or end in `.png` or `.pdf`"
  )
  expect_identical(phewasFlow:::.pf_file_sha256(result_rds), result_hash)

  config_as_pdf <- file.path(dirname(fixture$config), "publication.pdf")
  expect_true(file.copy(fixture$config, config_as_pdf))
  config_hash <- phewasFlow:::.pf_file_sha256(config_as_pdf)
  expect_error(
    phewasflow_cli(c(
      "plot", "--config", config_as_pdf, "--output",
      file.path(dirname(config_as_pdf), "publication.png")
    )),
    "collides with a protected analysis path.*publication.pdf"
  )
  expect_identical(
    phewasFlow:::.pf_file_sha256(config_as_pdf), config_hash
  )
  expect_error(
    phewasflow_cli(c(
      "plot", "--config", fixture$config, "--output",
      file.path(fixture$output, ".phewasflow", "publication")
    )),
    "cannot be inside the internal `.phewasflow` directory"
  )
  expect_error(
    phewasflow_cli(c(
      "plot", "--config", fixture$config, "--output",
      file.path(fixture$output, "phewas-volcano.png")
    )),
    "collides with a protected analysis path.*phewas-volcano.png"
  )
  expect_error(
    combine_phewas_shards(
      fixture$config,
      bundle_directory = file.path(dirname(fixture$config), "analysis.rds"),
      overwrite = TRUE
    ),
    "directory collides with an input"
  )
})

test_that("failed replacement leaves the previous result pair unchanged", {
  directory <- withr::local_tempdir()
  bundle <- file.path(directory, "bundle")
  result <- make_bundle_result()
  write_phewas_bundle(result, bundle)
  rds_path <- file.path(bundle, "results.rds")
  tsv_path <- file.path(bundle, "results.tsv")
  old_rds <- readBin(rds_path, what = "raw", n = file.info(rds_path)$size)
  old_tsv <- readBin(tsv_path, what = "raw", n = file.info(tsv_path)$size)

  changed <- data.table::copy(result)
  changed$estimate[[1L]] <- changed$estimate[[1L]] + 1
  attr(changed, "spec") <- attr(result, "spec", exact = TRUE)
  attr(changed, "run_metadata") <- attr(result, "run_metadata", exact = TRUE)
  testthat::local_mocked_bindings(
    .pf_write_rds = function(...) stop("deliberate staged RDS failure"),
    .package = "phewasFlow"
  )
  expect_error(
    write_phewas_bundle(changed, bundle, overwrite = TRUE),
    "deliberate staged RDS failure"
  )
  expect_identical(
    readBin(rds_path, what = "raw", n = file.info(rds_path)$size), old_rds
  )
  expect_identical(
    readBin(tsv_path, what = "raw", n = file.info(tsv_path)$size), old_tsv
  )
  expect_setequal(list.files(bundle), c("results.rds", "results.tsv"))
})

test_that("writing a bundle protects a preexisting TSV and directories", {
  directory <- withr::local_tempdir()
  result <- make_bundle_result()
  tsv_path <- file.path(directory, "results.tsv")
  writeLines("existing table", tsv_path)
  expect_error(write_phewas_bundle(result, directory), "Results already exist")
  expect_identical(readLines(tsv_path), "existing table")
  expect_false(file.exists(file.path(directory, "results.rds")))

  unlink(tsv_path)
  dir.create(tsv_path)
  retained_path <- file.path(tsv_path, "retain.txt")
  writeLines("keep this file", retained_path)
  expect_error(
    write_phewas_bundle(result, directory, overwrite = TRUE),
    "Cannot replace a directory"
  )
  expect_identical(readLines(retained_path), "keep this file")
  expect_false(file.exists(file.path(directory, "results.rds")))
})

test_that("result readers reject malformed objects", {
  directory <- withr::local_tempdir()
  dir.create(file.path(directory, "bad"))
  saveRDS(
    data.frame(phenotype = "p1", status = "ok"),
    file.path(directory, "bad", "results.rds")
  )
  expect_error(read_phewas_bundle(file.path(directory, "bad")), "phewas_result")
  expect_error(read_phewas_bundle(character()), "one non-empty path")
})

test_that("a failed atomic write leaves no published or temporary file", {
  directory <- withr::local_tempdir()
  path <- file.path(directory, "result.rds")
  expect_error(
    phewasFlow:::.pf_atomic_write(path, function(temporary) {
      saveRDS(list(value = 1), temporary)
      stop("deliberate writer failure")
    }),
    "deliberate writer failure"
  )
  expect_false(file.exists(path))
  expect_length(list.files(directory, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("a failed atomic replacement restores the previous file", {
  directory <- withr::local_tempdir()
  path <- file.path(directory, "result.rds")
  saveRDS(list(value = "old"), path)
  old_bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  rename_calls <- 0L
  testthat::local_mocked_bindings(
    .pf_file_rename = function(from, to) {
      rename_calls <<- rename_calls + 1L
      if (rename_calls %in% c(1L, 3L)) {
        return(FALSE)
      }
      base::file.rename(from, to)
    },
    .package = "phewasFlow"
  )
  expect_error(
    phewasFlow:::.pf_atomic_write(path, function(temporary) {
      saveRDS(list(value = "new"), temporary)
    }),
    "Could not publish file atomically"
  )
  expect_identical(
    readBin(path, what = "raw", n = file.info(path)$size), old_bytes
  )
  expect_setequal(list.files(directory), "result.rds")
})
