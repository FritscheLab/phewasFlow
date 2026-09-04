# Deterministic file-based sharding -----------------------------------------

.pf_as_config <- function(config) {
  if (is.character(config)) {
    return(read_phewasflow_config(config, check_data = FALSE))
  }
  if (!inherits(config, "phewasflow_config")) {
    return(validate_phewasflow_config(config, check_data = FALSE))
  }
  config
}

.pf_manifest_path <- function(config, path = NULL) {
  .pf_default(path, file.path(config$output$directory, ".phewasflow", "manifest.tsv"))
}

.pf_shard_directory <- function(config) {
  file.path(config$output$directory, ".phewasflow", "shards")
}

.pf_shard_path <- function(config, shard_id) {
  file.path(.pf_shard_directory(config), sprintf("shard-%04d.rds", as.integer(shard_id)))
}

.pf_manifest_data <- function(config) {
  phenotypes <- .pf_parse_metadata_lists(
    .pf_read_phenotype_metadata(config$inputs$phenotype_metadata)
  )
  spec <- validate_phewas_spec(.pf_spec_from_config(config, phenotypes))
  phenotypes <- sort(as.character(spec$phenotypes$phenotype), method = "radix")
  if (!length(phenotypes)) {
    stop("The specification contains no phenotypes.", call. = FALSE)
  }
  if (anyNA(phenotypes) || any(!nzchar(phenotypes)) || anyDuplicated(phenotypes)) {
    stop("Phenotypes must be unique, non-missing strings.", call. = FALSE)
  }
  n_shards <- config$execution$shards
  if (n_shards > length(phenotypes)) {
    stop(
      "`execution.shards` cannot exceed the number of phenotypes (",
      length(phenotypes), ").",
      call. = FALSE
    )
  }
  identity <- .pf_analysis_identity(config, spec, phenotypes)
  shard_id <- ((seq_along(phenotypes) - 1L) %% n_shards) + 1L
  shard_fingerprints <- vapply(seq_len(n_shards), function(index) {
    .pf_object_sha256(list(
      run_fingerprint = identity$run_fingerprint,
      shard_id = index,
      phenotypes = phenotypes[shard_id == index]
    ))
  }, character(1L))
  manifest <- data.frame(
    manifest_version = "2",
    run_fingerprint = identity$run_fingerprint,
    analysis_fingerprint = identity$analysis_fingerprint,
    inputs_fingerprint = identity$inputs_fingerprint,
    package_version = identity$package_version,
    shard_id = shard_id,
    target_order = seq_along(phenotypes),
    phenotype = phenotypes,
    shard_fingerprint = unname(shard_fingerprints[shard_id]),
    stringsAsFactors = FALSE
  )
  attr(manifest, "identity") <- identity
  attr(manifest, "spec") <- spec
  manifest
}

#' Create a deterministic shard manifest
#'
#' Phenotypes are sorted and assigned round-robin to the configured number
#' of shards. The manifest fingerprints the validated specification, input
#' files, package version, and exact target membership.
#'
#' @param config A configuration path or normalized config object.
#' @param path Output TSV path. Defaults to `.phewasflow/manifest.tsv` in the
#'   configured output directory.
#' @param write Whether to write the manifest.
#'
#' @return A manifest data frame.
#' @export
create_phewas_manifest <- function(config, path = NULL, write = TRUE) {
  config <- .pf_as_config(config)
  manifest <- .pf_manifest_data(config)
  path <- .pf_manifest_path(config, path)
  if (isTRUE(write)) {
    default_manifest <- .pf_manifest_path(config)
    .pf_assert_safe_write_target(
      config, path, purpose = "Manifest", allow = default_manifest
    )
    .pf_write_tsv(manifest, path)
  }
  attr(manifest, "manifest_path") <- path
  manifest
}

.pf_read_manifest <- function(config, path = NULL) {
  path <- .pf_manifest_path(config, path)
  if (!file.exists(path)) {
    stop("Manifest does not exist: ", path, call. = FALSE)
  }
  # Phenotype identifiers are text, even when they look like numbers or NA.
  manifest <- data.table::fread(
    path, sep = "\t", data.table = FALSE, colClasses = "character",
    na.strings = NULL
  )
  if (!"manifest_version" %in% names(manifest)) {
    stop("Manifest is missing column `manifest_version`.", call. = FALSE)
  }
  versions <- unique(as.character(manifest$manifest_version))
  if (identical(versions, "1")) {
    stop(
      "Manifest version 1 is stale and cannot be reused. Recreate the ",
      "manifest with this version of phewasFlow.",
      call. = FALSE
    )
  }
  if (!identical(versions, "2")) {
    stop(
      "Unsupported manifest version: ", paste(versions, collapse = ", "), ".",
      call. = FALSE
    )
  }
  required <- c(
    "manifest_version", "run_fingerprint", "analysis_fingerprint",
    "inputs_fingerprint", "package_version", "shard_id", "target_order",
    "phenotype", "shard_fingerprint"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("Manifest is missing column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (anyDuplicated(as.character(manifest$phenotype))) {
    stop("Manifest contains duplicate phenotypes.", call. = FALSE)
  }
  if (anyDuplicated(as.integer(manifest$target_order))) {
    stop("Manifest contains duplicate target-order values.", call. = FALSE)
  }

  expected <- .pf_manifest_data(config)
  identity <- attr(expected, "identity", exact = TRUE)
  expected_spec <- attr(expected, "spec", exact = TRUE)
  manifest <- manifest[order(as.integer(manifest$target_order)), required, drop = FALSE]
  expected <- expected[order(expected$target_order), required, drop = FALSE]
  rownames(manifest) <- NULL
  rownames(expected) <- NULL
  for (field in required) {
    observed <- if (field %in% c("shard_id", "target_order")) {
      as.integer(manifest[[field]])
    } else {
      as.character(manifest[[field]])
    }
    wanted <- if (field %in% c("shard_id", "target_order")) {
      as.integer(expected[[field]])
    } else {
      as.character(expected[[field]])
    }
    if (!identical(observed, wanted)) {
      stop(
        "Manifest does not match the current configuration or input files (field `",
        field, "`). Recreate it before running shards.",
        call. = FALSE
      )
    }
  }
  attr(manifest, "manifest_path") <- path
  attr(manifest, "identity") <- identity
  attr(manifest, "spec") <- expected_spec
  manifest
}

.pf_validate_shard_artifact <- function(artifact, manifest_rows, path,
                                        expected_spec = NULL) {
  if (!is.list(artifact)) {
    stop("Invalid shard artifact: ", path, call. = FALSE)
  }
  artifact_version <- as.character(artifact$artifact_version)
  if (identical(artifact_version, "1")) {
    stop(
      "Shard artifact version 1 is stale and cannot be reused: ", path,
      ". Re-run this shard with this version of phewasFlow.",
      call. = FALSE
    )
  }
  if (!identical(artifact_version, "2")) {
    stop("Invalid shard artifact version: ", path, call. = FALSE)
  }
  expected_id <- unique(as.integer(manifest_rows$shard_id))
  expected_targets <- as.character(manifest_rows$phenotype)
  expected_shard_fingerprint <- unique(as.character(manifest_rows$shard_fingerprint))
  expected_run_fingerprint <- unique(as.character(manifest_rows$run_fingerprint))
  if (length(expected_id) != 1L || length(expected_shard_fingerprint) != 1L ||
      length(expected_run_fingerprint) != 1L) {
    stop("Manifest contains mixed shard metadata.", call. = FALSE)
  }
  if (!identical(as.integer(artifact$shard_id), expected_id) ||
      !identical(as.character(artifact$run_fingerprint), expected_run_fingerprint) ||
      !identical(as.character(artifact$shard_fingerprint), expected_shard_fingerprint)) {
    stop("Shard artifact is stale or belongs to another analysis: ", path, call. = FALSE)
  }
  if (!identical(as.character(artifact$phenotypes), expected_targets)) {
    stop("Shard artifact target membership is invalid: ", path, call. = FALSE)
  }
  if (!is.list(artifact$run_metadata) || is.null(names(artifact$run_metadata))) {
    stop("Shard artifact has invalid run metadata: ", path, call. = FALSE)
  }
  tryCatch(
    .pf_validate_result_object(
      artifact$results,
      label = paste0("Shard result `", path, "`"),
      require_fingerprint = FALSE
    ),
    error = function(error) {
      stop(conditionMessage(error), call. = FALSE)
    }
  )
  if (!is.null(expected_spec)) {
    observed_spec <- attr(artifact$results, "spec", exact = TRUE)
    if (!identical(.pf_object_sha256(unclass(observed_spec)),
                   .pf_object_sha256(unclass(expected_spec)))) {
      stop("Shard result specification is inconsistent: ", path, call. = FALSE)
    }
    if (any(as.character(artifact$results$analysis_id) !=
            expected_spec$analysis_id)) {
      stop("Shard result analysis ID is inconsistent: ", path, call. = FALSE)
    }
  }
  observed_targets <- as.character(artifact$results$phenotype)
  if (anyDuplicated(observed_targets) || !setequal(observed_targets, expected_targets)) {
    stop("Shard result coverage is not exact: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

.pf_read_shard_artifact <- function(path, manifest_rows, expected_spec = NULL) {
  artifact <- tryCatch(
    readRDS(path),
    error = function(error) {
      stop("Could not read shard artifact `", path, "`: ", conditionMessage(error), call. = FALSE)
    }
  )
  .pf_validate_shard_artifact(
    artifact, manifest_rows, path, expected_spec = expected_spec
  )
  artifact
}

.pf_validate_shard_set <- function(artifacts) {
  if (!length(artifacts)) {
    stop("No shard artifacts were supplied.", call. = FALSE)
  }
  reference_names <- names(artifacts[[1L]]$results)
  reference_types <- vapply(artifacts[[1L]]$results, typeof, character(1L))
  for (index in seq_along(artifacts)[-1L]) {
    observed_names <- names(artifacts[[index]]$results)
    observed_types <- vapply(artifacts[[index]]$results, typeof, character(1L))
    if (!identical(observed_names, reference_names) ||
        !identical(observed_types, reference_types)) {
      stop("Shard result schemas are inconsistent.", call. = FALSE)
    }
  }

  consistent_fields <- c(
    "backend", "workers", "n_base", "exclusion_active", "exclusion_column",
    "exclusion_values", "n_before_exclusion", "n_excluded",
    "n_after_exclusion", "transformation", "transformation_mean",
    "transformation_sd", "corrections_applied"
  )
  metadata_signature <- function(artifact) {
    metadata <- artifact$run_metadata
    values <- lapply(consistent_fields, function(field) metadata[[field]])
    names(values) <- consistent_fields
    .pf_object_sha256(values)
  }
  signatures <- vapply(artifacts, metadata_signature, character(1L))
  if (length(unique(signatures)) != 1L) {
    stop("Shard run metadata are inconsistent.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Run one manifest shard
#'
#' Results are published as an atomic RDS artifact. A valid existing artifact
#' is reused unless `overwrite` is true.
#'
#' @param config A configuration path or normalized config object.
#' @param shard_id Positive shard number, commonly `SLURM_ARRAY_TASK_ID`.
#' @param manifest_path Optional manifest TSV path.
#' @param overwrite Refit and replace a valid existing artifact.
#'
#' @return The shard artifact path, invisibly.
#' @export
run_phewas_shard <- function(config, shard_id, manifest_path = NULL, overwrite = FALSE) {
  config <- .pf_as_config(config)
  shard_id <- .pf_positive_integer(as.numeric(shard_id), "shard_id")
  manifest <- .pf_read_manifest(config, manifest_path)
  manifest_rows <- manifest[as.integer(manifest$shard_id) == shard_id, , drop = FALSE]
  if (!nrow(manifest_rows)) {
    stop("Shard ", shard_id, " is not present in the manifest.", call. = FALSE)
  }
  manifest_rows <- manifest_rows[order(as.integer(manifest_rows$target_order)), , drop = FALSE]
  artifact_path <- .pf_shard_path(config, shard_id)
  if (file.exists(artifact_path) && !isTRUE(overwrite)) {
    .pf_read_shard_artifact(
      artifact_path, manifest_rows,
      expected_spec = attr(manifest, "spec", exact = TRUE)
    )
    message("Reusing valid shard artifact: ", artifact_path)
    return(invisible(normalizePath(artifact_path, winslash = "/", mustWork = TRUE)))
  }

  loaded <- .pf_load_config_inputs(config)
  spec <- validate_phewas_spec(.pf_spec_from_config(config, loaded$phenotypes), data = loaded$data)
  phenotypes <- as.character(manifest_rows$phenotype)
  results <- run_phewas(
    data = loaded$data,
    spec = spec,
    backend = config$execution$backend,
    workers = config$execution$workers,
    phenotype_ids = phenotypes,
    adjust = FALSE
  )
  if (!is.data.frame(results) || !"phenotype" %in% names(results)) {
    stop("`run_phewas()` did not return valid result rows.", call. = FALSE)
  }
  order_index <- match(phenotypes, as.character(results$phenotype))
  if (anyNA(order_index) || anyDuplicated(as.character(results$phenotype)) ||
      nrow(results) != length(phenotypes)) {
    stop("`run_phewas()` did not return exactly one row per shard target.", call. = FALSE)
  }
  result_spec <- attr(results, "spec", exact = TRUE)
  result_run_metadata <- attr(results, "run_metadata", exact = TRUE)
  results <- as.data.frame(data.table::copy(results), stringsAsFactors = FALSE)
  results <- results[order_index, , drop = FALSE]
  result_run_metadata$phenotypes <- as.character(results$phenotype)
  results <- data.table::as.data.table(results)
  class(results) <- unique(c("phewas_result", class(results)))
  attr(results, "spec") <- result_spec
  attr(results, "run_metadata") <- result_run_metadata

  identity <- attr(manifest, "identity", exact = TRUE)
  artifact <- list(
    artifact_version = "2",
    shard_id = shard_id,
    run_fingerprint = unique(as.character(manifest_rows$run_fingerprint)),
    analysis_fingerprint = unique(as.character(manifest_rows$analysis_fingerprint)),
    shard_fingerprint = unique(as.character(manifest_rows$shard_fingerprint)),
    package_version = unique(as.character(manifest_rows$package_version)),
    input_hashes = identity$input_hashes,
    phenotypes = phenotypes,
    created_at = .pf_timestamp(),
    run_metadata = attr(results, "run_metadata", exact = TRUE),
    results = results
  )
  .pf_write_rds(artifact, artifact_path)
  invisible(normalizePath(artifact_path, winslash = "/", mustWork = TRUE))
}

.pf_results_match_shards <- function(existing, combined, spec, run_metadata) {
  existing_metadata <- attr(existing, "run_metadata", exact = TRUE)
  expected_metadata <- run_metadata
  existing_metadata$created_at <- NULL
  expected_metadata$created_at <- NULL
  same_columns <- identical(names(existing), names(combined)) &&
    nrow(existing) == nrow(combined) &&
    all(vapply(names(existing), function(column) {
      identical(existing[[column]], combined[[column]])
    }, logical(1L)))
  same_columns && identical(
    .pf_object_sha256(unclass(attr(existing, "spec", exact = TRUE))),
    .pf_object_sha256(unclass(spec))
  ) && identical(
    .pf_object_sha256(existing_metadata),
    .pf_object_sha256(expected_metadata)
  )
}

#' Combine and globally correct completed PheWAS shards
#'
#' Every expected shard must match the manifest and provide its assigned
#' phenotype results. Multiple-testing correction is then applied once across
#' all successful associations.
#'
#' @param config A configuration path or normalized config object.
#' @param manifest_path Optional manifest TSV path.
#' @param bundle_directory Destination result directory. Defaults to the
#'   configured output directory.
#' @param overwrite Replace any existing results in the destination.
#'
#' @return The combined `phewas_result` object.
#' @export
combine_phewas_shards <- function(config, manifest_path = NULL,
                                  bundle_directory = NULL, overwrite = FALSE) {
  config <- .pf_as_config(config)
  manifest <- .pf_read_manifest(config, manifest_path)
  bundle_directory <- .pf_default(
    bundle_directory,
    config$output$directory
  )
  bundle_directory <- .pf_scalar_character(
    bundle_directory, "bundle_directory"
  )
  result_paths <- file.path(
    bundle_directory, c("results.rds", "results.tsv")
  )
  default_result_paths <- file.path(
    config$output$directory, c("results.rds", "results.tsv")
  )
  .pf_assert_safe_write_target(
    config,
    result_paths,
    purpose = "Result",
    allow = default_result_paths,
    directories = bundle_directory
  )
  run_fingerprint <- unique(as.character(manifest$run_fingerprint))
  if (length(run_fingerprint) != 1L) {
    stop("Manifest contains mixed analysis fingerprints.", call. = FALSE)
  }
  shard_ids <- sort(unique(as.integer(manifest$shard_id)))
  expected_names <- sprintf("shard-%04d.rds", shard_ids)
  shard_directory <- .pf_shard_directory(config)
  actual_names <- if (dir.exists(shard_directory)) {
    list.files(shard_directory, pattern = "^shard-[0-9]+\\.rds$")
  } else {
    character()
  }
  missing <- setdiff(expected_names, actual_names)
  unexpected <- setdiff(actual_names, expected_names)
  if (length(missing) || length(unexpected)) {
    details <- c(
      if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(unexpected)) paste0("unexpected: ", paste(unexpected, collapse = ", "))
    )
    stop("Shard set is not exact (", paste(details, collapse = "; "), ").", call. = FALSE)
  }

  artifacts <- lapply(shard_ids, function(index) {
    rows <- manifest[as.integer(manifest$shard_id) == index, , drop = FALSE]
    rows <- rows[order(as.integer(rows$target_order)), , drop = FALSE]
    .pf_read_shard_artifact(
      .pf_shard_path(config, index), rows,
      expected_spec = attr(manifest, "spec", exact = TRUE)
    )
  })
  artifact_runs <- unique(vapply(artifacts, function(x) as.character(x$run_fingerprint), character(1L)))
  if (!identical(artifact_runs, run_fingerprint)) {
    stop("Shard artifacts contain mixed analysis specifications.", call. = FALSE)
  }
  .pf_validate_shard_set(artifacts)

  combined <- data.table::rbindlist(
    lapply(artifacts, `[[`, "results"),
    use.names = TRUE,
    fill = FALSE
  )
  if (anyDuplicated(as.character(combined$phenotype))) {
    stop("Combined shards contain duplicate phenotypes.", call. = FALSE)
  }
  expected_targets <- as.character(manifest$phenotype[order(as.integer(manifest$target_order))])
  position <- match(expected_targets, as.character(combined$phenotype))
  if (anyNA(position) || nrow(combined) != length(expected_targets)) {
    stop("Combined shards do not provide exact target coverage.", call. = FALSE)
  }
  combined <- as.data.frame(combined, stringsAsFactors = FALSE)
  combined <- combined[position, , drop = FALSE]
  combined <- .apply_multiple_testing(combined, config$analysis$fdr_threshold)
  combined$testing_family_size <- sum(combined$status == "ok", na.rm = TRUE)
  combined <- data.table::as.data.table(combined)
  class(combined) <- unique(c("phewas_result", class(combined)))

  spec <- attr(manifest, "spec", exact = TRUE)
  identity <- attr(manifest, "identity", exact = TRUE)
  shard_run_metadata <- .pf_default(artifacts[[1L]]$run_metadata, list())
  run_metadata <- list(
    run_fingerprint = run_fingerprint,
    analysis_fingerprint = unique(as.character(manifest$analysis_fingerprint)),
    inputs_fingerprint = unique(as.character(manifest$inputs_fingerprint)),
    input_hashes = identity$input_hashes,
    package_version = unique(as.character(manifest$package_version)),
    shard_count = length(shard_ids),
    phenotype_count = length(expected_targets),
    phenotypes = expected_targets,
    correction_universe = sum(combined$status == "ok", na.rm = TRUE),
    corrections_applied = TRUE,
    transformation = spec$anchor_transform,
    transformation_mean = shard_run_metadata$transformation_mean,
    transformation_sd = shard_run_metadata$transformation_sd,
    n_base = shard_run_metadata$n_base,
    exclusion_active = shard_run_metadata$exclusion_active,
    exclusion_column = shard_run_metadata$exclusion_column,
    exclusion_values = shard_run_metadata$exclusion_values,
    n_before_exclusion = shard_run_metadata$n_before_exclusion,
    n_excluded = shard_run_metadata$n_excluded,
    n_after_exclusion = shard_run_metadata$n_after_exclusion,
    created_at = .pf_timestamp()
  )
  attr(combined, "spec") <- spec
  attr(combined, "run_metadata") <- run_metadata
  results_rds <- file.path(bundle_directory, "results.rds")
  if (!isTRUE(overwrite) && file.exists(results_rds)) {
    existing <- read_phewas_bundle(bundle_directory)
    if (!.pf_results_match_shards(existing, combined, spec, run_metadata)) {
      stop(
        "Existing results do not match the validated shard set. ",
        "Set `overwrite = TRUE` to replace them: ", bundle_directory,
        call. = FALSE
      )
    }
    # results.rds is canonical. Recreate the portable TSV on every safe reuse
    # so a missing or edited table cannot drift from it.
    .pf_write_tsv(existing, file.path(bundle_directory, "results.tsv"))
    message("Reusing combined results: ", bundle_directory)
    attr(existing, "bundle_directory") <- normalizePath(
      bundle_directory, winslash = "/", mustWork = TRUE
    )
    return(existing)
  }
  write_phewas_bundle(
    combined,
    directory = bundle_directory,
    run_metadata = run_metadata,
    overwrite = overwrite
  )
  attr(combined, "bundle_directory") <- normalizePath(
    bundle_directory,
    winslash = "/",
    mustWork = TRUE
  )
  combined
}
