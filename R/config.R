# YAML configuration ---------------------------------------------------------

.pf_default <- function(x, default) {
  if (is.null(x)) default else x
}

.pf_scalar_character <- function(x, name, choices = NULL, allow_null = FALSE) {
  if (is.null(x) && isTRUE(allow_null)) {
    return(NULL)
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be one non-empty string.", call. = FALSE)
  }
  if (!is.null(choices) && !x %in% choices) {
    stop(
      "`", name, "` must be one of: ", paste(choices, collapse = ", "), ".",
      call. = FALSE
    )
  }
  x
}

.pf_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop("`", name, "` must be a positive integer.", call. = FALSE)
  }
  as.integer(x)
}

.pf_resolve_path <- function(path, base_directory, must_exist = FALSE) {
  path <- .pf_scalar_character(path, "path")
  expanded <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", expanded)) {
    expanded <- file.path(base_directory, expanded)
  }
  if (isTRUE(must_exist) && !file.exists(expanded)) {
    stop("Configured file does not exist: ", expanded, call. = FALSE)
  }
  normalized <- normalizePath(expanded, winslash = "/", mustWork = must_exist)
  normalized
}

.pf_path_key <- function(path) {
  normalized <- vapply(as.character(path), function(value) {
    candidate <- path.expand(value)
    suffix <- character()
    while (!file.exists(candidate)) {
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      suffix <- c(basename(candidate), suffix)
      candidate <- parent
    }
    candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    for (component in suffix) {
      if (!nzchar(component) || identical(component, ".")) next
      candidate <- if (identical(component, "..")) {
        dirname(candidate)
      } else {
        file.path(candidate, component)
      }
    }
    candidate
  }, character(1L), USE.NAMES = FALSE)
  if (.Platform$OS.type == "windows") tolower(normalized) else normalized
}

.pf_path_is_within <- function(path, directory) {
  path <- .pf_path_key(path)
  directory <- sub("/+$", "", .pf_path_key(directory))
  identical(path, directory) || startsWith(path, paste0(directory, "/"))
}

.pf_plot_output_paths <- function(output, name = "plot.output") {
  output <- .pf_scalar_character(output, name)
  filename <- basename(output)
  extension <- tolower(tools::file_ext(filename))
  if (endsWith(filename, ".") ||
      (nzchar(extension) && !extension %in% c("png", "pdf"))) {
    stop(
      "`", name,
      "` must be a filename stem or end in `.png` or `.pdf`.",
      call. = FALSE
    )
  }
  stem <- if (nzchar(extension)) {
    substr(output, 1L, nchar(output) - nchar(extension) - 1L)
  } else {
    output
  }
  c(
    manhattan_png = paste0(stem, ".png"),
    manhattan_pdf = paste0(stem, ".pdf"),
    volcano_png = paste0(stem, "-volcano.png"),
    volcano_pdf = paste0(stem, "-volcano.pdf")
  )
}

.pf_config_write_paths <- function(config) {
  c(
    results_rds = file.path(config$output$directory, "results.rds"),
    results_tsv = file.path(config$output$directory, "results.tsv"),
    manifest = file.path(
      config$output$directory, ".phewasflow", "manifest.tsv"
    ),
    stats::setNames(
      file.path(
        config$output$directory, ".phewasflow", "shards",
        sprintf("shard-%04d.rds", seq_len(config$execution$shards))
      ),
      sprintf("shard_%04d", seq_len(config$execution$shards))
    ),
    .pf_plot_output_paths(config$plot$output)
  )
}

.pf_assert_config_path_safety <- function(config) {
  if (file.exists(config$output$directory) &&
      !dir.exists(config$output$directory)) {
    stop("`output.directory` is an existing file, not a directory.",
         call. = FALSE)
  }
  source_path <- attr(config, "config_path", exact = TRUE)
  read_paths <- c(unlist(config$inputs, use.names = FALSE), source_path)
  read_paths <- read_paths[!is.na(read_paths) & nzchar(read_paths)]
  write_paths <- .pf_config_write_paths(config)
  write_keys <- .pf_path_key(write_paths)
  if (anyDuplicated(write_keys)) {
    duplicated_names <- names(write_paths)[
      duplicated(write_keys) | duplicated(write_keys, fromLast = TRUE)
    ]
    stop(
      "Configured write paths collide: ",
      paste(duplicated_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  overlap <- which(write_keys %in% .pf_path_key(read_paths))
  if (length(overlap)) {
    stop(
      "Configured write path would overwrite an input or configuration file: ",
      write_paths[[overlap[[1L]]]],
      call. = FALSE
    )
  }
  state_directory <- file.path(config$output$directory, ".phewasflow")
  plot_paths <- .pf_plot_output_paths(config$plot$output)
  if (any(vapply(
    plot_paths, .pf_path_is_within, logical(1L), directory = state_directory
  ))) {
    stop("`plot.output` cannot be inside the internal `.phewasflow` directory.",
         call. = FALSE)
  }
  invisible(config)
}

.pf_assert_safe_write_target <- function(config, paths, purpose,
                                         allow = character(),
                                         directories = character()) {
  paths <- as.character(paths)
  protected_reads <- c(
    unlist(config$inputs, use.names = FALSE),
    attr(config, "config_path", exact = TRUE)
  )
  protected_reads <- protected_reads[
    !is.na(protected_reads) & nzchar(protected_reads)
  ]
  protected_writes <- .pf_config_write_paths(config)
  allow_keys <- .pf_path_key(allow)
  protected_writes <- protected_writes[
    !.pf_path_key(protected_writes) %in% allow_keys
  ]
  protected <- c(protected_reads, protected_writes)
  collision <- which(.pf_path_key(paths) %in% .pf_path_key(protected))
  if (length(collision)) {
    stop(
      purpose, " path collides with a protected analysis path: ",
      paths[[collision[[1L]]]],
      call. = FALSE
    )
  }
  if (length(directories)) {
    directory_collision <- which(
      .pf_path_key(directories) %in% .pf_path_key(protected_reads)
    )
    if (length(directory_collision)) {
      stop(
        purpose, " directory collides with an input or configuration file: ",
        directories[[directory_collision[[1L]]]],
        call. = FALSE
      )
    }
  }
  invisible(paths)
}

.pf_normalize_character_vector <- function(x, name) {
  if (is.null(x)) {
    return(character())
  }
  x <- unlist(x, recursive = TRUE, use.names = FALSE)
  if (!length(x)) {
    return(character())
  }
  if (!is.character(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
    stop("`", name, "` must contain unique, non-empty strings.", call. = FALSE)
  }
  x
}

.pf_reject_unknown_fields <- function(x, allowed, name) {
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) {
    stop("Unknown `", name, "` field(s): ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
}

.pf_normalize_config <- function(config, config_path = NULL) {
  if (!is.list(config) || is.data.frame(config)) {
    stop("The YAML document must contain a named mapping.", call. = FALSE)
  }
  if (is.null(names(config))) {
    stop("The YAML document must contain named top-level fields.", call. = FALSE)
  }
  .pf_reject_unknown_fields(
    config, c("schema_version", "inputs", "analysis", "output", "execution", "plot"),
    "top-level"
  )
  source_path <- .pf_default(config_path, attr(config, "config_path", exact = TRUE))
  if (is.null(source_path)) {
    base_directory <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  } else {
    source_path <- normalizePath(source_path, winslash = "/", mustWork = TRUE)
    base_directory <- dirname(source_path)
  }

  schema_version <- as.character(.pf_default(config$schema_version, ""))
  if (length(schema_version) != 1L || !schema_version %in% c("1", "1.0")) {
    stop("`schema_version` must be 1.", call. = FALSE)
  }
  if (!is.list(config$inputs)) {
    stop("`inputs` must be a mapping.", call. = FALSE)
  }
  input_names <- names(config$inputs)
  required_inputs <- c("data", "phenotype_metadata")
  missing_inputs <- setdiff(required_inputs, input_names)
  if (length(missing_inputs)) {
    stop(
      "Missing required input", if (length(missing_inputs) > 1L) "s" else "", ": ",
      paste(missing_inputs, collapse = ", "), ".",
      call. = FALSE
    )
  }
  unknown_inputs <- setdiff(input_names, required_inputs)
  if (length(unknown_inputs)) {
    stop("Unknown `inputs` field(s): ", paste(unknown_inputs, collapse = ", "), ".", call. = FALSE)
  }
  inputs <- list(
    data = .pf_resolve_path(config$inputs$data, base_directory, must_exist = TRUE),
    phenotype_metadata = .pf_resolve_path(
      config$inputs$phenotype_metadata,
      base_directory,
      must_exist = TRUE
    )
  )

  if (!is.list(config$analysis)) {
    stop("`analysis` must be a mapping.", call. = FALSE)
  }
  analysis <- config$analysis
  .pf_reject_unknown_fields(
    analysis,
    c(
      "analysis_id", "anchor", "direction", "covariates", "anchor_type",
      "anchor_reference", "anchor_levels", "anchor_scores", "outcome_type",
      "anchor_transform", "offset", "exclusion_column", "exclusion_values",
      "version", "fdr_threshold", "eligibility"
    ),
    "analysis"
  )
  analysis$analysis_id <- .pf_scalar_character(analysis$analysis_id, "analysis.analysis_id")
  analysis$anchor <- .pf_scalar_character(analysis$anchor, "analysis.anchor")
  analysis$direction <- .pf_scalar_character(
    analysis$direction,
    "analysis.direction",
    c("phenotypes_as_outcomes", "phenotypes_as_predictors")
  )
  analysis$covariates <- .pf_normalize_character_vector(
    .pf_default(analysis$covariates, character()),
    "analysis.covariates"
  )
  analysis$anchor_type <- .pf_scalar_character(
    analysis$anchor_type,
    "analysis.anchor_type",
    c("numeric", "binary", "count", "ordinal"),
    allow_null = TRUE
  )
  if (!is.null(analysis$anchor_reference)) {
    analysis$anchor_reference <- as.character(analysis$anchor_reference)
  }
  analysis$anchor_reference <- .pf_scalar_character(
    analysis$anchor_reference,
    "analysis.anchor_reference",
    allow_null = TRUE
  )
  anchor_levels <- unlist(
    .pf_default(analysis$anchor_levels, character()),
    recursive = TRUE,
    use.names = FALSE
  )
  if (length(anchor_levels)) {
    anchor_levels <- as.character(anchor_levels)
  }
  analysis$anchor_levels <- .pf_normalize_character_vector(
    anchor_levels,
    "analysis.anchor_levels"
  )
  anchor_scores <- .pf_default(analysis$anchor_scores, numeric())
  anchor_scores <- unlist(anchor_scores, recursive = TRUE, use.names = FALSE)
  if (!length(anchor_scores)) {
    anchor_scores <- numeric()
  }
  if (!is.numeric(anchor_scores) || anyNA(anchor_scores) || any(!is.finite(anchor_scores))) {
    stop("`analysis.anchor_scores` must contain finite numbers.", call. = FALSE)
  }
  analysis$anchor_scores <- as.numeric(anchor_scores)
  analysis$outcome_type <- .pf_scalar_character(
    analysis$outcome_type,
    "analysis.outcome_type",
    c("binary", "continuous", "count", "ordinal"),
    allow_null = TRUE
  )
  if (identical(analysis$direction, "phenotypes_as_predictors") && is.null(analysis$outcome_type)) {
    stop(
      "`analysis.outcome_type` is required when phenotypes are predictors.",
      call. = FALSE
    )
  }
  analysis$anchor_transform <- .pf_scalar_character(
    .pf_default(analysis$anchor_transform, "none"),
    "analysis.anchor_transform",
    c("none", "zscore")
  )
  analysis$offset <- .pf_scalar_character(
    analysis$offset,
    "analysis.offset",
    allow_null = TRUE
  )
  analysis$exclusion_column <- .pf_scalar_character(
    analysis$exclusion_column,
    "analysis.exclusion_column",
    allow_null = TRUE
  )
  exclusion_values_declared <- !is.null(analysis$exclusion_values)
  if (is.null(analysis$exclusion_column)) {
    if (exclusion_values_declared) {
      stop(
        "`analysis.exclusion_values` is only allowed with `analysis.exclusion_column`.",
        call. = FALSE
      )
    }
    analysis$exclusion_values <- NULL
  } else {
    exclusion_values <- unlist(
      analysis$exclusion_values,
      recursive = TRUE,
      use.names = FALSE
    )
    if (!length(exclusion_values) || !is.atomic(exclusion_values) ||
        anyNA(exclusion_values) || anyDuplicated(exclusion_values)) {
      stop(
        "`analysis.exclusion_values` must be a nonempty vector of unique, nonmissing values.",
        call. = FALSE
      )
    }
    analysis$exclusion_values <- exclusion_values
  }
  analysis$version <- .pf_scalar_character(
    as.character(.pf_default(analysis$version, "1.0")),
    "analysis.version"
  )
  if (is.null(analysis$fdr_threshold)) {
    stop("`analysis.fdr_threshold` must be explicitly specified.", call. = FALSE)
  }
  if (!is.numeric(analysis$fdr_threshold) || length(analysis$fdr_threshold) != 1L ||
      is.na(analysis$fdr_threshold) || analysis$fdr_threshold <= 0 ||
      analysis$fdr_threshold >= 1) {
    stop("`analysis.fdr_threshold` must be a number strictly between 0 and 1.", call. = FALSE)
  }
  if (is.null(analysis$eligibility)) {
    stop(
      "`analysis.eligibility` must explicitly specify min_n, min_cases, ",
      "min_controls, and min_outcome_levels.",
      call. = FALSE
    )
  }
  if (!is.list(analysis$eligibility)) {
    stop("`analysis.eligibility` must be a mapping.", call. = FALSE)
  }
  eligibility_fields <- c("min_n", "min_cases", "min_controls", "min_outcome_levels")
  missing_eligibility <- setdiff(eligibility_fields, names(analysis$eligibility))
  if (length(missing_eligibility)) {
    stop(
      "Missing eligibility field(s): ", paste(missing_eligibility, collapse = ", "), ".",
      call. = FALSE
    )
  }
  for (field in eligibility_fields) {
    analysis$eligibility[[field]] <- .pf_positive_integer(
      analysis$eligibility[[field]],
      paste0("analysis.eligibility.", field)
    )
  }
  if (analysis$eligibility$min_outcome_levels < 3L) {
    stop("`analysis.eligibility.min_outcome_levels` must be at least 3.", call. = FALSE)
  }
  unknown_eligibility <- setdiff(names(analysis$eligibility), eligibility_fields)
  if (length(unknown_eligibility)) {
    stop(
      "Unknown eligibility field(s): ", paste(unknown_eligibility, collapse = ", "), ".",
      call. = FALSE
    )
  }

  if (!is.list(config$output)) {
    stop("`output` must be a mapping with a `directory` field.", call. = FALSE)
  }
  .pf_reject_unknown_fields(config$output, "directory", "output")
  output <- list(
    directory = .pf_resolve_path(config$output$directory, base_directory, must_exist = FALSE)
  )

  execution <- .pf_default(config$execution, list())
  if (!is.list(execution)) {
    stop("`execution` must be a mapping.", call. = FALSE)
  }
  .pf_reject_unknown_fields(execution, c("shards", "backend", "workers"),
                            "execution")
  execution$shards <- .pf_positive_integer(.pf_default(execution$shards, 1L), "execution.shards")
  execution$backend <- .pf_scalar_character(
    .pf_default(execution$backend, "sequential"),
    "execution.backend",
    c("sequential", "psock")
  )
  execution$workers <- .pf_positive_integer(.pf_default(execution$workers, 1L), "execution.workers")

  plot <- .pf_default(config$plot, list())
  if (!is.list(plot)) {
    stop("`plot` must be a mapping.", call. = FALSE)
  }
  .pf_reject_unknown_fields(
    plot,
    c(
      "output", "significance", "group_display", "label", "highlight",
      "width", "height", "dpi"
    ),
    "plot"
  )
  plot$output <- .pf_resolve_path(
    .pf_default(plot$output, file.path(output$directory, "phewas.png")),
    base_directory,
    must_exist = FALSE
  )
  .pf_plot_output_paths(plot$output)
  plot$significance <- .pf_scalar_character(
    .pf_default(plot$significance, "bonferroni"),
    "plot.significance",
    c("bh", "bonferroni", "none")
  )
  plot$group_display <- .pf_scalar_character(
    .pf_default(plot$group_display, "legend"),
    "plot.group_display",
    c("auto", "x_axis", "legend")
  )
  plot$label <- .pf_normalize_character_vector(.pf_default(plot$label, character()), "plot.label")
  plot$highlight <- .pf_normalize_character_vector(
    .pf_default(plot$highlight, character()),
    "plot.highlight"
  )
  for (field in c("width", "height", "dpi")) {
    default <- switch(field, width = 10, height = 6, dpi = 300)
    plot[[field]] <- .pf_default(plot[[field]], default)
    if (!is.numeric(plot[[field]]) || length(plot[[field]]) != 1L ||
        !is.finite(plot[[field]]) || plot[[field]] <= 0) {
      stop("`plot.", field, "` must be a finite positive number.", call. = FALSE)
    }
  }

  normalized <- list(
    schema_version = "1",
    inputs = inputs,
    analysis = analysis,
    output = output,
    execution = execution,
    plot = plot
  )
  class(normalized) <- c("phewasflow_config", "list")
  attr(normalized, "config_path") <- source_path
  attr(normalized, "base_directory") <- base_directory
  .pf_assert_config_path_safety(normalized)
  normalized
}

.pf_parse_metadata_lists <- function(metadata) {
  parse_delimited <- function(x, numeric = FALSE) {
    if (is.list(x)) {
      return(x)
    }
    lapply(x, function(value) {
      if (length(value) == 0L || is.na(value) || !nzchar(trimws(as.character(value)))) {
        return(NULL)
      }
      pieces <- trimws(strsplit(as.character(value), "|", fixed = TRUE)[[1L]])
      if (isTRUE(numeric)) {
        converted <- suppressWarnings(as.numeric(pieces))
        if (anyNA(converted)) {
          stop("Ordinal `scores` must be numeric values separated by `|`.", call. = FALSE)
        }
        return(converted)
      }
      pieces
    })
  }
  if ("levels" %in% names(metadata)) {
    metadata$levels <- parse_delimited(metadata$levels)
  }
  if ("scores" %in% names(metadata)) {
    metadata$scores <- parse_delimited(metadata$scores, numeric = TRUE)
  }
  metadata
}

.pf_load_config_inputs <- function(config) {
  list(
    data = .pf_read_table(config$inputs$data, "Analysis data"),
    phenotypes = .pf_parse_metadata_lists(
      .pf_read_phenotype_metadata(config$inputs$phenotype_metadata)
    )
  )
}

#' Read a phewasFlow YAML configuration
#'
#' Relative paths in the file are resolved relative to the YAML file, not the
#' current working directory.
#'
#' @param path Path to a YAML configuration file.
#' @param check_data Also load the analysis table and validate the specification
#'   against it. Set to `FALSE` for a quicker structural check.
#'
#' @return A normalized `phewasflow_config` list.
#' @export
read_phewasflow_config <- function(path, check_data = TRUE) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  config <- yaml::read_yaml(path)
  validate_phewasflow_config(config, config_path = path, check_data = check_data)
}

#' Validate a phewasFlow configuration
#'
#' @param config A list returned by `yaml::read_yaml()` or a normalized config.
#' @param config_path Optional source YAML path used to resolve relative paths.
#' @param check_data Whether to validate the specification against the analysis
#'   data in addition to checking metadata and paths.
#'
#' @return A normalized `phewasflow_config` list.
#' @export
validate_phewasflow_config <- function(config, config_path = NULL, check_data = TRUE) {
  normalized <- .pf_normalize_config(config, config_path = config_path)
  phenotypes <- .pf_parse_metadata_lists(
    .pf_read_phenotype_metadata(normalized$inputs$phenotype_metadata)
  )
  spec <- .pf_spec_from_config(normalized, phenotypes)
  if (isTRUE(check_data)) {
    data <- .pf_read_table(normalized$inputs$data, "Analysis data")
    validate_phewas_spec(spec, data = data)
  }
  normalized
}

.pf_spec_from_config <- function(config, phenotypes) {
  analysis <- config$analysis
  phewas_spec(
    analysis_id = analysis$analysis_id,
    anchor = analysis$anchor,
    direction = analysis$direction,
    phenotypes = phenotypes,
    covariates = analysis$covariates,
    outcome_type = analysis$outcome_type,
    anchor_type = analysis$anchor_type,
    anchor_reference = analysis$anchor_reference,
    anchor_levels = analysis$anchor_levels,
    anchor_scores = analysis$anchor_scores,
    anchor_transform = analysis$anchor_transform,
    offset = analysis$offset,
    exclusion_column = analysis$exclusion_column,
    exclusion_values = analysis$exclusion_values,
    eligibility = analysis$eligibility,
    fdr_threshold = analysis$fdr_threshold,
    version = analysis$version
  )
}

#' Build a PheWAS specification from a configuration
#'
#' @param config A configuration path or a normalized config object.
#' @param data Optional analysis data used for validation.
#'
#' @return A validated `phewas_spec` object.
#' @export
config_to_phewas_spec <- function(config, data = NULL) {
  if (is.character(config)) {
    config <- read_phewasflow_config(config, check_data = FALSE)
  } else if (!inherits(config, "phewasflow_config")) {
    config <- validate_phewasflow_config(config, check_data = FALSE)
  }
  phenotypes <- .pf_parse_metadata_lists(
    .pf_read_phenotype_metadata(config$inputs$phenotype_metadata)
  )
  validate_phewas_spec(.pf_spec_from_config(config, phenotypes), data = data)
}

.pf_package_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("phewasFlow")),
    error = function(error) "0.1.0"
  )
}

.pf_canonicalize <- function(x) {
  if (is.factor(x)) {
    return(as.character(x))
  }
  if (inherits(x, "POSIXt")) {
    return(format(x, "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"))
  }
  if (is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    if ("phenotype" %in% names(x)) {
      x <- x[order(as.character(x$phenotype), method = "radix"), , drop = FALSE]
    }
    x <- x[, sort(names(x)), drop = FALSE]
    rownames(x) <- NULL
    return(lapply(x, .pf_canonicalize))
  }
  if (is.list(x)) {
    if (!is.null(names(x))) {
      x <- x[sort(names(x))]
    }
    return(lapply(x, .pf_canonicalize))
  }
  if (!is.null(names(x))) {
    x <- x[order(names(x), method = "radix")]
  }
  unname(x)
}

.pf_object_sha256 <- function(x) {
  serialized <- jsonlite::toJSON(
    .pf_canonicalize(x),
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = NA,
    dataframe = "columns"
  )
  unname(digest::digest(serialized, algo = "sha256", serialize = FALSE))
}

.pf_input_hashes <- function(config) {
  paths <- config$inputs[sort(names(config$inputs))]
  stats::setNames(lapply(paths, .pf_file_sha256), names(paths))
}

.pf_analysis_identity <- function(config, spec, phenotypes) {
  phenotypes <- sort(as.character(phenotypes), method = "radix")
  package_version <- .pf_package_version()
  input_hashes <- .pf_input_hashes(config)
  analysis_fingerprint <- .pf_object_sha256(unclass(spec))
  inputs_fingerprint <- .pf_object_sha256(input_hashes)
  run_fingerprint <- .pf_object_sha256(list(
    format_version = "2",
    package_version = package_version,
    analysis_fingerprint = analysis_fingerprint,
    input_hashes = input_hashes,
    phenotypes = phenotypes
  ))
  list(
    package_version = package_version,
    input_hashes = input_hashes,
    inputs_fingerprint = inputs_fingerprint,
    analysis_fingerprint = analysis_fingerprint,
    run_fingerprint = run_fingerprint
  )
}
