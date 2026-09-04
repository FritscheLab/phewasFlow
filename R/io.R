# File I/O -------------------------------------------------------------------

.pf_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

.pf_require_parent <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
    stop("Could not create directory: ", parent, call. = FALSE)
  }
  invisible(parent)
}

.pf_file_rename <- function(from, to) {
  file.rename(from, to)
}

.pf_atomic_write <- function(path, writer) {
  stopifnot(is.function(writer), length(path) == 1L, !is.na(path))
  if (dir.exists(path)) {
    stop("Cannot replace a directory with a file: ", path, call. = FALSE)
  }
  .pf_require_parent(path)
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  writer(temporary)
  if (!file.exists(temporary)) {
    stop("Atomic writer did not create its temporary file.", call. = FALSE)
  }
  if (!.pf_file_rename(temporary, path)) {
    if (!file.exists(path)) {
      stop("Could not publish file atomically: ", path, call. = FALSE)
    }
    # Windows cannot replace an existing file with file.rename(). Move the
    # previous file aside, then restore it if publication still fails.
    backup <- tempfile(
      pattern = paste0(".", basename(path), ".backup."),
      tmpdir = dirname(path)
    )
    if (!.pf_file_rename(path, backup)) {
      stop("Could not preserve the existing file before replacement: ", path,
           call. = FALSE)
    }
    if (!.pf_file_rename(temporary, path)) {
      restored <- .pf_file_rename(backup, path)
      if (!restored) {
        stop(
          "Could not publish or restore `", path,
          "`; the previous file remains at `", backup, "`.",
          call. = FALSE
        )
      }
      stop("Could not publish file atomically: ", path, call. = FALSE)
    }
    unlink(backup, force = TRUE)
  }
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

.pf_publish_file_set <- function(staged_paths, destination_paths) {
  if (length(staged_paths) != length(destination_paths) ||
      !length(staged_paths) || any(!file.exists(staged_paths))) {
    stop("Every staged file must exist and have one destination.", call. = FALSE)
  }
  destination_paths <- as.character(destination_paths)
  if (anyDuplicated(destination_paths)) {
    stop("Published file destinations must be unique.", call. = FALSE)
  }
  if (any(dir.exists(destination_paths))) {
    stop("Cannot replace a directory with a result file.", call. = FALSE)
  }
  invisible(lapply(destination_paths, .pf_require_parent))

  backups <- rep(NA_character_, length(destination_paths))
  published <- rep(FALSE, length(destination_paths))
  committed <- FALSE
  rollback <- function() {
    for (index in rev(which(published))) {
      if (file.exists(destination_paths[[index]])) {
        unlink(destination_paths[[index]], force = TRUE)
      }
    }
    for (index in which(!is.na(backups))) {
      if (file.exists(backups[[index]]) &&
          !.pf_file_rename(backups[[index]], destination_paths[[index]])) {
        warning(
          "Could not restore `", destination_paths[[index]],
          "`; the previous file remains at `", backups[[index]], "`.",
          call. = FALSE
        )
      }
    }
  }
  on.exit(if (!committed) rollback(), add = TRUE)

  for (index in seq_along(destination_paths)) {
    destination <- destination_paths[[index]]
    if (file.exists(destination)) {
      backup <- tempfile(
        pattern = paste0(".", basename(destination), ".backup."),
        tmpdir = dirname(destination)
      )
      if (!.pf_file_rename(destination, backup)) {
        stop("Could not preserve existing result file: ", destination,
             call. = FALSE)
      }
      backups[[index]] <- backup
    }
  }
  for (index in seq_along(destination_paths)) {
    if (!.pf_file_rename(staged_paths[[index]], destination_paths[[index]])) {
      stop("Could not publish result file: ", destination_paths[[index]],
           call. = FALSE)
    }
    published[[index]] <- TRUE
  }
  committed <- TRUE
  unlink(backups[!is.na(backups)], force = TRUE)
  invisible(destination_paths)
}

.pf_write_rds <- function(object, path) {
  .pf_atomic_write(path, function(temporary) {
    saveRDS(object, temporary, version = 3L)
  })
}

.pf_write_tsv <- function(x, path) {
  .pf_atomic_write(path, function(temporary) {
    data.table::fwrite(
      as.data.frame(x, stringsAsFactors = FALSE),
      temporary,
      sep = "\t",
      quote = TRUE,
      na = "NA"
    )
  })
}

.pf_file_sha256 <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }
  unname(digest::digest(file = path, algo = "sha256", serialize = FALSE))
}

.pf_read_table <- function(path, label = "input", character_columns = character()) {
  if (!file.exists(path)) {
    stop(label, " file does not exist: ", path, call. = FALSE)
  }
  read_delimited <- function(separator = NULL) {
    arguments <- list(input = path, data.table = FALSE)
    if (!is.null(separator)) {
      arguments$sep <- separator
    }
    header <- do.call(data.table::fread, c(arguments, list(nrows = 0L)))
    forced <- intersect(character_columns, names(header))
    if (length(forced)) {
      arguments$colClasses <- list(character = forced)
    }
    do.call(data.table::fread, arguments)
  }
  extension <- tolower(tools::file_ext(path))
  value <- switch(
    extension,
    rds = readRDS(path),
    csv = read_delimited(","),
    tsv = read_delimited("\t"),
    tab = read_delimited("\t"),
    txt = read_delimited(),
    stop(
      label,
      " must be an .rds, .csv, .tsv, .tab, or .txt file: ",
      path,
      call. = FALSE
    )
  )
  if (!is.data.frame(value)) {
    stop(label, " must contain a data frame.", call. = FALSE)
  }
  for (column in intersect(character_columns, names(value))) {
    value[[column]] <- as.character(value[[column]])
  }
  as.data.frame(value, stringsAsFactors = FALSE)
}

.pf_read_phenotype_metadata <- function(path) {
  .pf_read_table(
    path,
    label = "Phenotype metadata",
    character_columns = c(
      "phenotype", "phenotype_id", "phecode", "phenotype_column", "column"
    )
  )
}

.pf_normalize_legacy_result <- function(result, label = "PheWAS result") {
  if (!is.data.frame(result)) {
    return(result)
  }
  aliases <- list(
    phenotype = c("phenotype_id", "phecode", "phenotype_column", "column"),
    description = "label",
    group = "category",
    groupnum = "category_order"
  )
  specification <- attr(result, "spec", exact = TRUE)
  run_metadata <- attr(result, "run_metadata", exact = TRUE)
  out <- as.data.frame(
    if (inherits(result, "data.table")) data.table::copy(result) else result,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  fields_agree <- function(x, y, numeric = FALSE) {
    if (isTRUE(numeric)) {
      x <- suppressWarnings(as.numeric(x))
      y <- suppressWarnings(as.numeric(y))
    } else {
      x <- as.character(x)
      y <- as.character(y)
    }
    length(x) == length(y) && all(
      (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
    )
  }

  for (target in names(aliases)) {
    candidates <- intersect(c(target, aliases[[target]]), names(out))
    if (!length(candidates)) {
      next
    }
    reference <- candidates[[1L]]
    for (candidate in candidates[-1L]) {
      if (!fields_agree(
        out[[reference]], out[[candidate]], numeric = identical(target, "groupnum")
      )) {
        stop(
          label, " cannot be safely normalized because columns `", reference,
          "` and `", candidate, "` disagree for canonical `", target, "`.",
          call. = FALSE
        )
      }
    }
    if (!target %in% names(out)) {
      names(out)[match(reference, names(out))] <- target
    }
    for (alias in intersect(aliases[[target]], names(out))) {
      out[[alias]] <- NULL
    }
  }

  if (inherits(result, "phewas_result")) {
    class(out) <- c("phewas_result", "data.frame")
  }
  if (inherits(specification, "phewas_spec")) {
    specification <- tryCatch(
      validate_phewas_spec(specification),
      error = function(error) {
        stop(
          label, " has legacy specification metadata that cannot be safely ",
          "normalized: ", conditionMessage(error),
          call. = FALSE
        )
      }
    )
  }
  if (is.list(run_metadata) && "phenotype_ids" %in% names(run_metadata)) {
    if ("phenotypes" %in% names(run_metadata) && !fields_agree(
      run_metadata$phenotypes, run_metadata$phenotype_ids
    )) {
      stop(
        label, " cannot be safely normalized because run metadata fields ",
        "`phenotypes` and legacy `phenotype_ids` disagree.",
        call. = FALSE
      )
    }
    if (!"phenotypes" %in% names(run_metadata)) {
      run_metadata$phenotypes <- run_metadata$phenotype_ids
    }
    run_metadata$phenotype_ids <- NULL
  }
  attr(out, "spec") <- specification
  attr(out, "run_metadata") <- run_metadata
  out
}

.pf_validate_result_object <- function(result, label = "PheWAS result",
                                       require_fingerprint = TRUE) {
  if (!is.data.frame(result) || !inherits(result, "phewas_result")) {
    stop(label, " must be a `phewas_result` data frame.", call. = FALSE)
  }
  required <- c(
    "analysis_id", "spec_version", "phenotype", "description", "group",
    "groupnum", "color", "direction", "response", "predictor", "status",
    "reason_code", "log_p", "testing_family_size"
  )
  missing <- setdiff(required, names(result))
  if (length(missing)) {
    stop(
      label, " is missing required column(s): ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  legacy <- intersect(
    c(
      "phenotype_id", "phecode", "phenotype_column", "label", "category",
      "category_order", "column"
    ),
    names(result)
  )
  if (length(legacy)) {
    stop(
      label, " contains legacy result column(s): ", paste(legacy, collapse = ", "),
      ". Canonical results must use `phenotype`, `description`, `group`, ",
      "`groupnum`, and `color`.",
      call. = FALSE
    )
  }
  if (!nrow(result)) {
    stop(label, " contains no association rows.", call. = FALSE)
  }
  phenotypes <- as.character(result$phenotype)
  if (anyNA(phenotypes) || any(!nzchar(phenotypes)) ||
      anyDuplicated(phenotypes)) {
    stop(label, " has invalid or duplicated phenotypes.", call. = FALSE)
  }
  for (field in c("description", "group", "color")) {
    values <- as.character(result[[field]])
    if (anyNA(values) || any(!nzchar(values))) {
      stop(label, " has invalid `", field, "` values.", call. = FALSE)
    }
  }
  groupnum <- suppressWarnings(as.numeric(result$groupnum))
  if (length(groupnum) != nrow(result) || any(!is.finite(groupnum))) {
    stop(label, " has invalid `groupnum` values.", call. = FALSE)
  }
  valid_color <- vapply(as.character(result$color), function(value) {
    tryCatch(
      {
        grDevices::col2rgb(value)
        TRUE
      },
      error = function(error) FALSE
    )
  }, logical(1L))
  if (any(!valid_color)) {
    stop(label, " has invalid `color` values.", call. = FALSE)
  }
  analysis_ids <- unique(as.character(result$analysis_id))
  if (length(analysis_ids) != 1L || is.na(analysis_ids) ||
      !nzchar(analysis_ids)) {
    stop(label, " must contain one non-empty analysis ID.", call. = FALSE)
  }
  statuses <- as.character(result$status)
  if (anyNA(statuses) || any(!statuses %in% c("ok", "skipped", "error"))) {
    stop(label, " contains an invalid association status.", call. = FALSE)
  }
  directions <- unique(as.character(result$direction))
  if (length(directions) != 1L || is.na(directions) ||
      !directions %in% c("phenotypes_as_outcomes", "phenotypes_as_predictors")) {
    stop(label, " contains an invalid scan direction.", call. = FALSE)
  }

  specification <- attr(result, "spec", exact = TRUE)
  if (!inherits(specification, "phewas_spec") ||
      !identical(as.character(specification$analysis_id), analysis_ids)) {
    stop(label, " has missing or inconsistent specification metadata.",
         call. = FALSE)
  }
  canonical_specification <- tryCatch(
    validate_phewas_spec(specification),
    error = function(error) {
      stop(
        label, " has invalid specification metadata: ", conditionMessage(error),
        call. = FALSE
      )
    }
  )
  if (!identical(specification, canonical_specification)) {
    stop(
      label, " contains noncanonical specification metadata. New results must ",
      "not retain legacy phenotype aliases.",
      call. = FALSE
    )
  }
  run_metadata <- attr(result, "run_metadata", exact = TRUE)
  if (!is.list(run_metadata) || is.null(names(run_metadata))) {
    stop(label, " has missing or invalid run metadata.", call. = FALSE)
  }
  if ("phenotype_ids" %in% names(run_metadata)) {
    stop(
      label, " contains legacy run metadata field `phenotype_ids`; use ",
      "`phenotypes`.", call. = FALSE
    )
  }
  if ("phenotypes" %in% names(run_metadata) &&
      !identical(as.character(run_metadata$phenotypes), phenotypes)) {
    stop(label, " has run metadata that disagrees with `phenotype` rows.",
         call. = FALSE)
  }
  if (isTRUE(require_fingerprint)) {
    fingerprint <- run_metadata$run_fingerprint
    if (!is.character(fingerprint) || length(fingerprint) != 1L ||
        is.na(fingerprint) || !nzchar(fingerprint)) {
      stop(label, " has no valid run fingerprint.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Write PheWAS results
#'
#' Writes an RDS file that preserves the result attributes and a TSV file for
#' use outside R. Individual-level input data and fitted model objects are not
#' included.
#'
#' @param result A `phewas_result` object.
#' @param directory Destination directory.
#' @param run_metadata A named list of non-identifying run metadata.
#' @param overwrite Whether existing result files may be replaced.
#'
#' @return The normalized result directory, invisibly.
#' @export
write_phewas_bundle <- function(result, directory, run_metadata = list(),
                                overwrite = FALSE) {
  .pf_validate_result_object(
    result, label = "`result`", require_fingerprint = FALSE
  )
  if (!is.list(run_metadata) || (length(run_metadata) && is.null(names(run_metadata)))) {
    stop("`run_metadata` must be a named list.", call. = FALSE)
  }
  if (!length(run_metadata)) {
    inherited_metadata <- attr(result, "run_metadata", exact = TRUE)
    if (is.list(inherited_metadata)) {
      run_metadata <- inherited_metadata
    }
  }
  if (!is.character(directory) || length(directory) != 1L ||
      is.na(directory) || !nzchar(directory)) {
    stop("`directory` must be one non-empty path.", call. = FALSE)
  }

  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
    stop("Could not create result directory: ", directory, call. = FALSE)
  }
  results_rds <- file.path(directory, "results.rds")
  results_tsv <- file.path(directory, "results.tsv")
  if (any(file.exists(c(results_rds, results_tsv))) && !isTRUE(overwrite)) {
    stop(
      "Results already exist. Set `overwrite = TRUE` to replace them: ",
      directory,
      call. = FALSE
    )
  }

  run_metadata$run_fingerprint <- .pf_default(
    run_metadata$run_fingerprint,
    .pf_object_sha256(as.data.frame(result, stringsAsFactors = FALSE))
  )
  run_metadata$created_at <- .pf_default(run_metadata$created_at, .pf_timestamp())
  run_metadata$package_version <- .pf_default(run_metadata$package_version, .pf_package_version())

  # Serialize only documented result data and metadata.
  portable_result <- as.data.frame(
    data.table::copy(result),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  attributes(portable_result) <- attributes(portable_result)[
    c("names", "row.names", "class")
  ]
  class(portable_result) <- c("phewas_result", "data.frame")
  attr(portable_result, "spec") <- attr(result, "spec", exact = TRUE)
  attr(portable_result, "run_metadata") <- run_metadata
  .pf_validate_result_object(portable_result, label = "`result`")

  staging <- tempfile(pattern = ".phewas-results.", tmpdir = directory)
  if (!dir.create(staging)) {
    stop("Could not create result staging directory: ", staging, call. = FALSE)
  }
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  staged_tsv <- file.path(staging, "results.tsv")
  staged_rds <- file.path(staging, "results.rds")
  .pf_write_tsv(portable_result, staged_tsv)
  .pf_write_rds(portable_result, staged_rds)
  # Publish the human-readable table first and the canonical RDS last. If
  # either replacement fails, both previous files are restored.
  .pf_publish_file_set(
    c(staged_tsv, staged_rds),
    c(results_tsv, results_rds)
  )
  invisible(normalizePath(directory, winslash = "/", mustWork = TRUE))
}

#' Read PheWAS results
#'
#' @param directory Result directory created by [write_phewas_bundle()].
#'
#' @return A `phewas_result` object.
#' @export
read_phewas_bundle <- function(directory) {
  if (!is.character(directory) || length(directory) != 1L ||
      is.na(directory) || !nzchar(directory)) {
    stop("`directory` must be one non-empty path.", call. = FALSE)
  }
  path <- file.path(directory, "results.rds")
  if (!file.exists(path)) {
    stop("Result file does not exist: ", path, call. = FALSE)
  }
  result <- tryCatch(
    readRDS(path),
    error = function(error) {
      stop("Could not read result file `", path, "`: ", conditionMessage(error),
           call. = FALSE)
    }
  )
  label <- paste0("Result file `", path, "`")
  result <- .pf_normalize_legacy_result(result, label = label)
  .pf_validate_result_object(result, label = label)
  specification <- attr(result, "spec", exact = TRUE)
  run_metadata <- attr(result, "run_metadata", exact = TRUE)
  result <- data.table::as.data.table(result)
  class(result) <- unique(c("phewas_result", class(result)))
  attr(result, "spec") <- specification
  attr(result, "run_metadata") <- run_metadata
  result
}
