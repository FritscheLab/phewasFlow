`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.stop_bad_argument <- function(argument, message) {
  stop(sprintf("`%s` %s", argument, message), call. = FALSE)
}

.scalar_character <- function(x, argument, allow_empty = FALSE) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      (!allow_empty && !nzchar(x))) {
    .stop_bad_argument(argument, "must be one nonmissing character string.")
  }
  x
}

.plain_data_frame <- function(x, argument) {
  if (!is.data.frame(x)) {
    .stop_bad_argument(argument, "must be a data frame.")
  }
  # data.table uses reference semantics. Copy it before converting so every
  # core routine can freely work with a private, ordinary data frame.
  if (data.table::is.data.table(x)) {
    x <- data.table::copy(x)
  }
  out <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  if (anyDuplicated(names(out))) {
    duplicates <- unique(names(out)[duplicated(names(out))])
    .stop_bad_argument(
      argument,
      sprintf("has duplicated column names: %s.", paste(duplicates, collapse = ", "))
    )
  }
  out
}

.validate_id_table <- function(x, id, argument) {
  x <- .plain_data_frame(x, argument)
  if (!id %in% names(x)) {
    .stop_bad_argument(argument, sprintf("does not contain ID column `%s`.", id))
  }

  ids <- x[[id]]
  if (!is.atomic(ids) || is.list(ids)) {
    .stop_bad_argument(argument, "must use an atomic ID column.")
  }
  missing_id <- is.na(ids)
  if (is.character(ids) || is.factor(ids)) {
    missing_id <- missing_id | !nzchar(as.character(ids))
  }
  if (any(missing_id)) {
    .stop_bad_argument(argument, "contains missing or blank IDs.")
  }
  if (anyDuplicated(ids)) {
    .stop_bad_argument(argument, "must have one row per ID.")
  }
  x
}

#' Assemble analysis-ready PheWAS data
#'
#' Strictly left-join phenotype, anchor, and covariate tables to a declared
#' sample spine. Every input must contain one unique, nonmissing ID per row.
#' Non-ID column collisions are errors, and the input objects are never
#' modified by reference.
#'
#' @param sample Data frame defining the people and row order in the analysis.
#' @param phenotypes Optional data frame of analysis-ready phenotypes.
#' @param anchors Optional data frame of fixed predictors or outcomes.
#' @param covariates Optional data frame of adjustment variables.
#' @param id Name of the ID column shared by every supplied table.
#'
#' @return A `data.table` in sample-spine order. Its `join_report` attribute
#'   reports matches and unmatched rows for every input table.
#'
#' @details Input data tables are copied before any work is done, so this
#'   function never changes caller-owned objects by reference. Use
#'   [data.table::copy()] before subsequently modifying the returned table with
#'   `:=`, or `as.data.frame()` when base-data-frame semantics are preferred.
#' @import data.table
#' @export
assemble_phewas_data <- function(sample, phenotypes = NULL, anchors = NULL,
                                 covariates = NULL, id) {
  id <- .scalar_character(id, "id")
  sample <- .validate_id_table(sample, id, "sample")

  supplied <- list(
    phenotypes = phenotypes,
    anchors = anchors,
    covariates = covariates
  )
  supplied <- supplied[!vapply(supplied, is.null, logical(1))]

  out <- sample
  report <- list(data.frame(
    table = "sample",
    input_rows = nrow(sample),
    matched_input_rows = nrow(sample),
    unmatched_input_rows = 0L,
    unmatched_sample_rows = 0L,
    stringsAsFactors = FALSE
  ))

  for (table_name in names(supplied)) {
    table <- .validate_id_table(supplied[[table_name]], id, table_name)
    same_storage <- identical(typeof(table[[id]]), typeof(sample[[id]]))
    same_class <- identical(class(table[[id]]), class(sample[[id]]))
    same_factor_levels <- !is.factor(sample[[id]]) ||
      identical(levels(table[[id]]), levels(sample[[id]]))
    if (!same_storage || !same_class || !same_factor_levels) {
      stop(
        sprintf("ID column `%s` has a different type in `%s` and `sample`.",
                id, table_name),
        call. = FALSE
      )
    }
    new_columns <- setdiff(names(table), id)
    collisions <- intersect(new_columns, names(out))
    if (length(collisions)) {
      stop(
        sprintf(
          "Column collision while joining `%s`: %s.",
          table_name,
          paste(collisions, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    location <- match(out[[id]], table[[id]])
    if (length(new_columns)) {
      additions <- table[location, new_columns, drop = FALSE]
      row.names(additions) <- NULL
      out[new_columns] <- additions
    }

    input_matched <- table[[id]] %in% sample[[id]]
    report[[length(report) + 1L]] <- data.frame(
      table = table_name,
      input_rows = nrow(table),
      matched_input_rows = sum(input_matched),
      unmatched_input_rows = sum(!input_matched),
      unmatched_sample_rows = sum(is.na(location)),
      stringsAsFactors = FALSE
    )
  }

  out <- data.table::as.data.table(out)
  data.table::setattr(out, "id_column", id)
  data.table::setattr(
    out, "join_report", data.table::as.data.table(do.call(rbind, report))
  )
  out
}
