.metadata_value <- function(row, column) {
  if (!column %in% names(row)) {
    return(NULL)
  }
  value <- row[[column]]
  if (is.list(value) && length(value) == 1L) {
    value <- value[[1L]]
  }
  value
}

.is_missing_metadata <- function(value) {
  is.null(value) || length(value) == 0L ||
    (length(value) == 1L && (is.na(value) ||
      (is.character(value) && !nzchar(trimws(value)))))
}

.metadata_sequence <- function(value, numeric = FALSE) {
  if (is.list(value) && length(value) == 1L) {
    value <- value[[1L]]
  }
  if (.is_missing_metadata(value)) {
    return(if (numeric) numeric() else character())
  }

  if (is.character(value) && length(value) == 1L) {
    text <- trimws(value)
    if (startsWith(text, "[") && endsWith(text, "]") &&
        requireNamespace("jsonlite", quietly = TRUE)) {
      value <- jsonlite::fromJSON(text)
    } else if (grepl("|", text, fixed = TRUE)) {
      value <- strsplit(text, "|", fixed = TRUE)[[1L]]
    }
  }

  if (numeric) {
    if (is.factor(value)) value <- as.character(value)
    suppressWarnings(out <- as.numeric(value))
    if (length(out) != length(value) || any(!is.finite(out))) {
      stop("Ordinal scores must be finite numeric values.", call. = FALSE)
    }
    return(out)
  }
  as.character(value)
}

.metadata_scalar <- function(value) {
  if (is.list(value) && length(value) == 1L) {
    value <- value[[1L]]
  }
  if (.is_missing_metadata(value)) {
    return(NULL)
  }
  if (length(value) != 1L) {
    stop("Metadata fields that are not list-valued must be scalar.", call. = FALSE)
  }
  value
}

.validate_ordinal_declaration <- function(levels, scores = NULL, label,
                                          minimum_levels = 2L,
                                          scores_allowed = TRUE) {
  levels <- .metadata_sequence(levels)
  if (length(levels) < minimum_levels || anyNA(levels) || any(!nzchar(levels)) ||
      anyDuplicated(levels)) {
    stop(
      sprintf(
        "%s must declare at least %d unique, nonmissing ordered levels.",
        label, minimum_levels
      ),
      call. = FALSE
    )
  }

  if (!scores_allowed && !.is_missing_metadata(scores)) {
    stop(sprintf("%s must not declare trend scores.", label), call. = FALSE)
  }
  if (scores_allowed) {
    parsed_scores <- .metadata_sequence(scores, numeric = TRUE)
    if (!length(parsed_scores)) {
      parsed_scores <- seq_along(levels)
    }
    if (length(parsed_scores) != length(levels) || any(diff(parsed_scores) <= 0)) {
      stop(
        sprintf("%s scores must be strictly increasing and match its levels.", label),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

.metadata_alias_values_equal <- function(x, y, numeric = FALSE) {
  if (!is.atomic(x) || !is.atomic(y) || length(x) != length(y)) {
    return(FALSE)
  }

  x_missing <- is.na(x)
  y_missing <- is.na(y)
  equal <- x_missing & y_missing
  observed <- !x_missing & !y_missing
  if (numeric && any(observed)) {
    x_text <- as.character(x[observed])
    y_text <- as.character(y[observed])
    suppressWarnings(x_numeric <- as.numeric(x_text))
    suppressWarnings(y_numeric <- as.numeric(y_text))
    both_numeric <- is.finite(x_numeric) & is.finite(y_numeric)
    observed_equal <- (!both_numeric & x_text == y_text) |
      (both_numeric & x_numeric == y_numeric)
    equal[observed] <- observed_equal
  } else if (any(observed)) {
    equal[observed] <- as.character(x[observed]) == as.character(y[observed])
  }
  all(equal)
}

.normalise_phenotype_metadata <- function(phenotypes, direction) {
  phenotypes <- .plain_data_frame(phenotypes, "phenotypes")
  adopt_alias <- function(target, aliases, numeric = FALSE) {
    sources <- aliases[aliases %in% names(phenotypes)]
    if (!target %in% names(phenotypes) && length(sources)) {
      phenotypes[[target]] <<- phenotypes[[sources[[1L]]]]
    }
    if (target %in% names(phenotypes)) {
      for (source in sources) {
        if (!.metadata_alias_values_equal(
          phenotypes[[target]], phenotypes[[source]], numeric = numeric
        )) {
          stop(
            sprintf(
              "Phenotype metadata `%s` and legacy `%s` columns disagree.",
              target, source
            ),
            call. = FALSE
          )
        }
      }
    }
    phenotypes[sources] <<- NULL
    invisible(NULL)
  }

  # Accept the previous metadata vocabulary only when every supplied name has
  # the same value. Silent precedence would make the association key or plot
  # annotations ambiguous.
  adopt_alias(
    "phenotype", c("phenotype_id", "phecode", "phenotype_column", "column")
  )
  adopt_alias("description", "label")
  adopt_alias("group", "category")
  adopt_alias("groupnum", "category_order", numeric = TRUE)

  required <- c("phenotype", "description", "group", "groupnum", "color")
  if (identical(direction, "phenotypes_as_outcomes")) {
    required <- c(required, "outcome_type")
  }
  missing_columns <- setdiff(required, names(phenotypes))
  if (length(missing_columns)) {
    stop(
      sprintf(
        "`phenotypes` is missing required columns: %s.",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!nrow(phenotypes)) {
    .stop_bad_argument("phenotypes", "must contain at least one row.")
  }

  if (!"variable_type" %in% names(phenotypes)) {
    phenotypes$variable_type <- rep("binary", nrow(phenotypes))
  } else if (is.atomic(phenotypes$variable_type) &&
             length(phenotypes$variable_type) == nrow(phenotypes)) {
    variable_type <- as.character(phenotypes$variable_type)
    missing_variable_type <- vapply(
      variable_type, .is_missing_metadata, logical(1L)
    )
    variable_type[missing_variable_type] <- "binary"
    phenotypes$variable_type <- variable_type
  }

  scalar_text <- c(
    "phenotype", "description", "group", "color", "variable_type"
  )
  for (column in scalar_text) {
    values <- phenotypes[[column]]
    if (!is.atomic(values) || anyNA(values) || any(!nzchar(as.character(values)))) {
      stop(sprintf("Phenotype metadata `%s` must be nonmissing.", column),
           call. = FALSE)
    }
    phenotypes[[column]] <- as.character(values)
  }

  groupnum <- phenotypes$groupnum
  if (is.factor(groupnum)) groupnum <- as.character(groupnum)
  suppressWarnings(groupnum <- as.numeric(groupnum))
  if (any(!is.finite(groupnum))) {
    stop("Phenotype metadata `groupnum` must be finite numeric values.",
         call. = FALSE)
  }
  phenotypes$groupnum <- groupnum

  valid_color <- vapply(phenotypes$color, function(value) {
    tryCatch(
      {
        grDevices::col2rgb(value)
        TRUE
      },
      error = function(error) FALSE
    )
  }, logical(1L))
  if (any(!valid_color)) {
    stop(
      sprintf(
        "Phenotype metadata `color` contains invalid R color value(s): %s.",
        paste(unique(phenotypes$color[!valid_color]), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  group_rows <- split(seq_len(nrow(phenotypes)), phenotypes$group)
  inconsistent_groupnum <- names(group_rows)[vapply(
    group_rows,
    function(rows) length(unique(phenotypes$groupnum[rows])) != 1L,
    logical(1L)
  )]
  if (length(inconsistent_groupnum)) {
    stop(
      sprintf(
        "Each `group` must have one `groupnum`; inconsistent group(s): %s.",
        paste(inconsistent_groupnum, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  inconsistent_color <- names(group_rows)[vapply(
    group_rows,
    function(rows) length(unique(phenotypes$color[rows])) != 1L,
    logical(1L)
  )]
  if (length(inconsistent_color)) {
    stop(
      sprintf(
        "Each `group` must have one `color`; inconsistent group(s): %s.",
        paste(inconsistent_color, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  group_numbers <- vapply(
    group_rows, function(rows) phenotypes$groupnum[rows[[1L]]], numeric(1L)
  )
  duplicated_group_numbers <- unique(group_numbers[
    duplicated(group_numbers) | duplicated(group_numbers, fromLast = TRUE)
  ])
  if (length(duplicated_group_numbers)) {
    stop(
      sprintf(
        "`groupnum` must uniquely order groups; duplicated value(s): %s.",
        paste(duplicated_group_numbers, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(phenotypes$phenotype)) {
    stop("`phenotype` values must be unique within a specification.",
         call. = FALSE)
  }

  optional <- c("reference", "levels", "scores", "offset")
  for (column in setdiff(optional, names(phenotypes))) {
    phenotypes[[column]] <- rep(NA_character_, nrow(phenotypes))
  }

  allowed_variable_types <- c("numeric", "binary", "ordinal", "count")
  if (any(!phenotypes$variable_type %in% allowed_variable_types)) {
    stop(
      "Phenotype `variable_type` must be one of numeric, binary, ordinal, or count.",
      call. = FALSE
    )
  }

  binary_rows <- phenotypes$variable_type == "binary"
  reference <- phenotypes$reference
  if (is.atomic(reference) && length(reference) == nrow(phenotypes)) {
    reference <- as.character(reference)
    missing_reference <- vapply(
      reference, .is_missing_metadata, logical(1L)
    )
    reference[binary_rows & missing_reference] <- "0"
    phenotypes$reference <- reference
  } else if (is.list(reference) && length(reference) == nrow(phenotypes)) {
    for (i in which(binary_rows)) {
      if (.is_missing_metadata(reference[[i]])) {
        reference[[i]] <- "0"
      }
    }
    phenotypes$reference <- reference
  }

  if (identical(direction, "phenotypes_as_outcomes")) {
    phenotypes$outcome_type <- as.character(phenotypes$outcome_type)
    allowed_outcomes <- c("binary", "continuous", "count", "ordinal")
    if (anyNA(phenotypes$outcome_type) ||
        any(!phenotypes$outcome_type %in% allowed_outcomes)) {
      stop(
        "Phenotype `outcome_type` must be one of binary, continuous, count, or ordinal.",
        call. = FALSE
      )
    }
    required_variable_type <- c(
      binary = "binary", continuous = "numeric", count = "count",
      ordinal = "ordinal"
    )[phenotypes$outcome_type]
    if (any(phenotypes$variable_type != unname(required_variable_type))) {
      stop("Phenotype `variable_type` and `outcome_type` declarations disagree.",
           call. = FALSE)
    }
  } else {
    if (any(!phenotypes$variable_type %in% c("numeric", "binary", "ordinal"))) {
      stop(
        "Scanned predictors must have variable_type numeric, binary, or ordinal.",
        call. = FALSE
      )
    }
  }

  for (i in seq_len(nrow(phenotypes))) {
    row <- phenotypes[i, , drop = FALSE]
    variable_type <- row$variable_type[[1L]]
    outcome_type <- if ("outcome_type" %in% names(row)) {
      row$outcome_type[[1L]]
    } else {
      NULL
    }
    is_binary <- identical(variable_type, "binary")
    is_ordinal <- identical(variable_type, "ordinal")

    if (is_binary && .is_missing_metadata(.metadata_value(row, "reference"))) {
      stop(sprintf("Phenotype `%s` must declare a binary reference level.",
                   row$phenotype[[1L]]), call. = FALSE)
    }
    if (is_binary) {
      .metadata_scalar(.metadata_value(row, "reference"))
    }
    if (is_ordinal) {
      predictor_ordinal <- identical(direction, "phenotypes_as_predictors")
      .validate_ordinal_declaration(
        .metadata_value(row, "levels"),
        .metadata_value(row, "scores"),
        sprintf("Phenotype `%s`", row$phenotype[[1L]]),
        minimum_levels = if (identical(outcome_type, "ordinal")) 3L else 2L,
        scores_allowed = predictor_ordinal
      )
    }

    offset <- .metadata_scalar(.metadata_value(row, "offset"))
    if (!is.null(offset) &&
        (!is.character(offset) || is.na(offset) || !nzchar(offset))) {
      stop(sprintf("Phenotype `%s` has an invalid offset column name.",
                   row$phenotype[[1L]]), call. = FALSE)
    }
    if (!is.null(offset) && !identical(outcome_type, "count")) {
      stop(sprintf("Only count outcome `%s` may declare an offset.",
                   row$phenotype[[1L]]), call. = FALSE)
    }
  }

  canonical_order <- c(
    "phenotype", "description", "group", "groupnum", "color",
    "variable_type", "outcome_type", "reference", "levels", "scores",
    "offset"
  )
  phenotypes <- phenotypes[c(
    intersect(canonical_order, names(phenotypes)),
    setdiff(names(phenotypes), canonical_order)
  )]
  row.names(phenotypes) <- NULL
  phenotypes
}

.validate_exclusion <- function(exclusion_column, exclusion_values, data = NULL) {
  if (is.null(exclusion_column)) {
    if (!is.null(exclusion_values) && length(exclusion_values)) {
      stop("`exclusion_values` requires an `exclusion_column`.", call. = FALSE)
    }
    return(list(column = NULL, values = NULL))
  }
  exclusion_column <- .scalar_character(exclusion_column, "exclusion_column")
  if (is.factor(exclusion_values)) {
    exclusion_values <- as.character(exclusion_values)
  }
  supported <- is.character(exclusion_values) || is.numeric(exclusion_values) ||
    is.logical(exclusion_values)
  if (!supported || !length(exclusion_values) || anyNA(exclusion_values) ||
      anyDuplicated(exclusion_values)) {
    stop(
      "`exclusion_values` must be unique, nonmissing character, numeric, or logical values.",
      call. = FALSE
    )
  }
  if (is.character(exclusion_values) && any(!nzchar(exclusion_values))) {
    stop("Character `exclusion_values` must not be blank.", call. = FALSE)
  }
  if (is.numeric(exclusion_values) && any(!is.finite(exclusion_values))) {
    stop("Numeric `exclusion_values` must be finite.", call. = FALSE)
  }
  exclusion_values <- unname(exclusion_values)

  if (!is.null(data)) {
    if (!exclusion_column %in% names(data)) {
      stop(sprintf("Analysis data is missing exclusion column `%s`.",
                   exclusion_column), call. = FALSE)
    }
    marker <- data[[exclusion_column]]
    compatible <- if (is.factor(marker) || is.character(marker)) {
      is.character(exclusion_values)
    } else if (is.numeric(marker)) {
      is.numeric(exclusion_values)
    } else if (is.logical(marker)) {
      is.logical(exclusion_values)
    } else {
      FALSE
    }
    if (!compatible) {
      stop(
        sprintf(
          "Exclusion values are not type-compatible with column `%s`.",
          exclusion_column
        ),
        call. = FALSE
      )
    }
  }
  list(column = exclusion_column, values = exclusion_values)
}

.exclusion_mask <- function(data, exclusion_column, exclusion_values) {
  if (is.null(exclusion_column)) {
    return(rep(FALSE, nrow(data)))
  }
  declaration <- .validate_exclusion(exclusion_column, exclusion_values, data)
  marker <- data[[declaration$column]]
  if (is.factor(marker)) {
    marker <- as.character(marker)
  }
  !is.na(marker) & marker %in% declaration$values
}

.validate_eligibility <- function(eligibility) {
  if (!is.list(eligibility) || is.null(names(eligibility)) ||
      anyNA(names(eligibility)) || any(!nzchar(names(eligibility))) ||
      anyDuplicated(names(eligibility))) {
    .stop_bad_argument("eligibility", "must be an explicitly named list.")
  }
  required <- c("min_n", "min_cases", "min_controls", "min_outcome_levels")
  missing <- setdiff(required, names(eligibility))
  if (length(missing)) {
    .stop_bad_argument(
      "eligibility",
      sprintf("must explicitly define %s.", paste(missing, collapse = ", "))
    )
  }
  unknown <- setdiff(names(eligibility), required)
  if (length(unknown)) {
    .stop_bad_argument(
      "eligibility",
      sprintf("contains unknown thresholds: %s.", paste(unknown, collapse = ", "))
    )
  }
  out <- eligibility[required]
  for (name in required) {
    value <- out[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value > .Machine$integer.max || value < -.Machine$integer.max ||
        value != trunc(value)) {
      stop(sprintf("Eligibility threshold `%s` must be one integer.", name),
           call. = FALSE)
    }
    out[[name]] <- as.integer(value)
  }
  if (out$min_n < 2L || out$min_cases < 1L || out$min_controls < 1L ||
      out$min_outcome_levels < 3L) {
    stop(
      "Eligibility requires min_n >= 2, min_cases/min_controls >= 1, and min_outcome_levels >= 3.",
      call. = FALSE
    )
  }
  out
}

#' Define a PheWAS analysis specification
#'
#' A specification describes exactly one fixed anchor and one family of
#' phenotype tests. Metadata may use native list-columns for `levels` and
#' `scores`, or `|`-delimited strings when it comes from TSV/YAML files.
#' The canonical metadata vocabulary is `phenotype`, `description`, `group`,
#' `groupnum`, and `color`; `phenotype` names the analysis-data column. Legacy
#' `phenotype_id`/`phecode`/`phenotype_column`/`column`, `label`, `category`,
#' and `category_order` fields are accepted only when they agree with any
#' canonical fields supplied in the same table, and are not retained in the
#' validated specification.
#'
#' @param analysis_id Stable, unique name for the analysis.
#' @param anchor Name of the fixed predictor (phenotypes as outcomes) or fixed
#'   outcome (phenotypes as predictors).
#' @param direction Either `"phenotypes_as_outcomes"` or
#'   `"phenotypes_as_predictors"`.
#' @param phenotypes Data frame with phenotype metadata. Required columns are
#'   `phenotype`, `description`, `group`, `groupnum`, and `color`. Missing or
#'   blank `variable_type` values default to `"binary"`; binary phenotypes with
#'   a missing or blank `reference` default to `"0"`. Outcome scans additionally
#'   require `outcome_type`.
#' @param covariates Character vector of adjustment-variable columns.
#' @param outcome_type Fixed-anchor outcome type for phenotype-predictor scans.
#' @param anchor_type Fixed-anchor type. Required for phenotype-outcome scans.
#' @param anchor_transform Either `"none"` or `"zscore"`.
#' @param anchor_reference Explicit reference for a binary anchor.
#' @param anchor_levels Explicit ordered levels for an ordinal anchor.
#' @param anchor_scores Optional strictly increasing ordinal-trend scores.
#' @param offset Optional log-offset column for a fixed count outcome.
#' @param exclusion_column Optional cohort marker column. Matching rows are
#'   removed before base-cohort completeness and anchor transformation.
#' @param exclusion_values Explicit marker values to remove. Missing marker
#'   values are retained.
#' @param eligibility Explicit named list containing `min_n`, `min_cases`,
#'   `min_controls`, and `min_outcome_levels`. `min_cases` and `min_controls`
#'   set the minimum non-reference and reference counts for each scanned binary
#'   phenotype, whether it is modeled as an outcome or predictor. They also
#'   protect event and nonevent counts when the fixed outcome is binary.
#' @param fdr_threshold Threshold used to flag BH-adjusted results.
#' @param version Specification schema version.
#'
#' @return A validated object of class `phewas_spec`.
#' @export
phewas_spec <- function(analysis_id, anchor,
                        direction = c("phenotypes_as_outcomes",
                                      "phenotypes_as_predictors"),
                        phenotypes, covariates = character(),
                        outcome_type = NULL, anchor_type = NULL,
                        anchor_transform = c("none", "zscore"),
                        anchor_reference = NULL, anchor_levels = NULL,
                        anchor_scores = NULL, offset = NULL,
                        exclusion_column = NULL, exclusion_values = NULL,
                        eligibility,
                        fdr_threshold, version = "1.0") {
  if (missing(fdr_threshold)) {
    .stop_bad_argument("fdr_threshold", "must be explicitly supplied.")
  }
  analysis_id <- .scalar_character(analysis_id, "analysis_id")
  anchor <- .scalar_character(anchor, "anchor")
  direction <- match.arg(direction)
  anchor_transform <- match.arg(anchor_transform)
  version <- .scalar_character(version, "version")
  if (!identical(version, "1.0")) {
    stop(sprintf("Unsupported PheWAS specification version `%s`.", version),
         call. = FALSE)
  }

  if (!is.character(covariates) || anyNA(covariates) ||
      any(!nzchar(covariates)) || anyDuplicated(covariates)) {
    .stop_bad_argument("covariates", "must contain unique, nonmissing column names.")
  }
  if (anchor %in% covariates) {
    stop("The anchor cannot also be a covariate.", call. = FALSE)
  }
  if (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1L ||
      is.na(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold >= 1) {
    .stop_bad_argument("fdr_threshold", "must be one number strictly between 0 and 1.")
  }

  phenotypes <- .normalise_phenotype_metadata(phenotypes, direction)
  if (anchor %in% phenotypes$phenotype) {
    stop("The anchor cannot also be a scanned phenotype column.", call. = FALSE)
  }
  overlap <- intersect(covariates, phenotypes$phenotype)
  if (length(overlap)) {
    stop(sprintf("Covariates cannot also be scanned phenotypes: %s.",
                 paste(overlap, collapse = ", ")), call. = FALSE)
  }
  exclusion <- .validate_exclusion(exclusion_column, exclusion_values)

  allowed_outcomes <- c("binary", "continuous", "count", "ordinal")
  allowed_predictors <- c("numeric", "binary", "ordinal")
  if (identical(direction, "phenotypes_as_outcomes")) {
    if (is.null(anchor_type) || length(anchor_type) != 1L ||
        !anchor_type %in% allowed_predictors) {
      .stop_bad_argument(
        "anchor_type",
        "must be explicitly set to numeric, binary, or ordinal for an outcome scan."
      )
    }
    if (!is.null(outcome_type)) {
      stop("`outcome_type` is row-specific for a phenotype-outcome scan.",
           call. = FALSE)
    }
    if (!is.null(offset)) {
      stop("Use phenotype metadata `offset` for a phenotype-outcome scan.",
           call. = FALSE)
    }
    if (identical(anchor_type, "binary") && .is_missing_metadata(anchor_reference)) {
      stop("A binary anchor must declare `anchor_reference`.", call. = FALSE)
    }
    if (identical(anchor_type, "binary")) {
      .metadata_scalar(anchor_reference)
    }
    if (identical(anchor_type, "ordinal")) {
      .validate_ordinal_declaration(anchor_levels, anchor_scores,
                                    "The ordinal anchor", 2L, TRUE)
    }
    if (identical(anchor_transform, "zscore") &&
        !identical(anchor_type, "numeric")) {
      stop("Only a numeric predictor anchor can be z-standardized.", call. = FALSE)
    }
  } else {
    if (is.null(outcome_type) || length(outcome_type) != 1L ||
        !outcome_type %in% allowed_outcomes) {
      .stop_bad_argument(
        "outcome_type",
        "must explicitly identify the fixed outcome as binary, continuous, count, or ordinal."
      )
    }
    expected_anchor_type <- c(
      binary = "binary", continuous = "numeric", count = "count",
      ordinal = "ordinal"
    )[[outcome_type]]
    if (is.null(anchor_type)) {
      anchor_type <- expected_anchor_type
    } else if (length(anchor_type) != 1L ||
               !identical(anchor_type, expected_anchor_type)) {
      stop("`anchor_type` disagrees with the declared fixed `outcome_type`.",
           call. = FALSE)
    }
    if (identical(outcome_type, "binary") && .is_missing_metadata(anchor_reference)) {
      stop("A binary fixed outcome must declare `anchor_reference`.", call. = FALSE)
    }
    if (identical(outcome_type, "binary")) {
      .metadata_scalar(anchor_reference)
    }
    if (identical(outcome_type, "ordinal")) {
      .validate_ordinal_declaration(anchor_levels, anchor_scores,
                                    "The ordinal fixed outcome", 3L, FALSE)
    }
    if (!is.null(offset) && !identical(outcome_type, "count")) {
      stop("Only a count fixed outcome may declare `offset`.", call. = FALSE)
    }
    if (!is.null(offset) &&
        (!is.character(offset) || length(offset) != 1L || is.na(offset) ||
         !nzchar(offset))) {
      stop("`offset` must be one nonmissing column name.", call. = FALSE)
    }
    if (identical(anchor_transform, "zscore") &&
        !identical(outcome_type, "continuous")) {
      stop("Only a continuous fixed outcome can be z-standardized.", call. = FALSE)
    }
  }

  spec <- list(
    version = version,
    analysis_id = analysis_id,
    anchor = anchor,
    direction = direction,
    phenotypes = phenotypes,
    covariates = covariates,
    outcome_type = outcome_type,
    anchor_type = anchor_type,
    anchor_transform = anchor_transform,
    anchor_reference = anchor_reference,
    anchor_levels = anchor_levels,
    anchor_scores = anchor_scores,
    offset = offset,
    exclusion_column = exclusion$column,
    exclusion_values = exclusion$values,
    eligibility = .validate_eligibility(eligibility),
    fdr_threshold = as.numeric(fdr_threshold)
  )
  class(spec) <- "phewas_spec"
  spec
}

#' Validate a PheWAS specification
#'
#' @param spec A `phewas_spec` object.
#' @param data Optional analysis data used to verify all declared columns.
#'
#' @return `spec`, validated and safe to use.
#' @export
validate_phewas_spec <- function(spec, data = NULL) {
  if (!inherits(spec, "phewas_spec") || !is.list(spec)) {
    .stop_bad_argument("spec", "must be created by phewas_spec().")
  }
  required_fields <- c(
    "version", "analysis_id", "anchor", "direction", "phenotypes",
    "covariates", "anchor_type", "anchor_transform", "eligibility",
    "fdr_threshold", "exclusion_column", "exclusion_values"
  )
  if (length(setdiff(required_fields, names(spec)))) {
    stop("The specification is incomplete or uses an unsupported schema.",
         call. = FALSE)
  }
  if (!identical(spec$version, "1.0")) {
    stop(sprintf("Unsupported PheWAS specification version `%s`.", spec$version),
         call. = FALSE)
  }

  # Specifications are mutable lists. Reapply every constructor invariant so
  # edits cannot silently change the direction, encoding, or requested model.
  arguments <- intersect(names(formals(phewas_spec)), names(spec))
  validated <- do.call(phewas_spec, unclass(spec[arguments]))
  spec[names(validated)] <- validated

  if (!is.null(data)) {
    data <- .plain_data_frame(data, "data")
    .validate_exclusion(spec$exclusion_column, spec$exclusion_values, data)
    offset_columns <- character()
    if (identical(spec$direction, "phenotypes_as_outcomes")) {
      for (i in seq_len(nrow(spec$phenotypes))) {
        value <- .metadata_scalar(.metadata_value(
          spec$phenotypes[i, , drop = FALSE], "offset"
        ))
        if (!is.null(value)) {
          offset_columns <- c(offset_columns, as.character(value))
        }
      }
    } else if (!is.null(spec$offset)) {
      offset_columns <- as.character(spec$offset)
    }
    declared <- unique(c(
      spec$anchor, spec$covariates, spec$phenotypes$phenotype, offset_columns
    ))
    missing <- setdiff(declared, names(data))
    if (length(missing)) {
      stop(sprintf("Analysis data is missing declared columns: %s.",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
  }
  spec
}

#' @export
print.phewas_spec <- function(x, ...) {
  cat(sprintf(
    "<phewas_spec %s> %s: %d phenotypes, %s\n",
    x$version, x$analysis_id, nrow(x$phenotypes), x$direction
  ))
  invisible(x)
}

.association_variables <- function(spec, phenotype_row) {
  if (!is.data.frame(phenotype_row) || nrow(phenotype_row) != 1L) {
    stop("`phenotype_row` must be exactly one metadata row.", call. = FALSE)
  }
  forward <- identical(spec$direction, "phenotypes_as_outcomes")
  offset <- if (forward) {
    .metadata_scalar(.metadata_value(phenotype_row, "offset"))
  } else {
    spec$offset
  }
  list(
    response = if (forward) phenotype_row$phenotype[[1L]] else spec$anchor,
    predictor = if (forward) spec$anchor else phenotype_row$phenotype[[1L]],
    covariates = spec$covariates,
    offset = if (is.null(offset)) NULL else as.character(offset),
    complete_case_variables = unique(c(
      if (forward) phenotype_row$phenotype[[1L]] else spec$anchor,
      if (forward) spec$anchor else phenotype_row$phenotype[[1L]],
      spec$covariates,
      if (is.null(offset)) character() else as.character(offset)
    )),
    outcome_type = if (forward) {
      phenotype_row$outcome_type[[1L]]
    } else {
      spec$outcome_type
    },
    predictor_type = if (forward) spec$anchor_type else phenotype_row$variable_type[[1L]],
    outcome_reference = if (forward) {
      .metadata_scalar(.metadata_value(phenotype_row, "reference"))
    } else {
      spec$anchor_reference
    },
    outcome_levels = if (forward) {
      .metadata_sequence(.metadata_value(phenotype_row, "levels"))
    } else {
      .metadata_sequence(spec$anchor_levels)
    },
    predictor_reference = if (forward) {
      spec$anchor_reference
    } else {
      .metadata_scalar(.metadata_value(phenotype_row, "reference"))
    },
    predictor_levels = if (forward) {
      .metadata_sequence(spec$anchor_levels)
    } else {
      .metadata_sequence(.metadata_value(phenotype_row, "levels"))
    },
    predictor_scores = if (forward) {
      .metadata_sequence(spec$anchor_scores, numeric = TRUE)
    } else {
      .metadata_sequence(.metadata_value(phenotype_row, "scores"), numeric = TRUE)
    }
  )
}
