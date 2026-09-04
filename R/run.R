.normalise_row_index <- function(row_index, n) {
  if (is.null(row_index)) {
    return(seq_len(n))
  }
  if (is.logical(row_index)) {
    if (length(row_index) != n || anyNA(row_index)) {
      stop("A logical `row_index` must be nonmissing and match nrow(data).",
           call. = FALSE)
    }
    return(which(row_index))
  }
  if (!is.numeric(row_index) || anyNA(row_index) ||
      any(row_index != as.integer(row_index)) || any(row_index < 1L) ||
      any(row_index > n) || anyDuplicated(row_index)) {
    stop("A numeric `row_index` must contain unique valid row numbers.",
         call. = FALSE)
  }
  as.integer(row_index)
}

.prepare_analysis_data <- function(data, spec, row_index = NULL) {
  if (isTRUE(attr(data, "phewas_prepared")) && is.null(row_index)) {
    return(data)
  }
  spec <- validate_phewas_spec(spec, data)
  id_column <- attr(data, "id_column")
  data <- .plain_data_frame(data, "data")
  selected <- .normalise_row_index(row_index, nrow(data))
  data <- data[selected, , drop = FALSE]

  n_before_exclusion <- nrow(data)
  excluded <- .exclusion_mask(
    data, spec$exclusion_column, spec$exclusion_values
  )
  n_excluded <- sum(excluded)
  if (n_excluded) {
    data <- data[!excluded, , drop = FALSE]
    selected <- selected[!excluded]
  }
  n_after_exclusion <- nrow(data)

  base_variables <- unique(c(
    spec$anchor,
    spec$covariates,
    if (identical(spec$direction, "phenotypes_as_predictors") &&
        !is.null(spec$offset)) spec$offset else character()
  ))
  base_complete <- stats::complete.cases(data[, base_variables, drop = FALSE])
  source_index <- selected[base_complete]
  data <- data[base_complete, , drop = FALSE]
  row.names(data) <- NULL

  transform_mean <- NA_real_
  transform_sd <- NA_real_
  if (identical(spec$anchor_transform, "zscore")) {
    anchor <- data[[spec$anchor]]
    if (!is.numeric(anchor) || any(!is.finite(anchor))) {
      stop("A z-standardized anchor must be finite and numeric on the base cohort.",
           call. = FALSE)
    }
    if (length(anchor)) {
      transform_mean <- mean(anchor)
    }
    if (length(anchor) >= 2L) {
      transform_sd <- stats::sd(anchor)
      if (!is.finite(transform_sd) || transform_sd <= 0) {
        stop("A z-standardized anchor must have positive finite variation.",
             call. = FALSE)
      }
      data[[spec$anchor]] <- (anchor - transform_mean) / transform_sd
    }
  }

  sample_keys <- source_index
  if (!is.null(id_column) && id_column %in% names(data)) {
    sample_keys <- as.character(data[[id_column]])
  }
  attr(data, "phewas_prepared") <- TRUE
  attr(data, "source_index") <- source_index
  attr(data, "sample_keys") <- sample_keys
  attr(data, "n_base") <- nrow(data)
  attr(data, "n_before_exclusion") <- n_before_exclusion
  attr(data, "n_excluded") <- n_excluded
  attr(data, "n_after_exclusion") <- n_after_exclusion
  attr(data, "transformation_mean") <- transform_mean
  attr(data, "transformation_sd") <- transform_sd
  data
}

.subset_prepared_data <- function(data, row_index) {
  selected <- .normalise_row_index(row_index, nrow(data))
  source_index <- attr(data, "source_index")[selected]
  sample_keys <- attr(data, "sample_keys")[selected]
  n_base <- attr(data, "n_base")
  transform_mean <- attr(data, "transformation_mean")
  transform_sd <- attr(data, "transformation_sd")
  n_before_exclusion <- attr(data, "n_before_exclusion")
  n_excluded <- attr(data, "n_excluded")
  n_after_exclusion <- attr(data, "n_after_exclusion")
  data <- data[selected, , drop = FALSE]
  row.names(data) <- NULL
  attr(data, "phewas_prepared") <- TRUE
  attr(data, "source_index") <- source_index
  attr(data, "sample_keys") <- sample_keys
  attr(data, "n_base") <- n_base
  attr(data, "n_before_exclusion") <- n_before_exclusion
  attr(data, "n_excluded") <- n_excluded
  attr(data, "n_after_exclusion") <- n_after_exclusion
  attr(data, "transformation_mean") <- transform_mean
  attr(data, "transformation_sd") <- transform_sd
  data
}

.collapse_metadata <- function(value, numeric = FALSE) {
  values <- .metadata_sequence(value, numeric = numeric)
  if (!length(values)) NA_character_ else paste(values, collapse = "|")
}

.association_engine <- function(outcome_type) {
  unname(c(
    binary = "logistf", continuous = "lm", count = "glm.nb", ordinal = "clm"
  )[[outcome_type]])
}

.association_effect_measure <- function(outcome_type) {
  unname(c(
    binary = "odds_ratio",
    continuous = "mean_difference",
    count = "incidence_rate_ratio",
    ordinal = "common_odds_ratio"
  )[[outcome_type]])
}

.scalar_text_or_na <- function(value) {
  value <- .metadata_scalar(value)
  if (is.null(value)) NA_character_ else as.character(value)
}

.empty_result_row <- function(spec, phenotype_row, variables, prepared_data) {
  outcome_levels <- .collapse_metadata(variables$outcome_levels)
  predictor_levels <- .collapse_metadata(variables$predictor_levels)
  predictor_scores <- variables$predictor_scores
  if (identical(variables$predictor_type, "ordinal") &&
      !length(predictor_scores)) {
    predictor_scores <- seq_along(variables$predictor_levels)
  }

  data.frame(
    analysis_id = spec$analysis_id,
    spec_version = spec$version,
    phenotype = as.character(phenotype_row$phenotype[[1L]]),
    description = as.character(phenotype_row$description[[1L]]),
    group = as.character(phenotype_row$group[[1L]]),
    groupnum = as.numeric(phenotype_row$groupnum[[1L]]),
    color = as.character(phenotype_row$color[[1L]]),
    direction = spec$direction,
    response = variables$response,
    predictor = variables$predictor,
    engine = .association_engine(variables$outcome_type),
    outcome_type = variables$outcome_type,
    predictor_type = variables$predictor_type,
    transformation = spec$anchor_transform,
    transformation_mean = as.numeric(attr(prepared_data, "transformation_mean")),
    transformation_sd = as.numeric(attr(prepared_data, "transformation_sd")),
    transformation_source_n = as.integer(attr(prepared_data, "n_base")),
    reference = .scalar_text_or_na(variables$outcome_reference),
    outcome_reference = .scalar_text_or_na(variables$outcome_reference),
    predictor_reference = .scalar_text_or_na(variables$predictor_reference),
    levels = if (identical(spec$direction, "phenotypes_as_outcomes")) {
      outcome_levels
    } else {
      predictor_levels
    },
    outcome_levels = outcome_levels,
    predictor_levels = predictor_levels,
    scores = if (!length(predictor_scores)) NA_character_ else
      paste(predictor_scores, collapse = "|"),
    offset = if (is.null(variables$offset)) NA_character_ else variables$offset,
    effect_measure = .association_effect_measure(variables$outcome_type),
    formula = .display_formula(variables),
    effective_formula = NA_character_,
    covariates = if (!length(variables$covariates)) NA_character_ else
      paste(variables$covariates, collapse = "|"),
    effective_covariates = NA_character_,
    dropped_invariant_covariates = NA_character_,
    exclusion_active = !is.null(spec$exclusion_column),
    exclusion_column = if (is.null(spec$exclusion_column)) NA_character_ else
      spec$exclusion_column,
    exclusion_values = if (is.null(spec$exclusion_values)) NA_character_ else
      paste(spec$exclusion_values, collapse = "|"),
    n_before_exclusion = as.integer(attr(prepared_data, "n_before_exclusion")),
    n_excluded = as.integer(attr(prepared_data, "n_excluded")),
    n_after_exclusion = as.integer(attr(prepared_data, "n_after_exclusion")),
    status = "error",
    reason_code = "not_run",
    message = NA_character_,
    n_base = as.integer(attr(prepared_data, "n_base")),
    n_complete = NA_integer_,
    n_event = NA_integer_,
    n_nonevent = NA_integer_,
    outcome_counts = NA_character_,
    predictor_counts = NA_character_,
    estimate = NA_real_,
    std_error = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    statistic = NA_real_,
    test_method = NA_character_,
    df = NA_real_,
    native_effect = NA_real_,
    native_conf_low = NA_real_,
    native_conf_high = NA_real_,
    log_p = NA_real_,
    p_value = NA_real_,
    p_value_text = NA_character_,
    neg_log10_p = NA_real_,
    q_value = NA_real_,
    q_value_text = NA_character_,
    neg_log10_q = NA_real_,
    bonferroni_p = NA_real_,
    bonferroni_p_text = NA_character_,
    neg_log10_bonferroni_p = NA_real_,
    fdr_significant = NA,
    bonferroni_significant = NA,
    theta = NA_real_,
    ordinal_levels = NA_character_,
    converged = NA,
    iterations = NA_integer_,
    warnings = NA_character_,
    sample_fingerprint = NA_character_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.format_counts <- function(x) {
  counts <- table(x, useNA = "no")
  if (!length(counts)) {
    return(NA_character_)
  }
  paste(paste0(names(counts), "=", as.integer(counts)), collapse = "|")
}

.set_result_problem <- function(row, status, reason_code, message,
                                warning_messages = character()) {
  row$status <- status
  row$reason_code <- reason_code
  row$message <- message
  if (length(warning_messages)) {
    row$warnings <- paste(unique(warning_messages), collapse = " | ")
  }
  row
}

.fit_one_association <- function(data, spec, phenotype_row, row_index = NULL) {
  if (!is.data.frame(phenotype_row) || nrow(phenotype_row) != 1L) {
    stop("`phenotype_row` must be exactly one metadata row.", call. = FALSE)
  }
  if (isTRUE(attr(data, "phewas_prepared")) && !is.null(row_index)) {
    data <- .subset_prepared_data(data, row_index)
  } else if (!isTRUE(attr(data, "phewas_prepared"))) {
    data <- .prepare_analysis_data(data, spec, row_index = row_index)
  }
  variables <- .association_variables(spec, phenotype_row)
  result <- .empty_result_row(spec, phenotype_row, variables, data)

  complete <- stats::complete.cases(
    data[, variables$complete_case_variables, drop = FALSE]
  )
  complete_positions <- which(complete)
  association_data <- data[complete_positions, , drop = FALSE]
  source_index <- attr(data, "source_index")[complete_positions]
  sample_keys <- attr(data, "sample_keys")[complete_positions]
  result$n_complete <- nrow(association_data)
  result$sample_fingerprint <- digest::digest(sample_keys, algo = "sha256",
                                               serialize = TRUE)

  if (nrow(association_data) < spec$eligibility$min_n) {
    result <- .set_result_problem(
      result, "skipped", "insufficient_n",
      sprintf("Only %d complete cases; min_n is %d.", nrow(association_data),
              spec$eligibility$min_n)
    )
    attr(result, "fit") <- NULL
    attr(result, "sample_index") <- source_index
    return(result)
  }

  model_data <- tryCatch(
    .prepare_model_data(association_data, variables),
    phewas_skip = function(error) error,
    error = function(error) error
  )
  if (inherits(model_data, "error")) {
    is_skip <- inherits(model_data, "phewas_skip")
    result <- .set_result_problem(
      result,
      if (is_skip) "skipped" else "error",
      if (is_skip) model_data$reason_code else "data_error",
      conditionMessage(model_data)
    )
    attr(result, "fit") <- NULL
    attr(result, "sample_index") <- source_index
    return(result)
  }

  effective_covariates <- attr(model_data, "effective_covariates")
  dropped_invariant_covariates <- attr(
    model_data, "dropped_invariant_covariates"
  )
  result$effective_formula <- .display_formula(
    variables, covariates = effective_covariates
  )
  result$effective_covariates <- if (!length(effective_covariates)) {
    NA_character_
  } else {
    paste(effective_covariates, collapse = "|")
  }
  result$dropped_invariant_covariates <- if (
    !length(dropped_invariant_covariates)
  ) {
    NA_character_
  } else {
    paste(dropped_invariant_covariates, collapse = "|")
  }

  if (variables$outcome_type %in% c("binary", "ordinal")) {
    result$outcome_counts <- .format_counts(
      association_data[[variables$response]]
    )
  }
  if (variables$predictor_type %in% c("binary", "ordinal")) {
    result$predictor_counts <- .format_counts(
      association_data[[variables$predictor]]
    )
  }
  if (identical(variables$outcome_type, "binary")) {
    result$n_event <- sum(model_data$.response == 1L)
    result$n_nonevent <- sum(model_data$.response == 0L)
    if (result$n_event < spec$eligibility$min_cases) {
      result <- .set_result_problem(
        result, "skipped", "insufficient_cases",
        sprintf("Only %d cases; min_cases is %d.", result$n_event,
                spec$eligibility$min_cases)
      )
    } else if (result$n_nonevent < spec$eligibility$min_controls) {
      result <- .set_result_problem(
        result, "skipped", "insufficient_controls",
        sprintf("Only %d controls; min_controls is %d.", result$n_nonevent,
                spec$eligibility$min_controls)
      )
    }
  } else if (identical(variables$outcome_type, "ordinal")) {
    observed_levels <- sum(table(model_data$.response) > 0)
    if (observed_levels < spec$eligibility$min_outcome_levels) {
      result <- .set_result_problem(
        result, "skipped", "insufficient_outcome_levels",
        sprintf("Only %d ordinal outcome levels; the declared minimum is %d.",
                observed_levels, spec$eligibility$min_outcome_levels)
      )
    }
  } else if (length(unique(model_data$.response)) < 2L) {
    result <- .set_result_problem(
      result, "skipped", "constant_outcome",
      "The outcome is constant on complete cases."
    )
  }
  if (identical(result$reason_code, "not_run") &&
      identical(spec$direction, "phenotypes_as_predictors") &&
      identical(variables$predictor_type, "binary")) {
    nonreference_n <- sum(model_data$.predictor == 1L)
    reference_n <- sum(model_data$.predictor == 0L)
    if (nonreference_n < spec$eligibility$min_cases) {
      result <- .set_result_problem(
        result, "skipped", "insufficient_cases",
        sprintf(
          paste0(
            "Only %d non-reference observations for scanned binary ",
            "phenotype `%s`; min_cases is %d."
          ),
          nonreference_n, phenotype_row$phenotype[[1L]],
          spec$eligibility$min_cases
        )
      )
    } else if (reference_n < spec$eligibility$min_controls) {
      result <- .set_result_problem(
        result, "skipped", "insufficient_controls",
        sprintf(
          paste0(
            "Only %d reference observations for scanned binary phenotype ",
            "`%s`; min_controls is %d."
          ),
          reference_n, phenotype_row$phenotype[[1L]],
          spec$eligibility$min_controls
        )
      )
    }
  }
  if (identical(result$reason_code, "not_run") &&
      length(unique(model_data$.predictor)) < 2L) {
    result <- .set_result_problem(
      result, "skipped", "constant_predictor",
      "The predictor is constant on complete cases."
    )
  }
  if (!identical(result$reason_code, "not_run")) {
    attr(result, "fit") <- NULL
    attr(result, "sample_index") <- source_index
    return(result)
  }

  formula <- .model_formula(length(effective_covariates),
                            !is.null(variables$offset))
  warning_messages <- character()
  engine_result <- tryCatch(
    withCallingHandlers(
      .fit_engine(variables$outcome_type, formula, model_data),
      warning = function(warning) {
        warning_messages <<- c(warning_messages, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    phewas_skip = function(error) error,
    error = function(error) error
  )

  if (inherits(engine_result, "error")) {
    is_skip <- inherits(engine_result, "phewas_skip")
    message <- conditionMessage(engine_result)
    reason_code <- if (is_skip) {
      engine_result$reason_code
    } else if (grepl("convergence", message, ignore.case = TRUE)) {
      "nonconvergence"
    } else if (grepl("Package `", message, fixed = TRUE)) {
      "dependency_missing"
    } else {
      "model_error"
    }
    result <- .set_result_problem(
      result, if (is_skip) "skipped" else "error", reason_code, message,
      warning_messages
    )
    attr(result, "fit") <- NULL
    attr(result, "sample_index") <- source_index
    return(result)
  }
  if (!isTRUE(engine_result$converged)) {
    result <- .set_result_problem(
      result, "error", "nonconvergence", "The model did not converge.",
      warning_messages
    )
    result$converged <- FALSE
    attr(result, "fit") <- engine_result$fit
    attr(result, "sample_index") <- source_index
    return(result)
  }

  multiplicative <- variables$outcome_type %in% c("binary", "count", "ordinal")
  result$status <- "ok"
  result$reason_code <- "ok"
  result$message <- NA_character_
  result$estimate <- engine_result$estimate
  result$std_error <- engine_result$std_error
  result$conf_low <- engine_result$conf_low
  result$conf_high <- engine_result$conf_high
  result$statistic <- engine_result$statistic
  result$test_method <- engine_result$test_method
  result$df <- engine_result$df
  result$native_effect <- if (multiplicative) exp(engine_result$estimate) else
    engine_result$estimate
  result$native_conf_low <- if (multiplicative) exp(engine_result$conf_low) else
    engine_result$conf_low
  result$native_conf_high <- if (multiplicative) exp(engine_result$conf_high) else
    engine_result$conf_high
  result$log_p <- engine_result$log_p
  result$p_value <- .probability_from_log(engine_result$log_p)
  result$p_value_text <- .format_log_probability(engine_result$log_p)
  result$neg_log10_p <- .neg_log10_from_log(engine_result$log_p)
  result$theta <- engine_result$theta
  result$ordinal_levels <- engine_result$ordinal_levels
  result$converged <- TRUE
  result$iterations <- engine_result$iterations
  if (length(warning_messages)) {
    result$warnings <- paste(unique(warning_messages), collapse = " | ")
  }

  attr(result, "fit") <- engine_result$fit
  attr(result, "sample_index") <- source_index
  result
}

.psock_fit_one <- function(row, data, analysis_spec) {
  .fit_one_association(data, analysis_spec, row)
}

#' Run one PheWAS testing family
#'
#' Fits every selected phenotype in deterministic metadata order. Parallel
#' execution uses a cross-platform PSOCK cluster and yields rows in the same
#' order as sequential execution.
#'
#' @param data Analysis-ready data containing every column declared in `spec`.
#' @param spec A validated `phewas_spec`.
#' @param backend Either `"sequential"` or `"psock"`.
#' @param workers Positive number of PSOCK workers.
#' @param phenotype_ids Optional explicit subset of `phenotype` values. This
#'   argument name is retained for API compatibility; result and metadata
#'   columns use the canonical name `phenotype`.
#' @param adjust Apply BH and Bonferroni correction across all successful rows.
#'
#' @return A `data.table` subclass `phewas_result`, with `spec` and
#'   `run_metadata` attributes and one row per requested phenotype. `formula`
#'   and `covariates` record the requested model. `effective_formula`,
#'   `effective_covariates`, and `dropped_invariant_covariates` record the
#'   outcome-specific model after omitting covariates that are constant on its
#'   complete cases.
#'
#' @details Results support normal `data.table` operations. Use
#'   [data.table::copy()] before modifying a result with `:=`, or
#'   `as.data.frame()` for base-data-frame subsetting semantics. `log_p` is the
#'   authoritative natural-log p-value; `neg_log10_p`, the printable
#'   `p_value_text`, and all multiple-testing fields are derived without first
#'   rounding through a numeric p-value. Thus `p_value` may underflow to zero
#'   while the log-scale and text fields remain valid.
#'
#' @examples
#' example <- phewas_example_data(n = 200)
#' spec <- phewas_spec(
#'   analysis_id = "example",
#'   anchor = "pgs",
#'   direction = "phenotypes_as_outcomes",
#'   phenotypes = example$forward_metadata,
#'   covariates = "age",
#'   anchor_type = "numeric",
#'   anchor_transform = "zscore",
#'   eligibility = list(
#'     min_n = 100L, min_cases = 10L, min_controls = 10L,
#'     min_outcome_levels = 3L
#'   ),
#'   fdr_threshold = 0.05
#' )
#' result <- run_phewas(example$data, spec)
#' as.data.frame(result)[c("phenotype", "description", "engine", "estimate",
#'   "status")]
#' @export
run_phewas <- function(data, spec, backend = c("sequential", "psock"),
                       workers = 1L, phenotype_ids = NULL, adjust = TRUE) {
  backend <- match.arg(backend)
  if (!is.numeric(workers) || length(workers) != 1L || !is.finite(workers) ||
      workers < 1L || workers > .Machine$integer.max || workers != trunc(workers)) {
    stop("`workers` must be one positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)
  if (!is.logical(adjust) || length(adjust) != 1L || is.na(adjust)) {
    stop("`adjust` must be TRUE or FALSE.", call. = FALSE)
  }
  spec <- validate_phewas_spec(spec, data)

  metadata <- spec$phenotypes
  if (!is.null(phenotype_ids)) {
    if (!is.character(phenotype_ids) || anyNA(phenotype_ids) ||
        anyDuplicated(phenotype_ids)) {
      stop("`phenotype_ids` must contain unique, nonmissing IDs.", call. = FALSE)
    }
    missing_ids <- setdiff(phenotype_ids, metadata$phenotype)
    if (length(missing_ids)) {
      stop(sprintf("Unknown phenotype IDs: %s.", paste(missing_ids, collapse = ", ")),
           call. = FALSE)
    }
    metadata <- metadata[metadata$phenotype %in% phenotype_ids, , drop = FALSE]
  }
  if (!nrow(metadata)) {
    stop("No phenotypes were selected.", call. = FALSE)
  }

  started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  prepared_data <- .prepare_analysis_data(data, spec)
  rows <- lapply(seq_len(nrow(metadata)), function(i) metadata[i, , drop = FALSE])

  if (identical(backend, "sequential") || workers == 1L) {
    fitted <- lapply(rows, function(row) {
      .fit_one_association(prepared_data, spec, row)
    })
    effective_backend <- "sequential"
    effective_workers <- 1L
  } else {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    fitted <- parallel::parLapply(
      cluster,
      rows,
      .psock_fit_one,
      data = prepared_data,
      analysis_spec = spec
    )
    effective_backend <- "psock"
    effective_workers <- workers
  }

  # Strip per-fit attributes before row binding: public results never retain
  # individual-level row indices or fitted model objects.
  fitted <- lapply(fitted, function(row) {
    attr(row, "fit") <- NULL
    attr(row, "sample_index") <- NULL
    class(row) <- "data.frame"
    row
  })
  results <- do.call(rbind, fitted)
  row.names(results) <- NULL
  if (anyDuplicated(paste(results$analysis_id, results$phenotype, sep = "\r"))) {
    stop("The run produced duplicate analysis_id + phenotype keys.",
         call. = FALSE)
  }
  if (adjust) {
    results <- .apply_multiple_testing(results, spec$fdr_threshold)
  }
  results$testing_family_size <- sum(results$status == "ok")

  results <- data.table::as.data.table(results)
  data.table::setattr(
    results, "class", c("phewas_result", "data.table", "data.frame")
  )
  attr(results, "spec") <- spec
  attr(results, "run_metadata") <- list(
    started_at = started_at,
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    backend = effective_backend,
    workers = effective_workers,
    phenotypes = results$phenotype,
    n_base = attr(prepared_data, "n_base"),
    exclusion_active = !is.null(spec$exclusion_column),
    exclusion_column = spec$exclusion_column,
    exclusion_values = spec$exclusion_values,
    n_before_exclusion = attr(prepared_data, "n_before_exclusion"),
    n_excluded = attr(prepared_data, "n_excluded"),
    n_after_exclusion = attr(prepared_data, "n_after_exclusion"),
    transformation = spec$anchor_transform,
    transformation_mean = attr(prepared_data, "transformation_mean"),
    transformation_sd = attr(prepared_data, "transformation_sd"),
    testing_family_size = sum(results$status == "ok"),
    corrections_applied = adjust
  )
  results
}

#' @export
print.phewas_result <- function(x, ...) {
  full_result_markers <- c(
    "analysis_id", "spec_version", "phenotype", "direction", "response",
    "predictor", "status", "reason_code", "log_p", "testing_family_size"
  )
  if (!all(full_result_markers %in% names(x))) {
    # Column subsets retain the result class; print them as ordinary tables.
    selected <- data.table::copy(x)
    data.table::setattr(selected, "class", c("data.table", "data.frame"))
    print(selected, ...)
    return(invisible(x))
  }
  counts <- table(factor(x$status, levels = c("ok", "skipped", "error")))
  cat(sprintf(
    "<phewas_result> %d associations: %d ok, %d skipped, %d error\n",
    nrow(x), counts[["ok"]], counts[["skipped"]], counts[["error"]]
  ))
  shown <- c(
    "phenotype", "description", "outcome_type", "status", "estimate",
    "native_effect", "p_value_text", "neg_log10_p", "q_value",
    "reason_code"
  )
  shown <- intersect(shown, names(x))
  print.data.frame(as.data.frame(x)[shown], ...)
  invisible(x)
}
