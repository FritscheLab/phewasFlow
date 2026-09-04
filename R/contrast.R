.pf_assert_nested_specs <- function(base_spec, adjusted_spec, data) {
  base_spec <- validate_phewas_spec(base_spec, data)
  adjusted_spec <- validate_phewas_spec(adjusted_spec, data)
  if (!all(base_spec$covariates %in% adjusted_spec$covariates)) {
    stop("The base covariates must be a subset of the adjusted covariates.",
         call. = FALSE)
  }
  if (!length(setdiff(adjusted_spec$covariates, base_spec$covariates))) {
    stop("The adjusted specification must add at least one covariate.",
         call. = FALSE)
  }
  fields <- c(
    "version", "anchor", "direction", "phenotypes", "outcome_type",
    "anchor_type", "anchor_transform", "anchor_reference", "anchor_levels",
    "anchor_scores", "offset", "eligibility", "fdr_threshold",
    "exclusion_column", "exclusion_values"
  )
  different <- fields[!vapply(fields, function(field) {
    isTRUE(all.equal(base_spec[[field]], adjusted_spec[[field]],
                     check.attributes = FALSE))
  }, logical(1))]
  if (length(different)) {
    stop("Nested specifications differ outside analysis_id/covariates: ",
         paste(different, collapse = ", "), ".", call. = FALSE)
  }
  list(base_spec = base_spec, adjusted_spec = adjusted_spec)
}

.pf_contrast_value <- function(row, name, default = NA) {
  if (name %in% names(row)) row[[name]][[1L]] else default
}

.pf_fit_for_contrast <- function(data, spec, phenotype_row) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(
      .fit_one_association(data, spec, phenotype_row),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(value, "error")) {
    return(list(
      row = data.frame(), fit = NULL, sample_index = integer(),
      status = if (inherits(value, "phewas_skip")) "skipped" else "error",
      reason_code = value$reason_code %||% if (inherits(value, "phewas_skip")) {
        "ineligible"
      } else {
        "fit_error"
      },
      message = conditionMessage(value), warnings = warnings
    ))
  }
  list(
    row = as.data.frame(value, stringsAsFactors = FALSE),
    fit = attr(value, "fit", exact = TRUE),
    sample_index = attr(value, "sample_index", exact = TRUE) %||% integer(),
    status = as.character(.pf_contrast_value(value, "status", "ok")),
    reason_code = as.character(.pf_contrast_value(value, "reason_code", NA_character_)),
    message = as.character(.pf_contrast_value(value, "message", NA_character_)),
    warnings = unique(c(warnings, as.character(.pf_contrast_value(
      value, "warnings", NA_character_
    ))))
  )
}

.pf_subset_prepared <- function(data, row_index) {
  source_index <- attr(data, "source_index", exact = TRUE)
  sample_keys <- attr(data, "sample_keys", exact = TRUE)
  preserved <- attributes(data)[c(
    "phewas_prepared", "n_base", "n_before_exclusion", "n_excluded",
    "n_after_exclusion", "transformation_mean", "transformation_sd"
  )]
  out <- data[row_index, , drop = FALSE]
  row.names(out) <- NULL
  for (name in names(preserved)) attr(out, name) <- preserved[[name]]
  attr(out, "phewas_prepared") <- TRUE
  attr(out, "source_index") <- source_index[row_index]
  attr(out, "sample_keys") <- sample_keys[row_index]
  out
}

.pf_lm_contrast_r2 <- function(base_fit, adjusted_fit) {
  out <- c(delta_r_squared = NA_real_, partial_r_squared = NA_real_)
  if (!inherits(base_fit, "lm") || !inherits(adjusted_fit, "lm")) return(out)
  base_r2 <- summary(base_fit)$r.squared
  adjusted_r2 <- summary(adjusted_fit)$r.squared
  base_rss <- sum(stats::residuals(base_fit)^2)
  adjusted_rss <- sum(stats::residuals(adjusted_fit)^2)
  out[["delta_r_squared"]] <- adjusted_r2 - base_r2
  if (is.finite(base_rss) && base_rss > 0) {
    out[["partial_r_squared"]] <- (base_rss - adjusted_rss) / base_rss
  }
  out
}

.pf_empty_contrast_fit <- function(status, reason_code, message) {
  list(
    row = data.frame(), fit = NULL, sample_index = integer(), status = status,
    reason_code = reason_code, message = message, warnings = character()
  )
}

.pf_build_contrast_row <- function(phenotype_row, base, adjusted,
                                   base_spec, adjusted_spec, shared_n) {
  base_row <- base$row
  adjusted_row <- adjusted$row
  base_value <- function(name, default = NA) {
    if (nrow(base_row)) .pf_contrast_value(base_row, name, default) else default
  }
  adjusted_value <- function(name, default = NA) {
    if (nrow(adjusted_row)) .pf_contrast_value(adjusted_row, name, default) else default
  }
  both_ok <- identical(base$status, "ok") && identical(adjusted$status, "ok")
  status <- if (both_ok) "ok" else if ("error" %in% c(base$status, adjusted$status)) {
    "error"
  } else {
    "skipped"
  }
  reason <- if (both_ok) {
    NA_character_
  } else {
    paste(c(
      if (!identical(base$status, "ok")) paste0("base_", base$reason_code),
      if (!identical(adjusted$status, "ok")) paste0("adjusted_", adjusted$reason_code)
    ), collapse = ";")
  }
  message <- paste(stats::na.omit(c(
    if (!identical(base$status, "ok")) paste0("Base: ", base$message),
    if (!identical(adjusted$status, "ok")) paste0("Adjusted: ", adjusted$message)
  )), collapse = " | ")
  if (!nzchar(message)) message <- NA_character_
  samples_match <- identical(base$sample_index, adjusted$sample_index)
  if (both_ok && !samples_match) {
    status <- "error"
    reason <- "sample_mismatch"
    message <- "Base and adjusted models did not use identical analysis rows."
    both_ok <- FALSE
  }

  outcome_type <- as.character(adjusted_value(
    "outcome_type", base_value("outcome_type", NA_character_)
  ))
  r2 <- if (both_ok && identical(outcome_type, "continuous")) {
    .pf_lm_contrast_r2(base$fit, adjusted$fit)
  } else {
    c(delta_r_squared = NA_real_, partial_r_squared = NA_real_)
  }
  base_estimate <- as.numeric(base_value("estimate", NA_real_))
  adjusted_estimate <- as.numeric(adjusted_value("estimate", NA_real_))
  base_native <- as.numeric(base_value("native_effect", NA_real_))
  adjusted_native <- as.numeric(adjusted_value("native_effect", NA_real_))
  base_ok <- identical(base$status, "ok")
  adjusted_ok <- identical(adjusted$status, "ok")

  data.frame(
    analysis_id = paste0(base_spec$analysis_id, "__vs__", adjusted_spec$analysis_id),
    base_analysis_id = base_spec$analysis_id,
    adjusted_analysis_id = adjusted_spec$analysis_id,
    phenotype = as.character(phenotype_row$phenotype[[1L]]),
    description = as.character(phenotype_row$description[[1L]]),
    group = as.character(phenotype_row$group[[1L]]),
    groupnum = as.numeric(phenotype_row$groupnum[[1L]]),
    color = as.character(phenotype_row$color[[1L]]),
    direction = base_spec$direction,
    response = as.character(adjusted_value("response", base_value("response", NA_character_))),
    predictor = as.character(adjusted_value("predictor", base_value("predictor", NA_character_))),
    outcome_type = outcome_type,
    effect_measure = as.character(adjusted_value(
      "effect_measure", base_value("effect_measure", NA_character_)
    )),
    transformation = base_spec$anchor_transform,
    transformation_mean = as.numeric(adjusted_value("transformation_mean", NA_real_)),
    transformation_sd = as.numeric(adjusted_value("transformation_sd", NA_real_)),
    transformation_source_n = as.integer(adjusted_value("transformation_source_n", NA_integer_)),
    reference = as.character(adjusted_value("reference", base_value("reference", NA_character_))),
    outcome_reference = as.character(adjusted_value(
      "outcome_reference", base_value("outcome_reference", NA_character_)
    )),
    predictor_reference = as.character(adjusted_value(
      "predictor_reference", base_value("predictor_reference", NA_character_)
    )),
    offset = as.character(adjusted_value("offset", base_value("offset", NA_character_))),
    exclusion_active = as.logical(adjusted_value(
      "exclusion_active", base_value("exclusion_active", FALSE)
    )),
    exclusion_column = as.character(adjusted_value(
      "exclusion_column", base_value("exclusion_column", NA_character_)
    )),
    exclusion_values = as.character(adjusted_value(
      "exclusion_values", base_value("exclusion_values", NA_character_)
    )),
    n_before_exclusion = as.integer(adjusted_value(
      "n_before_exclusion", base_value("n_before_exclusion", NA_integer_)
    )),
    n_excluded = as.integer(adjusted_value(
      "n_excluded", base_value("n_excluded", NA_integer_)
    )),
    n_after_exclusion = as.integer(adjusted_value(
      "n_after_exclusion", base_value("n_after_exclusion", NA_integer_)
    )),
    status = status,
    reason_code = reason,
    message = message,
    base_status = base$status,
    base_reason_code = base$reason_code,
    base_message = base$message,
    adjusted_status = adjusted$status,
    adjusted_reason_code = adjusted$reason_code,
    adjusted_message = adjusted$message,
    n_complete = as.integer(if (both_ok) {
      adjusted_value("n_complete", shared_n)
    } else shared_n),
    n_base = as.integer(adjusted_value("n_base", base_value("n_base", NA_integer_))),
    base_covariates = paste(base_spec$covariates, collapse = "|"),
    adjusted_covariates = paste(adjusted_spec$covariates, collapse = "|"),
    added_covariates = paste(
      setdiff(adjusted_spec$covariates, base_spec$covariates), collapse = "|"
    ),
    base_formula = as.character(base_value("formula", NA_character_)),
    adjusted_formula = as.character(adjusted_value("formula", NA_character_)),
    base_estimate = if (base_ok) base_estimate else NA_real_,
    base_std_error = if (base_ok) as.numeric(base_value("std_error", NA_real_)) else NA_real_,
    base_conf_low = if (base_ok) as.numeric(base_value("conf_low", NA_real_)) else NA_real_,
    base_conf_high = if (base_ok) as.numeric(base_value("conf_high", NA_real_)) else NA_real_,
    base_native_effect = if (base_ok) base_native else NA_real_,
    base_native_conf_low = if (base_ok) as.numeric(base_value("native_conf_low", NA_real_)) else NA_real_,
    base_native_conf_high = if (base_ok) as.numeric(base_value("native_conf_high", NA_real_)) else NA_real_,
    base_log_p = if (base_ok) as.numeric(base_value("log_p", NA_real_)) else NA_real_,
    base_p_value = if (base_ok) as.numeric(base_value("p_value", NA_real_)) else NA_real_,
    adjusted_estimate = if (adjusted_ok) adjusted_estimate else NA_real_,
    adjusted_std_error = if (adjusted_ok) as.numeric(adjusted_value("std_error", NA_real_)) else NA_real_,
    adjusted_conf_low = if (adjusted_ok) as.numeric(adjusted_value("conf_low", NA_real_)) else NA_real_,
    adjusted_conf_high = if (adjusted_ok) as.numeric(adjusted_value("conf_high", NA_real_)) else NA_real_,
    adjusted_native_effect = if (adjusted_ok) adjusted_native else NA_real_,
    adjusted_native_conf_low = if (adjusted_ok) as.numeric(adjusted_value("native_conf_low", NA_real_)) else NA_real_,
    adjusted_native_conf_high = if (adjusted_ok) as.numeric(adjusted_value("native_conf_high", NA_real_)) else NA_real_,
    adjusted_log_p = if (adjusted_ok) as.numeric(adjusted_value("log_p", NA_real_)) else NA_real_,
    adjusted_p_value = if (adjusted_ok) as.numeric(adjusted_value("p_value", NA_real_)) else NA_real_,
    estimate_change = if (both_ok) adjusted_estimate - base_estimate else NA_real_,
    native_effect_change = if (both_ok) adjusted_native - base_native else NA_real_,
    native_effect_ratio = if (both_ok && is.finite(base_native) && base_native != 0) {
      adjusted_native / base_native
    } else NA_real_,
    delta_r_squared = unname(r2[["delta_r_squared"]]),
    partial_r_squared = unname(r2[["partial_r_squared"]]),
    sample_fingerprint = as.character(adjusted_value(
      "sample_fingerprint", base_value("sample_fingerprint", NA_character_)
    )),
    warnings = paste(unique(stats::na.omit(c(base$warnings, adjusted$warnings))),
                     collapse = " | "),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.pf_adjust_contrast_side <- function(results, side) {
  log_name <- paste0(side, "_log_p")
  prefix <- paste0(side, "_")
  n <- nrow(results)
  results[[paste0(prefix, "q_value")]] <- rep(NA_real_, n)
  results[[paste0(prefix, "q_value_text")]] <- rep(NA_character_, n)
  results[[paste0(prefix, "neg_log10_q")]] <- rep(NA_real_, n)
  results[[paste0(prefix, "bonferroni_p")]] <- rep(NA_real_, n)
  results[[paste0(prefix, "neg_log10_bonferroni_p")]] <- rep(NA_real_, n)
  ok <- which(results$status == "ok" & !is.na(results[[log_name]]))
  if (!length(ok)) return(results)
  log_q <- .adjust_log_bh(results[[log_name]][ok])
  log_bonf <- pmin(0, results[[log_name]][ok] + log(length(ok)))
  results[[paste0(prefix, "q_value")]][ok] <- .probability_from_log(log_q)
  results[[paste0(prefix, "q_value_text")]][ok] <- .format_log_probability(log_q)
  results[[paste0(prefix, "neg_log10_q")]][ok] <- .neg_log10_from_log(log_q)
  results[[paste0(prefix, "bonferroni_p")]][ok] <- .probability_from_log(log_bonf)
  results[[paste0(prefix, "neg_log10_bonferroni_p")]][ok] <-
    .neg_log10_from_log(log_bonf)
  results
}

#' Compare nested PheWAS specifications on identical samples
#'
#' Each phenotype is fitted twice on the exact same complete-case rows. The
#' adjusted specification may differ from the base specification only by its
#' analysis ID and by adding covariates. BH and Bonferroni values are then
#' recomputed over the mutually successful testing universe.
#'
#' @param data Analysis-ready individual-level data.
#' @param base_spec,adjusted_spec Nested `phewas_spec` objects.
#' @param phenotype_ids Optional subset of phenotype values.
#'
#' @return A `phewas_contrast` whose `associations` element contains paired
#'   estimates, changes, and shared-universe corrections.
#' @export
run_phewas_contrast <- function(data, base_spec, adjusted_spec,
                                phenotype_ids = NULL) {
  data <- .plain_data_frame(data, "data")
  specs <- .pf_assert_nested_specs(base_spec, adjusted_spec, data)
  base_spec <- specs$base_spec
  adjusted_spec <- specs$adjusted_spec
  ids <- as.character(adjusted_spec$phenotypes$phenotype)
  if (!is.null(phenotype_ids)) {
    phenotype_ids <- unique(as.character(phenotype_ids))
    unknown <- setdiff(phenotype_ids, ids)
    if (length(unknown)) {
      stop("Unknown phenotype value(s): ", paste(unknown, collapse = ", "),
           ".", call. = FALSE)
    }
    ids <- ids[ids %in% phenotype_ids]
  }
  prepared <- .prepare_analysis_data(data, adjusted_spec)
  rows <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    phenotype_row <- adjusted_spec$phenotypes[
      adjusted_spec$phenotypes$phenotype == ids[[i]], , drop = FALSE
    ]
    base_variables <- .association_variables(base_spec, phenotype_row)
    adjusted_variables <- .association_variables(adjusted_spec, phenotype_row)
    complete_variables <- unique(c(
      base_variables$complete_case_variables,
      adjusted_variables$complete_case_variables
    ))
    shared <- which(stats::complete.cases(prepared[, complete_variables, drop = FALSE]))
    if (!length(shared)) {
      base <- .pf_empty_contrast_fit("skipped", "no_complete_cases",
                                     "No shared complete-case rows.")
      adjusted <- base
    } else {
      shared_data <- .pf_subset_prepared(prepared, shared)
      base <- .pf_fit_for_contrast(shared_data, base_spec, phenotype_row)
      adjusted <- .pf_fit_for_contrast(shared_data, adjusted_spec, phenotype_row)
    }
    rows[[i]] <- .pf_build_contrast_row(
      phenotype_row, base, adjusted, base_spec, adjusted_spec, length(shared)
    )
  }
  associations <- if (length(rows)) do.call(rbind, rows) else data.frame()
  rownames(associations) <- NULL
  if (nrow(associations)) {
    associations <- .pf_adjust_contrast_side(associations, "base")
    associations <- .pf_adjust_contrast_side(associations, "adjusted")
  }
  diagnostics <- if (nrow(associations)) {
    associations[c("phenotype", "status", "reason_code", "message",
                   "n_complete", "sample_fingerprint", "warnings")]
  } else {
    data.frame()
  }
  associations <- data.table::as.data.table(associations)
  diagnostics <- data.table::as.data.table(diagnostics)
  out <- list(
    associations = associations,
    diagnostics = diagnostics,
    base_spec = base_spec,
    adjusted_spec = adjusted_spec,
    run_metadata = list(
      analysis_id = paste0(base_spec$analysis_id, "__vs__", adjusted_spec$analysis_id),
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      testing_universe = sum(associations$status == "ok")
    )
  )
  class(out) <- "phewas_contrast"
  out
}

#' @export
print.phewas_contrast <- function(x, ...) {
  cat("<phewas_contrast>", x$run_metadata$analysis_id, "\n")
  cat("  associations:", nrow(x$associations), "\n")
  cat("  mutually successful:", sum(x$associations$status == "ok"), "\n")
  invisible(x)
}
