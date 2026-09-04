# Result accessors ---------------------------------------------------------

.pf_associations <- function(x, arg = deparse(substitute(x))) {
  if (!inherits(x, "phewas_result") || !is.data.frame(x)) {
    stop(arg, " must be a phewas_result.", call. = FALSE)
  }
  out <- .pf_normalize_legacy_result(x, label = arg)
  out <- if (inherits(out, "data.table")) data.table::copy(out) else out
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  if (!all(c("phenotype", "status") %in% names(out))) {
    stop(arg, " associations must contain `phenotype` and `status`.",
         call. = FALSE)
  }
  phenotype <- trimws(as.character(out$phenotype))
  if (anyNA(phenotype) || any(!nzchar(phenotype))) {
    stop(arg, " contains a missing or empty phenotype.", call. = FALSE)
  }
  if (anyDuplicated(phenotype)) {
    stop(arg, " contains duplicated phenotype values.", call. = FALSE)
  }
  out$phenotype <- phenotype
  if (!"description" %in% names(out)) out$description <- phenotype
  out$description <- trimws(as.character(out$description))
  missing_description <- is.na(out$description) | !nzchar(out$description)
  out$description[missing_description] <- phenotype[missing_description]
  if (!"group" %in% names(out)) out$group <- "Uncategorized"
  out$group <- trimws(as.character(out$group))
  missing_group <- is.na(out$group) | !nzchar(out$group)
  out$group[missing_group] <- "Uncategorized"
  if (!"groupnum" %in% names(out)) out$groupnum <- NA_real_
  out$groupnum <- suppressWarnings(as.numeric(out$groupnum))
  if (!"color" %in% names(out)) out$color <- NA_character_
  out$color <- as.character(out$color)
  if ("category_color" %in% names(out)) {
    alias_color <- as.character(out$category_color)
    conflict <- !is.na(out$color) & nzchar(out$color) &
      !is.na(alias_color) & nzchar(alias_color) & out$color != alias_color
    if (any(conflict)) {
      stop(arg, " has conflicting legacy values for canonical `color` metadata.",
           call. = FALSE)
    }
    fill <- (is.na(out$color) | !nzchar(out$color)) &
      !is.na(alias_color) & nzchar(alias_color)
    out$color[fill] <- alias_color[fill]
    out$category_color <- NULL
  }
  front <- c(
    "phenotype", "description", "group", "groupnum", "color", "status"
  )
  out[, c(front, setdiff(names(out), front)), drop = FALSE]
}

.pf_result_spec <- function(x) {
  attr(x, "spec", exact = TRUE)
}

.pf_run_metadata <- function(x) {
  attr(x, "run_metadata", exact = TRUE)
}

.pf_column <- function(x, choices, default = NA) {
  hit <- choices[choices %in% names(x)]
  if (length(hit)) {
    x[[hit[[1L]]]]
  } else if (length(default) == nrow(x)) {
    default
  } else {
    rep(default, length.out = nrow(x))
  }
}

.pf_scalar_text <- function(x) {
  if (!length(x) || is.null(x) || all(is.na(x))) return(NA_character_)
  value <- unique(as.character(x[!is.na(x)]))
  if (length(value) == 1L) value else NA_character_
}

.pf_analysis_label <- function(x, fallback) {
  spec <- .pf_result_spec(x)
  meta <- .pf_run_metadata(x)
  value <- c(
    if (is.list(spec)) spec$analysis_id else NULL,
    if (is.list(meta)) meta$analysis_id else NULL
  )
  value <- value[!is.na(value) & nzchar(as.character(value))]
  if (length(value)) as.character(value[[1L]]) else fallback
}

.pf_equal_na <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
}

.pf_comparison_table <- function(result, arg) {
  x <- .pf_associations(result, arg)
  required <- c("estimate", "std_error", "p_value", "q_value", "direction")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(arg, " is missing comparison column(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (anyNA(x$direction) || any(!x$direction %in% c(
    "phenotypes_as_outcomes", "phenotypes_as_predictors"
  ))) {
    stop(arg, " contains invalid or missing direction metadata.", call. = FALSE)
  }

  data.frame(
    phenotype = x$phenotype,
    status = as.character(x$status),
    description = x$description,
    group = x$group,
    groupnum = x$groupnum,
    color = x$color,
    direction = as.character(x$direction),
    response = as.character(.pf_column(x, "response")),
    predictor = as.character(.pf_column(x, "predictor")),
    outcome_type = as.character(.pf_column(x, "outcome_type")),
    predictor_type = as.character(.pf_column(x, "predictor_type")),
    transformation = as.character(.pf_column(x, "transformation", "none")),
    transformation_mean = suppressWarnings(as.numeric(.pf_column(
      x, "transformation_mean", NA_real_
    ))),
    transformation_sd = suppressWarnings(as.numeric(.pf_column(
      x, "transformation_sd", NA_real_
    ))),
    reference = as.character(.pf_column(x, "reference")),
    outcome_reference = as.character(.pf_column(
      x, "outcome_reference", .pf_column(x, "reference")
    )),
    predictor_reference = as.character(.pf_column(x, "predictor_reference")),
    levels = as.character(.pf_column(x, "levels")),
    outcome_levels = as.character(.pf_column(x, "outcome_levels")),
    predictor_levels = as.character(.pf_column(x, "predictor_levels")),
    scores = as.character(.pf_column(x, "scores")),
    effect_measure = as.character(.pf_column(x, "effect_measure")),
    estimate = suppressWarnings(as.numeric(x$estimate)),
    std_error = suppressWarnings(as.numeric(x$std_error)),
    p_value = suppressWarnings(as.numeric(x$p_value)),
    q_value = suppressWarnings(as.numeric(x$q_value)),
    neg_log10_p = suppressWarnings(as.numeric(.pf_column(x, "neg_log10_p"))),
    native_effect = suppressWarnings(as.numeric(.pf_column(x, "native_effect"))),
    native_conf_low = suppressWarnings(as.numeric(.pf_column(x, "native_conf_low"))),
    native_conf_high = suppressWarnings(as.numeric(.pf_column(x, "native_conf_high"))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

.pf_fdr_threshold <- function(x, y, supplied) {
  if (!is.null(supplied)) {
    if (length(supplied) != 1L || !is.finite(supplied) || supplied <= 0 || supplied >= 1) {
      stop("fdr_threshold must be one number strictly between 0 and 1.",
           call. = FALSE)
    }
    return(as.numeric(supplied))
  }
  get_one <- function(z) {
    spec <- .pf_result_spec(z)
    if (is.list(spec)) spec$fdr_threshold else NULL
  }
  values <- unlist(list(get_one(x), get_one(y)), use.names = FALSE)
  values <- unique(as.numeric(values[!is.na(values)]))
  if (length(values) != 1L || !is.finite(values)) {
    stop("Supply fdr_threshold when the two results do not record the same threshold.",
         call. = FALSE)
  }
  values
}

#' Compare two PheWAS results
#'
#' Compares aggregate association results only. Correlations use phenotypes
#' successfully fitted in both runs, while discovery overlap uses every
#' successful phenotype in either run. Native-scale statistics are reported
#' only when predictor, transformation, reference, and effect-scale metadata
#' agree for every matched phenotype.
#'
#' @param x,y Two `phewas_result` objects.
#' @param fdr_threshold FDR cutoff. By default the common cutoff recorded in
#'   the result specifications is used.
#' @param x_name,y_name Optional labels for the two runs.
#'
#' @return A `phewas_comparison` containing `matched`, `metrics`, and
#'   `group_overlap` data frames.
#' @export
compare_phewas <- function(x, y, fdr_threshold = NULL, x_name = NULL,
                           y_name = NULL) {
  left <- .pf_comparison_table(x, "x")
  right <- .pf_comparison_table(y, "y")
  cutoff <- .pf_fdr_threshold(x, y, fdr_threshold)
  x_name <- x_name %||% .pf_analysis_label(x, "x")
  y_name <- y_name %||% .pf_analysis_label(y, "y")

  overlap <- merge(left, right, by = "phenotype", suffixes = c("_x", "_y"),
                   sort = TRUE)
  if (nrow(overlap)) {
    same_direction <- .pf_equal_na(overlap$direction_x, overlap$direction_y)
    reference_ok <- rep(TRUE, nrow(overlap))
    if (any(same_direction)) {
      reference_ok[same_direction] <-
        .pf_equal_na(overlap$outcome_reference_x[same_direction],
                     overlap$outcome_reference_y[same_direction]) &
        .pf_equal_na(overlap$predictor_reference_x[same_direction],
                     overlap$predictor_reference_y[same_direction])
    }
    x_forward <- overlap$direction_x == "phenotypes_as_outcomes"
    cross <- !same_direction
    if (any(cross & x_forward)) {
      take <- cross & x_forward
      reference_ok[take] <- .pf_equal_na(
        overlap$outcome_reference_x[take], overlap$predictor_reference_y[take]
      )
    }
    if (any(cross & !x_forward)) {
      take <- cross & !x_forward
      reference_ok[take] <- .pf_equal_na(
        overlap$predictor_reference_x[take], overlap$outcome_reference_y[take]
      )
    }
    if (any(!reference_ok)) {
      bad <- overlap$phenotype[!reference_ok]
      stop("Cannot compare incompatible reference levels for phenotype(s): ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    }
    levels_ok <- rep(TRUE, nrow(overlap))
    if (any(same_direction)) {
      levels_ok[same_direction] <-
        .pf_equal_na(overlap$outcome_levels_x[same_direction],
                     overlap$outcome_levels_y[same_direction]) &
        .pf_equal_na(overlap$predictor_levels_x[same_direction],
                     overlap$predictor_levels_y[same_direction])
    }
    if (any(cross & x_forward)) {
      take <- cross & x_forward
      levels_ok[take] <- .pf_equal_na(
        overlap$outcome_levels_x[take], overlap$predictor_levels_y[take]
      )
    }
    if (any(cross & !x_forward)) {
      take <- cross & !x_forward
      levels_ok[take] <- .pf_equal_na(
        overlap$predictor_levels_x[take], overlap$outcome_levels_y[take]
      )
    }
    if (any(!levels_ok)) {
      bad <- overlap$phenotype[!levels_ok]
      stop("Cannot compare incompatible ordinal orientation for phenotype(s): ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    }
    group_ok <- .pf_equal_na(overlap$group_x, overlap$group_y)
    if (any(!group_ok)) {
      bad <- overlap$phenotype[!group_ok]
      stop("Conflicting group metadata for phenotype(s): ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    }
  }

  matched <- overlap[
    overlap$status_x == "ok" & overlap$status_y == "ok",
    , drop = FALSE
  ]
  matched$z_x <- matched$estimate_x / matched$std_error_x
  matched$z_y <- matched$estimate_y / matched$std_error_y
  matched$significant_x <- !is.na(matched$q_value_x) & matched$q_value_x <= cutoff
  matched$significant_y <- !is.na(matched$q_value_y) & matched$q_value_y <= cutoff

  native_fields <- c(
    "direction", "response", "outcome_type", "predictor", "predictor_type",
    "transformation", "transformation_mean", "transformation_sd",
    "outcome_reference", "predictor_reference", "outcome_levels",
    "predictor_levels", "scores", "effect_measure"
  )
  native_ok <- nrow(matched) > 0L
  native_reason <- character()
  if (nrow(matched)) {
    for (field in native_fields) {
      same <- .pf_equal_na(matched[[paste0(field, "_x")]],
                           matched[[paste0(field, "_y")]])
      if (any(!same)) native_reason <- c(native_reason, field)
    }
    if (length(unique(matched$effect_measure_x[!is.na(matched$effect_measure_x)])) != 1L) {
      native_reason <- c(native_reason, "mixed effect measures")
    }
    native_ok <- !length(native_reason) &&
      all(is.finite(matched$native_effect_x)) &&
      all(is.finite(matched$native_effect_y))
  } else {
    native_reason <- "no mutually successful phenotypes"
  }

  safe_cor <- function(a, b, method = "pearson") {
    keep <- is.finite(a) & is.finite(b)
    if (sum(keep) < 2L || stats::sd(a[keep]) == 0 || stats::sd(b[keep]) == 0) {
      return(NA_real_)
    }
    stats::cor(a[keep], b[keep], method = method)
  }
  sign_keep <- is.finite(matched$estimate_x) & is.finite(matched$estimate_y) &
    matched$estimate_x != 0 & matched$estimate_y != 0
  native_slope <- native_intercept <- native_rmse <- native_cor <- NA_real_
  if (native_ok && nrow(matched) >= 2L) {
    native_x <- matched$native_effect_x
    native_y <- matched$native_effect_y
    native_cor <- safe_cor(native_x, native_y)
    model <- stats::lm(native_y ~ native_x)
    native_intercept <- unname(stats::coef(model)[[1L]])
    native_slope <- unname(stats::coef(model)[[2L]])
    native_rmse <- sqrt(mean(
      (matched$native_effect_y - matched$native_effect_x)^2,
      na.rm = TRUE
    ))
  }

  sig_x <- left$phenotype[left$status == "ok" & !is.na(left$q_value) &
                            left$q_value <= cutoff]
  sig_y <- right$phenotype[right$status == "ok" & !is.na(right$q_value) &
                             right$q_value <= cutoff]
  intersection <- intersect(sig_x, sig_y)
  union_ids <- union(sig_x, sig_y)
  metrics <- data.frame(
    n_x = sum(left$status == "ok"),
    n_y = sum(right$status == "ok"),
    n_mutually_successful = nrow(matched),
    fdr_x = length(sig_x),
    fdr_y = length(sig_y),
    fdr_intersection = length(intersection),
    fdr_union = length(union_ids),
    fdr_jaccard = if (length(union_ids)) length(intersection) / length(union_ids) else NA_real_,
    sign_concordance = if (any(sign_keep)) {
      mean(sign(matched$estimate_x[sign_keep]) == sign(matched$estimate_y[sign_keep]))
    } else NA_real_,
    standardized_estimate_correlation = safe_cor(matched$z_x, matched$z_y),
    standardized_estimate_rank_correlation = safe_cor(matched$z_x, matched$z_y, "spearman"),
    p_value_rank_correlation = safe_cor(
      matched$neg_log10_p_x, matched$neg_log10_p_y, "spearman"
    ),
    native_compatible = native_ok,
    native_correlation = native_cor,
    native_intercept = native_intercept,
    native_slope = native_slope,
    native_rmse = native_rmse,
    stringsAsFactors = FALSE
  )

  all_ids <- sort(unique(c(left$phenotype, right$phenotype)))
  groups <- data.frame(phenotype = all_ids, stringsAsFactors = FALSE)
  groups <- merge(groups, left[c("phenotype", "group", "status", "q_value")],
                  by = "phenotype", all.x = TRUE, sort = TRUE)
  names(groups)[2:4] <- c("group_x", "status_x", "q_value_x")
  groups <- merge(groups, right[c("phenotype", "group", "status", "q_value")],
                  by = "phenotype", all.x = TRUE, sort = TRUE)
  names(groups)[5:7] <- c("group_y", "status_y", "q_value_y")
  groups$group <- ifelse(!is.na(groups$group_x), groups$group_x, groups$group_y)
  groups$group[is.na(groups$group) | !nzchar(groups$group)] <- "Uncategorized"
  groups$sig_x <- groups$status_x == "ok" & !is.na(groups$q_value_x) &
    groups$q_value_x <= cutoff
  groups$sig_y <- groups$status_y == "ok" & !is.na(groups$q_value_y) &
    groups$q_value_y <= cutoff
  groups$sig_x[is.na(groups$sig_x)] <- FALSE
  groups$sig_y[is.na(groups$sig_y)] <- FALSE
  split_group <- split(groups, groups$group, drop = TRUE)
  group_overlap <- do.call(rbind, lapply(split_group, function(z) {
    both <- sum(z$sig_x & z$sig_y)
    either <- sum(z$sig_x | z$sig_y)
    data.frame(
      group = z$group[[1L]],
      tested_x = sum(z$status_x == "ok", na.rm = TRUE),
      tested_y = sum(z$status_y == "ok", na.rm = TRUE),
      significant_x = sum(z$sig_x),
      significant_y = sum(z$sig_y),
      both = both,
      only_x = sum(z$sig_x & !z$sig_y),
      only_y = sum(!z$sig_x & z$sig_y),
      either = either,
      jaccard = if (either) both / either else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(group_overlap) <- NULL
  group_overlap <- group_overlap[order(group_overlap$group), , drop = FALSE]

  out <- list(
    matched = data.table::as.data.table(matched),
    metrics = data.table::as.data.table(metrics),
    group_overlap = data.table::as.data.table(group_overlap),
    fdr_threshold = cutoff,
    run_labels = c(x = x_name, y = y_name),
    native_compatible = native_ok,
    native_incompatibility = unique(native_reason)
  )
  class(out) <- "phewas_comparison"
  out
}

#' @export
print.phewas_comparison <- function(x, ...) {
  cat("<phewas_comparison>", x$run_labels[["x"]], "vs", x$run_labels[["y"]], "\n")
  cat("  mutually successful:", x$metrics$n_mutually_successful, "\n")
  cat("  FDR intersection/union:", x$metrics$fdr_intersection, "/",
      x$metrics$fdr_union, "\n")
  cat("  native effects compatible:", x$native_compatible, "\n")
  invisible(x)
}
