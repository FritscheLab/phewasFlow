.stop_association <- function(reason_code, message) {
  condition <- structure(
    list(message = message, call = NULL, reason_code = reason_code),
    class = c("phewas_skip", "error", "condition")
  )
  stop(condition)
}

.encode_binary <- function(x, reference, role) {
  reference <- .metadata_scalar(reference)
  if (is.null(reference)) {
    .stop_association("missing_reference",
                      sprintf("The binary %s has no declared reference.", role))
  }
  values <- as.character(x)
  observed <- unique(values)
  if (length(observed) > 2L) {
    .stop_association(
      "invalid_binary_levels",
      sprintf("The binary %s has more than two observed levels.", role)
    )
  }
  reference <- as.character(reference)
  if (length(observed) == 2L && !reference %in% observed) {
    .stop_association(
      "reference_not_observed",
      sprintf("The declared reference `%s` is not observed in the %s.",
              reference, role)
    )
  }
  out <- as.integer(values != reference)
  other <- setdiff(observed, reference)
  attr(out, "reference") <- reference
  attr(out, "comparison") <- other
  out
}

.encode_ordinal <- function(x, levels, scores = NULL, role,
                            as_predictor = FALSE) {
  levels <- .metadata_sequence(levels)
  values <- as.character(x)
  unknown <- setdiff(unique(values), levels)
  if (length(unknown)) {
    .stop_association(
      "undeclared_ordinal_level",
      sprintf("The %s contains undeclared levels: %s.", role,
              paste(unknown, collapse = ", "))
    )
  }
  observed <- levels[levels %in% values]
  if (as_predictor && length(observed) < 2L) {
    .stop_association("constant_predictor",
                      sprintf("The ordinal %s has fewer than two observed levels.", role))
  }
  if (!as_predictor) {
    return(ordered(values, levels = observed))
  }

  scores <- .metadata_sequence(scores, numeric = TRUE)
  if (!length(scores)) {
    scores <- seq_along(levels)
  }
  unname(scores[match(values, levels)])
}

.prepare_outcome <- function(x, outcome_type, reference, levels) {
  if (identical(outcome_type, "binary")) {
    return(.encode_binary(x, reference, "outcome"))
  }
  if (identical(outcome_type, "continuous")) {
    if (!is.numeric(x) || any(!is.finite(x))) {
      .stop_association("invalid_continuous_outcome",
                        "A continuous outcome must contain finite numeric values.")
    }
    return(as.numeric(x))
  }
  if (identical(outcome_type, "count")) {
    if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0) ||
        any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
      .stop_association(
        "invalid_count_outcome",
        "A count outcome must contain finite, nonnegative integer values."
      )
    }
    return(as.numeric(x))
  }
  if (identical(outcome_type, "ordinal")) {
    return(.encode_ordinal(x, levels, role = "outcome", as_predictor = FALSE))
  }
  stop(sprintf("Unsupported outcome type `%s`.", outcome_type), call. = FALSE)
}

.prepare_predictor <- function(x, predictor_type, reference, levels, scores) {
  if (identical(predictor_type, "numeric")) {
    if (!is.numeric(x) || any(!is.finite(x))) {
      .stop_association("invalid_numeric_predictor",
                        "A numeric predictor must contain finite numeric values.")
    }
    return(as.numeric(x))
  }
  if (identical(predictor_type, "binary")) {
    return(.encode_binary(x, reference, "predictor"))
  }
  if (identical(predictor_type, "ordinal")) {
    return(.encode_ordinal(x, levels, scores, "predictor", TRUE))
  }
  stop(sprintf("Unsupported predictor type `%s`.", predictor_type), call. = FALSE)
}

.prepare_model_data <- function(data, variables) {
  model_data <- data.frame(
    .response = .prepare_outcome(
      data[[variables$response]], variables$outcome_type,
      variables$outcome_reference, variables$outcome_levels
    ),
    .predictor = .prepare_predictor(
      data[[variables$predictor]], variables$predictor_type,
      variables$predictor_reference, variables$predictor_levels,
      variables$predictor_scores
    ),
    check.names = FALSE
  )

  effective_covariates <- character()
  dropped_invariant_covariates <- character()
  if (length(variables$covariates)) {
    for (covariate in variables$covariates) {
      value <- data[[covariate]]
      if (!(is.numeric(value) || is.integer(value) || is.logical(value) ||
            is.factor(value))) {
        .stop_association(
          "invalid_covariate_type",
          sprintf(
            "Covariate `%s` must be numeric, logical, or a factor with explicit levels.",
            covariate
          )
        )
      }
      if (is.numeric(value) && any(!is.finite(value))) {
        .stop_association(
          "invalid_covariate_values",
          sprintf("Covariate `%s` contains non-finite values.", covariate)
        )
      }
      if (length(unique(value)) < 2L) {
        dropped_invariant_covariates <- c(
          dropped_invariant_covariates, covariate
        )
        next
      }
      if (is.factor(value)) {
        value <- droplevels(value)
      }
      effective_covariates <- c(effective_covariates, covariate)
      model_data[[paste0(
        ".covariate_", length(effective_covariates)
      )]] <- value
    }
  }

  if (!is.null(variables$offset)) {
    offset <- data[[variables$offset]]
    if (!is.numeric(offset) || any(!is.finite(offset))) {
      .stop_association("invalid_offset",
                        "A declared log-offset must contain finite numeric values.")
    }
    model_data$.offset <- as.numeric(offset)
  }
  attr(model_data, "effective_covariates") <- effective_covariates
  attr(model_data, "dropped_invariant_covariates") <-
    dropped_invariant_covariates
  model_data
}

.model_formula <- function(covariate_count, has_offset) {
  terms <- c(".predictor", if (covariate_count) {
    paste0(".covariate_", seq_len(covariate_count))
  })
  if (has_offset) {
    terms <- c(terms, "offset(.offset)")
  }
  stats::as.formula(paste(".response ~", paste(terms, collapse = " + ")))
}

.quote_formula_name <- function(x) {
  paste0("`", gsub("`", "\\`", x, fixed = TRUE), "`")
}

.display_formula <- function(variables, covariates = variables$covariates) {
  terms <- c(.quote_formula_name(variables$predictor),
             vapply(covariates, .quote_formula_name, character(1)))
  if (!is.null(variables$offset)) {
    terms <- c(terms, sprintf("offset(%s)",
                              .quote_formula_name(variables$offset)))
  }
  paste(.quote_formula_name(variables$response), "~",
        paste(terms, collapse = " + "))
}

.check_model_matrix <- function(formula, model_data) {
  terms <- stats::delete.response(stats::terms(formula))
  matrix <- stats::model.matrix(terms, data = model_data)
  if (ncol(matrix) < 2L || qr(matrix)$rank < ncol(matrix)) {
    .stop_association(
      "rank_deficient",
      "The association design matrix is rank deficient on complete cases."
    )
  }
  invisible(TRUE)
}

.fit_logistf_engine <- function(formula, model_data) {
  if (!requireNamespace("logistf", quietly = TRUE)) {
    stop("Package `logistf` is required for binary outcomes.", call. = FALSE)
  }
  maximum_iterations <- 1000L
  fit <- logistf::logistf(
    formula,
    data = model_data,
    pl = TRUE,
    control = logistf::logistf.control(maxit = maximum_iterations),
    plcontrol = logistf::logistpl.control(maxit = maximum_iterations)
  )
  coefficient_index <- match(".predictor", names(fit$coefficients))
  if (is.na(coefficient_index)) {
    stop("logistf did not return the requested predictor coefficient.",
         call. = FALSE)
  }
  test <- logistf::logistftest(fit, test = coefficient_index)
  likelihood_ratio <- max(0, 2 * unname(test$loglik[["full"]] -
                                         test$loglik[["null"]]))
  log_p <- stats::pchisq(likelihood_ratio, df = test$df,
                         lower.tail = FALSE, log.p = TRUE)

  estimate <- unname(fit$coefficients[[coefficient_index]])
  std_error <- sqrt(unname(fit$var[coefficient_index, coefficient_index]))
  full_iterations <- unname(fit$iter[["full"]])
  profile_iterations <- fit$pl.iter[coefficient_index, c("Lower", "Upper")]
  converged <- is.finite(estimate) && is.finite(std_error) &&
    all(is.finite(c(fit$ci.lower[[coefficient_index]],
                    fit$ci.upper[[coefficient_index]]))) &&
    is.finite(full_iterations) && full_iterations < maximum_iterations &&
    all(is.finite(profile_iterations)) &&
    all(profile_iterations < maximum_iterations)
  if (!converged) {
    stop("logistf failed its model or profile-likelihood convergence checks.",
         call. = FALSE)
  }

  list(
    fit = fit,
    engine = "logistf",
    estimate = estimate,
    std_error = std_error,
    conf_low = unname(fit$ci.lower[[coefficient_index]]),
    conf_high = unname(fit$ci.upper[[coefficient_index]]),
    statistic = sign(estimate) * sqrt(likelihood_ratio),
    df = unname(test$df),
    log_p = min(0, unname(log_p)),
    test_method = "penalized_likelihood_ratio",
    theta = NA_real_,
    converged = TRUE,
    iterations = as.integer(max(c(full_iterations, profile_iterations))),
    ordinal_levels = NA_character_
  )
}

.fit_lm_engine <- function(formula, model_data) {
  fit <- stats::lm(formula, data = model_data)
  summary_fit <- summary(fit)
  coefficients <- summary_fit$coefficients
  if (!".predictor" %in% row.names(coefficients)) {
    stop("lm did not return the requested predictor coefficient.", call. = FALSE)
  }
  estimate <- unname(coefficients[".predictor", "Estimate"])
  std_error <- unname(coefficients[".predictor", "Std. Error"])
  statistic <- unname(coefficients[".predictor", "t value"])
  df <- unname(fit$df.residual)
  critical <- stats::qt(0.975, df = df)
  list(
    fit = fit,
    engine = "lm",
    estimate = estimate,
    std_error = std_error,
    conf_low = estimate - critical * std_error,
    conf_high = estimate + critical * std_error,
    statistic = statistic,
    df = df,
    log_p = .log_two_sided_t(statistic, df),
    test_method = "t_test",
    theta = NA_real_,
    converged = isTRUE(all(is.finite(c(estimate, std_error, statistic)))),
    iterations = NA_integer_,
    ordinal_levels = NA_character_
  )
}

.fit_negative_binomial_engine <- function(formula, model_data) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package `MASS` is required for count outcomes.", call. = FALSE)
  }
  fit <- MASS::glm.nb(formula, data = model_data, link = log)
  coefficients <- summary(fit)$coefficients
  if (!".predictor" %in% row.names(coefficients)) {
    stop("glm.nb did not return the requested predictor coefficient.",
         call. = FALSE)
  }
  estimate <- unname(coefficients[".predictor", "Estimate"])
  std_error <- unname(coefficients[".predictor", "Std. Error"])
  statistic <- unname(coefficients[".predictor", "z value"])
  critical <- stats::qnorm(0.975)
  list(
    fit = fit,
    engine = "glm.nb",
    estimate = estimate,
    std_error = std_error,
    conf_low = estimate - critical * std_error,
    conf_high = estimate + critical * std_error,
    statistic = statistic,
    df = NA_real_,
    log_p = .log_two_sided_normal(statistic),
    test_method = "wald_z",
    theta = unname(fit$theta),
    # glm.nb's flag covers the GLM coefficients only; failed dispersion or
    # alternating iterations are recorded separately in th.warn.
    converged = isTRUE(fit$converged) && is.null(fit$th.warn) &&
      is.finite(fit$SE.theta) && fit$theta > 0 &&
      all(is.finite(c(estimate, std_error, statistic, fit$theta))),
    iterations = as.integer(fit$iter),
    ordinal_levels = NA_character_
  )
}

.fit_ordinal_engine <- function(formula, model_data) {
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    stop("Package `ordinal` is required for ordinal outcomes.", call. = FALSE)
  }
  fit <- ordinal::clm(formula, data = model_data, link = "logit")
  coefficients <- coef(summary(fit))
  if (!".predictor" %in% row.names(coefficients)) {
    stop("clm did not return the requested predictor coefficient.", call. = FALSE)
  }
  estimate <- unname(coefficients[".predictor", "Estimate"])
  std_error <- unname(coefficients[".predictor", "Std. Error"])
  statistic <- unname(coefficients[".predictor", "z value"])
  critical <- stats::qnorm(0.975)
  converged <- identical(fit$convergence$code, 0L) &&
    all(is.finite(c(estimate, std_error, statistic)))
  list(
    fit = fit,
    engine = "clm",
    estimate = estimate,
    std_error = std_error,
    conf_low = estimate - critical * std_error,
    conf_high = estimate + critical * std_error,
    statistic = statistic,
    df = NA_real_,
    log_p = .log_two_sided_normal(statistic),
    test_method = "wald_z",
    theta = NA_real_,
    converged = converged,
    iterations = as.integer(max(fit$niter)),
    ordinal_levels = paste(levels(model_data$.response), collapse = "|")
  )
}

.fit_engine <- function(outcome_type, formula, model_data) {
  .check_model_matrix(formula, model_data)
  switch(
    outcome_type,
    binary = .fit_logistf_engine(formula, model_data),
    continuous = .fit_lm_engine(formula, model_data),
    count = .fit_negative_binomial_engine(formula, model_data),
    ordinal = .fit_ordinal_engine(formula, model_data),
    stop(sprintf("Unsupported outcome type `%s`.", outcome_type), call. = FALSE)
  )
}
