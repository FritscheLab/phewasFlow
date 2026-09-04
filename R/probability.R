.log_two_sided_normal <- function(statistic) {
  if (is.na(statistic)) {
    return(NA_real_)
  }
  min(0, log(2) + stats::pnorm(abs(statistic), lower.tail = FALSE,
                               log.p = TRUE))
}

.log_two_sided_t <- function(statistic, df) {
  if (is.na(statistic) || is.na(df) || df <= 0) {
    return(NA_real_)
  }
  min(0, log(2) + stats::pt(abs(statistic), df = df, lower.tail = FALSE,
                            log.p = TRUE))
}

.probability_from_log <- function(log_probability) {
  ifelse(is.na(log_probability), NA_real_, exp(pmin(0, log_probability)))
}

.neg_log10_from_log <- function(log_probability) {
  ifelse(is.na(log_probability), NA_real_, -log_probability / log(10))
}

.format_log_probability <- function(log_probability, digits = 6L) {
  vapply(log_probability, function(value) {
    if (is.na(value)) {
      return(NA_character_)
    }
    if (is.infinite(value) && value < 0) {
      return("0")
    }
    if (!is.finite(value) || value > 0) {
      return(NA_character_)
    }

    probability <- exp(value)
    if (probability > 0) {
      return(formatC(probability, format = "g", digits = digits))
    }

    exponent <- floor(value / log(10))
    mantissa <- exp(value - exponent * log(10))
    paste0(
      sprintf(paste0("%.", max(0L, digits - 1L), "f"), mantissa),
      "e", sprintf("%+.0f", exponent)
    )
  }, character(1))
}

.adjust_log_bh <- function(log_p) {
  m <- length(log_p)
  if (!m) {
    return(numeric())
  }
  ordering <- order(log_p, decreasing = TRUE)
  reverse_order <- order(ordering)
  ranks_from_largest <- m:1L
  candidates <- log_p[ordering] + log(m) - log(ranks_from_largest)
  pmin(0, cummin(candidates))[reverse_order]
}

.apply_multiple_testing <- function(results, fdr_threshold = 0.05) {
  if (!is.data.frame(results)) {
    stop("`results` must be a data frame.", call. = FALSE)
  }
  if (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1L ||
      is.na(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold >= 1) {
    stop("`fdr_threshold` must be strictly between zero and one.", call. = FALSE)
  }
  if (!all(c("status", "log_p") %in% names(results))) {
    stop("Results must contain `status` and `log_p` columns.", call. = FALSE)
  }

  out <- as.data.frame(results, stringsAsFactors = FALSE, check.names = FALSE)
  n <- nrow(out)
  out$q_value <- rep(NA_real_, n)
  out$q_value_text <- rep(NA_character_, n)
  out$neg_log10_q <- rep(NA_real_, n)
  out$bonferroni_p <- rep(NA_real_, n)
  out$bonferroni_p_text <- rep(NA_character_, n)
  out$neg_log10_bonferroni_p <- rep(NA_real_, n)
  out$fdr_significant <- rep(NA, n)
  out$bonferroni_significant <- rep(NA, n)

  eligible <- which(out$status == "ok" & !is.na(out$log_p))
  if (!length(eligible)) {
    return(out)
  }
  if (any(out$log_p[eligible] > 0)) {
    stop("Successful results contain an invalid positive log p-value.",
         call. = FALSE)
  }

  log_p <- out$log_p[eligible]
  log_q <- .adjust_log_bh(log_p)
  log_bonferroni <- pmin(0, log_p + log(length(log_p)))

  out$q_value[eligible] <- .probability_from_log(log_q)
  out$q_value_text[eligible] <- .format_log_probability(log_q)
  out$neg_log10_q[eligible] <- .neg_log10_from_log(log_q)
  out$bonferroni_p[eligible] <- .probability_from_log(log_bonferroni)
  out$bonferroni_p_text[eligible] <- .format_log_probability(log_bonferroni)
  out$neg_log10_bonferroni_p[eligible] <-
    .neg_log10_from_log(log_bonferroni)
  out$fdr_significant[eligible] <- log_q <= log(fdr_threshold)
  out$bonferroni_significant[eligible] <-
    log_bonferroni <= log(fdr_threshold)
  out
}
