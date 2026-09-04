#' Generate a small synthetic PheWAS data set
#'
#' The returned data contain no real participants. They are designed for
#' examples and tests of every supported outcome family and both scan
#' directions.
#'
#' @param n Number of synthetic participants.
#' @param seed Random-number seed.
#'
#' @return A list of `data.table`s containing `sample`, `phenotypes`,
#'   `anchors`, `covariates`, the joined `data`, and forward/reverse phenotype
#'   metadata.
#'
#' @examples
#' example <- phewas_example_data(n = 200)
#' names(example)
#' head(example$data)
#' @export
phewas_example_data <- function(n = 400L, seed = 2026L) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 100L || n > .Machine$integer.max || n != trunc(n)) {
    stop("`n` must be one integer of at least 100.", call. = FALSE)
  }
  n <- as.integer(n)
  set.seed(seed)

  id <- sprintf("person_%04d", seq_len(n))
  age <- round(stats::runif(n, 30, 80), 1)
  sex <- factor(sample(c("female", "male"), n, replace = TRUE))
  pgs <- as.numeric(scale(stats::rnorm(n)))

  pheno_binary <- stats::rbinom(
    n, 1, stats::plogis(-1.1 + 0.65 * pgs + 0.012 * (age - 50))
  )
  pheno_continuous <- 0.75 * pgs + 0.018 * (age - 50) + stats::rnorm(n)
  count_mean <- exp(0.5 + 0.30 * pgs + 0.01 * (age - 50))
  pheno_count <- stats::rnbinom(n, mu = count_mean, size = 1.2)
  ordinal_latent <- 0.55 * pgs + 0.015 * (age - 50) + stats::rnorm(n)
  pheno_ordinal <- cut(
    ordinal_latent,
    breaks = c(-Inf, -0.45, 0.45, Inf),
    labels = c("low", "medium", "high"),
    ordered_result = TRUE
  )

  endpoint_continuous <-
    0.9 * pheno_binary + 0.35 * pheno_continuous +
    0.015 * (age - 50) + stats::rnorm(n)
  endpoint_binary <- stats::rbinom(
    n, 1, stats::plogis(-1 + 0.8 * pheno_binary + 0.2 * pheno_continuous)
  )
  endpoint_count <- stats::rnbinom(
    n,
    mu = exp(0.4 + 0.35 * pheno_binary + 0.12 * pheno_continuous),
    size = 1.5
  )
  endpoint_ordinal <- cut(
    0.7 * pheno_binary + 0.25 * pheno_continuous + stats::rnorm(n),
    breaks = c(-Inf, -0.25, 0.75, Inf),
    labels = c("mild", "moderate", "severe"),
    ordered_result = TRUE
  )

  sample_table <- data.frame(id = id, stringsAsFactors = FALSE)
  phenotype_table <- data.frame(
    id = id,
    pheno_binary = pheno_binary,
    pheno_continuous = pheno_continuous,
    pheno_count = pheno_count,
    pheno_ordinal = pheno_ordinal,
    stringsAsFactors = FALSE
  )
  anchor_table <- data.frame(
    id = id,
    pgs = pgs,
    endpoint_continuous = endpoint_continuous,
    endpoint_binary = endpoint_binary,
    endpoint_count = endpoint_count,
    endpoint_ordinal = endpoint_ordinal,
    stringsAsFactors = FALSE
  )
  covariate_table <- data.frame(
    id = id,
    age = age,
    sex = sex,
    stringsAsFactors = FALSE
  )
  joined <- cbind(
    sample_table,
    phenotype_table[-1L],
    anchor_table[-1L],
    covariate_table[-1L]
  )

  forward_metadata <- data.frame(
    phenotype = c(
      "pheno_binary", "pheno_continuous", "pheno_count", "pheno_ordinal"
    ),
    description = c(
      "Binary phenotype", "Continuous phenotype", "Count phenotype",
      "Ordinal phenotype"
    ),
    group = c("Diagnoses", "Measurements", "Utilization", "Symptoms"),
    groupnum = 1:4,
    color = c("#3B82F6", "#10B981", "#F59E0B", "#8B5CF6"),
    variable_type = c("binary", "numeric", "count", "ordinal"),
    outcome_type = c("binary", "continuous", "count", "ordinal"),
    reference = c("0", NA, NA, NA),
    stringsAsFactors = FALSE
  )
  forward_metadata$levels <- I(list(
    c("0", "1"), character(), character(), c("low", "medium", "high")
  ))
  forward_metadata$scores <- I(list(
    numeric(), numeric(), numeric(), numeric()
  ))
  forward_metadata$offset <- NA_character_

  reverse_metadata <- forward_metadata[c(1L, 2L, 4L), setdiff(
    names(forward_metadata), c("outcome_type", "offset")
  )]
  reverse_metadata$scores[[3L]] <- c(1, 2, 3)

  list(
    sample = data.table::as.data.table(sample_table),
    phenotypes = data.table::as.data.table(phenotype_table),
    anchors = data.table::as.data.table(anchor_table),
    covariates = data.table::as.data.table(covariate_table),
    data = data.table::as.data.table(joined),
    forward_metadata = data.table::as.data.table(forward_metadata),
    reverse_metadata = data.table::as.data.table(reverse_metadata)
  )
}
