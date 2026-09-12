# Runnable with Rscript, or source this file to inspect the returned objects.
# All participant records and associations in this example are synthetic.

run_phewas_bridge <- function() {
  if (!requireNamespace("PheWAS", quietly = TRUE)) {
    stop("Install the optional PheWAS package; see this example's README.",
         call. = FALSE)
  }
  if (!requireNamespace("phewasFlow", quietly = TRUE)) {
    stop("Install phewasFlow before running this example.", call. = FALSE)
  }
  # Upstream addPhecodeInfo() looks up its pheinfo dataset on the search path.
  suppressPackageStartupMessages(library("PheWAS", character.only = TRUE))
  set.seed(20260912)
  n <- 400L
  sample <- data.frame(id = sprintf("S%04d", seq_len(n)))
  anchors <- data.frame(id = sample$id, prs = stats::rnorm(n))
  covariates <- data.frame(
    id = sample$id, age = stats::runif(n, 40, 80),
    sex = sample(c("Female", "Male"), n, replace = TRUE)
  )

  # Two distinct dates can establish a case; a single date remains uncertain.
  # These simulated ICD9CM codes exercise PheWAS mapping and rollup behavior.
  codes <- c("250.00", "401.9", "493.90")
  records <- lapply(seq_along(codes), function(j) {
    case <- stats::runif(n) < stats::plogis(-0.7 + 0.25 * anchors$prs)
    count <- ifelse(case, 2L, 0L)
    count[seq_len(8L)] <- 1L
    index <- rep(seq_len(n), count)
    data.frame(
      id = sample$id[index], vocabulary_id = "ICD9CM", code = codes[[j]],
      date = as.Date("2020-01-01") + sequence(count)
    )
  })
  records <- do.call(rbind, records)

  phenotypes <- PheWAS::createPhenotypes(
    records, min.code.count = 2L,
    id.sex = covariates[c("id", "sex")],
    full.population.ids = sample$id
  )
  phenotype_ids <- setdiff(names(phenotypes), "id")
  stopifnot(length(phenotype_ids) > 0L)
  stopifnot(all(vapply(phenotypes[phenotype_ids], is.logical, logical(1))))

  metadata <- as.data.frame(PheWAS::addPhecodeInfo(
    data.frame(phenotype = phenotype_ids),
    groupnums = TRUE, groupcolors = TRUE
  ))
  # addPhecodeInfo() uses an inner join. Never silently lose unknown codes.
  stopifnot(!anyDuplicated(metadata$phenotype))
  stopifnot(setequal(metadata$phenotype, phenotype_ids))
  metadata <- metadata[match(phenotype_ids, metadata$phenotype),
                       c("phenotype", "description", "group", "groupnum", "color")]
  metadata$variable_type <- "binary"
  metadata$outcome_type <- "binary"
  metadata$reference <- "FALSE"

  covariates$sex <- factor(covariates$sex, levels = c("Female", "Male"))
  analysis <- phewasFlow::assemble_phewas_data(
    sample, phenotypes = phenotypes, anchors = anchors,
    covariates = covariates, id = "id"
  )
  spec <- phewasFlow::phewas_spec(
    analysis_id = "synthetic_phewas_bridge", anchor = "prs",
    direction = "phenotypes_as_outcomes", phenotypes = metadata,
    covariates = c("age", "sex"), anchor_type = "numeric",
    anchor_transform = "zscore",
    eligibility = list(
      min_n = 100L, min_cases = 10L, min_controls = 10L,
      min_outcome_levels = 3L
    ),
    fdr_threshold = 0.05
  )
  result <- phewasFlow::run_phewas(analysis, spec)
  list(phenotypes = phenotypes, metadata = metadata, analysis = analysis,
       spec = spec, result = result)
}

if (sys.nframe() == 0L) {
  example <- run_phewas_bridge()
  print(example$result)
}
